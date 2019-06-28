# screen
# qrsh -pe smp 20 -l mem_requested=20G
# cd /share/ScratchGeneral/briglo/scRNA/chrcha
# module load briglo/R/3.6.0
# R
# a lot of these things should be simpler with all the functions that ive written, but heres what i acutally did


library(Seurat)
library(reticulate)
#library(scGPS)
library(locfit) #required but not dependency
library(e1071)
use_python("/share/ClusterShare/software/contrib/briglo/miniconda3/envs/magic/bin/python")

source("source/seuratFunctions.r")
source("source/scGPSfunctions.r")
load("data/cell.cyclegenes.rdata")

##make a list of seurat objects
snam<-dir("data")
snam<-snam[-grep("FB|fnam",snam)]

id<- lapply(snam, function(x) {
    ddir<-paste0(getwd(),"/data/",x,"/outs/filtered_feature_bc_matrix/")
#system(paste0("gunzip ",ddir,"*"))
rd <- Read10X(ddir)
rd<-rd[grepl("hg19",rownames(rd)),]
rownames(rd)<-gsub("hg19_","",rownames(rd))
return(CreateSeuratObject(counts=rd,project=gsub('^M',"",x), min.cells=3,min.features=200))
})
#save(id,file="r_objects/190531_allRawSeurat.rdata")

#the fancy new integration method
kill<-c(1,17,22) #retaining just the 4 by 4 experiment
anchors <- FindIntegrationAnchors(object.list = id[-kill], dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30) #this takes ages
integrated@active.assay="RNA"
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
integrated<-subset(integrated, subset = percent.mt < 10)
integrated<-CellCycleScoring(object = integrated, s.features = toupper(s.genes), g2m.features = toupper(g2m.genes), set.ident = FALSE)

#the usual Seurat thing
integrated@active.assay="integrated"
integrated<-SCTransform(integrated, vars.to.regress = c("percent.mt","nFeature_RNA","S.Score","G2M.Score"), verbose = FALSE)
integrated<-RunPCA(integrated, verbose = FALSE)
ElbowPlot(integrated,ndims=50) ; dev.off()
integrated<-RunUMAP(integrated, dims = 1:30, verbose = T)
#DimPlot(integrated, reduction = "umap", pt.size = 0.1,group.by="Phase") ; dev.off()
integrated<-FindNeighbors(integrated, dims = 1:20, verbose = FALSE)
integrated<-FindClusters(integrated, verbose = FALSE)
integrated<-RunTSNE(integrated, dims.use=1:17, do.fast=T)
#save(integrated,file="r_objects/190606_4by4_IntegratedSeurat.rdata")

####this is to map between an old manual merge that used experiment to annotate barcodes that i used for magic
# load("190504_magicALLgenes.rdata")
# comid<-data.frame(do.call(rbind,lapply(strsplit(rownames(integrated@meta.data), "_"), function(X) c(X[1],X[2]))))
# comid$X3<-gsub("_","\\.",integrated@meta.data$orig.ident)
# rownames(comid)<-paste0(comid$X3,"_",comid$X1)

# magid<-data.frame(do.call(rbind,lapply(strsplit(rownames(magicAll), "_"), function(X) c(gsub("\\.","_",X[1]),X[2]))),row.names=rownames(m
# agicAll))

# m<-match(rownames(comid),rownames(magid))
# n<-match(rownames(integrated@assays$SCT@data),colnames(magicAll))

# nd<-t(magicAll[m,])
# colnames(nd)<-rownames(integrated@meta.data)
# nd<-nd[n,]
# rownames(nd)<-rownames(integrated@assays$SCT@data)
# nd[is.na(nd)]=0
# nd<-Matrix(nd,sparse=T)
#####################################

###redo SCGPS slice by locaton
x<-data.frame(do.call(rbind, lapply(strsplit(integrated@meta.data$orig.ident,"_"), function(x) x[1:2])))
colnames(x)<-c("id","location")
x$id<-gsub("^M","",x$id)
x$location<-gsub("LN","Lymph",x$location)
integrated@meta.data<-data.frame(cbind(integrated@meta.data,x))
cu<-split(rownames(integrated@meta.data), integrated@meta.data$location)
comb<-lapply(cu, function(X) {
tmp<-subset(integrated,cells=X)
tmp<-RunUMAP(tmp, dims = 1:30, verbose = T)
tmp<-FindNeighbors(tmp, dims = 1:30, verbose = FALSE)
 for (res in seq(0,0.4,.05)) tmp<-FindClusters(tmp, verbose = FALSE,resolution=res,force.recalc=T)
#tmp<-RunTSNE(tmp, dims.use=1:17, do.fast=T)
return(tmp)
})
save(comb,file="BYLOCATION.rdata")

#making scGPS objects
#manually assigning cluster based on clustertree vis of comb object, basically before grps start switching
comb$Liver<-SetIdent(comb$Liver,value=comb$Liver@meta.data$SCT_snn_res.0.15)
comb$Lung<-SetIdent(comb$Lung,value=comb$Lung@meta.data$SCT_snn_res.0.25)
comb$Lymph<-SetIdent(comb$Lymph,value=comb$Lymph@meta.data$SCT_snn_res.0.25)
comb$PT<-SetIdent(comb$PT,value=comb$PT@meta.data$SCT_snn_res.0.15)
#this bit is simpler now because of function seurat2scGPS (i hope)
prim<-new_scGPS_object(ExpressionMatrix=as(comb$PT@assays$RNA@counts,'matrix'),  CellMetadata = data.frame(X=comb$PT@active.ident,row.names=rownames(comb$PT@meta.data)),GeneMetadata=data.frame(rownames(comb$PT@assays$RNA@counts),row.names=rownames(comb$PT@assays$RNA@counts)))
mets<-lapply(comb[-4], function(y){
return(new_scGPS_object(ExpressionMatrix=as(y@assays$RNA@counts,'matrix'),  CellMetadata = data.frame(X=y@active.ident,row.names=rownames(y@meta.data)),GeneMetadata=data.frame(rownames(y@assays$RNA@counts),row.names=rownames(y@assays$RNA@counts))))
})
save(prim,mets,file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190614_scGPSinputData.rdata")

#making gene list for primary tumour
cluster_prim <- as.numeric(as.character(colData(prim)[,1]))+1
allID <- unique(cluster_prim)

    DEgenes<-find_markers(expression_matrix=assay(prim)[comb$PT@assays$SCT@var.features,], 
    cluster = as.numeric(colData(prim)[,"X"]), 
    selected_cluster=unique(as.numeric(colData(prim)[,"X"]) ))

save(DEgenes,file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190613_DEgenesByPTcluster")

##making gene list for lung mets #not done as yet
cluster_lung <- as.numeric(as.character(colData(mets$Lung)[,1]))+1
allID <- unique(cluster_lung)
DEgenes<-find_markers(expression_matrix=assay(mets$Lung)[comb$Lung@assays$SCT@var.features,], 
    cluster = as.numeric(colData(mets$Lung)[,"X"]), 
    selected_cluster=allID) 
save(DEgenes,file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190613_DEgenesByLungCluster.rdata")

# ###########in bash
# module load briglo/R/3.6.0
# qsub -V -cwd -pe smp 20 -l mem_requested=20G -N scGPS -b y -j y -m ae -M b.gloss@garvan.org.au Rscript source/190613_runscGPS.r #the long way, worked though but took 24 hours
################


##analysing outputs
load("r_objects/190614_scGPS_byLocation.rdata")
lpt<-mklasso(scGPS_location)
 pp3<-plotLasso(lpt)
save(lpt,pp3,file="r_objects/190614_scGPS_networkDiag.rdata")

#plot scgps


#integrating the location based data into master dataset
md<-lapply(comb, function(X) return(data.frame(cbind(X@reductions$umap@cell.embeddings,X@reductions$tsne@cell.embeddings,X@reductions$pca@cell.embeddings[,1:3],X@meta.data[,c(12,13,17:19)]))))
names(md)<-NULL
mmd<-data.frame(do.call(rbind,md))
colnames(mmd)<-paste0("LocSpec_",colnames(mmd))
load("r_objects/190606_4by4_IntegratedSeurat.rdata")
for (res in c(seq(0,0.4,.05),.5,1)) integrated<-FindClusters(integrated, verbose = FALSE,resolution=res,force.recalc=T)
integrated@meta.data<-data.frame(cbind(integrated@meta.data,integrated@reductions$umap@cell.embeddings,integrated@reductions$tsne@cell.embeddings,integrated@reductions$pca@cell.embeddings[,1:3], mmd[rownames(integrated@meta.data),]))

#adding magic and PHATE and PCA
colnames(data_MAGIC$result)<-paste0("magicEMT_",colnames(data_MAGIC$result))
colnames(data_MAGIC_PCA$result)<-paste0("magicEMT_",colnames(data_MAGIC_PCA$result))
x<-data.frame(integrated@meta.data[,c("orig.ident","tSNE_1")])
x$barcode<-unlist(lapply(strsplit(rownames(x), "_"), function(x) x[1]))
x$newRN<-paste0(gsub("_",".",x$orig.ident),"_",x$barcode)
m<-match(x$newRN, gsub("^M","",rownames(data_MAGIC$result)))
integrated@meta.data<-data.frame(cbind(integrated@meta.data,data_MAGIC$result[m,],data_MAGIC_PCA$result[m,1:3],data_PHATE$embedding[m,]))

md<-integrated@meta.data#(file="r_objects/190606_4by4_IntegratedSeurat_metadata.rdata")
save(md,file="r_objects/190606_4by4_IntegratedSeurat_metadata.rdata")

save(integrated,file="r_objects/190606_4by4_IntegratedSeurat.rdata")



#AR ANALYSIS originally done for r_objects/mergeSeurat_AllSample1000Random_varNemtNbatchTSNEandMonocleAndSCORE.rdata")
#make distinctions based on AR and ZEB1 results
integrated@meta.data$isDP<-ifelse( integrated@meta.data$LocSpec_id=="43978", 
ifelse(integrated@meta.data$magicEMT_AR>0.0035 & integrated@meta.data$magicEMT_ZEB1>.65, "DP", ifelse(integrated@meta.data$magicEMT_AR>0.0035,"ARonly","neither")),ifelse(integrated@meta.data$LocSpec_id=="43979",
ifelse(integrated@meta.data$magicEMT_AR<0.005 & integrated@meta.data$magicEMT_ZEB1>.8,"superZeb", ifelse(integrated@meta.data$magicEMT_AR>0.005 & integrated@meta.data$magicEMT_ZEB1>.65,"DP","neither")),
ifelse(integrated@meta.data$magicEMT_AR>0.0035 , "DP","neither")))
ggplot(integrated@meta.data,aes(x=magicEMT_ZEB1,y=magicEMT_AR,color=isDP)) + geom_point(alpha=.5) + facet_wrap(~LocSpec_id) ; dev.off()

##useful genelist
library(cmapR)
library(dplyr)
x<-parse.gmx("data/Androgen_response_genesets.gmx")
allAR<-unique(unlist(lapply(x,function(X) return(X$entry))))

integrated<-SetIdent(integrated,value=integrated@meta.data$isDP)
m<-FindAllMarkers(integrated)
m %>% filter(p_val_adj<0.05) -> comb
save(m,comb,file="newAllMarkers.rdata")


hasEntrez<-getEntrez(integrated,'all')
scomb<-split(comb,comb$cluster)
sscomb<-lapply(scomb,function(x) splitDirection(x,getEntrezObj=hasEntrez))
oi<-lapply(sscomb,reactomeClusts)


