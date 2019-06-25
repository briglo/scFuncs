screen
qrsh -pe smp 12 -l h_vmem=10G
source ~/.profile
module load briglo/R/3.4.2
cd /share/ScratchGeneral/briglo/scRNA/chrcha
R

library(monocle)


#Tp remove non-relevant clusters - not part of the  differentiating cells

pd <- AnnotatedDataFrame(alldat@meta.data)

gene.names <- as.array(unlist(ALLs15@raw.data@Dimnames[1]))
fd <- data.frame(gene.names)
rownames(fd) <- gene.names
colnames(fd) <- c("gene_short_name")
fd <- AnnotatedDataFrame(fd)

# To create the monocle object
HSMM <- newCellDataSet(alldat@raw.data[, rownames(pd)], phenoData = pd, featureData = fd, lowerDetectionLimit=1, expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM)

ordering_genes <- subset(disp_table,mean_expression >= 0.05 &dispersion_empirical >= 1 * dispersion_fit)$gene_id

#Now use these genes to order cells and plot the cell trajectory
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2)
HSMM <- orderCells(HSMM, reverse=FALSE)
plot_cell_trajectory(HSMM, color_by="orig.ident")
dev.off()
 save(HSMM,file="r_objects/190309_VarGeneMonocle.rdata")

## reorganising baset on EMT genes
library(msigdbr)
h_df = msigdbr(species = "Homo sapiens", category = "H")
h_df %>% filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% select(gene_symbol) -> emtgenes
ordering_genes <- emtgenes$gene_symbol
HSMM_emt <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM_emt)
HSMM_emt <- reduceDimension(HSMM_emt, max_components = 2, method = 'DDRTree')
HSMM_emt <- orderCells(HSMM_emt)
plot_cell_trajectory(HSMM_emt, color_by = "orig.ident")
save(HSMM_emt,file="r_objects/190309_EMTGeneMonocle.rdata")


###done april 2019 for three
# d1
screen
qrsh -pe smp 16 
source .profile 
module load briglo/R/3.4.2
cd /share/ScratchGeneral/briglo/scRNA/chrcha
# R  
library(SingleCellExperiment)
library(monocle)
library(Seurat)
fnam<-dir("r_objects/",pattern="id_")
id<-gsub("id_|.rdata","",fnam)
load("r_objects/190322_clean_norm_human_sce_SCORE_GPS_minus43978.rdata")

for(n in id[1]) {
    load(paste0("r_objects/id_",n,".rdata"))
SCOREclust<-data.frame(do.call(rbind,lapply(GPSdat[grep(n,names(GPSdat))], function(x){
colnames(x)<-colData(x)$Barcode
return(data.frame(row.names=gsub("-[12]","",rownames(colData(x))),SCOREclust=as.numeric(colData(x)$Cluster)))
})))
rownames(SCOREclust)<-gsub("\\.",'fuck',rownames(SCOREclust))
rownames(SCOREclust)<-gsub("_",'\\.',rownames(SCOREclust))
rownames(SCOREclust)<-gsub("fuck",'_',rownames(SCOREclust))
SCOREclust<-SCOREclust[rownames(alldat@meta.data),]
SCOREclust[is.na(SCOREclust)]=0
alldat@meta.data<-cbind(alldat@meta.data,SCOREclust)

#Tp remove non-relevant clusters - not part of the  differentiating cells
###subsetted ALLs15 by epicells
clusters <- alldat@meta.data[,c("varRes.0.3")]
cell.names <-row.names(alldat@meta.data)
clusters <- droplevels(clusters)

pd <- data.frame(clusters)
rownames(pd) <- cell.names
colnames(pd) <- c("SEURAT")
pd <- AnnotatedDataFrame(pd)

gene.names <- as.array(unlist(alldat@raw.data@Dimnames[1]))
fd <- data.frame(gene.names)
rownames(fd) <- gene.names
colnames(fd) <- c("gene_short_name")
fd <- AnnotatedDataFrame(fd)

# To create the monocle object
HSMM <- newCellDataSet(as(as.matrix(alldat@raw.data[, cell.names]), "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit=1, expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
#estimateDispersions failed locally, moved to cluster

HSMM <- estimateDispersions(HSMM) 
disp_table <- dispersionTable(HSMM)
table(alldat@meta.data$orig.ident)

ordering_genes <- subset(disp_table,mean_expression >= 0.05 &dispersion_empirical >= 1 * dispersion_fit)$gene_id
library(msigdbr)
h_df = msigdbr(species = "Homo sapiens", category = "H")
h_df %>% filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% select(gene_symbol) -> emtgenes

genList<-list("vargen"=ordering_genes,'emtgen'=emtgenes$gene_symbol)
source("source/seuratFunctions.r")
mon_run<-lapply(genList, function(x) {
    tmp<-plotTraj(monoclObj=HSMM,genevec=x)
    return(data.frame(DDRTREE=t(monocle::reducedDimS(tmp)),as(phenoData(tmp)[,3:4],"data.frame")))
    })

save(mon_run,file="r_objects/allSampleMerge_monocleVarEMT.rdata")

#Now use these genes to order cells and plot the cell trajectory
HSMM <- setOrderingFilter(HSMM, ordering_genes)
#plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2)
HSMM <- orderCells(HSMM, reverse=FALSE) #this takes forever probably coz lots of cells? Failed for 73 need to redp with more than 8 cores lol
# plot_cell_trajectory(HSMM, color_by="SEURAT")
# plot_cell_trajectory(HSMM, color_by="SCORE")
# plot_cell_trajectory(HSMM, color_by="State") + facet_wrap (~orig.ident)
# dev.off()

# save(HSMM,file="r_objects/43975_HSMM_monocleDiffgen.rdata")

kdat<-data.frame(DDRTREE=t(monocle::reducedDimS(HSMM)),as(phenoData(HSMM)[,3:4],"data.frame"))
kdat<-kdat[rownames(alldat@meta.data),]
trimdat@meta.data<-cbind(trimdat@meta.data,kdat)
save(trimdat,file=)
save(alldat, file=paste0("r_objects/id_",n,".rdata"))
rm(alldat,HSSM)
gc()
}

## reorganising baset on EMT genes
library(msigdbr)
h_df = msigdbr(species = "Homo sapiens", category = "H")
h_df %>% filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% select(gene_symbol) -> emtgenes
ordering_genes <- emtgenes$gene_symbol
HSMM <- setOrderingFilter(HSMM, ordering_genes)
#plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

kdat<-data.frame(DDRTREE_emt=t(monocle::reducedDimS(HSMM)),PSEUDO_emt=as(phenoData(HSMM)[,3:4],"data.frame"))
kdat<-kdat[rownames(alldat@meta.data),]
alldat@meta.data<-cbind(alldat@meta.data,kdat)
save(alldat, file="r_objects/mergeSeurat_AllSample1000Random_emtmon.rdata")
###up to here apr 16


###do combined object
load("r_objects/mergeSeurat_AllSample1000Random_varNemtTSNE.rdata") 
library(Seurat)
library(monocle)
clusters <- trimdat@meta.data[,c("varRes.0.3")]
cell.names <-row.names(trimdat@meta.data)
clusters <- droplevels(clusters)

pd <- data.frame(clusters)
rownames(pd) <- cell.names
colnames(pd) <- c("SEURAT")
pd <- AnnotatedDataFrame(pd)

gene.names <- as.array(unlist(trimdat@raw.data@Dimnames[1]))
fd <- data.frame(gene.names)
rownames(fd) <- gene.names
colnames(fd) <- c("gene_short_name")
fd <- AnnotatedDataFrame(fd)

# To create the monocle object
HSMM <- newCellDataSet(as(as.matrix(trimdat@raw.data[, cell.names]), "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit=1, expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM) 


save(HSMM, file='r_objects/190429_subset1000allMonocle.rdata')


disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,mean_expression >= 0.05 &dispersion_empirical >= 1 * dispersion_fit)$gene_id
library(msigdbr)
h_df = msigdbr(species = "Homo sapiens", category = "H")
h_df %>% filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% select(gene_symbol) -> emtgenes

HSMM_var <- setOrderingFilter(HSMM, ordering_genes)
HSMM_var <- reduceDimension(HSMM_var, max_components=2)
HSMM_var <- orderCells(HSMM_var, reverse=FALSE)
save(HSMM_var, file='r_objects/190429_subset1000allMonocle_var.rdata')

HSMM_emt <- setOrderingFilter(HSMM, emtgenes$gene_symbol)
HSMM_emt <- reduceDimension(HSMM_emt, max_components=2)
HSMM_emt <- orderCells(HSMM_emt, reverse=FALSE)
save(HSMM_emt, file='r_objects/190429_subset1000allMonocle_var.rdata')
