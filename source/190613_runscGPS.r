###scgps
library(scGPS)
library(locfit) #required but not dependency
library(e1071) #required but not dependency
library(Seurat)

load("/share/ScratchGeneral/briglo/scRNA/chrcha/BYLOCATION.rdata") #assigning cluster based on clustertree
comb$Liver<-SetIdent(comb$Liver,value=comb$Liver@meta.data$SCT_snn_res.0.15)
comb$Lung<-SetIdent(comb$Lung,value=comb$Lung@meta.data$SCT_snn_res.0.25)
comb$Lymph<-SetIdent(comb$Lymph,value=comb$Lymph@meta.data$SCT_snn_res.0.25)
comb$PT<-SetIdent(comb$PT,value=comb$PT@meta.data$SCT_snn_res.0.15)

load("/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190613_DEgenesByPTcluster")

# for (i in 1:4) UMAPPlot(comb[[i]]) 
# dev.off()
#devtools::install_github("IMB-Computational-Genomics-Lab/scGPS")
#require RcppParallel which jose got me to being able to compile by typing "scl enable devtoolset-2 bash" into terminal before firing up R to install the package
prim<-new_scGPS_object(ExpressionMatrix=as(comb$PT@assays$RNA@counts,'matrix'),  CellMetadata = data.frame(X=comb$PT@active.ident,row.names=rownames(comb$PT@meta.data)),GeneMetadata=data.frame(rownames(comb$PT@assays$RNA@counts),row.names=rownames(comb$PT@assays$RNA@counts)))
cluster_prim <- as.numeric(as.character(colData(prim)[,1]))+1
allID <- unique(cluster_prim)

#     DEgenes<-find_markers(expression_matrix=assay(prim)[comb$PT@assays$SCT@var.features,], 
#     cluster = as.numeric(colData(prim)[,"X"]), 
#     selected_cluster=unique(as.numeric(colData(prim)[,"X"]) ))

# save(DEgenes,file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190613_DEgenesByPTcluster")

scGPS_location<-lapply(comb[-4], function(y){
met<-new_scGPS_object(ExpressionMatrix=as(y@assays$RNA@counts,'matrix'),  CellMetadata = data.frame(X=y@active.ident,row.names=rownames(y@meta.data)),GeneMetadata=data.frame(rownames(y@assays$RNA@counts),row.names=rownames(y@assays$RNA@counts)))
cluster_met <- as.numeric(as.character(colData(met)[,1]))+1
lapply(as.numeric(allID), function(c_selectID){
    message("doing scGPS for ",allID)
    genes=DEgenes[[grep(c_selectID,names(DEgenes))]]$id[DEgenes[[grep(c_selectID,names(DEgenes))]]$padj<0.01 & abs(DEgenes[[grep(c_selectID,names(DEgenes))]]$AdjustedLogFC)>.5]
return(bootstrap(nboots = 3, mixedpop1 = prim, 
    mixedpop2 = met, genes = genes, c_selectID  = c_selectID,
    listData = list(), cluster_mixedpop1 = cluster_prim, 
    cluster_mixedpop2 = cluster_met, trainset_ratio = 0.7))
})
})

 save(scGPS_location, file="/share/ScratchGeneral/briglo/scRNA/chrcha/r_objects/190614_scGPS_byLocation.rdata")
