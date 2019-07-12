#' make2sampReactomePipe
#'
#' wrapper for seurat to reactome but specifying 
#'
#' @param mkEntObj an output from  mkEnt or an element of splotDirection
#' @param seuratObj the Seurat object
#' @param groupBy_1 the metadata to group by (like cell type)
#' @param groupID_1 a specific value of groupBy, required
#' @param groupBy_2 the metadata to group by (like a result from annotateByCuts)
#' @param groupID_2 two values of groupBy_2, required
#' @param splitBy the metadata to split by (like experiment)
#' @param splitID a subset of splitBy if required, default NULL (all)
#'
#' @return list with markers and a list of compareClusterObject for plotting
#'
#' @examples
#' enrichment<-make2sampReactomePipe(integrated,groupBy_1="orig_final.ident",groupID_1="Cancer_2",groupBy_2="PFN1_GZMA_GZMB",groupID_2=c("hi_Low_Low","Low_Low_Low"),splitBy="orig.ident",splitID=NULL)
#' for i in 1:length(enrichment$CP_result)) dotplot(enrichment$CP_result[[i]]) + ggtitle(names(enrichment$CP_result)[i])
#'
#' @export	 
make2sampReactomePipe<-function(integrated,groupBy_1,groupID_1,groupBy_2,groupID_2,splitBy,splitID=NULL){
    require(Seurat)
    require(clusterProfiler)
    require(dplyr)
    message('making markers, might take a while')
seuratObj<-SetIdent(seuratObj, value=paste0(seuratObj@meta.data[,groupBy_1],"_",seuratObj@meta.data[,groupBy_2]))
if(is.null(splitID)) {
    cellvec<-split(rownames(seuratObj@meta.data),integrated@meta.data[,splitBy])
    mm<-lapply(cellvec, function(X){
        so<-subset(seuratObj,cells=X)
mar<-FindMarkers(so,ident.1=paste0(groupID_1,"_",groupID_2[1]),ident.2=paste0(groupID_1,"_",groupID_2[2]))
return(mar)
    })
    
    } else {
        cell.cull<-rownames(seuratObj@meta.data)[seuratObj@meta.data[,splitBy] %in% splitID]
so<-subset(seuratObj, cells=cell.cull)
cellvec<-split(rownames(so@meta.data),so@meta.data[,splitBy])
    mm<-lapply(cellvec, function(X){
        so<-subset(seuratObj,cells=X)
mar<-FindMarkers(so,ident.1=paste0(groupID_1,"_",groupID_2[1]),ident.2=paste0(groupID_1,"_",groupID_2[2]))
return(mar)
    })
    }
   hasEntrez<-getEntrez(seuratObj,'all')
return(list('hasEntrez'=hasEntrez, 'markerLists'=mm))
}


#' makeReactomeForMarkers
#'
#' wrapper for output of FindClusters to reactome but will crap out if no significance
#'
#' @param findMarkersOut a typical output from a FindMarkers seurat operation
#' @param hasEntrezObj an output from getEntrez on the seurat object
#'
#' @return list with markers and a list of compareClusterObject for plotting
#'
#' @examples
#' enrichment<-make2sampReactomePipe(integrated,groupBy_1="orig_final.ident",groupID_1="Cancer_2",groupBy_2="PFN1_GZMA_GZMB",groupID_2=c("hi_Low_Low","Low_Low_Low"),splitBy="orig.ident",splitID=NULL)
#' out<-makeReactomeForMarkers(enrichment$markerLists[[1]],enrichment$hasEntrez)
#'
#' @export

makeReactomeForMarkers<-function(findMarkersOut,hasEntrezObj){
    findMarkersOut$gene<-rownames(findMarkersOut)
    findMarkersOut %>% filter(p_val_adj<0.05) -> comb
    scomb<-split(comb$gene,ifelse(comb$avg_logFC>0,"up","down"))
    en<-mkEnt(scomb,hasEntrezObj)
    oi<-reactomeClusts(en)
return(oi)
    }


library(Seurat)
seuratObj=integrated
groupBy_1="orig_final.ident"
groupID_1="Cancer_2"
groupBy_2="PFN1_GZMA_GZMB"
groupID_2=c("hi_Low_Low","Low_Low_Low")
splitBy="orig.ident"
splitID=NULL