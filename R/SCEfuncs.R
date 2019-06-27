#' trimNscent
#'
#' takes a singleCellExperiment object and runs landSCENT
#'
#' @param sce_obj a SingleCellExperiment object
#' @param id a name for the sample
#'
#' @return invisible list of plot objects
#'
#' @examples
#' NULL
#'
#' @export	 

trimNscent<-function(sce_obj,id){
    message('libraries')
    require(LandSCENT)
require(SingleCellExperiment)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(scater)
message('tweaking expression values for SCENT')
sizeFactors(sce_obj) <- librarySizeFactors(sce_obj)
sce_obj<- normalize(sce_obj, log_exprs_offset = 1.1)
rowData(sce_obj)$anno.v <- mapIds(org.Hs.eg.db, keys = gsub('hg19_',"",rowData(sce_obj)$Symbol), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
sce_obj<-sce_obj[!is.na(rowData(sce_obj)$anno.v ),]
rownames(sce_obj)<-rowData(sce_obj)$anno.v
colnames(sce_obj)<-colData(sce_obj)$Barcode
message('doing SCENT')
return(DoLandSCENT(exp.m = assay(sce_obj, i = "logcounts"), ppiA.m = net17Jan16.m,
                             mc.cores = 2, coordinates = NULL,
                             PLOT = FALSE, PDF = FALSE))

}


#' SCEpipe
#'
#' Annotate mouse and human, plus mito, trim, calculate QC for a list of SCE object
#'
#' @param SCElist a list of SingleCellExperiment objects from a species mixing experiment
#' @param countsProp cut off of minimum proportion of counts in human, default 0.8
#' @param mitoProp cut off of maximum  proportion of mito counts in human, default 0.1
#'
#' @return a list of trimmed and cleaned
#'
#' @examples
#' NULL
#'
#' @export	 
SCEpipe<- function(SCElist,countsProp=.8,mitoProp=.1){
return(lapply(SCElist,function(x){
mh<- ifelse(grepl("mm10_",rowData(x)$ID),"mouse","human")
mg<-ifelse(grepl("_MT-",rowData(x)$Symbol,ignore.case=T),"mito","gene")
rowData(x)$qcplot<-paste0(mh,"_",mg)
x<-calculateQCMetrics(x,feature_controls= split(1:nrow(rowData(x)),rowData(x)$qcplot)[-1])
x<-x[grep("hg19_",rowData(x)$Symbol), colData(x)$total_counts_endogenous/colData(x)$total_counts>countsProp & colData(x)$total_counts_human_mito/(colData(x)$total_counts_endogenous+colData(x)$total_counts_human_mito)<mitoProp]
X<-calculateQCMetrics(clearSpikes(x))
X <- X[,!isOutlier(X$total_counts, nmads=3,type="lower", log=TRUE)]
X <- X[nexprs(X, byrow=TRUE) >= 3,]
sizeFactors(X) <- librarySizeFactors(X)
return(normalize(X))
}))
}
