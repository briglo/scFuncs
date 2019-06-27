#SCENT
#i spent 3/4/19 figuring out that mclapply wouldnt run nicely on the cluster (at least the way imprelemted in SCENT, so asking for a lit of Vmem and only asking R to use a couple of cores seems to be stable.... will see if that works out)
#COMMAND
# cd /share/ScratchGeneral/briglo/scRNA/chrcha
# module load samcli/R/3.5.1
#  qsub -pe smp 16 -V -cwd -N tryScent -l h_vmem=10G -b y -j y Rscript ./source/SCENT.r
###individual steps if this doesnt work
# Integration.l <- DoIntegPPI(exp.m = tss, ppiA.m = net17Jan16.m)
# str(Integration.l)
# SR.o <- CompSRana(Integration.l, local = TRUE, mc.cores = 8)
# InferPotency.o <- InferPotency(SR.o)
# InferLandmark.o <- InferLandmark(InferPotency.o, pheno.v = phenoExample.v,
#                                  reduceMethod = "PCA", clusterMethod = "PAM", k_pam = 2)
#mapping info to PPI
#on cluster 25/3/19 test on one kills everything
# screen
# qrsh -pe smp 12 -l h_vmem=10G
# source ~/.profile
# module load samcli/R/3.5.1
# cd /share/ScratchGeneral/briglo/scRNA/chrcha
# R
args = commandArgs(trailingOnly=TRUE)

setwd("/share/ScratchGeneral/briglo/scRNA/chrcha")
Sys.setenv(MC_CORES=12)
library(parallel)

library("LandSCENT")
#data(net13Jun12.m)
data(net17Jan16.m)


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
# message("integrate")
# Integration.l <- DoIntegPPI(exp.m = assay(sce_obj, i = "logcounts"), ppiA.m = net17Jan16.m)
# message("Thebigbugger")
# SR.o <- CompSRana(Integration.l, local = TRUE, mc.cores = 2)
# message("dunno")
# InferPotency.o <- InferPotency(SR.o)
# return(list('integration'=Integration.l, SRana=SR.o, state=InferPotency.o ))

# InferPotency.o$potencyState
# InferLandmark.o <- InferLandmark(InferPotency.o, pheno.v = phenoExample.v,
#                                  reduceMethod = "PCA", clusterMethod = "PAM",
#                                  k_pam = 2)

scent<-DoLandSCENT(exp.m = assay(sce_obj, i = "logcounts"), ppiA.m = net17Jan16.m,
                             mc.cores = 2, coordinates = NULL,
                             PLOT = FALSE, PDF = FALSE)
save(scent,file=paste0("r_objects/",id,"_scent.rdata")) 
rm(scent)
gc()
}

#Test
# load('r_objects/190322_clean_norm_human_sce_SCORE_GPS_minus43978.rdata')
# #tmp<-trimNscent(GPSdat[[12]]) #if this works then all i have to do is wait :S it did!
# #save(tmp,file='r_objects/SCENT_43979.rdata')
# sink(paste0("scentlog_",args[1],'to',args[2],".txt"))
# for (i in args[1]:args[2]) trimNscent(sce_obj=GPSdat[[i]],id=names(GPSdat)[i])
# sink()

#TestAGAIN
load('r_objects/190322_clean_norm_human_sce.rdata')
GPSdat<-cleanhudat[grep("78|79",names(cleanhudat))]
rm(cleanhudat)
sink(paste0("scentlog_again_",args[1],'to',args[2],".txt"))
for (i in args[1]:args[2]) trimNscent(sce_obj=GPSdat[[i]],id=names(GPSdat)[i])
sink()


#might have to do on combined datasets




#####plott
# load("r_objects/id_43975.rdata")
# fnam<-grep("43975",dir("r_objects",pattern='scent'), value=T)
# comscent<-lapply(fnam, function(x){
# id<-gsub('_scent.rdata','',x)
# load(paste0("r_objects/",x))
# return(data.frame(row.names=paste0(gsub("_","\\.",id),"_",gsub("-[1-9]","",colnames(scent[[1]]))),scentPotency=scent$SR,scentPotencyState=scent$potencyState))
# }
# )

# allscent<-data.frame(do.call(rbind,comscent))
# tallscent<-allscent[rownames(allscent) %in% rownames(alldat@meta.data),]
# tallscent<-tallscent[rownames(alldat@meta.data),]
# alldat@meta.data<-data.frame(cbind(alldat@meta.data,tallscent))
# ggplot(alldat@meta.data,aes(x=orig.ident,fill=as.factor(scentPotencyState))) + geom_bar(position='fill') #low potency cells are missing from distant mets


# outs<-list('integration'=Integration.l, SRana=SR.o)
# save(outs,file='r_objects/190403_iHopeTHisworks.rdata')


# #THIS IS REALLY INEFFICIENT and explains why its a dick to RAM... (CompSRana<-function(Integration.l) {
#  Integration.l <- CompMaxSR(Integration.l)
#     maxSR <- Integration.l$maxSR
#     idx.l <- as.list(seq_len(ncol(Integration.l$expMC)))
#     out.l <- mclapply(idx.l, CompSRanaPRL, exp.m = Integration.l$expMC, 
#         adj.m = Integration.l$adjMC, local = TRUE, maxSR = maxSR, mc.cores=1)
#     SR.v <- sapply(out.l, function(v) return(v[[1]]))
#     invP.v <- sapply(out.l, function(v) return(v[[2]]))
#     S.v <- sapply(out.l, function(v) return(v[[3]]))
#     NS.v <- sapply(out.l, function(v) return(v[[4]]))
#     Integration.l$SR <- SR.v
#     Integration.l$inv <- invP.v
#     Integration.l$s <- S.v
#     Integration.l$ns <- NS.v



#     CompMaxSR <- function(Integration.l){
    
#     adj.m <- Integration.l$adjMC
    
#     # find right eigenvector of adjacency matrix
#     fa <- function(x,extra=NULL) {
#         as.vector(adj.m %*% x)
#     }
#     ap.o <- igraph::arpack(fa,options=list(n=nrow(adj.m),nev=1,which="LM"), sym=TRUE)
#     v <- ap.o$vectors
#     lambda <- ap.o$values
    
#     # maximum entropy
#     MaxSR <- log(lambda)
    
#     Integration.l$maxSR <- MaxSR
    
#     return(Integration.l)
# }

# CompSRanaPRL <- function(idx,
#                          exp.m,
#                          adj.m,
#                          local=TRUE,
#                          maxSR=NULL)
# {
    
#     # compute outgoing flux around each node
#     exp.v <- exp.m[,idx];
#     sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
#     invP.v <- exp.v*sumexp.v;
#     nf <- sum(invP.v);
#     invP.v <- invP.v/nf;
#     p.m <- t(t(adj.m)*exp.v)/sumexp.v;
#     S.v <- apply(p.m,1,CompS);
#     SR <- sum(invP.v*S.v);
#     # if provided then normalise relative to maxSR
#     if(is.null(maxSR)==FALSE){
#         SR <- SR/maxSR;
#     }
#     if(local){
#         NS.v <- apply(p.m,1,CompNS);
#     }
#     else {
#         NS.v <- NULL;
#     }
#     return(list(sr=SR,inv=invP.v,s=S.v,ns=NS.v));
# }

# CompNS <- function(p.v){
    
#     tmp.idx <- which(p.v>0);
#     if(length(tmp.idx)>1){
#         NLS <- -sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
#     }
#     else {
#         # one degree nodes have zero entropy, avoid singularity.
#         NLS <- 0;
#     }
#     return(NLS);
# }

# CompS <- function(p.v){
    
#     tmp.idx <- which(p.v>0);
#     LS <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
#     return(LS);
# }
d1
screen
qrsh -pe smp 12
cd /share/ScratchGeneral/briglo/scRNA/chrcha
module load briglo/R/3.6.0
R
fnam<-dir("r_objects/", pattern="scent")[-c(1,13,14)] (removing weird ones)
scentdat<-lapply(fnam, function (x){
load(paste0("r_objects/", x))
return(data.frame(scent$InferLandmark.l$cl,scent.SR=scent$SR,scent.pca=scent$coordinates,scent.ps=scent$potencyState))

})
names(scentdat)<-gsub("_",".",gsub("id_|_scent.rdata","",fnam))
save(scentdat,file="r_objects/190513_scentSummaryData.rdata")


#this one worked too
scentdat<-vector(mode='list',length=length(fnam))
for (i in 1: length(fnam)) {
    load(paste0("r_objects/", fnam[i]))
scentdat[[i]]<-data.frame(scent$InferLandmark.l$cl,scent.SR=scent$SR,scent.pca=scent$coordinates,scent.ps=scent$potencyState)
rm(scent)
}


load("r_objects/190513_scentSummaryData.rdata")
load("r_objects/mergeSeurat_AllSample1000Random_varNemtNbatchTSNEandMonocleAndSCORE_R3.6.rdata")

for (i in 1:length(scentdat)) rownames(scentdat[[i]])<- paste0(names(scentdat)[i],"_",gsub("-[1-9]","",rownames(scentdat[[i]])))

names(scentdat)<-NULL
ascent<-data.frame(do.call(rbind,scentdat))

#check whats what 
table(integrated@meta.data$orig.ident,rownames(integrated@meta.data) %in% rownames(ascent))


scent<-ascent[rownames(alldat@meta.data),]
alldat@meta.data<-data.frame(cbind(alldat@meta.data,scent))
ggplot(integrated@meta.data,aes(x=location,fill=as.factor(scent.ps))) + geom_bar(position='fill') + facet_wrap(~id)


di<-data.frame(do.call(rbind,strsplit(rownames(ascent),"_")))
didi<-data.frame(do.call(rbind,lapply(strsplit(as.character(di$X1),"\\."), function(x) x[1:2])))
ddat<-data.frame(cbind(ascent,didi))
ddat$X1<-gsub("M","",ddat$X1)
ddat$X2<-gsub("LN","Lymph",ddat$X2)
pdf("190515_tryplots.pdf")
ggplot(ddat,aes(x=X2,fill=as.factor(ps))) + geom_bar(position='fill') + facet_wrap(~X1)