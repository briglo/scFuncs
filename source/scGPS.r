qrsh -pe smp 16
cd /share/ScratchGeneral/briglo/scRNA/chrcha
module load briglo/R/3.5.1
##reading in raw data
R
source("source/read10X.r")
library(Matrix)
library(SingleCellExperiment)
fnam<-dir("data/")

dat<-lapply(fnam, function(X){
ddir<-paste0(pwd(),"/data/",X,"/outs/filtered_feature_bc_matrix/")
x<-read10xCounts(ddir)
isSpike(x, "Mouse") <-grepl("mm10_",rownames(x))
sizeFactors(x) <- Matrix::colSums(assay(x))
return(x)
})

save(dat,file="r_objects/allSingleCell.rdata")

#coz i cant compile any mother fucking thing, the rest is local
library(scater)
library(Rtsne)

load("r_objects/allSingleCell.rdata")
#Annotate mouse and human, plus mito, calculate QC
qdat2<-lapply(dat,function(x){
mh<- ifelse(grepl("mm10_",rowData(x)$ID),"mouse","human")
    mg<-ifelse(grepl("_MT-",rowData(x)$Symbol,ignore.case=T),"mito","gene")
rowData(x)$qcplot<-paste0(mh,"_",mg)
return(calculateQCMetrics(x,feature_controls= split(1:nrow(rowData(x)),rowData(x)$qcplot)[-1]))
}

#plot mean counts by gene anno
pps<-lapply(qdat2,function(x) {
   return(ggplot(as(rowData(x),'data.frame'),aes(x=log(mean_counts),fill=qcplot)) + geom_density(alpha=.5))
})
plot_grid(plotlist=pps,labels=names(pps))
ggsave("../plots/190321_countDist_humanMouse.pdf")


#trim by hard cutoff (more than 80% genes of interest (human) and less than 10% mito, ONY HUMAN GENES)
tqhudat<-lapply(qdat2, function(x){
 return(x[grep("hg19_",rowData(x)$Symbol), colData(x)$total_counts_endogenous/colData(x)$total_counts>.8 & colData(x)$total_counts_human_mito/(colData(x)$total_counts_endogenous+colData(x)$total_counts_human_mito)<.1])  
})

#plot post trim spreads
pps<-lapply(tqhudat,function(x) {
   return(ggplot(as(rowData(x),'data.frame'),aes(x=log(mean_counts),fill=qcplot)) + geom_density(alpha=.5))
})
plot_grid(plotlist=pps,labels=names(pps))
ggsave("../plots/190321_countDist_postSampleHumanGeneTrim.pdf")

#Clear and recalculate QC
 rqtqhudat<-lapply(tqhudat,function(X) {
 X<-clearSpikes(X) #coz removed mouse genes
 return(calculateQCMetrics(X))
 })

#clean object not normalized (coz it crashed lol)
cleanhudat<-lapply(rqtqhudat,function(X,minCellExp=3) {
    X <- X[,!isOutlier(X$total_counts, nmads=3,type="lower", log=TRUE)]
    X <- X[nexprs(X, byrow=TRUE) >= 3,]
    sizeFactors(X) <- librarySizeFactors(X)
    return(X)
    # return(normalize(X))
})

 save(cleanhudat,file="~/Desktop/briglo_backup/scRNA/chrcha/r_objects/190322_clean_human_sce.rdata")

#clean object normalized
lapply(cleanhudat,normalize)
save(cleanhudat,file="~/Desktop/briglo_backup/scRNA/chrcha/r_objects/190322_clean_norm_human_sce.rdata")

#this uses too much cpu to do locally... but im doing it anyweay becasye rccp didnt fucking compile anywhere on the cluster
SCOREdat<-lapply(cleanhudat,function(X) {
    require(DropletUtils)
require(scater)
require(Matrix)
require(scGPS)
assay(X) <- as(assay(X),"matrix")
return(CORE_scGPS_bagging(X, remove_outlier = c(0), PCA=FALSE , bagging_run = 20, subsample_proportion = .8))
})
save(SCOREdat,file="~/Desktop/briglo_backup/scRNA/chrcha/r_objects/190322_clean_norm_human_sce_SCORE.rdata")

pdf(file="plots/190327_COREplot.pdf")
p<-lapply(SCOREdat, function(x) print(plot_CORE(x$tree, x$Cluster, color_branch = c("#208eb7", "#6ce9d3", "#1c5e39", "#8fca40", "#154975", "#b1c8eb"))))
dev.off()

pdf(file="plots/190327_COREplotOPT.pdf")
lapply(SCOREdat, function(x) print(plot_optimal_CORE(original_tree= x$tree, optimal_cluster = unlist(x$Cluster[x$optimal_index]), shift = -100)))
dev.off()

#trimming to samples that have the tissues we want to compare
sid<-split(names(SCOREdat),gsub("M","",unlist(lapply(strsplit(names(SCOREdat),"_"), function(x) x[1]))))
tsid<-sid[c(2:5)]
tsid[['43978']]<-tsid$samp_43978[c(2,4,6,7)] #consider chosing other replicates?

#this might need to be redone, clusters looks fucked
GPSdat<-vector(mode='list',length=length(unlist(tsid[c(1,2,4)])))
for (x in unlist(tsid[c(1,2,4)])) {
    print(x) 
   tmp<-cleanhudat[[x]]
colData(tmp)$Cluster=unlist(SCOREdat[[x]]$Cluster[SCOREdat[[x]]$optimal_index])
GPSdat[[as.character(x)]]<-tmp
}


GPSdat<-GPSdat[-c(1:12)]
save(GPSdat,file="~/Desktop/briglo_backup/scRNA/chrcha/r_objects/190322_clean_norm_human_sce_SCORE_GPS_minus43978.rdata")

vgene<-lapply(names(GPSdat), function(x){
     print(x)
     fit <- trendVar(GPSdat[[x]], parametric=TRUE,use.spikes=F)
     decomp <- decomposeVar(GPSdat[[x]], fit)
     return( rownames(decomp)[decomp$FDR<0.01 & decomp$mean>=1])
})
names(vgene)<-names(GPSdat)

DEgenes<-lapply(names(GPSdat), function(x){
    print(x)
      genes=vgene[[x]]
 return(findMarkers_scGPS(expression_matrix=as(assay(GPSdat[[x]])[vgene[[x]],],'matrix'), 
    cluster = as.numeric(colData(GPSdat[[x]])[,"Cluster"]), 
    selected_cluster=unique(as.numeric(colData(GPSdat[[x]])[,"Cluster"])) ))
})
names(DEgenes)<-names(GPSdat)
save(DEgenes,vgene, file="r_objects/190327_DEfromVargenes.rdata")



#notes should work on a sparse matrix
#up next actual scGPS
#need to do a matrix

sid<-split(names(SCOREdat),gsub("M","",unlist(lapply(strsplit(names(SCOREdat),"_"), function(x) x[1]))))
 tsid<-sid[c(2:5)]
tsid[['43978']]<-tsid[['43978']][c(2,4,6,7)]


for (id in c(43975,43979)) {
x<-GPSdat[grep(id, names(GPSdat))]
save(x, file=paste0('r_objects/id_',id[1],'_sce.rdata'))
}

##########################scGPS
#1) to make my computer stop spazzing out i split 190322_clean_norm_human_sce_SCORE_GPS_minus43978.rdata by one sample for now
setwd("/Users/briglo/Desktop/briglo_backup/scRNA/chrcha")
load('r_objects/tmp.rdata') #split by sample name 73 done, 75 in progress, 

load("r_objects/190327_DEfromVargenes.rdata")
library(scGPS)

#make scGPS data with a stupid matrix
 prim<-NewscGPS(ExpressionMatrix=as(assay(x[[4]]),'matrix'),  CellMetadata = data.frame(X=colData(x[[4]])$Cluster,row.names=colData(x[[4]])$Barcode),GeneMetadata=data.frame(rowData(x[[4]])[,1:2],row.names=x[[4]]$ID))

#genes<-DEgenes[[grep("PT",names(x),value=T)]]$id# use all var genes from all samples FIRST RESULT
genes<-DEgenes[[grep("M43975",names(DEgenes))]] #USED MARKER GENES FOR PRIMARY april 1st run (need to discuss), its faster, only little diff, need to look at overlap.

cluster_prim <- colData(prim)[,1] #describes what cell is what cluster
 allID <- unique(prim$X)# not sure but theoretically want to iterate over all source and all target, TBD
 mets<-names(x)[1:3] 

#this is the big bastard, turns out the scGPS object is very specific
sink(file="scgpsRunPTClustgenes_43975.txt")

allSCGPS_43975_PTCgen<-lapply(mets, function(y){
met<-NewscGPS(ExpressionMatrix=as(assay(x[[y]]),'matrix'),  CellMetadata = data.frame(X=colData(x[[y]])$Cluster,row.names=colData(x[[y]])$Barcode),GeneMetadata=data.frame(rowData(x[[y]])[,1:2],row.names=x[[y]]$ID))
cluster_met <- colData(met)[,1]
lapply(allID, function(c_selectID){
return(bootstrap_scGPS(nboots = 3, mixedpop1 = prim, 
    mixedpop2 = met, genes = unique(genes[[grep(c_selectID,names(genes))]]$id), c_selectID  = c_selectID,
    listData = list(), cluster_mixedpop1 = cluster_prim, 
    cluster_mixedpop2 = cluster_met, trainset_ratio = 0.7))
})
})
sink()
 save(allSCGPS_43975_PTCgen, file="r_objects/190329_scGPS_43975_PTmarkers.rdata")

 lpt<-mkasso(allSCGPS_43975_PTCgen)
 pp3<-plotLasso(lpt)
save(pp2,file="r_objects/190329_scGPS_43975_PTmarkerGenes_networkDiag.rdata")


####

combined <- rbind(LASSO_C1S1,LASSO_C2S1,LASSO_C3S1)
nboots = 3

combined_D3obj <-list(Nodes=combined[,(nboots+3):(nboots+4)],
    Links=combined[,c((nboots+2):(nboots+1),ncol(combined))]) 
combined <- combined[is.na(combined$Value) != TRUE,]


library(networkD3)

Node_source <- as.vector(sort(unique(combined_D3obj$Links$Source)))
Node_target <- as.vector(sort(unique(combined_D3obj$Links$Target)))
Node_all <-unique(c(Node_source, Node_target))

#assign IDs for Source (start from 0)
Source <-combined_D3obj$Links$Source
Target <- combined_D3obj$Links$Target

for(i in 1:length(Node_all)){
  Source[Source==Node_all[i]] <-i-1
  Target[Target==Node_all[i]] <-i-1
}

combined_D3obj$Links$Source <- as.numeric(Source)
combined_D3obj$Links$Target <- as.numeric(Target)
combined_D3obj$Links$LinkColor <- combined$NodeGroup

#prepare node info 
node_df <-data.frame(Node=Node_all)
node_df$id <-as.numeric(c(0, 1:(length(Node_all)-1)))

suppressMessages(library(dplyr))
n <- length(unique(node_df$Node))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
Color = getPalette(n)
node_df$color <- Color
suppressMessages(library(networkD3))
p3<-sankeyNetwork(Links =combined_D3obj$Links, Nodes = node_df,
    Value = "Value", NodeGroup ="color", LinkGroup = "LinkColor",
    NodeID="Node", Source="Source", Target="Target", fontSize = 22)

save(p1,p2,p3, file="r_objects/190329_scGPS_43973_networkDiag.rdata")

################ plotting enrichment
plotDElist<- function(DElist,...){
ggs<-unique(as.character(ano[rownames(ano) %in% gsub("hg19_","",DElist$id[1:100]),'Symbol']))
uni<-unique(as.character(ano[rownames(ano) %in% gsub("hg19_","",genes),'Symbol']))
dotplot(annotate_scGPS(ggs,pvalueCutoff=0.05, gene_symbol=TRUE,universe=uni), showCategory=10, font.size = 6,...)
}

pdf(file="plots/190329_enrichmentplots.pdf")
et<-lapply(1:length(DEgenes), function(i){
    lapply(1:length(i), function(j) plotDElist(DEgenes[[i]][[j]],title=paste(names(DEgenes)[i],names(DEgenes[[i]][j]))))
})
dev.off()

###had to modify it coz the native scGPS function doesnt allow for a background gene list
annotate_scGPS<-function (DEgeneList, pvalueCutoff = 0.05, gene_symbol = TRUE, 
    species = "human",universe=genes) 
{
    if (species == "human") {
        if (gene_symbol == TRUE) {
            convert_to_gene_ID = clusterProfiler::bitr(DEgeneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            universe_to_ID =clusterProfiler::bitr(universe,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            print("Original gene number in geneList")
            print(length(DEgeneList))
            print("Number of genes successfully converted")
            print(nrow(convert_to_gene_ID))
        }
        else {
            stop("The list must contain human gene symbols")
        }
    }
    else if (species == "mouse") {
        if (gene_symbol == TRUE) {
            convert_to_gene_ID = clusterProfiler::bitr(DEgeneList, 
                fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
            print("Original gene number in geneList")
            print(length(DEgeneList))
            print("Number of genes successfully converted")
            print(nrow(convert_to_gene_ID))
        }
        else {
            stop("The list must contain mouse gene symbols")
        }
    }
    Reactome_pathway_test <- ReactomePA::enrichPathway(gene = convert_to_gene_ID$ENTREZID, pvalueCutoff = 0.05,readable = TRUE, universe=universe_to_ID$ENTREZID)
    output_df <- as.data.frame(Reactome_pathway_test)
    return(Reactome_pathway_test)
}
