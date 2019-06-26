library(Seurat)
library(Rmagic)
library(ggplot2)
library(Matrix)
library(viridis)
library(reticulate)
library(Matrix)
library(phateR)
library(Seurat)
library(reticulate)
library(scGPS)
library(locfit) #required but not dependency
library(e1071)
use_python("/share/ClusterShare/software/contrib/briglo/miniconda3/envs/magic/bin/python")



#' buildMasterSeurat
#'
#' turns a list of Seurat objects into a master object
#'
#' @param seuratList the list of Seurat objects, usually from makeSeuratList
#' @param ndims dimensions to use in integration,SCT and UMAP, default 1:30
#' @param mt.cutoff numeric percentage cutoff for mito expression, defaullt 10
#' @param regress.cell.cycle logical, whether to regress out cell cycle along with mito and nFeature
#' @param res.vec vector of resolutions to find clusters at, defaults to c(seq(0,1,0.1),1.5,2)
#'
#' @return an integrated master seurat object
#'
#' @examples
#' NULL
#'
#' @export
buildMasterSeurat<-function(seuratList,ndims=1:30,mt.cutoff=10,regress.cell.cycle=TRUE,res.vec=c(seq(0,1,0.1),1.5,2)){
anchors <- FindIntegrationAnchors(object.list = seuratList, dims = ndims)
integrated <- IntegrateData(anchorset = anchors, dims = ndims) #this takes ages
integrated@active.assay="RNA"
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
integrated<-subset(integrated, subset = percent.mt < mt.cutoff)
integrated<-CellCycleScoring(object = integrated, s.features = toupper(s.genes), g2m.features = toupper(g2m.genes), set.ident = FALSE)
integrated@active.assay="integrated"
if(regress.cell.cycle) {
    integrated<-SCTransform(integrated, vars.to.regress = c("percent.mt","nFeature_RNA","S.Score","G2M.Score"), verbose = FALSE) } else {
        integrated<-SCTransform(integrated, vars.to.regress = c("percent.mt","nFeature_RNA"), verbose = FALSE)
    }
integrated<-RunPCA(integrated, verbose = FALSE)
integrated<-RunUMAP(integrated, dims = ndims, verbose = T)
integrated<-FindNeighbors(integrated, dims = ndims, verbose = FALSE)
for (res in resvec) integrated<-FindClusters(integrated,resolution=res, verbose = FALSE)
return(integrated)
}



#' seurat2scGPS
#'
#' turns a Seurat object into a scGPS object
#'
#' @param seuratObj the Seurat object
#'
#' @return a scGPS object
#'
#' @examples
#' NULL
#'
#' @export
seurat2scGPS<-function(seuratObj) {
return(new_scGPS_object(ExpressionMatrix=as(seuratObj@assays$RNA@counts,'matrix'),  CellMetadata = seuratObj@meta.data,GeneMetadata=data.frame(rownames(seuratObj@assays$RNA@counts),row.names=rownames(seuratObj@assays$RNA@counts))))
}

#' seurat2monocle
#'
#' turns a Seurat object into a monocle object
#'
#' @param seuratObj the Seurat object
#'
#' @return a monocle object
#'
#' @examples
#' NULL
#'
#' @export
seurat2monocle<-function(seuratObj) {
pd <- AnnotatedDataFrame(seuratObj@meta.data)
gene.names <- as.array(unlist(rownames(seuratObj@assays$RNA@counts)))
fd <- data.frame(gene.names)
rownames(fd) <- gene.names
colnames(fd) <- c("gene_short_name")
fd <- AnnotatedDataFrame(fd)
HSMM <- newCellDataSet(as(as.matrix(seuratObj@assays$RNA@counts[, colnames(seuratObj@meta.data)]), "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit=1, expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM) 
return(HSMM)
}

#' seurat2magic
#'
#' turns a Seurat object into a default magic/phate resuly
#'
#' @param seuratObj the Seurat object
#' @param geneList list of genes to calculate magic for default to a handfull of EMT genes
#'
#' @return a matrix of MAGIC expression
#'
#' @examples
#' NULL
#'
#' @export
seurat2magic<-function(seuratObj,geneList=c("ZEB1","SNAI1",'SNAI2',"VIM")){

dat<-t(integrated@assays$RNA@counts)
metadat<-integrated@meta.data
keep_cols <- colSums(dat > 0) > 10
dat <- dat[,keep_cols]
  keep_rows <- Matrix::rowSums(dat) > 1000 & Matrix::rowSums(dat) < 25000
  dat <- dat[keep_rows,]
dat <- library.size.normalize(dat)
dat <- sqrt(dat)

data_MAGIC_all <- magic(dat, k=15, genes=geneList)
data_MAGIC_PCA <- magic(dat, genes="pca_only", 
                         t=4, init=data_MAGIC)
data_PHATE <- phate(dat, knn=4, decay=100, t=20)
tmp<- data.frame(cbind(data_MAGIC$result,data_MAGIC_PCA$result[,1:3],data_PHATE$embedding))
colnames(tmp)<-paste0("MAGIC_",colnames(tmp))
return(tmp)
}


#' makeMetaScore
#'
#' adds a gene metascore to meta.data
#'
#' @param seuratObj the Seurat object
#' @param geneList list of genes to calculate magic for default to a handfull of EMT genes
#'
#' @return a dataframe of metadata +/- coords and metascore
#'
#' @examples
#' integrated<-plotMetaScore(integrated,sample(rownames(integrated@assays$SCT@scale.data),10))
#' ggplot(integrated@meta.data,aes(x=UMAP_1,y=UMAP_2,color=metascore)) + geom_point(alpha=.5) + scale_colour_gradientn(colours = jet.colors(7))
#'
#' @export
makeMetaScore<-function(seuratObj,geneList,reduction=NA){
    require(Seurat) # for playing with object
    require(ggplot2) # just to make sure can plot
    require(matlab) # for the jet color palette
    require(dplyr) # for the gene contribution transformation
    message(table(geneList %in% rownames(seuratObj@assays$integrated@scale.data))['TRUE']," out of ",length(geneList)," genes entered were used to generate score\n") # just shows how many are actually contributing to the score
    hm<-matrixStats::rowVars(seuratObj@assays$SCT@scale.data[rownames(seuratObj@assays$SCT@scale.data) %in% geneList,])
    names(hm)<-rownames(seuratObj@assays$SCT@scale.data[rownames(seuratObj@assays$SCT@scale.data) %in% geneList,])
    message("top 10 contributing genes (by percentage) contributing to signature")
    print((hm*100/sum(hm)) %>% sort(decreasing=T) %>% signif(2) %>% head(10) )
    
    seuratObj$metascore=colSums(seuratObj@assays$SCT@scale.data[rownames(seuratObj@assays$SCT@scale.data) %in% geneList,])
    return(seuratObj)
}



#' getEntrez
#'
#' use biomart to map gene symbols from Seurat object to entrez ids
#'
#' @param seuratObj the Seurat object
#' @param geneList list of genes to calculate magic for default variable genes from SCT
#'
#' @return a matrix of MAGIC expression
#'
#' @examples
#' hasEntrez<-getEntrez(integrated)
#'
#' @export	 

	 getEntrez<-function(seuratObj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
         human<- useMart(biomart='ensembl', dataset = "hsapiens_gene_ensembl")
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=integrated@assays$SCT@var.features,mart=human))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=rownames(integrated@assays$SCT$data),mart=human))
	 }
	 }#


#' mkEnt
#'
#' retrieves entrez ids from a getEntrez for a list of vectors of gene symbols
#'
#' @param listOids a list of gene symbol vectors
#' @param getEntrezObj lan output from getEnrez functin
#'
#' @return a list of entrez id vectors for enrichment
#'
#' @examples
#' entList<-mkEnt(split(markerList$gene,markerList$cluster),hasEntrez)
#'
#' @export	 

	 mkEnt<-function(listOids,getEntrezObj){
	     require(clusterProfiler)
	     require(ReactomePA)
	     return(lapply(listOids,function(x) as.character(na.omit(getEntrezObj$entrezgene[match(x,getEntrezObj$hgnc_symbol)]))))
	  } 


#' splitDirection
#'
#' turns a  Seurat FindMarkers object into a bunch of entrez IDs
#'
#' @param marker a Find[All]Markers object (may need marker$gene to be added)
#'
#' @return a matrix of MAGIC expression
#'
#' @examples
#' NULL
#'
#' @export	 
	 splitDirection<-function(marker) lapply(split(marker,marker$avg_logFC>0), function(x) mkEnt(split(x$gene,x$cluster))) # splits a FindAllMarkers object into cluster and up/down


#' reactomeClusts
#'
#' wrapper for compareCluster reactome for a list of entrez ID vectors
#'
#' @param mkEntObj an output from  mkEnt or an element of splotDirection
#'
#' @return compareClusterObject for plotting
#'
#' @examples
#' oi<-reactomeClusts(mkEntObj)
#' maybe dotplot(oi)
#' oi<-lapply(splitDirectionObj,reactomeClusts)
#' dotplot(oi[[1]])
#'
#' @export	 

	 reactomeClusts<-function(mkEntObj) {
	     require(clusterProfiler)
	     require(ReactomePA)
	     message("making object, takes a while")
	     return(compareCluster(geneCluster = mkEntObj, fun = "enrichPathway",organism='human',universe=as.character(na.omit(hasEntrez$entrezgene)),pvalueCutoff=0.05,readable=T)
	 )
	 } #finds enriched reactome clusters for a list of entrez ids requires a background getEntrez object called hasEntrez



########up to here

#' seurat2magic
#'
#' turns a Seurat object into a default magic/phate resuly
#'
#' @param seuratObj the Seurat object
#' @param geneList list of genes to calculate magic for default to a handfull of EMT genes
#'
#' @return a matrix of MAGIC expression
#'
#' @examples
#' NULL
#'
#' @export	 

preprocCellphone<-function(pval=0.01,varval=0.05){
message("you should have run this in an cellphone results directory")
message("pval cutoff=",pval)
message("variation value=",varval)

fnam<-c(pval='pvalues.txt',mean="means.txt")
dat<-lapply(fnam,function(x) read.table(x, header=T, stringsAsFactors = F, sep='\t', comment.char = ''))
dat<-lapply(dat,function(x){
rownames(x) = x$interacting_pair
x<-x[,-c(1:9)]
})
tdat<-lapply(dat,function(x) x[rowSums(dat$pval<pval)>0,colSums(dat$pval<pval)>0])
tmp<-tdat$mean
tmp[tdat$pval>pval]=0
tdat$plotdat<-tmp[matrixStats::rowVars(data.matrix(tmp))>varval,]
#print(pheatmap::pheatmap(data.matrix(tdat$plodat))) #this needs a big display object
return(invisible(tdat))
}

#' seurat2magic
#'
#' turns a Seurat object into a default magic/phate resuly
#'
#' @param seuratObj the Seurat object
#' @param geneList list of genes to calculate magic for default to a handfull of EMT genes
#'
#' @return a matrix of MAGIC expression
#'
#' @examples
#' NULL
#'
#' @export	 

intgraph<- function(scoremat=df,scoreCut=0.3, numberCut=0, numberSplit=35){
   require(ggraph)
   require(dplyr)
   require(igraph)
   require(cowplot)
   message("use like this: intgraph(scoreCut=0.3, numberCut=0, numberSplit=35)")
 tmp<-scoremat[scoremat[,paste0("meancut_",scoreCut)]>numberCut,c('source','target',paste0("meancut_",scoreCut))]
 colnames(tmp)[3]<-"score"
 an<-data.frame(ID=unique(c(scoremat$source,scoremat$target)))
 gr<-graph_from_data_frame(tmp,directed = T,vertices=an)
 deg=degree(gr, mode ="all")
 
 colname<-paste0("countAboveMean",numberCut)
  mypalette <- palette(brewer.pal(n=10,name='Paired'))[c(6,8,4,2,10,7,3,1,9)]
  V(gr)$color=mypalette
  plot(gr,arrow.size=.2,edve.curved=-.1,edge.width=E(gr)$score)#,vertex.color=V(gr)$color)
  # ggraph(gr, layout = 'linear')  +     geom_edge_arc(aes(width =score,alpha=score>numberSplit)) + geom_node_point()
   
   #aes(color=ID,shape=type,size=count)) + geom_node_label(aes(label=ID),label.size=.2,nudge_y=-1) + coord_flip() + ggtitle(paste0("scoreCut=",scoreCut, ", numberCut=",numberCut,", numberSplit=",numberSplit))
   #p2<-ggplot(tmp,aes(x=score)) + geom_bar()
   #ggdraw() + draw_plot(p1 + theme(legend.justification="bottom"), 0,0,1,1) + draw_plot(p2 + geom_vline(xintercept=numberSplit) ,0.75, 0.7, 0.2, 0.3)
   }
