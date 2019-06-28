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
library(matrixStats)
#library(scGPS)
library(locfit) #required but not dependency
library(e1071)
library(RColorBrewer)
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
    require(Seurat)
    data(cell.cyclegenes)
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
#' turns a Seurat object into a default magic/phate result
#'
#' @param seuratObj the Seurat object
#' @param geneList a vector of gene symbols
#' @param return.raw logical, return MAGIC object for init, default FALSE
#'
#' @return a matrix of MAGIC expression
#'
#' @examples
#' magicRand<-seurat2magic(seuratObj=integrated,geneList=sample(rownames(integrated@assays$SCT@scale.data),10),return.raw=F)
#'
#' @export
seurat2magic<-function(seuratObj,geneList=c("ZEB1","SNAI1",'SNAI2',"VIM"),return.raw=FALSE){

require(Seurat)
require(Rmagic)
require(Matrix)
message("building data object, this will be a while coz i dont optimise these things \n in fact im going to have a shower right now")
dat<-t(integrated@assays$RNA@counts)
metadat<-integrated@meta.data
keep_cols <- colSums(dat > 0) > 10
dat <- dat[,keep_cols]
message(sum(keep_cols),' out of ', length(keep_cols)," genes retained (",round(sum(keep_cols)/length(keep_cols),3),")")
  keep_rows <- Matrix::rowSums(dat) > 1000 & Matrix::rowSums(dat) < 25000
  dat <- dat[keep_rows,]
message(sum(keep_rows),' out of ', length(keep_rows)," cells retained (",round(sum(keep_rows)/length(keep_rows),3),")")
dat <- library.size.normalize(dat)
dat <- sqrt(dat)

message("running MAGIC")
data_MAGIC_all <- magic(dat, k=15, genes=geneList)
message("making MAGIC_PCA")
data_MAGIC_PCA <- magic(dat, genes="pca_only", 
                         t=4, init=data_MAGIC)
 message("running PHATE")                        
data_PHATE <- phate(dat, knn=4, decay=100, t=20)
tmp<- data.frame(cbind(data_MAGIC$result,data_MAGIC_PCA$result[,1:3],data_PHATE$embedding))
colnames(tmp)<-paste0("MAGIC_",colnames(tmp))
if(return.raw) return(data_MAGIC_all) else return(tmp)
}

#' subsetSeurat2cellPhone
#'
#' trims a seurat object to a smaller number of cells (ranked by number of genes expressed) per user defined cluster. writes out text files for cellphoneDB. Makes cellphoneDB a little quicker to run
#'
#' @param seuratObj the Seurat object
#' @param annoColumn name of metadata column you want
#' @param no.cells the number of cells per cluster you want
#' @param prefix file prefix for writing out data for cellphoneDB
#'
#' @return invisible(trimmedSeurat)
#'
#' @examples
#' ti<-subsetSeurat2cellPhone(seuratObj=integrated,annoColumn="SCT_snn_res.0.15",no.cells=50,prefix="small")
#' 
#' #on cluster
#' module load briglo/miniconda/3
#' source activate cellphone
#' qsub -V -cwd -b y -j y -pe smp 8 -N cpdb_1 "cellphonedb method statistical_analysis small_meta.txt small_counts.txt --project-name small --threshold 10 --threads 8"
#' EASY!!!
#'
#' @export
subsetSeurat2cellPhone<-function(seuratObj,annoColumn,no.cells,prefix) {
seuratObj@meta.data$newAnno<-paste0(prefix,"_",as.character(seuratObj@meta.data[,annoColumn]))
seuratObj<-SetIdent(seuratObj,value=seuratObj@meta.data$newAnno)
x<-seuratObj@meta.data[,c('nFeature_RNA','newAnno')]
sx<-split(x,x$newAnno)
osx<-lapply(sx, function(x) x[order(x$nFeature_RNA,decreasing=T),])
cells.use<-as.character(unlist(lapply(osx, function(x) head(rownames(x),no.cells))))
trimSeuratObject<-subset(seuratObj, cells=cells.use)
mappers<-getEntrez(trimSeuratObject,'all')
 
 da<-trimSeuratObject@assays$SCT@data[rownames(trimSeuratObject@assays$SCT@data) %in% mappers$hgnc_symbol[!is.na(mappers$ensembl_gene_id)],]
 m<-match(rownames(da),mappers$hgnc_symbol)
 rownames(da)<-mappers$ensembl_gene_id[m]
cd<-data.frame(Cell=rownames(trimSeuratObject@meta.data),cell_type=trimSeuratObject@meta.data$newAnno)
write.table(data.matrix(da),file=paste0(prefix,"_counts.txt"),quote=F,sep='\t')
 write.table(cd,file=paste0(prefix,"_meta.txt"),quote=F,sep='\t',row.names=F)
return(invisible(trimSeuratObject))
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
#' integrated<-makeMetaScore(integrated,sample(rownames(integrated@assays$SCT@scale.data),10))
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

#' plotSankey
#'
#' compares populations between metadata columns
#'
#' @param seuratObj the Seurat object
#' @param idvars character vector of medatadata column names
#'
#' @return a gawd awful javascript doohicky that i dont understand
#'
#' @examples
#' x<-plotSankey(integrated,c("orig.ident","SCT_snn_res.0.8"))
#'
#' @export
plotSankey<-function(seuratObj,idvar=c("varRes.0.3","emt_res.0.3")){
require(flipPlots)
message('try install_github("Displayr/flipPlots") if this doesnt work')
require(dplyr)
seuratObj@meta.data[,match(idvar,colnames(seuratObj@meta.data))] %>% arrange(.[,1]) %>% group_by_all() %>% summarise(COUNT = n()) ->> my.data
 #my.data<-as.factor(my.data[,1])

SankeyDiagram(my.data[, -grep("COUNT",colnames(my.data))],link.color = "Source",weights = my.data$COUNT,,max.categories = 1000)

}


#' getEntrez
#'
#' use biomart to map gene symbols from Seurat object to entrez ids
#'
#' @param seuratObj the Seurat object
#' @param geneList list of genes to calculate magic for default variable genes from SCT
#'
#' @return a dataframe of gene symbol and entrezID
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
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=rownames(integrated@assays$SCT@data),mart=human))
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
#' @param getEntrezObj
#'
#' @return a matrix of MAGIC expression
#'
#' @examples
#' NULL
#'
#' @export	 
	 splitDirection<-function(marker,getEntrezObj) lapply(split(marker,marker$avg_logFC>0), function(x) mkEnt(split(x$gene,x$cluster),getEntrezObj)) # splits a FindAllMarkers object into cluster and up/down


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


#' makeReactomePipe
#'
#' wrapper for seurat to reactome
#'
#' @param mkEntObj an output from  mkEnt or an element of splotDirection
#'
#' @return list with markers and a list of compareClusterObject for plotting
#'
#' @examples
#' enrichment<-makeReactomePipe(integrated,"SCT_snn_res.0.15")
#' for i in 1:length(enrichment$CP_result)) dotplot(enrichment$CP_result[[i]]) + ggtitle(names(enrichment$CP_result)[i])
#'
#' @export	 
makeReactomePipe<-function(seuratObj,metaDataColumn){
    require(Seurat)
    require(clusterProfiler)
    require(dplyr)
    message('making markers, might take a while')
seuratObj<-SetIdent(seuratObj,value=seuratObj@meta.data[,metaDataColumn])
m<-FindAllMarkers(seuratObj)
m %>% filter(p_val_adj<0.05) -> comb
message('making enrichment objects, might take a while again')
comb$cluster<-as.character(comb$cluster)
hasEntrez<-getEntrez(seuratObj,'all')
scomb<-splitDirection(comb,hasEntrez)
oi<-lapply(scomb, reactomeClusts)
return(list(
"markers"=m,
"CP_result"=oi,
"masterID"=hasEntrez))
}

#' preprocCellPhone
#'
#' takes CellPhoneDB output and turns it into something useable
#'
#' @param dir the directory containing cellPhoneDB output, defaults to .
#' @param pval p value cutoff for reporting an interaction, default 0.01
#' @param varval variance value cutoff for retaining interactions that are different between conditions, default 0
#'
#' @return a list containing raw data trimmed for 0 and one trimmed for p and var (plotdat)
#'
#' @examples
#' cd("PATH/TO/CELLPHONEDB/OUT/small")
#' x<-preprocCellphone(varval=0,pval=.05)
#' 
#' @export	 
preprocCellphone<-function(prefix,pval=0.01,varval=0){
    require(matrixStats)
message("you should run this in an cellphone results directory")
message("pval cutoff=",pval)
message("variation value=",varval)

fnam<-c(pval="pvalues.txt",mean='means.txt')
dat<-lapply(fnam,function(x) read.table(x, header=T, stringsAsFactors = F, sep='\t', comment.char = ''))
dat<-lapply(dat,function(x){
rownames(x) = x$interacting_pair
x<-x[,-c(1:9)]
})
tdat<-lapply(dat,function(x) x[rowSums(dat$pval<pval)>0,colSums(dat$pval<pval)>0])
tmp<-tdat$mean
tmp[tdat$pval>pval]=0
trimdat<-tmp[matrixStats::rowVars(data.matrix(tmp))>varval,]
target<-unlist(lapply(strsplit(colnames(trimdat),"_"), function(x) paste0(x[3],"_",x[4])))
source<-unlist(lapply(strsplit(colnames(trimdat),"_"), function(x) paste0(x[1],"_",x[2])))
df<-data.frame(target,source,stringsAsFactors=F)
con<-lapply(seq(0,1,.1), function(y) colSums(x>y))
names(con)<-paste0('countAboveMean',seq(0,1,.1))
tdat$countdat<-data.frame(source,target,do.call(cbind,con))
#print(pheatmap::pheatmap(data.matrix(tdat$plodat))) #this needs a big display object
return(invisible(tdat))
}


#' intgraph
#'
#' takes a preprocCellphone$countdat object and plots a graph
#'
#' @param cellphoneDB_data the countdat output from plotdat object
#' @param scoreCut mean score cutoff to use (seq(0,1,0.1)), default 0.3
#' @param numberCut the number of interactions required to be plotted, default 0
#' @param numberSplit the number of interactions to distinguish high from low, default 35
#'
#' @return invisible list of plot objects
#'
#' @examples
#' intgraph(x$countdat)
#'
#' @export	 

intgraph<- function(cellphoneDB_data,scoreCut=0.3, numberCut=0, numberSplit=35){
require(ggraph)
   require(dplyr)
   require(igraph)
   require(cowplot)
   message("use like this: intgraph(scoreCut=0.3, numberCut=0, numberSplit=35)")
 tmp<-cellphoneDB_data[cellphoneDB_data[,paste0("countAboveMean",scoreCut)]>numberCut,c('source','target',paste0("countAboveMean",scoreCut))]
 colnames(tmp)[3]<-"score"
 an<-data.frame(ID=unique(c(tmp$source,tmp$target)))
 gr<-graph_from_data_frame(tmp,directed = F,vertices=an)
 deg=degree(gr, mode ="all")
 
   p1<-ggraph(gr, layout = 'linear', circular=T) +     geom_edge_arc(aes(width =score,alpha=score>numberSplit,colour='black')) #+ geom_node_point(aes(color=ID,size=deg)) + geom_node_text(aes(label=ID), repel=F) + ggtitle(paste0("scoreCut=",scoreCut, ", numberCut=",numberCut,", numberSplit=",numberSplit))
   pdf("intgraph_1.pdf")
   p2<-ggplot(tmp,aes(x=score)) + geom_bar()
   plot(gr,arrow.size=.2,edge.curved=-.1,edge.width=log(E(gr)$score)) 
      dev.off()
      pdf("intgraph_2.pdf") 
   print(ggdraw() + draw_plot(p1 + theme(legend.justification="bottom"), 0,0,1,1) + draw_plot(p2 + geom_vline(xintercept=numberSplit) ,0.75, 0.7, 0.2, 0.3) )
   dev.off()
return(invisible(list(p1=p1,p2=p2,graph=gr)))
   }