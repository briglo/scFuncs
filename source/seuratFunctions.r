############functions#################
#QOL
sys<-function() system('open .')
#Modified input function coz its a shit
Read10X<- function (data.dir = NULL) {
    full.data <- list()
    for (i in seq_along(data.dir)) {
        run <- data.dir[i]
        if (!dir.exists(run)) {
            stop("Directory provided does not exist")
        }
        if (!grepl("\\/$", run)) {
            run <- paste(run, "/", sep = "")
        }
        barcode.loc <- paste0(run, "barcodes.tsv")
        gene.loc <- paste0(run, "features.tsv")
        matrix.loc <- paste0(run, "matrix.mtx")
        if (!file.exists(barcode.loc)) {
            stop("Barcode file missing")
        }
        if (!file.exists(gene.loc)) {
            stop("Gene name file missing")
        }
        if (!file.exists(matrix.loc)) {
            stop("Expression matrix file missing")
        }
        data <- readMM(file = matrix.loc)
        cell.names <- readLines(barcode.loc)
        gene.names <- readLines(gene.loc)
        if (all(grepl(pattern = "\\-1$", x = cell.names))) {
            cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                FUN = ExtractField, field = 1, delim = "-")))
        }
        rownames(x = data) <- make.unique(names = as.character(x = sapply(X = gene.names, 
            FUN = ExtractField, field = 2, delim = "\\t")))
        if (is.null(x = names(x = data.dir))) {
            if (i < 2) {
                colnames(x = data) <- cell.names
            }
            else {
                colnames(x = data) <- paste0(i, "_", cell.names)
            }
        }
        else {
            colnames(x = data) <- paste0(names(x = data.dir)[i], 
                "_", cell.names)
        }
        full.data <- append(x = full.data, values = data)
    }
    full.data <- do.call(cbind, full.data)
    return(full.data)
}

ExtractField <- function (string, field = 1, delim = "_") 
{
    fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), 
        split = ",")))
    if (length(x = fields) == 1) {
        return(strsplit(x = string, split = delim)[[1]][field])
    }
    return(paste(strsplit(x = string, split = delim)[[1]][fields], 
        collapse = delim))
}

testInput<-function(rawdata,minGenes=200,minCells=3,totExpr=1e4,GacceptHi=20000,NacceptHi=50000,mitcut=0.10){
    require(Seurat)
    require(Matrix)
    message(paste("mingenes=",minGenes))
    message(paste("mincells=",minCells))
    message(paste("total expr=",totExpr))
    
    print("setting up")
    
    elf5 <- CreateSeuratObject(raw.data=rawdata, min.genes = 1, min.cells=3, names.delim="-", project = "chrcha") 
    print("annotating mito genes")
    mito.genes <- grep("^MT-", rownames(elf5@data), value = T)
    percent.mito <- Matrix::colSums(elf5@raw.data[mito.genes, ])/Matrix::colSums(elf5@raw.data)
    elf5 <- AddMetaData(elf5, percent.mito, "percent.mito")    
	  print(table(elf5@ident))
elf5<-CellCycleScoring(object = elf5, s.genes = toupper(s.genes), g2m.genes = toupper(g2m.genes), set.ident = FALSE)

	# message("subsetting data")
   # message(paste("gene accept hi=",GacceptHi))
    #elf5 <- SubsetData(elf5, subset.name = "nGene", accept.high = GacceptHi)
	# print(table(elf5@ident))
	   #  message(paste("nUMI accept hi=",NacceptHi))
  #  elf5 <- SubsetData(elf5, subset.name = "nUMI", accept.high = NacceptHi)
	#  print(table(elf5@ident))
	#  message(paste("mitoCutoff=",mitcut))
	 elf5<-SubsetData(elf5, subset.name = "percent.mito", accept.high = mitcut)
	  print(table(elf5@ident))
	 elf5<-NormalizeData(elf5,normalization.method = "LogNormalize",scale.factor = 10000)
	 elf5<-ScaleData(object = elf5)
	 elf5 <- FindVariableGenes(elf5,  mean.function = ExpMean, dispersion.function = LogVMR,  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    print(paste(length(elf5@var.genes),"variable genes"))
    return(elf5)
}

mkPCA<-function(suraObj,pcComp=100,genelist=suraObj@var.genes){
    message("PCs to compare=",pcComp)
    require(Seurat)
    require(Matrix)
    message("making objects and PCAing, go grab a coffee")
    clus<-suraObj
    clus <- RunPCA(clus, pc.genes=clus@var.genes, pcs.compute = pcComp,do.print=F)
    return(clus)
	PCElbowPlot(clus,num.pc=pcComp)
    }
	 #as of april 2019, fond a few tweaks in applying a uniform pipeline in this function. (assigning var genes and using same dimvec)
	 addmetadata<-function(seuratObj,dimvec=1:20) {
	 	#message('cell cycle')
	 	seuratObj<-CellCycleScoring(object = seuratObj, s.genes = toupper(s.genes), g2m.genes = toupper(g2m.genes), set.ident = FALSE)
	 	message('cluster0')
	 	seuratObj <- FindClusters(seuratObj, dims.use = dimvec, genes.use =seuratObj@var.genes, resolution = 0,print.output = FALSE, save.SNN = TRUE,force.recalc =T)
	 	message('cluster tree')
	 	for (res in seq(0.1,1,0.1)) {
             message("res",res)
             seuratObj <- FindClusters(seuratObj, dims.use = dimvec, resolution = res, genes.use =seuratObj@var.genes, print.output = FALSE,reuse.SNN=T,force.recalc =T)
         }
	 	message('graphing lies')
	 	seuratObj<- RunTSNE(seuratObj, dims.use = dimvec, genes.use =seuratObj@var.genes, do.fast=T, check_duplicates = FALSE)
	 	return(seuratObj)
	 }


 addmetadata_cca<-function(seuratObj,alignmentDim=1:20) {
	 	#message('cell cycle')
	 	#seuratObj<-CellCycleScoring(object = seuratObj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
	 	message("aligning subspace in ", max(alignmentDim) ," dimensions. stay cool")
	 	seuratObj <- AlignSubspace(object = seuratObj, reduction.type = "cca", grouping.var = "isElfMouse", 
    dims.align = 1:20)
	message('graphing lies')
	 	seuratObj<- RunTSNE(seuratObj, dims.use = 1:20,reduction.type = "cca.aligned", do.fast=T,check_duplicates = FALSE)
		 message('cluster0')
	 	seuratObj <- FindClusters(seuratObj, dims.use = 1:20,reduction.type = "cca.aligned", resolution = 0,print.output = FALSE, save.SNN = TRUE)
	 	message('cluster tree')
	 	for (res in seq(0.1,2.8,0.1)) seuratObj <- FindClusters(seuratObj, resolution = res, reduction.type = "cca.aligned", print.output = FALSE)
return(seuratObj)
	 }

plotTraj<-function(monoclObj,genevec,filename=NULL,metadataDesired=NULL,...){
message('depending on your gene list, this might take ages \n you have been warned')
#pdf(file=paste0("plots/",filename,".pdf"))
tmp<-monoclObj
tmp <- setOrderingFilter(tmp, genevec)
#print(plot_ordering_genes(tmp))
#ncenters <- length(unlist(sapply(metadataDesired, function(x) unique(as(phenoData(HSMM),'data.frame')[,x]))))
tmp <- reduceDimension(tmp, max_components=2,method = 'DDRTree')#,auto_param_selection = F, nCenter = ncenters)
tmp <- orderCells(tmp, tmp <- orderCells(tmp, reverse=FALSE))
#sapply(c(metadataDesired,"Pseudotime"), function(x) print(plot_cell_trajectory(tmp, color_by=x)))
return(tmp)
#dev.off()
}

###########plotting!

compPlot<-function(seuratObj=epidat,ccaObj=epi_cca,var_method="pca",groupVar="isElfMouse") {
p1 <- DimPlot(object = seuratObj, reduction.use = var_method, group.by = groupVar, pt.size = 0.5, 
    do.return = TRUE,plot.title='seurat')
p2 <- DimPlot(object = ccaObj, reduction.use = var_method, group.by = groupVar, pt.size = 0.5, 
    do.return = TRUE,plot.title='cca')
plot_grid(p1, p2)
}	 

plotPie<-function(seuratObj=epidat,primeFactor="isElfMouse",res=0.7){
    x<-table(seuratObj@meta.data[,primeFactor],seuratObj@meta.data[,paste0('res.',as.character(res))])
    x<-x/rowSums(x)
     for(i in 1:ncol(x)) pie(x[,i],main=paste('res=',res,'clus=',colnames(x)[i]))
    plot(1)
return(x)
}

plot_spec_Pie<-function(seuratObj=epidat,primeFactor="isElfMouse",res=0.7,clusters){
    x<-table(seuratObj@meta.data[seuratObj@meta.data[,paste0('res.',as.character(res))] %in% clusters,primeFactor],seuratObj@meta.data[seuratObj@meta.data[,paste0('res.',as.character(res))] %in% clusters , paste0('res.',as.character(res))])
    x<-x/rowSums(x)
    for(i in 1:ncol(x)) pie(x[,i],main=paste('res=',res,'clus=',colnames(x)[i]))
    plot(1)
return(x)
}


	 plotstuff=function(seuratObj){
	 	PCElbowPlot(seuratObj,num.pc=100) 
	 	PCElbowPlot(seuratObj,num.pc=50) 
	 	#JackStrawPlot(alldat, PCs = 1:30)
	 	PCAPlot(seuratObj,group.by='res.0')
	 	PCAPlot(seuratObj,group.by='res.0.1')
	 	#PCAPlot(seuratObj,group.by='res.1')
	 	#PCAPlot(seuratObj,group.by='res.2.2')
	 	PCAPlot(seuratObj,group.by='orig.ident')
	 	FeaturePlot(seuratObj,c("percent.mito",'nGene','nUMI','CD44'),cols.use=c('grey','blue'))
	 	TSNEPlot(seuratObj,group.by='orig.ident')
	 	TSNEPlot(seuratObj,group.by='res.0')
	 	TSNEPlot(seuratObj,group.by='res.0.1')
	 	#TSNEPlot(seuratObj,group.by='res.1')
	 	#TSNEPlot(seuratObj,group.by='res.2.2')
	 }

plo<-function(x="gorsmy_dimRed_x",y="gorsmy_dimRed_y",contVar){
    #print(grep("_x|_y",colnames(epidat@meta.data),value=T))
    #print(grep("HALLMARK",colnames(epidat@meta.data),value=T))
   require(ggplot2)
ggplot(epidat@meta.data,aes_string(x=x , y=y, color=contVar)) + geom_point(alpha=.5) +facet_wrap(~isElfMouse) + scale_colour_gradient2(low = "blue", high = "red") + theme(legend.position="top",legend.text=element_text(size=5))
}

plotMetaScore<-function(seuratObj=alldat,genelist){
    require(Seurat) # for playing with object
    require(ggplot2) # just to make sure can plot
    require(matlab) # for the jet color palette
    require(dplyr) # for the gene contribution transformation
    message(table(genelist %in% rownames(seuratObj@scale.data))['TRUE']," out of ",length(genelist)," genes entered were used to generate score\n") # just shows how many are actually contributing to the score
    hm<-matrixStats::rowVars(seuratObj@scale.data[rownames(seuratObj@scale.data) %in% genelist,]) #
    names(hm)<-rownames(seuratObj@scale.data[rownames(seuratObj@scale.data) %in% genelist,])
    message("top contributing genes (by percentage) contributing to signature")
    print((hm*100/sum(hm)) %>% sort(decreasing=T) %>% signif(2) %>% head(10) )
    print(ggplot(data.frame(cbind(seuratObj@meta.data,seuratObj@dr$tsne@cell.embeddings),metascore=colSums(seuratObj@scale.data[rownames(seuratObj@scale.data) %in% genelist,])),aes(x=emt.emt_DDRTREE.1,y=emt.emt_DDRTREE.2,color=metascore)) + geom_point(alpha=.5) + scale_colour_gradientn(colours = jet.colors(7))) #generates the object on the fly and plots it
}


	 
	 ###########enriching
	 getEntrez<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=seuratobj@var.genes,mart=human))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=seuratobj@raw.data,mart=human))
	 }
	 }#pulls available data for all variable genesgenes in a seurat object takes a whilew

	 splMar<-function(marker) lapply(split(marker,marker$avg_logFC>0), function(x) mkEnt(split(x$gene,x$cluster))) # splits a FindAllMarkers object into cluster and up/down

	 mkEnt<-function(listOids){
	     require(clusterProfiler)
	     require(ReactomePA)
	     return(lapply(listOids,function(x) as.character(na.omit(hasEntrez$entrezgene[match(x,hasEntrez$hgnc_symbol)]))))
	  } # returns a trimmed list of entrez ids from a list of gene symbols

	 reactomeClusts<-function(mkEntObj) {
	     require(clusterProfiler)
	     require(ReactomePA)
	     message("making object, takes a while")
	     return(compareCluster(geneCluster = mkEntObj, fun = "enrichPathway",organism='human',universe=as.character(na.omit(hasEntrez$entrezgene)),pvalueCutoff=0.05,readable=T)
	 )
	 } #finds enriched reactome clusters for a list of entrez ids requires a background getEntrez object called hasEntrez

	 hallmarkClusts<-function(mkEntObj) {
	     require(clusterProfiler)
	     require(ReactomePA)
	     message("making object, takes a while")
	     return(compareCluster(geneCluster = mkEntObj, fun = "enricher", organism='human', universe=as.character(na.omit(hasEntrez$entrezgene)), pvalueCutoff=0.05, TERM2GENE=hallmark)
	 )
	 }#finds enriched hallmark clusters for a list of entrez ids requires a background getEntrez object called hasEntrez and a hallmark object from hallmark_mouse_TERM2GENE.rdata'


	 keggClusts<-function(mkEntObj) {
	     require(clusterProfiler)
	     require(ReactomePA)
	     message("making object, takes a while")
	     return(compareCluster(geneCluster = mkEntObj, fun = "enrichKEGG",organism='mouse',universe=as.character(na.omit(hasEntrez$entrezgene)),pvalueCutoff=0.05)
)
}


###cellphone
getHumanEns<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
		 human<- useMart(biomart='ensembl', dataset = "hsapiens_gene_ensembl")
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=seuratobj@var.genes,mart=human))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=rownames(seuratobj@data),mart=human))
	 }}



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
