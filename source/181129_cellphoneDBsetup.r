setting up for cellphone  

#in R
###EPI
load("r_objects/181005_epi_ano_monocle_tern_cca_POISSONgsva.rdata")
epiID<-data.frame(row.names=rownames(epidat@meta.data),res=paste0("epi_",epidat@meta.data$res.0.7))

###FIB
load("/Users/briglo/Dropbox/PYMT scRNAseq paper 2018/180725_fib_ano.rdata")
fibID<-data.frame(row.names=rownames(fibdat@meta.data),res=paste0("fib_",fibdat@meta.data$res.0.1))

clusterdat<-data.frame(rbind(epiID,fibID))
clusterdat$cellID=rownames(clusterdat)
clusterdat$cellType=as.character(clusterdat$res)
save(clusterdat,file="r_objects/181128_IDsForCellPhone.rdata")

#getData
tda<-SubsetData(ALLs15,cells.use=clusterdat$cellID)

getEntrez<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
		 human<- useMart(biomart='ensembl', dataset = "hsapiens_gene_ensembl")
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=toupper(seuratobj@var.genes),mart=human))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=toupper(rownames(seuratobj@data)),mart=human))
	 }
	 }#p
mappers<-getEntrez(tda,'all') #Elf5 isnt a var.gene(WTF)

 da<-tda@data[toupper(rownames(tda@data)) %in% mappers$hgnc_symbol[!is.na(mappers$ensembl_gene_id)],]
 m<-match(toupper(rownames(da)),mappers$hgnc_symbol)
 rownames(da)<-mappers$ensembl_gene_id[m]
 save(tda,da,mappers,clusterdat,file="r_objects/181129_cellPhoneIntermediates.rdata")

 cd<-data.frame(Cell=clusterdat[colnames(da),'cellID'],cell_type=clusterdat[colnames(da),'cellType'])

 write.table(data.matrix(da),file="fib_epi_counts.txt",quote=F,sep='\t')
 write.table(cd,file="fib_epi_meta.txt",quote=F,sep='\t',row.names=F)

###############on cluster
# module load briglo/miniconda/3
# source activate cellphone
# qsub -V -cwd -b y -j y -pe smp 8 -N cpdb_1 "cellphonedb method statistical_analysis fib_epi_meta.txt fib_epi_counts.txt --project-name fibEpi10 --threshold 10 --threads 15"
# EASY!!!
d#########################

This is lifted from DotPlot in Seurat
p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot,
        y = id)) + geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
        scale.func(range = c(0, dot.scale), limits = c(scale.min,
            scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())



		load("r_objects/181129_cellPhoneIntermediates.rdata")

		#to to : remap Ens IDS BACK to mouse symbols
		do Dotplot of neoclusters 
load("r_objects/181129_cellPhoneIntermediates.rdata")
tda@meta.data$cellphoneDB_id<-clusterdat[rownames(tda@meta.data),'cellType']
tda<-SetIdent(tda,ident.use=tda@meta.data$cellphoneDB_id)

cd("cellphoneDB/fibEpi/")
ms<-read.table("significant.txt",header=T,stringsAsFactors=F,sep='\t')

x<-strsplit(ms$interacting_pair,"_")
df<-tools::toTitleCase(tolower(unlist(lapply(x, function(x) x[1]))))
ip<-tools::toTitleCase(tolower(unlist(lapply(x, function(x) x[2]))))
tddf<-ddf[!(grepl(" ",ddf$df) | grepl(" ",ddf$ip)),]
pdf("out.pdf")
apply(tddf,1,function(x) print(DotPlot(tda,x)))
dev.off()