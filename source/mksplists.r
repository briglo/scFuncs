setwd('/share/ScratchGeneral/PYMT/briglo')
load("../ALLs15.rdata")
library(Seurat)
library(Matrix)
library(dplyr)
library(reshape2)

spl<-list(
res.0.15=list(c(0,1)),
res.0.2=list(c(3,9)),
res.0.25=list(c(3,7)),
res.0.3=list(c(3,7)),
res.0.35=list(c(0,2)),
res.0.6=list(c(11,12)),
res.0.7=list(c(1,4),c(12,13)),
res.0.9 =list(c(15,8)),
res.1.2 =list(c(1,8 ),c(2,5)),
res.1.4 =list(c(18,7)),
res.1.6 =list(c(15,6),c(18,21),c(11,17)),
res.1.8 =list(c(11,7),c( 21,8)),
res.2 =list(c(1,8 ),c(16,8 )),
res.2.5 =list(c(12,2),c( 26,3 ),c(0,1 ),c(21,7 ),c(19,25),c( 24,29 ))
)

treeSplits<-vector(mode='list',length=length(spl))
for (i in 1:length(spl)) {
ALLs15<-SetIdent(ALLs15,ident.use=ALLs15@meta.data[,names(spl)[i]])
treeSplits[[i]]<-lapply(spl[[i]], function(x) FindMarkers(ALLs15,ident.1=x[1],ident.2=x[2],min.pct = 0.25, thresh.use = 0.25))
}
names(treeSplits)<-names(spl)
save(treeSplits,spl,file="171221_markerSplits.rdata")

###done interactively
res<-paste0('res.',c(0.05,0.1))
primeSplits<-vector(mode='list',length=2)
names(primeSplits)<-res
for (i in 1:2) {
	ALLs15<-SetIdent(ALLs15,ident.use=ALLs15@meta.data[,res[i]])
	tmp<-FindAllMarkers(ALLs15,min.pct = 0.25, thresh.use = 0.25)
	tmp %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
	pdf(file=paste0("180103_",res[i],'_heatmap.pdf'),width=20,height=12.5)
	DoHeatmap(ALLs15,genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
	dev.off()
	primeSplits[[i]]<-tmp
}

tmp2<-ALLs15@meta.data
tmp2<-tmp2[,grep("res",colnames(tmp2))] 

geno_nums<-apply(tmp2,2, function(x) table(ALLs15@meta.data$genotype,x))
 geno_props<-lapply(geno_nums, function(x) x/rowSums(x))

save(geno_nums,geno_props,treeSplits,primeSplits,spl,file='180103_markerLists.rdata')

plo<-function(res,clust){
x<-geno_props[[res]]
x<-x[,match(clust,colnames(x))]
print(x)
apply(x,2,pie)
plot(NA)
} #used for rest using spl object


plopi<-function(res){
x<-geno_props[[paste0('res.',res)]]
par(mfrow=c(4,5))
apply(x,2,pie)
} #used for 0.05 and 0.1




# pdf(file='180103_heatmaps.pdf')
lapply(primeSplits, function(x){

x %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
DoHeatmap(ALLs15, use.raw = T,genes.use = top10$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)