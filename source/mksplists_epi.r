setwd('/share/ScratchGeneral/PYMT/briglo')
load("r_objects/new_epi_ano.rdata")
library(Seurat)
library(Matrix)
library(dplyr)
library(reshape2)
source('source/clustering_tree.R')
plotClusteringTreeSeurat(epidat)
 dev.copy2pdf(file="plots/180108_epiTrees.pdf")

##needs tp be redoone
spl<-list(
res.0.1=list(c(0,2),c(3,6)),
res.0.25=list(c(1,6)),
res.0.3=list(c(0,2)),
res.0.35=list(c(0,2)),
res.0.45=list(c(0,8)),
res.0.6=list(c(1,5)),
res.0.7=list(c(1,3)),
res.0.8 =list(c(1,3),c(0,6)),
res.1 =list(c(3,6 ),c(2,9)),
res.1.2 =list(c(2,4),c(0,12)),
res.1.6 =list(c(10,8),c(12,9),c(13,7)),
res.1.8 =list(c(10,9)),
res.2 =list(c(12,18 )),
res.2.2 =list(c(12,4),c( 14,20 ),c(0,1 ))
)

treeSplits<-vector(mode='list',length=length(spl))
for (i in 1:length(spl)) {
    print(names(spl)[i])
epidat<-SetIdent(epidat,ident.use=epidat@meta.data[,names(spl)[i]])
treeSplits[[i]]<-lapply(spl[[i]], function(x) FindMarkers(epidat,ident.1=x[1],ident.2=x[2],min.pct = 0.25, thresh.use = 0.25))
}
names(treeSplits)<-names(spl)

###done interactively
res<-paste0('res.',c(0.05))
primeSplits<-vector(mode='list',length=2)
names(primeSplits)<-res
for (i in 1) {
	epidat<-SetIdent(epidat,ident.use=epidat@meta.data[,res[i]])
	tmp<-FindAllMarkers(epidat,min.pct = 0.25, thresh.use = 0.25)
	tmp %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
	pdf(file=paste0("180103_",res[i],'_heatmap.pdf'),width=20,height=12.5)
	DoHeatmap(epidat,genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
	dev.off()
	primeSplits[[i]]<-tmp
}

tmp2<-epidat@meta.data
tmp2<-tmp2[,grep("res",colnames(tmp2))] 

geno_nums<-apply(tmp2,2, function(x) table(epidat@meta.data$isElfMouse,x))
 geno_props<-lapply(geno_nums, function(x) x/rowSums(x))

save(geno_nums,geno_props,treeSplits,primeSplits,spl,file='180108_EpiMarkerLists.rdata')
##pie
plo<-function(res,clust){
x<-geno_props[[res]]
x<-x[,match(clust,colnames(x))]
print(x)
apply(x,2,pie)
plot(1)
} #used for rest using spl object
par(mfrow=c(7,12))
pie(c(1,1))
plot(1)
apply(geno_props[['res.0.05']],2,pie)
plot(1)
plot(1)
for (i in 1:length(spl)) for(j in 1:length(spl[[i]])) plo(names(spl)[i],spl[[i]][[j]])
dev.copy2pdf(file="180110_epiPie.pdf")

# primary split heatmap and lists
lapply(primeSplits, function(x){
primeSplits[[i]] %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
DoHeatmap(epidat,genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
})
dev.copy2pdf(file='180110_epi_primeSplits_heatmap.pdf')

lapply(split(top10$gene,top10$cluster),function(x) data.frame(x))
#regexped to make lists

##splitlists
lapply(treeSplits, function(x) lapply(x, function(y){
y$gene<-rownames(y)
y %>% top_n(10, avg_logFC) -> top10
data.frame(top10$gene)
})
)

#want to make sankey of elf5 vs wildtype
##based on this https://plot.ly/r/sankey-diagram/
#cant get the mother fucker working
load("r_objects/new_epi_ano.rdata")
source("source/clustering_tree.R")
n<-getTreeNodes(tmp2[,c(1,19,26)])
e<-getTreeEdges(tmp2[,c(1,19,26)],n)


mynode = list(
      label = n$Node,
      color = rainbow(length(n$Node)),
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    )
mylink = list(
      source = e$FromCluster,
      target = e$ToCluster,
      value =  e$TransPropFrom,
      label =  NA#e$ToClus
    )
    
onode = list(
      label = json_data$data[[1]]$node$label,
      color = json_data$data[[1]]$node$color,
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    )

olink = list(
      source = json_data$data[[1]]$link$source,
      target = json_data$data[[1]]$link$target,
      value =  json_data$data[[1]]$link$value,
      label =  json_data$data[[1]]$link$label
    )

p <- plot_ly(
    type = "sankey",
    domain = c(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = "TWh",

    node=mynode,

    link=mylink
  ) %>% 
  layout(
    title = "will this fucking thing ever work",
    font = list(
      size = 10
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
)
p

#failure thus far BUT i want nodes to be nodes but i want groups to be elf5...