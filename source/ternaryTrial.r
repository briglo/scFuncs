
#for Seurat obj
library(Seurat)
load("r_objects/new_epi_ano.rdata")
load("../../extraData/GSE63310_RAW/180227_gorsmy_sigGenes.rdata")
allgenes<-unique(unlist(signatureGenes))[unique(unlist(signatureGenes)) %in% rownames(epidat@raw.data)]
expdat<-FetchData(epidat,allgenes)>0
signatureGenes<-lapply(signatureGenes, function(x) x[x %in% allgenes])
scores<-data.frame(do.call(cbind, lapply(signatureGenes, function(x) apply(expdat[,x],1,function(y) sum(y)/length(y)))))
clusterano<-FetchData(epidat,c('res.1.2',"isElfMouse","Phase"))
limma::plotDensities(scores[,1:3])
scores<-cbind(scores,clusterano[rownames(scores),])
scores$res.1.2<-as.numeric(as.character(scores$res.1.2))
dev.copy2pdf(file="180417_gorsmy_davgal_scorespread.pdf")

save(scores,signatureGenes,expdat,file='180417_gorsmy_res1.2_scores.rdata')




#this wont work if seurat is loaded... or if scMCA has tidyversed ggplot
# install.packages("ggplot2")
#  install.packages("ggtern")
# fixes it if you get stupidf environment error

#generic Plot tools
library(ggtern)
library(reshape2)

load("180227_gorsmy_res1.2_scores.rdata")
nscore=scores[c(rownames(scores)[scores$isElfMouse=="Elf5"],sample(rownames(scores)[scores$isElfMouse=="WT"],sum(scores$isElfMouse=="Elf5"))),]
ggtern(nscore,aes(LP,basal,ML,colour=factor(res.1.2))) + geom_point(alpha=1) 
ggsave('180417_gorsmy_ternres1.2_davgal.pdf')
#or
pdf(file='ternaryres.pdf')
ggtern(nscore,aes(LP,basal,ML,colour=factor(res.1))) + geom_point(alpha=.5) + facet_wrap(~isElfMouse)
dev.off()
print()

#or
mn<-melt(nscore,id.vars=c("isElfMouse","Phase","res.1"))
ggplot(mn,aes(x=value,group=isElfMouse,fill=isElfMouse)) + geom_density(alpha=.5)+ facet_grid(res.1~variable,scale='free')


ggplot(mn,aes(y=value,x=factor(res.1),fill=isElfMouse)) + geom_violin(alpha=.5) +facet_wrap(~variable,ncol=1)



#######################for single cells#################################
library(ggtern)
load("../../extraData/GSE63310_RAW/180227_gorsmy_sigGenes.rdata")
expdat<-t(data.frame(read.table("../../extraData/GSE95432_RAW/GSE95432_GeneCounts.txt",skip=1)[,-2],row.names=1))>0
sigGene_ent<-lapply(sigGene_ent, function(x) x[x %in% colnames(expdat)])
scores<-data.frame(do.call(cbind, lapply(sigGene_ent, function(x) apply(expdat[,x],1,function(y) sum(y)/length(y)))))
scores$res.1<-c(rep('basal','96'),rep('luminal',90))
limma::plotDensities(scores[,1:3])
dev.copy2pdf(file="180227_gorsmy_GSE95432_scorespread.pdf")


save(scores,sigGene_ent,expdat,file='GSE95432_180227_gorsmy_scores.rdata')

#this wont work if seurat is loaded
#generic Plot tools
library(ggtern)
library(reshape2)


load('GSE95432_180227_gorsmy_scores.rdata')
ggtern(scores,aes(LP,basal,ML,colour=res.1)) + geom_point(alpha=.5)
ggsave('GSE95432_180227_gorsmy_tern.pdf')


load("180227_gorsmy_scores.rdata")
ours<-scores
load('GSE95432_180227_gorsmy_scores.rdata')
ref<-scores
ref$isElfMouse<-ref$res.1
ref$res.1<-(-1)
comb<-data.frame(rbind(ref,ours[,1:5]))
mc<-melt(comb,id.vars=c('res.1','isElfMouse'))
ggplot(mc,aes(x=value,group=isElfMouse,fill=isElfMouse)) + geom_density(alpha=.5)+ facet_wrap(~variable,scale='free')
dev.copy2pdf(file="ScoreSpread_gorsmy.pdf")
#or
ggplot(mc,aes(y=value,x=factor(res.1),fill=factor(res.1))) + geom_violin(alpha=.5,scale='width') +facet_wrap(~variable,ncol=1)


