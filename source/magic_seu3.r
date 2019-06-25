module load briglo/miniconda/3 
source activate magic
module load briglo/R/3.6.0
R

library(Seurat)
library(Rmagic)
library(ggplot2)
library(Matrix)
library(viridis)
library(reticulate)
libary(Matrix)
library(phateR)

use_condaenv("magic")
load('r_objects/190516_CCAintegrateAllCells.rdata')

dat<-t(integrated@assays$RNA@counts)
metadat<-integrated@meta.data
#rm(alldat)
keep_cols <- colSums(dat > 0) > 10
dat <- dat[,keep_cols]

#   keep_rows <- Matrix::rowSums(dat) > 1000 & Matrix::rowSums(dat) < 5000
#   dat <- dat[keep_rows,]
# ggplot() +
#   geom_histogram(aes(x=Matrix::rowSums(data), bins=50)) +
#   geom_vline(xintercept = 1000, color='red')

  keep_rows <- Matrix::rowSums(dat) > 1000 & Matrix::rowSums(dat) < 25000
  dat <- dat[keep_rows,]
dat <- library.size.normalize(dat)
dat <- sqrt(dat)

data_MAGIC_all <- magic(dat, k=15, genes=c("ZEB1","AR","SNAI1",'SNAI2',"VIM","CD44","LCN2","ITGB4"))

data_MAGIC_all <-magic(dat, k=15, init=data_MAGIC)
magicAll<-data_MAGIC_all$result
save(magicAll, file="r_objects/190524_magicALLgenes.rdata")
data_MAGIC_PCA <- magic(dat, genes="pca_only", 
                         t=4, init=data_MAGIC)

data_PHATE <- phate(dat, knn=4, decay=100, t=20)

save(data_MAGIC,data_MAGIC_PCA,data_PHATE,file="r_objects/190524_magicPHATEall.rdata")

colnames(data_MAGIC$result)<-paste0("magicEMT_",colnames(data_MAGIC$result))
colnames(data_MAGIC_PCA$result)<-paste0("magicEMT_",colnames(data_MAGIC_PCA$result))
summ<-data.frame(data_MAGIC$result,data_MAGIC_PCA$result[,1:3],data_PHATE$embedding)
integrated@meta.data<-cbind(integrated@meta.data,summ[rownames(integrated@meta.data),])

save(integrated,file='r_objects/190516_CCAintegrateAllCells.rdata')


cors_AR<-apply(magicAll,2,function(x) cor(x,magicAll$AR,method='spearman'))

cors_CD44<-apply(magicAll,2,function(x) cor(x,magicAll$CD44,method='spearman'))

save(cors_AR,cors_CD44, file="190524_magicCors_ARcd44.rdata")

magicCor<-cor(magicAll,method='spearman')
save(magicCor, file="190524_magicCors.rdata")
