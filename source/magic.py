#base on https://nbviewer.jupyter.org/github/KrishnaswamyLab/MAGIC/blob/master/python/tutorial_notebooks/emt_tutorial.ipynb
#in R
x<-t(alldat@raw.data)
keepcells<-Matrix::rowSums(x)>1500
x<-x[keepcells,]
keepgenes<-Matrix::colSums(x)>0
x<-x[,keepgenes]
writeMM(x, file="comb100.mtx")
write.table(quote=F,col.names=F,row.names=F,x=rownames(x),file='barcodes.tsv')
write.table(quote=F,col.names=F,row.names=F,x=colnames(x),file='genes.tsv')

import magic
import scprep
import scipy.io
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import feather

raw_data = scprep.io.load_mtx('/Users/briglo/Desktop/comb100.mtx',gene_names='/Users/briglo/Desktop/genes.tsv', cell_names='/Users/briglo/Desktop/barcodes.tsv', sparse=True)

scprep.plot.plot_library_size(raw_data, cutoff=1500)
#emt_data = scprep.filter.filter_library_size(emt_data, cutoff=1500)

emt_data = scprep.normalize.library_size_normalize(raw_data)
emt_data = scprep.transform.sqrt(emt_data)

magic_op = magic.MAGIC()
emt_magic = magic_op.fit_transform(emt_data, genes=['VIM', 'CDH1', 'ZEB1','LCN2'])
feather.write_dataframe(emt_magic, "/Users/briglo/Desktop/emt_magic.feather")
feather.write_dataframe( pd.DataFrame(emt_data.index), "/Users/briglo/Desktop/emt_magic_barcode.feather")


ll_magic = magic_op.transform(genes="all_genes") 
feather.write_dataframe(all_magic, "/Users/briglo/Desktop/all_magic.feather")
feather.write_dataframe( pd.DataFrame(emt_data.index), "/Users/briglo/Desktop/all_magic_barcode.feather")

# back in R
allmagic<-read_feather("~/Desktop/emt_magic.feather")
bc<-read_feather("~/Desktop/emt_magic_barcode.feather")
mdat<-data.frame(row.names=bc$`0`,allmagic)
m<-match(rownames(alldat@meta.data), rownames(mdat))
# trimdat@meta.data$magic_AR<-mdat$AR[m]
# trimdat@meta.data$magic_SNAI1<-mdat$SNAI1[m]
# trimdat@meta.data$magic_SNAI2<-mdat$SNAI2[m]
# trimdat@meta.data$magic_ZEB1<-mdat$ZEB1[m] 
# trimdat@meta.data$magic_VIM<-mdat$VIM[m] 
# trimdat@meta.data$magic_CD44<-mdat$CD44[m] 
alldat@meta.data$magic_LCN2<-mdat$LCN2[m] 