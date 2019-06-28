# scFuncs
a repo of functions that i used alot playing with scRNA
*a product of three separate scRNA expression studies and what sorta remained useful over all three with a lot of software changes between 2017-2019*
**most of the shitty hardcoding is based on the nature of the data in this project, i did my best**

* Rationalised R functions in /R Recently updated for Seurat3,
* workflow employed at the end of my interaction with christine chaffer is in /extraScripts
* R objects that were particularly painful to get are in  /share/ClusterShare/thingamajigs/SCCG/briglo/chrcha_mets
* The actual working scripts used to make this repo are in /source

##a simple (if time comsuming) Seurat workflow
```R
setwd("PATH/TO/DATA")
datadirs<- paste0(dir(),"/outs/filtered_feature_bc_matrix")
bysample<-makeSeuratList(datadirs)
integrated<-buildMasterSeurat(bysample)

# look at population composition 
plotSankey(integrated,c("orig.ident","Phase"))

#look at reactome enrichment for clusters
markers<-makeReactomePipe(integrated,"SCT_snn_res.0.8")
for i in 1:length(markers$CP_result)) dotplot(markers$CP_result[[i]]) + ggtitle(names(markers$CP_result)[i])

# look at combined expression of some marker genes in UMAP space
integrated<-makeMetaScore(integrated,markers$markers$gene[1:50])
ggplot(integrated@meta.data,aes(x=UMAP_1,y=UMAP_2,color=metascore)) + geom_point(alpha=.5) + scale_colour_gradientn(colours = jet.colors(7))
```