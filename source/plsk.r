 
x<-do.call(rbind,lapply(strsplit(trimdat@meta.data$orig.ident,"\\."), function(x) x[1:2]))
colnames(x)<-c('id','location')
x[,"id"]<-gsub("M","",x[,"id"])
x[,"location"]<-gsub("LN","Lymph",x[,"location"])
trimdat@meta.data<-data.frame(cbind(trimdat@meta.data,x))



plsk<-function(seuratObj,idvar=c("varRes.0.3","emt_res.0.3")){
require(flipPlots)
message('try install_github("Displayr/flipPlots") if this doesnt work')
require(dplyr)
seuratObj@meta.data[,match(idvar,colnames(seuratObj@meta.data))] %>% arrange(.[,1]) %>% group_by_all() %>% summarise(COUNT = n()) ->> my.data
 #my.data<-as.factor(my.data[,1])

SankeyDiagram(my.data[, -grep("COUNT",colnames(my.data))],link.color = "Source",weights = my.data$COUNT,,max.categories = 1000)

}

plsk(seuratObj=trimdat,idvar=c("orig.ident","SCOREclust"))

SankeyDiagram(df[, -grep("COUNT",colnames(df))],link.color = "Source",weights = df$COUNT)
return(my.data)
}
um<-plsk(smd[[1]])



plsk<-function(seuratObj,idvar=c("varRes.0.3","emt_res.0.3")){
require(flipPlots)
message('try install_github("Displayr/flipPlots") if this doesnt work')
require(dplyr)
seuratObj[,match(idvar,colnames(seuratObj))] %>% arrange(.[,1]) %>% group_by_all() %>% summarise(COUNT = n()) ->> my.data
 #my.data<-as.factor(my.data[,1])

SankeyDiagram(my.data[, -grep("COUNT",colnames(my.data))],link.color = "Source",weights = my.data$COUNT,,max.categories = 1000)

}