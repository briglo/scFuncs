mklasso<-function(scPredObj) {
    lapply(scPredObj, function(x) {
    nPredSubpop=length(unique(grep("subpop",unlist(x[[1]]$ElasticNetPredict),value=T)))
    ngroup<-c("#7570b3","#1b9e77","#e7298a")
    return(lapply(1:length(x), function(y){
        reformat_LASSO(c_selectID=y,mp_selectID = 1,LSOLDA_dat=x[[y]],nPredSubpop =nPredSubpop,Nodes_group=ngroup[y],nboots=3)}))})}


plotLasso<- function(mklassoObj) {
    lapply(mklassoObj, function(lassoObj) {
    require(networkD3)
    require(dplyr)
    combined <- do.call(rbind,lassoObj)
    nboots = 3
    combined_D3obj <-list(Nodes=combined[,(nboots+3):(nboots+4)],
         Links=combined[,c((nboots+2):(nboots+1),ncol(combined))]) 
    combined <- combined[is.na(combined$Value) != TRUE,]

    Node_source <- as.vector(sort(unique(combined_D3obj$Links$Source)))
    Node_target <- as.vector(sort(unique(combined_D3obj$Links$Target)))
    Node_all <-unique(c(Node_source, Node_target))

    #assign IDs for Source (start from 0)
    Source <-combined_D3obj$Links$Source
    Target <- combined_D3obj$Links$Target

    for(i in 1:length(Node_all)){
    Source[Source==Node_all[i]] <-i-1
    Target[Target==Node_all[i]] <-i-1
    }
    combined_D3obj$Links$Source <- as.numeric(Source)
    combined_D3obj$Links$Target <- as.numeric(Target)
    combined_D3obj$Links$LinkColor <- combined$NodeGroup

    #prepare node info 
    node_df <-data.frame(Node=Node_all)
    node_df$id <-as.numeric(c(0, 1:(length(Node_all)-1)))

    n <- length(unique(node_df$Node))
    getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    Color = getPalette(n)
    node_df$color <- Color

return(sankeyNetwork(Links =combined_D3obj$Links, Nodes = node_df,
        Value = "Value", NodeGroup ="color", LinkGroup = "LinkColor",
        NodeID="Node", Source="Source", Target="Target", fontSize = 22))
    })
}
