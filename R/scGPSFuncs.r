#' mklasso
#'
#' turns a scpred object into a LASSO object
#'
#' @param scPredObj an output of scPred
#'
#' @return LASSO object that can be tabulated or plotted
#'
#' @examples
#' NULL
#'
#' @export
mklasso<-function(scPredObj) {
    require(scGPS)
    require(locfit) #required but not dependency
    require(e1071)
    lapply(scPredObj, function(x) {
    nPredSubpop=length(unique(grep("subpop",unlist(x[[1]]$ElasticNetPredict),value=T)))
    ngroup<-c("#7570b3","#1b9e77","#e7298a")
    return(lapply(1:length(x), function(y){
        reformat_LASSO(c_selectID=y,mp_selectID = 1,LSOLDA_dat=x[[y]],nPredSubpop =nPredSubpop,Nodes_group=ngroup[y],nboots=3)}))})}

#' plotLasso
#'
#' turns a mklasso object into a sankey diagram
#'
#' @param mklassoObj an output of mklasso
#'
#' @return LASSO object that can be tabulated or plotted
#'
#' @examples
#' NULL
#'
#' @export
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


#' annotate_scGPS
#'
#' had to modify it coz the native scGPS function doesnt allow for a background gene list
#'
#' @param DEgeneList the list of genes being interrogated
#' @param pvalueCutoff self explanatory default 0.05
#' @param gene_symbol logical
#' @param species character default "human"
#' @param universe my addiditon, a vector of symbols representing the background
#'
#' @return a reactome pathway analysis
#'
#' @examples
#' NULL
#'
#' @export
annotate_scGPS<-function (DEgeneList, pvalueCutoff = 0.05, gene_symbol = TRUE, 
    species = "human",universe) 
{
    if (species == "human") {
        if (gene_symbol == TRUE) {
            convert_to_gene_ID = clusterProfiler::bitr(DEgeneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            universe_to_ID =clusterProfiler::bitr(universe,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            print("Original gene number in geneList")
            print(length(DEgeneList))
            print("Number of genes successfully converted")
            print(nrow(convert_to_gene_ID))
        }
        else {
            stop("The list must contain human gene symbols")
        }
    }
    else if (species == "mouse") {
        if (gene_symbol == TRUE) {
            convert_to_gene_ID = clusterProfiler::bitr(DEgeneList, 
                fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
            print("Original gene number in geneList")
            print(length(DEgeneList))
            print("Number of genes successfully converted")
            print(nrow(convert_to_gene_ID))
        }
        else {
            stop("The list must contain mouse gene symbols")
        }
    }
    Reactome_pathway_test <- ReactomePA::enrichPathway(gene = convert_to_gene_ID$ENTREZID, pvalueCutoff = 0.05,readable = TRUE, universe=universe_to_ID$ENTREZID)
    output_df <- as.data.frame(Reactome_pathway_test)
    return(Reactome_pathway_test)
}
