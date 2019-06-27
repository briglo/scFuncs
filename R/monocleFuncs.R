#' plotTraj
#'
#' turns a monocle object into a monocle object ordered by a gene list
#'
#' @param monocleObj the monocle object 
#' @param geneVec the character vector of genes to assemble cells on
#'
#' @return a monocle object with Pseudo vectors and state in me
#'
#' @examples
#' NULL
#'
#' @export

plotTraj<-function(monocleObj,geneVec){
message('depending on your gene list, this might take ages \n you have been warned')
tmp<-monocleObj
tmp <- setOrderingFilter(tmp, geneVec)
tmp <- reduceDimension(tmp, max_components=2,method = 'DDRTree')
tmp <- orderCells(tmp, reverse=FALSE)
return(tmp)
}