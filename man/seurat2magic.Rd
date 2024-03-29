% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratFuncs.R
\name{seurat2magic}
\alias{seurat2magic}
\title{seurat2magic}
\usage{
seurat2magic(seuratObj, geneList = c("ZEB1", "SNAI1", "SNAI2", "VIM"),
  return.raw = FALSE)
}
\arguments{
\item{seuratObj}{the Seurat object}

\item{geneList}{a vector of gene symbols}

\item{return.raw}{logical, return MAGIC object for init, default FALSE}
}
\value{
a matrix of MAGIC expression
}
\description{
turns a Seurat object into a default magic/phate result
}
\examples{
magicRand<-seurat2magic(seuratObj=integrated,geneList=sample(rownames(integrated@assays$SCT@scale.data),10),return.raw=F)

}
