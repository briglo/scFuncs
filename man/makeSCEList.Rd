% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readInFuncs.R
\name{makeSCEList}
\alias{makeSCEList}
\title{makeSCEList}
\usage{
makeSCEList(data.dirs)
}
\arguments{
\item{data.dirs}{a vector of  paths to cellranger output}

\item{min.cells}{minimum cells expressing gene to retain (per data.dir), defaults to 3}

\item{min.features}{minimum genes expressed per cell to retain (per data.dir), defaults to 200}
}
\value{
a list of seurat objects
}
\description{
reads multiple 10X directories into a list of SingleCellExperiment Objects
}
\examples{
snam<-dir()
id <- makeSeuratList(snam)

}
