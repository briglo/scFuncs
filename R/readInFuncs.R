#' Modified read 10x function 
#'
#' Made when the on board function in seurat was broken
#'
#' @param data.dir a directory path containing cellranger output
#'
#' @return a matrix of counts 
#'
#' @examples
#' raw.data <- Read10X("PATH/TO/OUT")
#'
#' @export

Read10X<- function (data.dir = NULL) {
    full.data <- list()
    for (i in seq_along(data.dir)) {
        run <- data.dir[i]
        if (!dir.exists(run)) {
            stop("Directory provided does not exist")
        }
        if (!grepl("\\/$", run)) {
            run <- paste(run, "/", sep = "")
        }
        barcode.loc <- paste0(run, "barcodes.tsv")
        gene.loc <- paste0(run, "features.tsv")
        matrix.loc <- paste0(run, "matrix.mtx")
        if (!file.exists(barcode.loc)) {
            stop("Barcode file missing")
        }
        if (!file.exists(gene.loc)) {
            stop("Gene name file missing")
        }
        if (!file.exists(matrix.loc)) {
            stop("Expression matrix file missing")
        }
        data <- readMM(file = matrix.loc)
        cell.names <- readLines(barcode.loc)
        gene.names <- readLines(gene.loc)
        if (all(grepl(pattern = "\\-1$", x = cell.names))) {
            cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                FUN = ExtractField, field = 1, delim = "-")))
        }
        rownames(x = data) <- make.unique(names = as.character(x = sapply(X = gene.names, 
            FUN = ExtractField, field = 2, delim = "\\t")))
        if (is.null(x = names(x = data.dir))) {
            if (i < 2) {
                colnames(x = data) <- cell.names
            }
            else {
                colnames(x = data) <- paste0(i, "_", cell.names)
            }
        }
        else {
            colnames(x = data) <- paste0(names(x = data.dir)[i], 
                "_", cell.names)
        }
        full.data <- append(x = full.data, values = data)
    }
    full.data <- do.call(cbind, full.data)
    return(full.data)
}


#' helper function for Read10X 
#'
#' turns barcodes.csv into "useful" ids
#'
#' @param string i dont know, ou just need it
#'
#' @return a vector of ids
#'
#' @examples
#' NULL
#'
#' @export

ExtractField <- function (string, field = 1, delim = "_") 
{
    fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), 
        split = ",")))
    if (length(x = fields) == 1) {
        return(strsplit(x = string, split = delim)[[1]][field])
    }
    return(paste(strsplit(x = string, split = delim)[[1]][fields], 
        collapse = delim))
}


#' makeSeuratList 
#' 
#' reads multiple 10X directories into a list of Seurat Objects
#'
#' @param data.dirs a vector of  paths to cellranger output
#' @param min.cells minimum cells expressing gene to retain (per data.dir), defaults to 3
#' @param min.features minimum genes expressed per cell to retain (per data.dir), defaults to 200
#'
#' @return a list of seurat objects
#'
#' @examples
#' snam<-dir()
#' id <- makeSeuratList(snam)
#'
#' @export

makeSeuratList<-function(data.dirs=NULL,min.cells=3,min.features=200){
return(lapply(data.dirs, function(x) {
rd <- Read10X(x)
rd<-rd[grepl("hg19",rownames(rd)),]
rownames(rd)<-gsub("hg19_","",rownames(rd))
return(CreateSeuratObject(counts=rd,project=gsub('^M',"",x), min.cells=3,min.features=200))
}))
}
