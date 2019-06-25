 read10xCounts<-function (samples, col.names = FALSE, ...)
{
    nsets <- length(samples)
    full_data <- vector("list", nsets)
    gene_info_list <- vector("list", nsets)
    cell_info_list <- vector("list", nsets)
    for (i in seq_len(nsets)) {
        run <- samples[i]
        barcode.loc <- file.path(run, "barcodes.tsv")
        gene.loc <- file.path(run, "features.tsv")
        matrix.loc <- file.path(run, "matrix.mtx")
        data_mat <- read10xMatrix(matrix.loc, ...)
        cell.names <- readLines(barcode.loc)
        gene.info <- read.delim(gene.loc, header = FALSE, colClasses = "character",
            stringsAsFactors = FALSE, quote = "", comment.char = "")
        full_data[[i]] <- data_mat
        gene_info_list[[i]] <- gene.info
        cell_info_list[[i]] <- DataFrame(Sample = run, Barcode = cell.names)
    }
    if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
        stop("gene information differs between runs")
    }
    gene_info <- gene_info_list[[1]]
    colnames(gene_info) <- c("ID", "Symbol")
    full_data <- do.call(cbind, full_data)
    rownames(full_data) <- gene_info$ID
    cell_info <- do.call(rbind, cell_info_list)
    if (col.names && nsets == 1L) {
        colnames(full_data) <- cell_info$Barcode
    }
    SingleCellExperiment(list(counts = full_data), rowData = gene_info,
        colData = cell_info)
}

read10xMatrix<-function (file, hdf5.out = FALSE, chunk.size)
{
    if (is.character(file)) {
        fhandle <- file(file, open = "r")
        on.exit(close(fhandle))
    }
    else {
        fhandle <- file
    }
    if (!hdf5.out) {
        out <- readMM(file)
        out <- as(out, "dgCMatrix")
        return(out)
    }
    chunk.size <- as.integer(chunk.size)
    if (chunk.size <= 0L) {
        stop("chunk size should be a positive integer")
    }
    type <- readLines(fhandle, 1)
    if (!grepl("%%MatrixMarket matrix coordinate (real|integer) general",
        type)) {
        stop("expected numeric/integer matrix in MatrixMarket format")
    }
    if (strsplit(type, split = " ")[[1]][4] == "integer") {
        input.x <- integer()
    }
    else {
        input.x <- double()
    }
    dims <- scan(fhandle, what = integer(), nmax = 3, comment.char = "%",
        quiet = TRUE)
    nr <- dims[1]
    nc <- dims[2]
    nz <- dims[3]
    out <- .Call(cxx_load_tenx_to_hdf5, fhandle, chunk.size,
        input.x, nr, nc, nz)
    return(out)
}