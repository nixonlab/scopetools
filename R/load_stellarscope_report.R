#' Load count matrix from Stellarscope output files
#'
#' @param dirname Stellarscope output directory.
#' @param count_file_suffix File name suffix for count matrix in MTX format. Default is 'TE_counts.mtx'.
#' @param feature_file_suffix File name suffix for feature list. Default is 'features.tsv'.
#' @param barcode_file_suffix File name suffix for barcode list. Default is 'barcodes.tsv'.
#'
#' @return A sparse matrix of of class "dgTMatrix" with features as rows and cells as columns.
#' @export
#'
#' @examples
load_stellarscope_report <-
  function(dirname,
           count_file_suffix = 'TE_counts.mtx',
           feature_file_suffix = 'features.tsv',
           barcode_file_suffix = 'barcodes.tsv') {
    tecounts_mtx <- Sys.glob(file.path(dirname, paste0('*', count_file_suffix)))
    features_tsv <- Sys.glob(file.path(dirname, paste0('*', feature_file_suffix)))
    barcodes_tsv <- Sys.glob(file.path(dirname, paste0('*', barcode_file_suffix)))

    ret <- Matrix::readMM(file = tecounts_mtx)
    ret <- Matrix::t(ret)
    features <- read.delim(file = features_tsv,
                           header = FALSE,
                           stringsAsFactors = FALSE)
    barcodes <- read.delim(file = barcodes_tsv,
                           header = FALSE,
                           stringsAsFactors = FALSE)

    rownames(ret) <- features$V1
    colnames(ret) <- barcodes$V1

    ret
  }
