#' Load count matrix from Stellarscope output files
#'
#' @param stellarscope_dir Stellarscope output directory.
#' @param exp_tag Experiment tag used in stellarscope run. Default is "stellarscope".
#' @param count_file_suffix File name suffix for count matrix in MTX format.
#'                          Default is 'TE_counts.mtx'.
#' @param feature_file_suffix File name suffix for feature list. Default is
#'                            'features.tsv'.
#' @param barcode_file_suffix File name suffix for barcode list. Default is
#'                            'barcodes.tsv'.
#'
#' @return A sparse matrix of of class "dgTMatrix" with features as rows and
#'         cells as columns.
#' @export
#'
#' @examples
load_stellarscope_report <-
  function(stellarscope_dir,
           exp_tag = 'stellarscope',
           count_file_suffix = 'TE_counts.mtx',
           feature_file_suffix = 'features.tsv',
           barcode_file_suffix = 'barcodes.tsv') {

    prefix <- file.path(stellarscope_dir, sprintf('%s-', exp_tag))
    tecounts_mtx <- Sys.glob(paste0(prefix, '*', count_file_suffix))
    features_tsv <- Sys.glob(paste0(prefix, '*', feature_file_suffix))
    barcodes_tsv <- Sys.glob(paste0(prefix, '*', barcode_file_suffix))

    stopifnot(length(tecounts_mtx) == 1)
    stopifnot(length(features_tsv) == 1)
    stopifnot(length(barcodes_tsv) == 1)

    counts <- Matrix::readMM(file = tecounts_mtx)
    # counts <- Matrix::t(counts)
    features <- read.delim(file = features_tsv,
                           header = FALSE,
                           stringsAsFactors = FALSE)
    barcodes <- read.delim(file = barcodes_tsv,
                           header = FALSE,
                           stringsAsFactors = FALSE)
    if(all(dim(counts) == c(nrow(barcodes), nrow(features)))) {
        counts <- counts # print('input is barcode x feature (old way')
    } else if (all(dim(counts) == c(nrow(features), nrow(barcodes)))) {
        # print('input is feature x barcode (new way)')
        counts <- Matrix::t(counts)
    } else {
      msg <- paste0(
        sprintf("Matrix dimensions (%d, %d) ", nrow(counts), ncol(counts)),
        sprintf("do not match barcode (%d) ", nrow(barcodes)),
        sprintf("or feature (%d) dimensions.", nrow(features))
      )
      stop(msg)
    }
    rownames(counts) <- barcodes$V1
    colnames(counts) <- features$V1

    counts
  }
