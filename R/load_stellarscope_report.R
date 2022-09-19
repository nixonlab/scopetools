#' Load count matrix from Stellarscope output files
#'
#' @param stellarscope_dir Stellarscope output directory.
#' @param exp_tag Experiment tag used in stellarscope run.
#'    Default is "stellarscope".
#' @param count_filename File name pattern for count matrix (MTX format). May
#'    contain wildcards.
#'    Default is '*TE_counts.mtx'.
#' @param feature_filename File name pattern for feature list (tab-separated,
#'    features in first column). May contain wildcards.
#'    Default is '*features.tsv'.
#' @param barcode_filename File name pattern for barcode list (tab-separated,
#'    barcodes in first column). May contain wildcards.
#'    Default is '*barcodes.tsv'.
#'
#' @return A sparse matrix of of class "dgTMatrix" with features as rows and
#'         cells as columns.
#' @export
#'
#' @examples
load_stellarscope_report <-
    function(
        stellarscope_dir,
        exp_tag = 'stellarscope',
        count_filename = '*TE_counts.mtx',
        feature_filename = '*features.tsv',
        barcode_filename = '*barcodes.tsv'
    ) {
        if(!is.null(exp_tag)) {
            prefix <- file.path(stellarscope_dir, sprintf('%s-', exp_tag))
        } else {
            prefix <- file.path(stellarscope_dir, '')
        }

        tecounts_mtx <- Sys.glob(paste0(prefix, count_filename))
        if(length(tecounts_mtx) == 0)
            stop("No count file found.")
        if(length(tecounts_mtx) > 1)
            stop(paste(
                c("Multiple matching count files:", tecounts_mtx),
                collapse='\n\t'
            ))

        features_tsv <- Sys.glob(paste0(prefix, feature_filename))
        if(length(features_tsv) == 0)
            stop("No features file found.")
        if(length(features_tsv) > 1)
            stop(paste(
                c("Multiple matching features files:", features_tsv),
                collapse='\n\t'
            ))

        barcodes_tsv <- Sys.glob(paste0(prefix, barcode_filename))
        if(length(barcodes_tsv) == 0)
            stop("No features file found.")
        if(length(barcodes_tsv) > 1)
            stop(paste(
                c("Multiple matching features files:", barcodes_tsv),
                collapse='\n\t'
            ))

        stopifnot(length(tecounts_mtx) == 1)
        stopifnot(length(features_tsv) == 1)
        stopifnot(length(barcodes_tsv) == 1)

    load_labeled_mtx(tecounts_mtx, features_tsv, barcodes_tsv)
  }


#' Load a sparse matrix with row and column labels
#'
#' Sparse matrixes in the MTX format have integer indexes and do not have
#' row or column names included. The labels may be included in separate
#' files. This function takes three input files and returns a labeled sparse
#' matrix. Dimensions of the returned sparse matrix are match the provided row
#' and column names, transposing the matrix if needed.
#'
#' @param mtx_file Sparse matrix file in MTX format
#' @param rowname_file TSV file with row labels in first column.
#' @param colname_file TSV file with column labels in first column.
#'
#' @return A sparse matrix of of class "dgTMatrix" with row and columns labels.
#' @export
#'
#' @examples
load_labeled_mtx <-
    function(mtx_file, rowname_file, colname_file) {
        m <- Matrix::readMM(file = mtx_file)
        rnames <- read.delim(file = rowname_file,
                             header = FALSE,
                             stringsAsFactors = FALSE
                             )
        cnames <- read.delim(file = colname_file,
                             header = FALSE,
                             stringsAsFactors = FALSE
                             )

        if (all(dim(m) == c(nrow(rnames), nrow(cnames)))) {
            m <- m
        } else if (all(dim(m) == c(nrow(cnames), nrow(rnames)))) {
            # matrix is transposed relative to provided rownames and colnames
            m <- Matrix::t(m)
        } else {
            stop(paste0(
                sprintf("Matrix dimensions (%d, %d) ", nrow(m), ncol(m)),
                sprintf("do not match row (%d) ", nrow(rnames)),
                sprintf("or column (%d) dimensions.", nrow(cnames))
            ))
        }
        rownames(m) <- rnames$V1
        colnames(m) <- cnames$V1

        m
    }



