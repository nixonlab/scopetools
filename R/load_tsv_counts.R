#' Load counts from multiple tab-delimited reports
#'
#' This returns a dataframe where the rows are genes/loci and columns are
#' samples. If `files` is a named vector, the names are used as sample
#' names; otherwise sample names can be provided by `colnames`. If not
#' provided, unique strings will be extracted from the report filenames.
#'
#' @param files A character vector with paths to the reports.
#' @param colnames Output column names. If `files` is a named vector, the names
#'    will be used as column names. Otherwise finds a unique string
#'    from the report filenames.
#' @param all_locs A character vector with all genes/loci to be included
#'    in the output; typically this is a list of all loci in the annotation.
#'    If not provided, the loci will be the union of loci in all reports, in
#'    an arbitrary order.
#' @param header Whether the first line of the input files is a header.
#' @param gene_id_column Report column number containing the gene ID.
#' @param count_column Report column number containing the count.
#' @param exclude_features Gene/feature IDs to be excluded from dataframe.
#'
#' @return A dataframe where the rows are genes/loci and the columns are
#'    samples (when each report corresponds to one sample).
#' @export
#'
#' @examples
#' ddir <- system.file("extdata", package="scopetools")
#' tsv_files <- Sys.glob(file.path(ddir, '*.ReadsPerGene.out.tab'))
#'
#' gene_counts <- load_tsv_counts(tsv_files)
#'
load_tsv_counts <-
    function(files,
             colnames = names(files),
             all_locs = NULL,
             header = FALSE,
             gene_id_column = 1,
             count_column = 2,
             exclude_features = c()
    ) {
        if (is.null(colnames)) {
            # colnames not provided and `files` was not a named vector
            colnames <- unique_names_from_file_list(files)
        }

        count_list <- lapply(1:length(files), function(i) {
            report <-
                read.table(
                    files[i],
                    sep = '\t',
                    header = header,
                    stringsAsFactors = F
                )
            ret <- report[, count_column]
            names(ret) <- report[, 1]
            ret
        })
        if (is.null(all_locs)) {
            # locus list was not provided
            all_locs <- unique(do.call(c, lapply(count_list, names)))
        }
        all_locs <- setdiff(all_locs, exclude_features)

        count.df <- lapply(1:length(count_list), function(i) {
            cts <- count_list[[i]]
            ret <- dplyr::left_join(
                data.frame(gene_id = all_locs, stringsAsFactors = F),
                data.frame(gene_id = names(cts), count = cts),
                by = 'gene_id'
            )
            ret[is.na(ret)] <- 0
            # stopifnot(all(ret$gene_id == all_locs))
            ret$gene_id <- NULL
            names(ret) <- c(colnames[i])
            ret
        }) %>% dplyr::bind_cols()

        row.names(count.df) <- all_locs
        count.df

    }


#' Load counts from multiple STAR tab-delimited reports
#'
#' This returns a dataframe where the rows are genes/loci and columns are
#' samples. If `files` is a named vector, the names are used as sample
#' names; otherwise sample names can be provided by `colnames`. If not
#' provided, unique strings will be extracted from the report filenames.
#'
#' @param files A character vector with paths to the reports.
#' @param colnames Output column names. If `files` is a named vector, the names
#'    will be used as column names. Otherwise finds a unique string
#'    from the report filenames.
#' @param stranded Strandedness of data. Use 'unstranded' for unstranded data,
#'    'read1' if the first read aligns to the RNA, and 'read2' if the second
#'    read aligns to the RNA. dUTP-based methods should use 'read2'.
#'
#' @return A dataframe where the rows are genes/loci and the columns are
#'    samples (when each report corresponds to one sample).
#' @export
#'
#' @examples
#' ddir <- system.file("extdata", package="scopetools")
#' tsv_files <- Sys.glob(file.path(ddir, '*.ReadsPerGene.out.tab'))
#'
#' gene_counts <- load_star_counts(tsv_files)
load_star_counts <-
    function(files,
             colnames = names(files),
             stranded = c('unstranded', 'read1', 'read2')
    ) {
        star_exclude_features <-
            c('N_unmapped',
              'N_multimapping',
              'N_noFeature',
              'N_ambiguous')

        stranded_cols <- c('unstranded' = 2,
                           'read1' = 3,
                           'read2' = 4)
        stranded <- match.arg(stranded)

        load_tsv_counts(
            files,
            colnames,
            header = FALSE,
            gene_id_column = 1,
            count_column = stranded_cols[stranded],
            exclude_features = star_exclude_features
        )
    }
