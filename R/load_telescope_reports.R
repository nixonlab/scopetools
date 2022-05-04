#' Load counts from multiple Telescope reports
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
#' @param count_column The column in the Telescope report to use.
#'    Defaults is 'final_count'
#'
#' @return A dataframe where the rows are genes/loci and the columns are
#'    samples (when each report corresponds to one sample).
#' @export
#'
#' @examples
#' ddir <- system.file("extdata", package="scopetools")
#' t_files <- file.path(ddir, paste0(samples.encode$GEO_sample, '-telescope_report.tsv'))
#' names(t_files) <- samples.encode$Name
#'
#' # Provide all_locs to include all annotated loci (preferred)
#' herv_counts <- load_telescope_reports(t_files, all_locs=HERVrmsk.hg38.v2$locus)
#'
#' # Contstruct table with loci that appear in reports
#' \dontrun{
#' herv_counts <- load_telescope_reports(t_files)
#' }
#'
load_telescope_reports <-
  function(files,
           colnames = names(files),
           all_locs = NULL,
           count_column = 'final_count') {
    if (is.null(colnames)) {
      # colnames not provided and `files` was not a named vector
      colnames <- unique_names_from_file_list(files)
    }

    count_list <- lapply(1:length(files), function(i) {
      report <-
        read.table(
          files[i],
          sep = '\t',
          header = T,
          stringsAsFactors = F
        )
      ret <- report[, count_column]
      names(ret) <- report$transcript
      ret
    })

    if (is.null(all_locs)) {
      # locus list was not provided
      all_locs <- unique(do.call(c, sapply(count_list, names)))
    }

    count.df <- lapply(1:length(count_list), function(i) {
      cts <- count_list[[i]]
      ret <- dplyr::left_join(
        data.frame(transcript = all_locs, stringsAsFactors = F),
        data.frame(transcript = names(cts), count = cts),
        by = 'transcript'
      )
      ret[is.na(ret)] <- 0
      # stopifnot(all(ret$transcript == all_locs))
      ret$transcript <- NULL
      names(ret) <- c(colnames[i])
      ret
    }) %>% dplyr::bind_cols()

    row.names(count.df) <- all_locs
    count.df
  }
