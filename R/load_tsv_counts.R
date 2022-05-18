#' Load counts from multiple tab-delimited reports
#'
#' This returns a dataframe where the rows are genes/loci and columns are
#' samples. If `files` is a named vector, the names are used as sample
#' names; otherwise sample names can be provided by `colnames`. If not
#' provided, unique strings will be extracted from the report filenames.
#'
load_tsv_counts <- function(files, colnames = names(files), all_locs = NULL, count_column = 'final_count') {
  if (is.null(colnames)) {
    # colnames not provided and `files` was not a named vector
    colnames <- unique_names_from_file_list(files)
  }
  lapply(1:length(files), function(i) {

  })
}
