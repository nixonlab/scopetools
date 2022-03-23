#' Load counts from multiple Telescope reports
#'
#' This returns a dataframe where the rows are genes/loci with one column per
#' report. (Typically each report corresponds to one sample). The  and the columns are
#'
#'
#' @param files A character vector with paths to the reports. If a named vector, the
#'    names will be used as the column names.
#' @param count_column The column in the Telescope report to use.
#'    Defaults is 'final_count'
#' @param all_locs A character vector with the loci to include, possibly
#'    extracted from the annotation. If not provided, all loci appearing
#'    in the reports will be used in an arbitrary order.
#'
#' @return A dataframe where the rows are genes/loci and the columns are
#'    reports. Typically each report will correspond to one sample.
#' @export
#'
#' @examples
#' \dontrun{
#'    t_files <- file.path('samples', samples$sample_id, 'telescope.report.tsv')
#'    names(t_files) <- samples$sample_id
#'    herv_counts <- load_telescope_reports(t_files)
#'
#'    # If you have the annotation locs in a table, this will maintain
#'    # the ordering
#'    herv_counts <- load_telescope_reports(t_files, retro.annot$locus)
#'
#' }
load_telescope_reports <- function(files, count_column='final_count', all_locs=NULL) {

    if(is.null(names(files))) {
         # `files` was not a named vector
         # Create sample names by removing longest prefix and suffix from file names
        colnames <- gsub(paste0('^', Biobase::lcPrefix(files)), '', files)
        colnames <- gsub(paste0(Biobase::lcSuffix(colnames), '$'), '', colnames)
    } else {
      colnames <- names(files)
    }
    count_list <- lapply(1:length(files), function(i){
        report <- read.table(files[i], sep='\t', header=T, stringsAsFactors=F)
        ret <- report[,count_column]
        names(ret) <- report$transcript
        ret
    })

    if(is.null(all_locs)) {
        # locus list was not provided
        all_locs <- unique(do.call(c, sapply(count_list, names)))
    }

  count.df <- lapply(1:length(count_list), function(i){
    cts <- count_list[[i]]
    ret <- dplyr::left_join(
      data.frame(transcript=all_locs, stringsAsFactors = F),
      data.frame(transcript=names(cts), count=cts),
      by='transcript'
    )
    ret[is.na(ret)] <- 0
    # stopifnot(all(ret$transcript == all_locs))
    ret$transcript <- NULL
    names(ret) <- c(colnames[i])
    ret
  }) %>% bind_cols

  row.names(count.df) <- all_locs
  count.df
}
