#' Find unique names from list of files
#' @param files A character vector with filenames
#' @NoRd
unique_names_from_file_list <- function(files) {
    ret <- files

    fsplit <- strsplit(ret, '/')
    for(i in length(fsplit[[1]]):1) {
      tmp <- sapply(fsplit, '[[', i)
      if(!any(duplicated(tmp))) {
        ret <- tmp
        break
      }
    }
    # Remove suffix
    ret <- gsub(paste0(Biobase::lcSuffix(ret), '$'), '', ret)
    ret
}
