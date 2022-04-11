#' Code to prepare `samples.encode` dataset goes here
#'
#'
samples.encode <- readr::read_tsv(file.path(system.file("extdata", package="scopetools"), "encode.tsv"))

usethis::use_data(samples.encode, overwrite = TRUE)
