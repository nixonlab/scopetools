#' Code to prepare `HERVrmsk.hg38.v2` dataset goes here
#'
#'
HERVrmsk.hg38.v2 <- readr::read_tsv("https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/genes.tsv.gz", na=c('.'))

HERVrmsk.hg38.v2 <- HERVrmsk.hg38.v2 %>%
  tidyr::separate(locus, c("family"), sep='_', remove=F, extra='drop') %>%
  dplyr::mutate(
    te_class = factor('LTR', levels=c('LTR','LINE')),
  )

usethis::use_data(HERVrmsk.hg38.v2, overwrite = TRUE)
