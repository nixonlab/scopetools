#' Code to prepare `retro.hg38.v1` dataset goes here
#'
#'
retro.hg38.v1 <- readr::read_tsv("https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/retro.hg38.v1/genes.tsv.gz", na=c('.'))

retro.hg38.v1 <- retro.hg38.v1 %>%
    tidyr::separate(locus, c("family"), sep='_', remove=F, extra='drop') %>%
    dplyr::mutate(
      te_class = factor(ifelse(is.na(l1base_id), 'LTR', 'LINE'), levels=c('LTR','LINE')),
    )

usethis::use_data(retro.hg38.v1, overwrite = TRUE)
