#' Code to prepare `herv_families` dataset goes here
#'
#' Loads and annotates family information from the HERV_rmsk.hg38.v2 annotation
#' which only includes 60 HERV families. Includes supergroup and genus
#' classification from Vargiu et al. 2016 (10.1186/s12977-015-0232-y).
#' Future versions may include more families.
#'
#'
herv_families <- readr::read_tsv("https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/families.tsv",
                                 col_names = c('family', 'intModel', 'ltr_model'),
                                 comment='#')



# Supergroup families from Varigu et al.
supergroups <- list('HML'=c("HML1", "HML2", "HML3", "HML4", "HML5", "HML6", "HERVK11", "HERVK11D", "HERVKC4", "HERVK14C"),
                    'HERVERI'=c("HERVE", "HARLEQUIN", "HERV3", "HERVI", "HERV1", "HERVEA"),
                    'HERVW9'=c("HERVW", "HERV9", "HERV30", "MER41", "LTR19"),
                    'HERVIPADP'=c('HERVIP10F', 'HERVIP10FH', 'HERVP71A'),
                    'HERVHF'=c('HERVH','HERVH48', 'HERVFC1', 'HERVFC2', 'HERVFH19', 'HERVFH21', 'LTR46'),
                    'HERVFRDLIKE'=c('HERVFRD', 'PABLA', 'PABLB', 'PRIMA41'),
                    'HEPSI'=c('PRIMA4','MER4','MER4B','MER34B','MER101'),
                    'HUERSP'=c('HUERSP1','HUERSP2','HUERSP3','HUERSP3B',"LTR25",'MER61'),
                    'HSERVIII'=c('HERVL', 'HERVL18', 'HERVL32', 'HERVL40', 'HERVL66', 'HERVL74', 'ERVL', 'ERVLB4', 'ERVLE', 'LTR57','ERV316A3'),
                    'MLLV'=c('HERVS71'),
                    'ERV1'=c('HERV4', 'PRIMAX', 'LTR23')
)

herv_families$supergroup <- ''
for(sg in names(supergroups)){
    herv_families[herv_families$family %in% supergroups[[sg]],]$supergroup <- sg
}

herv_families$group <- c('HERVK', 'HERVK', 'HERVK', 'HERVK', 'HERVK', 'HERVK',
                         'HERVK', 'HERVK', 'HERVK', 'HERVK', 'HERVW', 'HERVW',
                         'HERVW', 'HERVE', 'HERVE', 'HERVF', 'HERVF', 'HERVF',
                         'HERVF', 'HERVF', 'HERVH', 'HERVH', 'HERVI', 'HERVI',
                         'HERVI', 'HERVL', 'HERVL', 'HERVL', 'HERVL', 'HERVL',
                         'HERVL', 'ERVL', 'ERVL', 'ERVL', 'ERV1', 'ERV1',
                         'HERVP', 'HERVS', 'PAB', 'PAB', 'MER4', 'MER4',
                         'MER4', 'MER4', 'MER4', 'MER4', 'PRIMA', 'PRIMA',
                         'PRIMA', 'MISC', 'MISC', 'MISC', 'MISC', 'MISC',
                         'ERV3', 'HARL', 'HUERS', 'HUERS', 'HUERS', 'HUERS')

herv_genus <- list('ERV1'=c('MLLV',"ERV1",'HUERSP','HEPSI','HERVFRDLIKE','HERVHF','HERVIPADP','HERVW9','HERVERI'),
                   'ERV2'=c('HML'),
                   'ERV3'=c('HSERVIII')
                   )
herv_families$genus <- ''
for(hg in names(herv_genus)){
    herv_families[herv_families$supergroup %in% herv_genus[[hg]],]$genus <- hg
}

usethis::use_data(herv_families, overwrite = TRUE)
