#' Filter cells from Seurat object
#'
#' This function expects that the Seurat object includes TE feature metadata.
#' Specifically, `feattype` should be set to "TE" for TE features and `te_class`
#' should indicate whether the TE is an "LTR" or "LINE" retroelement. Seurat
#' objects created with [scopetools::load_stellarscope_seurat()] will include
#' the necessary metadata.
#'
#' The function adds metadata for the mitochondrial percentage (`percent.mt`),
#' HERV percentage (`percent.HERV`), L1 percentage (`percent.L1`), and total TE
#' percentage (`percent.TE`). Outliers are detected based on the median absolute
#' deviation (MAD) implemented in [scuttle::isOutlier()]. Cells are excluded
#' based on extreme RNA count, extreme feature count, or high MT percentage.
#'
#' @param object Seurat object (with TE feature metadata)
#'
#' @return Seurat object with cells passing QC
#' @export
#'
#' @examples
#' # Seurat object
#' pbmc <- stellarscope_cell_qc(pbmc.orig)
#'
stellarscope_cell_qc <- function(object) {
  # Feature metadata
  fmeta <- object[['RNA']]@meta.features

  # Select mitochondrial features and calculate mt percentage
  mt_feats <- grepl('^MT-', fmeta$symbol)
  object[['percent.mt']] <-
    Seurat::PercentageFeatureSet(object, features = fmeta[mt_feats, 'id'])

  # Select HERV features and calculate HERV percentage
  herv_feats <- !is.na(fmeta$te_class) & fmeta$te_class == 'LTR'
  object[['percent.HERV']] <-
    Seurat::PercentageFeatureSet(object, features = fmeta[herv_feats, 'id'])

  # Select L1 features and calculate L1 percentage
  l1_feats <- !is.na(fmeta$te_class) & fmeta$te_class == 'LINE'
  object[['percent.L1']] <-
    Seurat::PercentageFeatureSet(object, features = fmeta[l1_feats, 'id'])

  # Select TE features and calculate TE percentage
  te_feats <- fmeta$feattype == 'TE'
  object[['percent.TE']] <-
    Seurat::PercentageFeatureSet(object, features = fmeta[te_feats, 'id'])

  # Determine which values in a numeric vector are outliers based on the
  # median absolute deviation (MAD).
  qc.ncount_rna <-
    scater::isOutlier(object$nCount_RNA, log = TRUE, type = "both")
  qc.nfeature_rna <-
    scater::isOutlier(object$nFeature_RNA, log = TRUE, type = "both")
  qc.percent_mt <-
    scater::isOutlier(object$percent.mt,  type = "higher")

  # (Use the outlier vectors directly)
  # The cell is excluded if isOutlier=TRUE for any of the three criteria:
  # extreme RNA count, extreme feature count, high MT percentage.
  # NOTE: nCount_RNA is included in the logical expression since Seurat will
  #       an error if no variables attached to the object are included. Should
  #       evaluate to TRUE unless nCount_RNA is zero, which should be excluded
  #       anyways.
  object <- subset(object,
                   subset = nCount_RNA &
                     !(qc.ncount_rna | qc.nfeature_rna | qc.percent_mt))
  object
}
