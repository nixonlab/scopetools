#' Load count matrix from Stellarscope output files
#'
#' This returns a `Seurat` object with raw counts from Stellarscope and
#' STARsolo (if provided). Only cell barcodes present in both are
#' retained. Feature metadata can also be included for TEs and genes.
#'
#' @param stellarscope_dir Stellarscope output directory
#' @param starsolo_dir STARsolo output directory, optional.
#' @param exp_tag Experiment tag used in stellarscope run.
#'    Default is "stellarscope".
#' @param project Project name for the `Seurat` object
#' @param min.cells Include features detected in at least this many cells.
#'    Will subset the counts matrix as well. To reintroduce excluded features,
#'    create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are
#'    detected.
#' @param remove.nofeat Remove "no feature" category from Stellarscope counts.
#'    Default is TRUE.
#' @param TE_count_filename File name pattern for TE count matrix (MTX format).
#'    May contain wildcards.Default is '*TE_counts.mtx'.
#' @param TE_feature_filename File name pattern for TE feature list
#'    (tab-separated, features in first column). May contain wildcards.
#'    Default is '*features.tsv'.
#' @param TE_barcode_filename File name pattern for TE barcode list
#'    (tab-separated, barcodes in first column). May contain wildcards.
#'    Default is '*barcodes.tsv'.
#' @param GENE_count_filename File name pattern for gene count matrix
#'    (MTX format). May contain wildcards. Default is '*matrix.mtx'.
#' @param GENE_feature_filename File name pattern for gene feature list
#'    (tab-separated, features in first column). May contain wildcards.
#'    Default is '*features.tsv'.
#' @param GENE_barcode_filename File name pattern for gene barcode list
#'    (tab-separated, barcodes in first column). May contain wildcards.
#'    Default is '*barcodes.tsv'.
#' @param TE_metadata Data frame with metadata related to the TE features
#' @param GE_metadata Data frame with metadata related to the gene features
#'
#' @return A `Seurat` object with counts from stellarscope and STARsolo.
#' @export
#'
#' @examples
load_stellarscope_seurat <-
  function(stellarscope_dir,
           starsolo_dir = NULL,
           exp_tag = 'stellarscope',
           project = 'StellarscopeProject',
           min.cells = 0,
           min.features= 0,
           remove.nofeat = TRUE,
           TE_count_filename = '*TE_counts.mtx',
           TE_feature_filename = '*features.tsv',
           TE_barcode_filename = '*barcodes.tsv',
           GENE_count_filename = '*matrix.mtx',
           GENE_feature_filename = '*features.tsv',
           GENE_barcode_filename = '*barcodes.tsv',
           TE_metadata = data.frame(retro.hg38.v1, row.names='locus'),
           GE_metadata = NULL
  ) {
    # Load TE count matrix
    counts.TE <- load_stellarscope_report(
        stellarscope_dir = stellarscope_dir,
        exp_tag = exp_tag,
        count_filename = TE_count_filename,
        feature_filename = TE_feature_filename,
        barcode_filename = TE_barcode_filename
    )

    # Remove no feature
    if(remove.nofeat) {
        rem <- which(rownames(counts.TE) == '__no_feature')
        if(length(rem)>0)
            counts.TE <- counts.TE[-rem, ]
    }

    # Create TE feature metadata
    meta.TE <- dplyr::left_join(
        data.frame(gene_id=rownames(counts.TE)),
        TE_metadata,
        by='gene_id'
    ) %>%
        dplyr::mutate(feattype='TE',
                      symbol=gene_id
        ) %>%
        dplyr::select(id=gene_id, feattype, symbol, te_class, te_family=family)

    stopifnot(all(meta.TE$id == rownames(counts.TE)))

    # Other (starsolo) matrix
    dogenes <- !is.null(starsolo_dir)

    if(!dogenes) {
      warning('Gene counts are not loaded\n')
      ret <- Seurat::CreateSeuratObject(
          counts.TE,
          project=project,
          min.cells=min.cells,
          min.features=min.features
      )
      ret[['RNA']] <- Seurat::AddMetaData(ret[['RNA']], meta.TE)
      return(ret)
    }

    counts.GE <- load_stellarscope_report(
        stellarscope_dir = starsolo_dir,
        exp_tag = NULL,
        count_filename = GENE_count_filename,
        feature_filename = GENE_feature_filename,
        barcode_filename = GENE_barcode_filename
    )

    if(is.null(GE_metadata)) {
        # if no df was provided, get symbols from starsolo output
        features_tsv <- Sys.glob(file.path(starsolo_dir, GENE_feature_filename))
        GE_metadata <- read.delim(file = features_tsv,
                                  header = FALSE,
                                  stringsAsFactors = FALSE
                                  )
        names(GE_metadata) <- c('id', 'symbol', 'assay')
    }

    meta.GE <- dplyr::left_join(
        data.frame(id=rownames(counts.GE)),
        GE_metadata,
        by='id'
    ) %>%
        dplyr::mutate(feattype='GE',
               te_class=NA,
               te_family=NA
        ) %>%
        dplyr::select(id, feattype, symbol, te_class, te_family)

    stopifnot(all(meta.GE$id == rownames(counts.GE)))

    # harmonize barcodes - found in both
    all_bc <- intersect(colnames(counts.TE), colnames(counts.GE))
    # harmonize barcodes - union
    # all_bc <- union(colnames(counts.TE), colnames(counts.GE))

    all_feat <- rbind(meta.GE, meta.TE)
    # replace NAs with blank, easier downstream
    all_feat[is.na(all_feat)] <- ''

    all_feat$orig_id <- all_feat$id
    all_feat$id <- gsub('_', '-', all_feat$id)
    rownames(all_feat) <- all_feat$id

    # TODO: this should be fixed for "union"
    counts <- rbind(counts.GE[,all_bc], counts.TE[,all_bc])

    stopifnot(dim(counts)[1] == nrow(all_feat))
    stopifnot(dim(counts)[2] == length(all_bc))

    ret <- Seurat::CreateSeuratObject(counts, project=project, min.cells=min.cells, min.features=min.features)
    ret[['RNA']] <- Seurat::AddMetaData(ret[['RNA']], all_feat)

    stopifnot(all(rownames(ret) == all_feat$id))
    stopifnot(all(colnames(ret) == all_bc))
    ret
  }
