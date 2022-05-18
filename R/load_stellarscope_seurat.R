#' Load count matrix from Stellarscope output files
#'
#' This returns a `Seurat` object with raw counts from Stellarscope and
#' STARsolo (if provided). Only cell barcodes present in both are
#' retained. Feature metadata can also be included for TEs and genes.
#'
#' @param stellarscope_dir Stellarscope output directory
#' @param starsolo_dir STARsolo output directory
#' @param project Project name for the `Seurat` object
#' @param min.cells Include features detected in at least this many cells.
#'    Will subset the counts matrix as well. To reintroduce excluded features,
#'    create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are
#'    detected.
#' @param TE_count_file File name for Stellarscope count matrix in MTX format.
#'    Default is a file matching '*TE_counts.mtx' in `stellarscope_dir`.
#' @param TE_feature_file File name for Stellarscope feature list in TSV
#'    format. Default is a file matching '*features.tsv' in `stellarscope_dir`.
#' @param TE_barcode_file File name for Stellarscope barcode list in TSV
#'    format. Default is a file matching '*barcodes.tsv' in `stellarscope_dir`.
#' @param GENE_count_file File name for STARsolo count matrix in MTX format.
#'    Default is a file matching '*matrix.mtx' in `starsolo_dir`.
#' @param GENE_feature_file File name for STARsolo feature list in TSV format.
#'    Default is a file matching '*features.tsv' in `starsolo_dir`.
#' @param GENE_barcode_file File name for STARsolo barcode list in TSV format.
#'    Default is a file matching '*barcodes.tsv' in `starsolo_dir`.
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
           project = 'StellarscopeProject',
           min.cells = 0,
           min.features= 0,
           TE_count_file = Sys.glob(file.path(stellarscope_dir, '*TE_counts.mtx')),
           TE_feature_file = Sys.glob(file.path(stellarscope_dir, '*features.tsv')),
           TE_barcode_file = Sys.glob(file.path(stellarscope_dir, '*barcodes.tsv')),
           GENE_count_file = Sys.glob(file.path(starsolo_dir, '*matrix.mtx')),
           GENE_feature_file = Sys.glob(file.path(starsolo_dir, '*features.tsv')),
           GENE_barcode_file = Sys.glob(file.path(starsolo_dir, '*barcodes.tsv')),
           TE_metadata = data.frame(retro.hg38.v1, row.names='locus'),
           GE_metadata = NULL
  ) {

    counts.TE <- Matrix::readMM(file = TE_count_file)
    counts.TE <- Matrix::t(counts.TE)
    features.TE <- read.delim(file = TE_feature_file,
                              header = FALSE,
                              stringsAsFactors = FALSE)
    barcodes.TE <- read.delim(file = TE_barcode_file,
                              header = FALSE,
                              stringsAsFactors = FALSE)

    if('dummy' %in% features.TE$V1) {
      rem <- which(features.TE$V1 == 'dummy')
      features.TE <- data.frame(V1=features.TE[-rem, 1])
      counts.TE <- counts.TE[-rem,]
    }
    stopifnot(nrow(features.TE) == nrow(counts.TE))
    stopifnot(nrow(barcodes.TE) == ncol(counts.TE))

    rownames(counts.TE) <- features.TE$V1
    colnames(counts.TE) <- barcodes.TE$V1

    features.TE <- data.frame(id=features.TE$V1, feattype='TE', symbol=features.TE$V1, te_class=TE_metadata[features.TE$V1, ]$te_class, te_family=TE_metadata[features.TE$V1, ]$family)

    dogenes <- TRUE
    if(length(GENE_count_file) == 0) {warning('Missing gene count matrix\n'); dogenes = dogenes & FALSE;}
    if(length(GENE_feature_file) == 0) {warning('Missing gene feature tsv\n'); dogenes = dogenes & FALSE;}
    if(length(GENE_barcode_file) == 0) {warning('Missing gene barcode tsv\n'); dogenes = dogenes & FALSE;}

    if(!dogenes) {
      warning('Gene counts are not loaded\n')
      ret <- Seurat::CreateSeuratObject(counts.TE, project=project, min.cells=min.cells, min.features=min.features)
      ret[['RNA']] <- Seurat::AddMetaData(ret[['RNA']], features.TE)
      return(ret)
    }

    counts.GE <- Matrix::readMM(file = GENE_count_file)
    features.GE <- read.delim(file = GENE_feature_file,
                              header = FALSE,
                              stringsAsFactors = FALSE)
    barcodes.GE <- read.delim(file = GENE_barcode_file,
                              header = FALSE,
                              stringsAsFactors = FALSE)

    rownames(counts.GE) <- features.GE$V1
    colnames(counts.GE) <- barcodes.GE$V1

    features.GE <- data.frame(id=features.GE$V1, feattype='GE', symbol=features.GE$V2, te_class=NA, te_family=NA)

    all_bc <- intersect(colnames(counts.TE), colnames(counts.GE))
    all_feat <- rbind(features.GE, features.TE)
    all_feat$orig_id <- all_feat$id
    all_feat$id <- gsub('_', '-', all_feat$id)
    rownames(all_feat) <- all_feat$id

    counts <- rbind(counts.GE[,all_bc], counts.TE[,all_bc])

    stopifnot(dim(counts)[1] == nrow(all_feat))
    stopifnot(dim(counts)[2] == length(all_bc))

    ret <- Seurat::CreateSeuratObject(counts, project=project, min.cells=min.cells, min.features=min.features)
    ret[['RNA']] <- Seurat::AddMetaData(ret[['RNA']], all_feat)

    stopifnot(all(rownames(ret) == all_feat$id))
    stopifnot(all(colnames(ret) == all_bc))
    ret
  }
