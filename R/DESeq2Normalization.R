#' DESeq2 normalization for the selected DEE2 data
#'
#' @param sExpr SummarizedExperiment object that is created via DEE2 R interface 
#'    or restored from downloaded files.
#' @param counts.cutoff threshold of filtering the 
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom BiocGenerics estimateSizeFactors counts

DESeq2Normalization <- function(sExpr, counts.cutoff = 5, excludeFail = TRUE) {
  sExpr <- .filterData(sExpr, counts.cutoff)

  # Create DESeq.ds from the summarizedExperiment object
  DESeq.ds <- DESeqDataSetFromMatrix(countData = (assay(sExpr, "counts") + 1),
    colData = colData(sExpr),
    rowData = rowData(sExpr),
    design = ~1)

  # DESeq2 Default normalization method
  DESeq.dsDefault <- estimateSizeFactors(DESeq.ds)
  norm.counts <- counts(DESeq.dsDefault, normalized = TRUE)
  log.norm.counts <- log2(norm.counts + 1)

  assay(sExpr, 'count') <- log.norm.counts

  sExpr
}