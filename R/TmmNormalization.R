#' TMM normalization for the selected DEE2 data
#'
#' This function allows you to express your love of cats.
#' @param sExpr SummarizedExperiment object that is created via DEE2 R interface 
#'    or restored from downloaded files.
#' @param counts.cutoff threshold of filtering the 
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment 
#' @importClassesFrom edgeR DGEList
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom edgeR DGEList filterByExpr calcNormFactors cpm
#' @references https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

TmmNormalization <- function(dee2Data) {
  sExpr <- dee2Data$sExpr

  dgeList <- edgeR::DGEList(counts = SummarizedExperiment::assay(sExpr))
  keep <- edgeR::filterByExpr(dgeList)
  dgeList <- dgeList[keep, , keep.lib.sizes=FALSE]
  dgeList <- edgeR::calcNormFactors(dgeList)
  log.norm.counts <- edgeR::cpm(dgeList, log=TRUE)

  geneSel <- rownames(log.norm.counts)
  normalizedSExpr <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=log.norm.counts),
    rowData=SummarizedExperiment::rowData(sExpr)[geneSel, ], 
    colData=SummarizedExperiment::colData(sExpr))
  dee2Data[['normalizedSExpr']] <- normalizedSExpr

  dee2Data
}