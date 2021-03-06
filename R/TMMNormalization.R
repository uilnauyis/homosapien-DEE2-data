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

TMMNormalization <- function(sExpr) {
  dgeList <- DGEList(counts = assay(sExpr, 'counts'))
  keep <- filterByExpr(dgeList)
  dgeList <- dgeList[keep, , keep.lib.sizes=FALSE]
  dgeList <- calcNormFactors(dgeList)
  log.norm.counts <- cpm(dgeList, log=TRUE)

  geneSel <- rownames(log.norm.counts)
  normalizedSExpr <- SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=log.norm.counts),
    rowData = rowData(sExpr)[geneSel, ], 
    colData = colData(sExpr))
  
  ## Return the normalized SummarizedExperiment object
  normalizedSExpr
}