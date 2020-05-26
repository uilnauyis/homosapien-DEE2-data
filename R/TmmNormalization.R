#' TMM normalization for the selected DEE2 data
#'
#' This function allows you to express your love of cats.
#' @param sExpr SummarizedExperiment object that is created via DEE2 R interface 
#'    or restored from downloaded files.
#' @param counts.cutoff threshold of filtering the 
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment 
#' @importClassesFrom edgeR DGEList
#' @importFrom SummarizedExperiment assay
#' @importFrom edgeR DGEList filterByExpr calcNormFactors cpm

TmmNormalization <- function(sExpr, counts.cutoff = 5, excludeFail = TRUE) {
  
  sExpr <- .filterData(sExpr, counts.cutoff)

  # https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  dgeList <- DGEList(counts = assay(sExpr))
  keep <- filterByExpr(dgeList)
  dgeList <- dgeList[keep, , keep.lib.sizes=FALSE]
  dgeList <- calcNormFactors(dgeList)
  norm.counts <- cpm(dgeList)
  log.norm.counts <- log2(norm.counts + 1)

  assay(sExpr, withDimnames=FALSE) <- log.norm.counts

  sExpr
}