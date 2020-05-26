################################################################################
## FUNCTION: .zscore
################################################################################
#' Calculates z-score. We use z-score as a measure of how close the standards
#' cluster to each other
#' @param sExpr raw sExpr data object 
#' @param counts.cutoff the threshold for filtering genes. For a specific gene,
#'    the total cound 
#' @return sExpr 
#' @importFrom stringr regex str_detect
#' @importFrom SummarizedExperiment colData

.filterData <- function(sExpr, counts.cutoff = 5, excludeFail = TRUE) {
  # remove samples that are marked as 'FAIL'
  if (excludeFail) {
    pData <- colData(sExpr)
    sExpr <- sExpr[, !str_detect(pData$QC_summary, regex('FAIL.*'))]
  }

  # remove genes smaller than the cutoff
  sExpr <- sExpr[rowSums(assay(sExpr)) > counts.cutoff, ]

  sExpr
}