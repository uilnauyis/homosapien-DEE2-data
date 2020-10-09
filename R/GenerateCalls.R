#' GenerateCalls
#'
#' @param normalizedSExpr SummarizedExperiment
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importMethodsFrom SummarizedExperiment assay assays
GenerateCalls <- function(normalizedSExpr) {
  calls <- SummarizedExperiment::assay(normalizedSExpr, 'counts')
  
  ## if the normalized count is greater than zero, set the conrisbonding call
  ## to 1, otherwise the call remains 0
  calls[calls > 0] <- 1
  
  SummarizedExperiment::assay(normalizedSExpr, "calls") <- calls
  
  normalizedSExpr
}