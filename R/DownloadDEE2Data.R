#' Download DEE2 Data.
#' 
#' By default, 'bulkDownloading' is set to false, as the current version of DEE2 
#' could not get consistent data with bulk downloading. 
#'
#' @export  
#' @param species
#' @param sRRAccessions
#' @param bulkDownloading
#' @importFrom getDEE2 getDEE2
#' @importFrom BiocGenerics cbind
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importMethodsFrom BiocGenerics cbind
#' @references https://github.com/markziemann/dee2
#' @examples 
#' ## Load the dependencies
#' library(DEE2HsapienData)
#' 
#' ## Create a list of SRR Accessions of which we want to retrieve the counts 
#' from DEE2 library
#' SRRvec = c("SRR1783836", "SRR1783837", "SRR1783838", "SRR1999221", "SRR2153338", 
#'            "SRR2153409", "SRR2153289")
#' 
#' ## Creat a 'temp' directory in current working directory if it is not created
#' ## already. In this tutorial, we use it as the workspace 
#' dee2Data <- DownloadDEE2Data('hsapiens', SRRvec)
DownloadDEE2Data <- function(species, sRRAccessions, 
                              bulkDownloading = FALSE) {
  if (bulkDownloading) {
    dee2Data <- getDEE2(species, sRRAccessions)
  } else {
    dee2Data <- do.call(cbind, lapply(sRRAccessions, getDEE2, species = species))
  }
  dee2Data
}