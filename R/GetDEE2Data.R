#' Get a 'SummarizedExperiment' object 'sExpr' from DEE2 library and set its 
#' gene info based on the data from biomaRt.
#'
#' @export  
#' @param sampleInfoPath A path to a spreadsheet file that stores 'category', 
#'      'general_cell_type', 'parental_cell_type' and 'SRR_accession' for each 
#'      sample.
#' @param species the species of the DEE2 data. Available species includes 
#'      'athaliana', 'celegans','dmelanogaster','drerio','ecoli', 'hsapiens',
#'      'mmusculus','rnorvegicus' and 'scerevisiae'.
#' @param counts either 'GeneCounts' or 'Tx2Gene'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData 
#' @importFrom S4Vectors DataFrame
#' @importFrom rio import
#' @importFrom getDEE2 getDEE2
#' @references https://github.com/markziemann/dee2/blob/master/AccessDEEfromR.md
#' @examples 
#' ## Load the dependencies
#' library(SummarizedExperiment)
#' 
#' ## Creat a 'temp' directory in current working directory if it is not created
#' ## already. In this tutorial, we use it as the workspace 
#' tempDirPath <- paste(getwd(), 'temp', sep = '/')
#' dir.create(tempDirPath)
#' 
#' ## get the 'list' object. 
#' dee2Data <- GetDEE2Data(
#'  system.file("extdata", "SelectedDEE2_v2.ods", package = "DEE2HsapienData"),
#'  paste(tempDirPath, 'DEE2Data', sep = '/') , 
#'  'hsapiens')
GetDEE2Data <- function(
    sampleInfoPath, 
    species, 
    counts = 'GeneCounts', 
    bulkDownloading = FALSE) {  
    ## 'counts' should be either 'Tx2Gene' or 'GeneCounts'
    stopifnot(counts == 'GeneCounts' ||
              counts == 'Tx2Gene')

    stopifnot(is.element(species, 
                         c('athaliana', 'celegans', 'dmelanogaster', 'drerio', 
                           'ecoli', 'hsapiens', 'mmusculus', 'rnorvegicus', 
                           'scerevisiae')))              

    ## get sample information from the file specified by 'sampleInfoPath'
    samplesInfoDFrame <- import(sampleInfoPath)
    
    ## extract targeted SRR accessions from the spreadsheet
    SRR_accessions <- unique((samplesInfoDFrame[['SRR_accession']]))

    dee2Data <- DownloadDEE2Data(species, SRR_accessions,
                                  bulkDownloading = bulkDownloading)

    dee2Data <- .excludeFailedSamples(dee2Data)

    dee2Data <- .constructSExprFromDEE2Data(dee2Data, samplesInfoDFrame)

    dee2Data
}

################################################################################
## FUNCTION: .excludeFailedSamples
################################################################################
#' .excludeFailedSamples function
#'
#' @param dee2Data a vector of dataframes that stores DEE2 data.
#' @return dee2Data 
#' @importFrom stringr str_detect regex
.excludeFailedSamples <- function(dee2Data) {
  ## find the list of samples that are marked as 'FAIL' and add this list to 
  ## selected dee2 Data
  colDat <- colData(dee2Data)

  fail <- str_detect(
      colDat$QC_SUMMARY,
      regex('FAIL.*'))

  ## exclude samples marked as 'FAIL' from dee2 data 
  sel <- !fail
  dee2Data <- dee2Data[, sel]

  dee2Data
}

################################################################################
## FUNCTION: .constructSExprFromDEE2Data
################################################################################
#' .constructSExprFromDEE2Data function
#' 
#' @param dee2Data a vector of dataframes that stores DEE2 data.
#' @param samplesInfoDFrame
#' @return dee2Data 
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#' @importFrom S4Vectors DataFrame
.constructSExprFromDEE2Data <- function(dee2Data, samplesInfoDFrame) {
    ## Exclude duplicated samples
    samplesInfoDFrame <- samplesInfoDFrame[
      !duplicated(samplesInfoDFrame[, "SRR_accession", drop = FALSE]), 
      ,
      drop = FALSE]
    
    ## Extract category', 'cell_type', 'general_cell_type', 'parental_cell_type', 
    ## and 'SRR_accession'. These columns are necessary for CellScore analysis.
    samplesInfoDFrame <- 
      samplesInfoDFrame[
        ,
        c('category', 'cell_type', 'general_cell_type', 
          'parental_cell_type', "SRR_accession"), drop = FALSE
      ]
    
    ## Get all SRR accessions and add them to the coldata 
    colDat <- colData(dee2Data)
    colDat$SRR_accession <- rownames(colDat)
    
    ## merge colDat and sampleMetadataDframe to get the complete
    ## sample metedata
    ## todo: verify the order (check left-join)
    mergedColData <- merge(colDat, 
                           samplesInfoDFrame, 
                           by = 'SRR_accession')
    
    rownames(mergedColData) <- mergedColData$SRR_accession
    SummarizedExperiment::colData(dee2Data) <- DataFrame(mergedColData)

    dee2Data
}