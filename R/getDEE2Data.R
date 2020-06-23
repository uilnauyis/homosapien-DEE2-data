#' Construct a 'list' object that includes seven 'data.frame' objects 
#' ( "GeneCounts", "TxCounts", "GeneInfo", "TxInfo", "QcMx", "MetadataSummary" 
#' and "MetadataFull"), and a 'SummarizedExperiment' object 'sExpr'. In the case 
#' when there are SRR accessions not present in DEE2 database or not matching
#' the 'species' parameter, an list object 
#' named 'absent' will be in the returned list containing those accessions. The 
#' 'sExpr' object in the returned list contains an assay containing the same 
#' data as in "GeneCounts" object. The column data of 'sExpr' is a 'data.frame' 
#' that combines the data in the spreadsheet indicated by 'sampleInfoPath' and 
#' in 'MetadataFull' object. The row data of 'sExpr' includes the same data as 
#' in 'GeneInfo' object.
#' 
#' for details of  "GeneCounts", "TxCounts", "GeneInfo", "TxInfo", "QcMx", 
#' "MetadataSummary" and "MetadataFull" dataframs, refer to 
#' https://github.com/markziemann/dee2/blob/master/AccessDEEfromR.md
#'
#' @export  
#' @param sampleInfoPath A path to a spreadsheet file that stores 'category', 
#'      'general_cell_type', 'parental_cell_type' and 'SRR_accession' for each 
#'      sample.
#' @param outFile By default, this parameter is NULL. If it is defined, then 
#'      files that contains DEE2 data will be downloaded to the current working 
#'      directory in a 'zip' file named as the value of this parameter.
#' @param species the species of the DEE2 data. Available species includes 
#'      'athaliana', 'celegans','dmelanogaster','drerio','ecoli', 'hsapiens',
#'      'mmusculus','rnorvegicus' and 'scerevisiae'.
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
#' dee2Data <- getDEE2Data(
#'  system.file("extdata", "SelectedDEE2_v2.ods", package = "DEE2HsapienData"),
#'  paste(tempDirPath, 'DEE2Data', sep = '/') , 
#'  'hsapiens')
getDEE2Data <- function(
    sampleInfoPath, 
    outFile,
    species = c('athaliana', 'celegans', 'dmelanogaster', 'drerio', 'ecoli',
        'hsapiens', 'mmusculus', 'rnorvegicus', 'scerevisiae')) {   
    ## get sample information from the file specified by 'sampleInfoPath'
    samplesInfoDFrame <- rio::import(sampleInfoPath)

    ## Download and assgin DEE2 data. If 'outFile' is specified, tsv files will
    ## be saved at path specified by 'outFile' parameter 
    dee2Data <- .downloadDEE2Data(species, samplesInfoDFrame, outFile)

    ## order all dataframes in 'dee2Data' by the same order of samples in 
    ## 'GeneCount' dataframe
    dee2Data <- .orderSamples(dee2Data) 

    dee2Data <- .excludeFailedSamples(dee2Data)

    dee2Data <- .constructSExprFromDEE2Data(dee2Data, samplesInfoDFrame)

    dee2Data
}

################################################################################
## FUNCTION: .downloadDEE2Data
################################################################################
#' .downloadDEE2Data function
#'
#' @param species 
#' @param samplesInfoDFrame
#' @param outFile
#' @importFrom getDEE2 getDEE2
#' @return sExpr 
.downloadDEE2Data <- function(species, samplesInfoDFrame, outFile) {
    SRR_accessions <- unique((samplesInfoDFrame[['SRR_accession']]))
    dee2Data <- getDEE2::getDEE2(species, SRR_accessions, outfile = outFile)
    dee2Data
}

################################################################################
## FUNCTION: .orderSamples
################################################################################
#' .orderSamples function
#'
#' @param dee2Data a vector of dataframes that stores DEE2 data.
#' @return dee2Data 
.orderSamples <- function(dee2Data) {
    orderSamples <- colnames(dee2Data$GeneCounts)
    dee2Data$GeneCounts <- dee2Data$GeneCounts[, orderSamples, drop = FALSE]
    dee2Data$TxCounts <- dee2Data$TxCounts[, orderSamples, drop = FALSE]
    dee2Data$QcMx <- dee2Data$QcMx[, orderSamples, drop = FALSE]
    dee2Data$MetadataSummary <- dee2Data$MetadataSummary[orderSamples, , drop = FALSE]
    dee2Data$MetadataFull <- dee2Data$MetadataFull[orderSamples, , drop = FALSE]

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
  fail <- stringr::str_detect(
      dee2Data$QcMx['QC_SUMMARY', , drop = FALSE], 
      stringr::regex('FAIL.*'))
  dee2Data[['fail']] <- colnames(dee2Data$QcMx[, fail, drop = FALSE])

  ## exclude samples marked as 'FAIL' from dee2 data 
  sel <- !fail
  dee2Data$GeneCounts <- dee2Data$GeneCounts[, sel, drop = FALSE]
  dee2Data$TxCounts <- dee2Data$TxCounts[, sel, drop = FALSE]
  dee2Data$QcMx <- dee2Data$QcMx[, sel, drop = FALSE]
  dee2Data$MetadataSummary <- dee2Data$MetadataSummary[sel, , drop = FALSE]
  dee2Data$MetadataFull <- dee2Data$MetadataFull[sel, , drop = FALSE]

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

    geneCountsDFrame <- dee2Data[["GeneCounts"]]
    geneInfoDFrame <- dee2Data[["GeneInfo"]]

    counts <- data.matrix(geneCountsDFrame)

    sampleMetadataDframe <- .constructSampleMetadata(samplesInfoDFrame, 
        dee2Data$MetadataFull)

    ## construct sExpr object for the selected DEE2 data. First, we need to get 
    ## the order of both column names and row names of the count matrix, so that
    ## we could use the order to correctly align the counts with sample metedata
    ## (coldata) and gene metadata 
    orderSamples <- colnames(counts)
    orderGenes <- rownames(counts)

    sExpr <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts))
    SummarizedExperiment::colData(sExpr) <- 
        S4Vectors::DataFrame(sampleMetadataDframe[orderSamples, , drop = FALSE])
    SummarizedExperiment::rowData(sExpr) <- 
        S4Vectors::DataFrame(geneInfoDFrame[orderGenes, , drop = FALSE])

    dee2Data[['sExpr']] <- sExpr

    dee2Data
}

################################################################################
## FUNCTION: .constructSampleMetadata
################################################################################
#' .constructSampleMetadata function
#' 
#' @param samplesInfoDFrame 
#' @param metaDataFull 
#' @return sampleMetadataDframe
#' @importFrom stringr str_detect
.constructSampleMetadata <- function(samplesInfoDFrame, metaDataFull) {
    ## Exclude duplicated samples
    samplesInfoDFrame <- samplesInfoDFrame[
        !duplicated(samplesInfoDFrame[, "SRR_accession", drop = FALSE]), ,drop = FALSE]

    selectedSamplesSRRAccessions <- samplesInfoDFrame$SRR_accession

    samplesInfoDFrame <- 
        samplesInfoDFrame[
            samplesInfoDFrame$SRR_accession %in% selectedSamplesSRRAccessions,
            c('category', 'cell_type', 'general_cell_type', 
                'parental_cell_type', "SRR_accession"), drop = FALSE
        ]

    ## Add "SRR_accession" column to sampleMetadataDframe for merging 
    ## "samplesInfoDFrame" into "sampleMetadataDframe" 
    sampleMetadataDframe <- 
        metaDataFull[rownames(metaDataFull) %in% selectedSamplesSRRAccessions, , drop = FALSE]
    sampleMetadataDframe[, "SRR_accession"] <- rownames(metaDataFull)
    
    ## merge sampleMetadataDframe and sampleMetadataDframe to get the complete
    ## sample metedata
    sampleMetadataDframe <- merge(samplesInfoDFrame, sampleMetadataDframe, 
        by = 'SRR_accession')
    rownames(sampleMetadataDframe) <- sampleMetadataDframe[, "SRR_accession", ]    
    sampleMetadataDframe[, "SRR_accession"] <- NULL

    sampleMetadataDframe
}