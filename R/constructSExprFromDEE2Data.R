#' Construct the SummarizedExperiment object for the raw counts data.
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
#' @importFrom rio import
#' @importFrom getDEE2 getDEE2
#' @references https://github.com/markziemann/dee2/blob/master/AccessDEEfromR.md

constructSExprFromDEE2Data <- function(
    sampleInfoPath, 
    outFile,
    species = c('athaliana', 'celegans', 'dmelanogaster', 'drerio', 'ecoli',
        'hsapiens', 'mmusculus', 'rnorvegicus', 'scerevisiae')) {   
    samplesInfoDFrame <- import(sampleInfoPath)

    dee2Data <- .downloadDEE2Data(species, samplesInfoDFrame, outFile)
    geneCountsDFrame <- dee2Data[["GeneCounts"]]
    geneInfoDFrame <- dee2Data[["GeneInfo"]]

    counts <- data.matrix(geneCountsDFrame)

    sampleMetadataDframe <- .constructSampleMetadata(sampleInfoPath, 
        dee2Data$MetadataFull)

    ## construct sExpr object for the selected DEE2 data
    orderSamples <- colnames(counts)
    orderGenes <- rownames(counts)

    sExpr <- SummarizedExperiment(assays=list(counts=counts))
    colData(sExpr) <- DataFrame(sampleMetadataDframe[orderSamples, ])
    rowData(sExpr) <- DataFrame(geneInfoDFrame[orderGenes, ])

    dee2Data[['sExpr']] <- sExpr

    dee2Data
}

.downloadDEE2Data <- function(species, samplesInfoDFrame, outFile) {
    SRR_accessions <- unique((samplesInfoDFrame[['SRR_accession']]))
    dee2Data <- getDEE2(species, SRR_accessions, outfile = outFile)
}

.constructSampleMetadata <- function(sampleInfoPath, metaDataFull) {
  
    ## get sample information from the file
    samplesInfoDFrame <- import(sampleInfoPath)

    ## Exclude duplicated samples
    samplesInfoDFrame <- samplesInfoDFrame[
        !duplicated(samplesInfoDFrame[, "SRR_accession"]), ]

    selectedSamplesSRRAccessions <- samplesInfoDFrame$SRR_accession

    samplesInfoDFrame <- 
        samplesInfoDFrame[
            samplesInfoDFrame$SRR_accession %in% selectedSamplesSRRAccessions,
            c('category', 'cell_type', 'general_cell_type', 
                'parental_cell_type', "SRR_accession")
        ]

    ## Add "SRR_accession" column to sampleMetadataDframe for merging 
    ## "samplesInfoDFrame" into "sampleMetadataDframe" 
    sampleMetadataDframe <- 
        metaDataFull[rownames(metaDataFull) %in% selectedSamplesSRRAccessions, ]
    sampleMetadataDframe[, "SRR_accession"] <- rownames(metaDataFull)
    
    ## merge sampleMetadataDframe and sampleMetadataDframe to get the complete
    ## sample metedata
    sampleMetadataDframe <- merge(samplesInfoDFrame, sampleMetadataDframe, 
        by = 'SRR_accession')
    rownames(sampleMetadataDframe) <- sampleMetadataDframe[, "SRR_accession"]    
    sampleMetadataDframe[, "SRR_accession"] <- NULL

    sampleMetadataDframe
}