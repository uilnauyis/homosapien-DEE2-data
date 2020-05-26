#' Construct the SummarizedExperiment object for the raw counts data.
#'
#' @param sExpr SummarizedExperiment object that is created via DEE2 R interface 
#'    or restored from downloaded files.
#' @param counts.cutoff threshold of filtering the 
#' @export  
#' @param sampleInfoPath a path to a .ods file that stores category, 
#'      general_cell_type, parental_cell_type and SRR_accession for each 
#'      sample
#' @param dee2DataDir directory to store files storing DEE2 data
#' @param species the species of the DEE2 data
#' @param overwrite a boolean parameter. If it is true, 'constructSExprFromDEE2Data'
#'      function will overwrite the files in the directory specified by parameter
#'      dee2DataDir. By default, the value is false.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#' @importFrom rio import
#' @importFrom getDEE2 getDEE2

constructSExprFromDEE2Data <- function(
    sampleInfoPath, 
    dee2DataDir, 
    species = c('athaliana', 'celegans','dmelanogaster','drerio','ecoli',
        'hsapiens','mmusculus','rnorvegicus','scerevisiae'), 
    overwrite = FALSE) {

    metaDataFull <- NULL
    geneCountsDFrame <- NULL
    geneInfoDFrame <- NULL
    sampleMetadataDframe <- NULL
    
    samplesInfoDFrame <- import(sampleInfoPath)

    sExprPath <- paste(dee2DataDir, "/selectedDEE2SExpr.rds", sep='')
    if (overwrite == FALSE && 
        file.exists(sExprPath)) {
        sExpr <- readRDS(sExprPath)
        return(sExpr)
    }

    # If the directory storing DEE2 data does not exist, or the directory is 
    # empty, or the user want to overwrite the existing data
    if (!dir.exists(dee2DataDir) || 
        length(list.files(dee2DataDir)) == 0 ||
        overwrite == TRUE) {
        dee2Data <- .downloadDEE2Data(species, samplesInfoDFrame, dee2DataDir)
        geneCountsDFrame <- dee2Data[["GeneCounts"]]
        geneInfoDFrame <- dee2Data[["GeneInfo"]]
    }

    ## get counts data from 'GeneCountMatrix.tsv' in the dee2 data folder if 
    ## it was not restored
    if (is.null(geneCountsDFrame)) {
        geneCountsDFrame <- read.table(
            file = paste(dee2DataDir, "/GeneCountMatrix.tsv", sep=''), 
            sep = '\t', 
            row.names=1, 
            header=TRUE)
    }

    counts <- data.matrix(geneCountsDFrame)

    ## get gene metadata
    if (is.null(geneInfoDFrame)) {
        geneInfoDFrame <- read.table(
            file = paste(dee2DataDir, "/GeneInfo.tsv", sep=''), 
            sep = '\t', 
            row.names=1, 
            header=TRUE)
    }

    ## construct metadata dataframe for the selected samples
    metaDataFull <- read.table(
        file = paste(dee2DataDir, "/MetadataFull.tsv", sep=''), 
        sep = '\t', 
        row.names=1, 
        header=TRUE, 
        comment.char = '' )

    sampleMetadataDframe <- .constructSampleMetadata(sampleInfoPath, metaDataFull)

    ## construct sExpr object for the selected DEE2 data
    orderSamples <- colnames(counts)
    orderGenes <- rownames(counts)

    sExpr <- SummarizedExperiment(assays=list(counts=counts))
    colData(sExpr) <- DataFrame(sampleMetadataDframe[orderSamples, ])
    rowData(sExpr) <- DataFrame(geneInfoDFrame[orderGenes, ])

    saveRDS(sExpr, sExprPath)

    sExpr
}


.downloadDEE2Data <- function(species, samplesInfoDFrame, dee2DataDir) {
    SRR_accessions <- unique((samplesInfoDFrame[['SRR_accession']]))

    # Create data folder if not exists
    dir.create(dee2DataDir)

    zipFilePath <- paste(dee2DataDir, '/data.zip', sep = '')
    dee2Data <- getDEE2(species, SRR_accessions, outfile = zipFilePath)

    unzip(zipFilePath, exdir = dee2DataDir)

    dee2Data
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

    ## Use SRR accessions as rownames and remove "SRR_accession" column for combining
    ## merging the 'category', 'general_cell_type' and 'parental_cell_type' with 
    ## other metadata
    rownames(samplesInfoDFrame) <- samplesInfoDFrame[, "SRR_accession"]

    ## Add "SRR_accession" column to sampleMetadataDframe 
    sampleMetadataDframe <- 
        metaDataFull[rownames(metaDataFull) %in% selectedSamplesSRRAccessions, ]
    sampleMetadataDframe[, "SRR_accession"] <- rownames(metaDataFull)

    rownames(sampleMetadataDframe) <- sampleMetadataDframe[, "SRR_accession"]
    ## merge sampleMetadataDframe and sampleMetadataDframe to get the complete
    ## sample metedata
    sampleMetadataDframe <- merge(samplesInfoDFrame, sampleMetadataDframe, 
        by = 'SRR_accession')
    rownames(sampleMetadataDframe) <- sampleMetadataDframe[, "SRR_accession"]    
    sampleMetadataDframe[, "SRR_accession"] <- NULL

    sampleMetadataDframe
}
