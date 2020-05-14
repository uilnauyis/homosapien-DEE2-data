library("SummarizedExperiment")
library("rio")
library("getDEE2")

constructSExprFromDEE2Data <- function() {
    unzip('./SelectedDEE2_v2.zip', exdir = 'SelectedDEE2_v2_csv')

    ## get counts
    geneCountsDFrame <- read.table(
        file = "./SelectedDEE2_v2_csv/GeneCountMatrix.tsv", 
        sep = '\t', 
        row.names=1, 
        header=TRUE)

    counts <- data.matrix(geneCountsDFrame)

    ## get rowData
    geneInfoDFrame <- read.table(
        file = "./SelectedDEE2_v2_csv/GeneInfo.tsv", 
        sep = '\t', 
        row.names=1, 
        header=TRUE)

    ## get colData
    samplesDFrame <- import("./SelectedDEE2_v2.ods")

    samplesDFrame <- samplesDFrame[!duplicated(samplesDFrame[, "SRR_accession"]), ]

    validSamplesSRRAccessions <- colnames(geneCountsDFrame)

    validSamplesDFrame <- 
        samplesDFrame[samplesDFrame$SRR_accession %in% validSamplesSRRAccessions,]

    rownames(validSamplesDFrame) <- validSamplesDFrame[, "SRR_accession"]
    validSamplesDFrame[, "validSamplesDFrame"] <- NULL
    validSamplesDFrame <- validSamplesDFrame[validSamplesSRRAccessions, ]

    ## construct sExpt object for the selected DEE2 data
    sExpr <- SummarizedExperiment(assays=list(counts=counts),
                         rowData=geneInfoDFrame, colData=validSamplesDFrame)

    saveRDS(sExpr, "./selectedDEE2SExpr.rds")

    sExpr
}

downloadDEE2Data <- function() {
    dFrame <- import('./SelectedDEE2_v2.ods')
    SRR_accessions <- unique((dFrame[['SRR_accession']]))

    dee2Data <- getDEE2('hsapiens', SRR_accessions, 
        outfile='./SelectedDEE2_v2')
}

downloadDEE2Data()
sExpr <- constructSExprFromDEE2Data()
