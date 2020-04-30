library("SummarizedExperiment")
library("readODS")

constructSExprFromDEE2Data <- function() {
    ## get counts
    geneCountsDFrame <- read.table(file = "./DEE2Data/GeneCountMatrix.tsv", sep = '\t', 
        row.names=1, header=TRUE)

    counts <- data.matrix(geneCountsDFrame)

    ## get rowData
    geneInfoDFrame <- read.table(file = "./DEE2Data/GeneInfo.tsv", sep = '\t', 
        row.names=1, header=TRUE)[, "GeneSymbol", drop=FALSE]

    ## get colData
    samplesDFrame <- read_ods(path = "./DEE2Data/SelectedDEE2_v2.ods", 
        sheet = 1, col_names = TRUE)

    samplesDFrame <- samplesDFrame[!duplicated(samplesDFrame[, "SRR_accession"]), ]

    validSamplesSRRAccessions <- colnames(geneCountsDFrame)

    validSamplesDFrame <- 
        samplesDFrame[samplesDFrame$SRR_accession %in% validSamplesSRRAccessions,]

    rownames(validSamplesDFrame) <- validSamplesDFrame[, "SRR_accession"]
    validSamplesDFrame[, "validSamplesDFrame"] <- NULL
    validSamplesDFrame <- validSamplesDFrame[validSamplesSRRAccessions, ]

    ## construct sExpt object for the selected DEE2 data
    SummarizedExperiment(assays=list(counts=counts),
                         rowData=geneInfoDFrame, colData=validSamplesDFrame)
}

sExpr <- constructSExprFromDEE2Data()
saveRDS(sExpr, "./rdsFiles/selectedDEE2SExpr.rds")