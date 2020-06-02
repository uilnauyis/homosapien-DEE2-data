################################################################################
## FUNCTION: .filterData
################################################################################
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

################################################################################
## FUNCTION: .transcriptLevelAnalysis
################################################################################
#' @param species
#' @param srrAccessions
#' @param outDir
#' @param txInfo 
#' @return txi.kallisto.tsv
#' @importFrom tximport tximport
#' @importFrom SummarizedExperiment SummarizedExperiment
.transcriptLevelAnalysis <- function(species, srrAccessions, txInfo, outDir = NULL) {
    files <- .downloadAbundanceTsv(species, srrAccessions, outDir)

    txInfo <- data.frame(TxID = rownames(txInfo), txInfo)[, 1:2]
    rownames(txInfo) <- c(1:nrow(txInfo))

    txi.kallisto.tsv <- tximport(files, 
        type = "kallisto", 
        tx2gene = txInfo, 
        ignoreAfterBar = TRUE)
    txi.kallisto.tsv
}

################################################################################
## FUNCTION: .downloadAbundanceTsv
################################################################################
#' @param species
#' @param srrAccessions
#' @param outDir 
#' @return abundancefiles 
.downloadAbundanceTsv <- function(species, srrAccessions, outDir = NULL) {
    abundancefiles <- lapply(srrAccessions, function(srrAccession) {
        abundanceUrl <- file.path('http://dee2.io/data', 
            species, srrAccession, paste(srrAccession, '.ke.tsv.gz', sep = ""))
        abundanceName <- NULL
        if(is.null(outDir)){
            outDir <- tempfile()
            abundanceName <- paste(outDir, '.abundance.tsv.gz', sep = "")
        } else {
            abundanceName <- file.path(outDir, 
                paste(srrAccession, '.abundance.tsv.gz', sep = ""))
        }
        download.file(abundanceUrl, destfile=abundanceName)
        abundanceDf <- read.table(gzfile(abundanceName)) 
        abundanceDf <- abundanceDf[2:nrow(abundanceDf), ]
        colnames(abundanceDf) <- c('target_id', 'length', 'eff_length', 
            'est_counts', 'tpm')
        write.table(abundanceDf, 
            file=gzfile(abundanceName), 
            row.names=FALSE, 
            quote=FALSE, 
            sep='\t')
        abundanceName
    })
    abundancefiles <- file.path(abundancefiles)
    names(abundancefiles) <- srrAccessions
    abundancefiles
}
