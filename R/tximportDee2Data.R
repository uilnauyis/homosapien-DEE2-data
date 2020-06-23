#' Download abundance.tsv files (transcript-level estimate) for specified species 
#' and SSR Accessions and get the gene-level estimate
#' 
#' @param species
#' @param srrAccessions
#' @param outDir
#' @param txInfo dee2Data$TxInfo
#' @export 
#' @return dee2Data
#' @importFrom tximport tximport
#' @import readr 
#' @references https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
#' @examples 
#' 
#' ## This operation may requires up to 8GB of memory
#' txi <- tximportDee2Data('hsapiens', dee2Data, )
tximportDee2Data <- function(species, dee2Data, outDir = NULL) {
    srrAccessions <- colnames(dee2Data$GeneCounts)
    txInfo <- dee2Data$TxInfo


    files <- .downloadAbundanceTsv(species, srrAccessions, outDir)

    # Subset txInfor object to retain the first two columns only, which are for 
    # transcript IDs and gene IDS
    txInfo <- data.frame(TxID = rownames(txInfo), txInfo)[, 1:2]
    rownames(txInfo) <- c(1:nrow(txInfo))

    txi.kallisto.tsv <- tximport::tximport(files, 
        type = "kallisto", 
        tx2gene = txInfo, 
        ignoreAfterBar = TRUE)

    ## Remove version numbers in Gene IDs
    rownames(txi.kallisto.tsv$counts) <- 
        lapply(
            rownames(txi.kallisto.tsv$counts), 
            sub, 
            pattern = "\\.\\d+$", 
            replacement = "")

    rownames(txi.kallisto.tsv$abundance) <- 
        lapply(
            rownames(txi.kallisto.tsv$abundance), 
            sub, 
            pattern = "\\.\\d+$", 
            replacement = "")

    rownames(txi.kallisto.tsv$length) <- 
        lapply(
            rownames(txi.kallisto.tsv$length), 
            sub, 
            pattern = "\\.\\d+$", 
            replacement = "")

    dee2Data[['txi']] <- txi.kallisto.tsv
    dee2Data
}