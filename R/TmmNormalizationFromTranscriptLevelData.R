#' TMM normalization based on gene-level count data summarized from 
#' transcript-level estimates
#'
#' @param species
#' @param srrAccessions
#' @param txInfo
#' @param outDir
#' @export 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom edgeR calcNormFactors DGEList scaleOffset filterByExpr
#' @importFrom csaw calculateCPM
#' @references https://bioconductor.org/packages/release/bioc/html/tximport.html
TmmNormalizationFromTranscriptLevelData <- function(dee2Data) {
    txi <- dee2Data$txi

    cts <- txi$counts
    normMat <- txi$length

    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normMat <- normMat/exp(rowMeans(log(normMat)))
    normCts <- cts/normMat

    # Computing effective library sizes from scaled counts, to account for
    # composition biases between samples.
    eff.lib <- calcNormFactors(normCts) * colSums(normCts)

    # Combining effective library sizes with the length factors, and calculating
    # offsets for a log-link GLM.
    normMat <- sweep(normMat, 2, eff.lib, "*")
    normMat <- log(normMat)

    # Creating a DGEList object for use in edgeR.
    y <- DGEList(cts)
    y <- scaleOffset(y, normMat)
    # filtering
    keep <- filterByExpr(y)
    y <- y[keep, ]

    sExpr <- SummarizedExperiment(
        assays = list(counts = y$counts, offset = y$offset))
    sExpr$totals <- y$samples$lib.size

    log.cpms <- calculateCPM(sExpr, use.offsets = TRUE, log = TRUE)

    normalizedSExpr <- sExpr
    assay(normalizedSExpr, 'counts', withDimnames=FALSE) <- log.cpms
    colData(normalizedSExpr) <- colData(dee2Data$sExpr)
    rowDataRaw <- rowData(dee2Data$sExpr)
    rowData(normalizedSExpr) <- 
        rowDataRaw[rownames(rowDataRaw) %in% rownames(normalizedSExpr), ]
    dee2Data[['normalizedSExpr']] <- normalizedSExpr
  
    dee2Data 
}