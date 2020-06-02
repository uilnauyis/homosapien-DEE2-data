#' edgeR normalization based on Transcript-level count data
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
edgeRTranscriptAnalysis <- function(species, srrAccessions, txInfo, outDir) {
    txi <- .transcriptLevelAnalysis(species, srrAccessions, txInfo, outDir = NULL)

    cts <- txi$counts
    normMat <- txi$length

    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normMat <- normMat/exp(rowMeans(log(normMat)))
    normCts <- cts/normMat

    # Computing effective library sizes from scaled counts, to account for
    # composition biases between samples.
    library(edgeR)
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
    library(csaw)
    cpms <- calculateCPM(sExpr, use.offsets = TRUE, log = FALSE)
}

