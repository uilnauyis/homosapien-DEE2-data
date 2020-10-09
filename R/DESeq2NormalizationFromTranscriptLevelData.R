#' DESeq2 normalization based on gene-level count data summarized from 
#' transcript-level estimates
#'
#' @param dee2Data SummarizedExperiment
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom BiocGenerics estimateSizeFactors counts
#' @references https://bioconductor.org/packages/release/bioc/html/tximport.html
#' @references http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#' @references https://uclouvain-cbio.github.io/BSS2019/rnaseq_gene_summerschool_belgium_2019.html
DESeq2NormalizationFromTranscriptLevelData <- function(dee2Data, counts.cutoff = 10) {
    ## DESeq2 performs independent filtering of lowly expressed genes internally,
    ## Thus we skip the filter step in DESeq2 normalization
    txi <- dee2Data$txi
    sExpr <- dee2Data$sExpr

    sExpr <- .filterData(sExpr, counts.cutoff)

    dESeq.ds <- DESeq2::DESeqDataSetFromTximport(
        txi = txi, 
        colData = SummarizedExperiment::colData(sExpr), 
        design=~1)

    # DESeq2 Default normalization method
    dESeq.dsDefault <- BiocGenerics::estimateSizeFactors(dESeq.ds)
    norm.counts <- BiocGenerics::counts(dESeq.dsDefault, normalized = TRUE)
    log.norm.counts <- log2(norm.counts + 1)    
    
    # reconstruct 'SExpr' parameter to include normalized data
    geneSel <- rownames(log.norm.counts)
    normalizedSExpr <- SummarizedExperiment::SummarizedExperiment(
        assays=list(counts=log.norm.counts),
        rowData=SummarizedExperiment::rowData(sExpr)[geneSel, ], 
        colData=SummarizedExperiment::colData(sExpr))

    ## Return the normalized SummarizedExperiment object
    normalizedSExpr
}

