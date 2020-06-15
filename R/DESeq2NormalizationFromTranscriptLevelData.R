#' Deseq2 normalization based on gene-level count data summarized from 
#' transcript-level estimates
#'
#' @param dee2Data
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom BiocGenerics estimateSizeFactors counts
#' @references https://bioconductor.org/packages/release/bioc/html/tximport.html
#' @references http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#' @references https://uclouvain-cbio.github.io/BSS2019/rnaseq_gene_summerschool_belgium_2019.html
Deseq2NormalizationFromTranscriptLevelData <- function(dee2Data) {
    ## DESeq2 performs independent filtering of lowly expressed genes internally,
    ## Thus we skip the filter step in DESeq2 normalization
    txi <- dee2Data$txi
    sExpr <- dee2Data$sExpr

    sExpr <- .filterData(sExpr)

    DESeq.ds <- DESeqDataSetFromTximport(
        txi = txi, 
        colData = colData(sExpr), 
        design=~1)

    # DESeq2 Default normalization method
    DESeq.dsDefault <- estimateSizeFactors(DESeq.ds)
    norm.counts <- counts(DESeq.dsDefault, normalized = TRUE)
    log.norm.counts <- log2(norm.counts + 1)    
    
    # reconstruct 'SExpr' parameter to include normalized data
    geneSel <- rownames(log.norm.counts)
    normalizedSExpr <- SummarizedExperiment(assays=list(counts=log.norm.counts),
        rowData=rowData(sExpr)[geneSel, ], colData=colData(sExpr))
    dee2Data[['normalizedSExpr']] <- normalizedSExpr
  
    dee2Data
}

