#' Deseq2 normalization for the selected DEE2 data
#' 'In some experiment, there might be gene-dependent dependencies which vary 
#' across samples. For instance, GC-content bias or length bias might vary across 
#' samples coming from different labs or processed at different times...'
#'
#' @param dee2Data
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom BiocGenerics estimateSizeFactors counts
#' @references http://www.sthda.com/english/wiki/rna-sequencing-data-analysis-counting-normalization-and-differential-expression#normalization-using-deseq2-size-factors
#' @references https://uclouvain-cbio.github.io/BSS2019/rnaseq_gene_summerschool_belgium_2019.html

Deseq2Normalization <- function(dee2Data, counts.cutoff = 10) {
  ## DESeq2 performs independent filtering of lowly expressed genes internally,
  ## Thus we skip the filter step in DESeq2 normalization
  sExpr <- dee2Data$sExpr

  sExpr <- .filterData(sExpr, counts.cutoff)

  # Create DESeq.ds from the summarizedExperiment object
  DESeq.ds <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(sExpr, "counts"),
    colData = SummarizedExperiment::colData(sExpr), 
    rowData = SummarizedExperiment::rowData(sExpr),
    design = ~1)

  # DESeq2 Default normalization method
  DESeq.dsDefault <- BiocGenerics::estimateSizeFactors(DESeq.ds)
  norm.counts <- BiocGenerics::counts(DESeq.dsDefault, normalized = TRUE)
  log.norm.counts <- log2(norm.counts + 1)

  # 
  normalizedSExpr <- sExpr
  SummarizedExperiment::assay(normalizedSExpr, 'counts') <- log.norm.counts
  dee2Data[['normalizedSExpr']] <- normalizedSExpr

  dee2Data
}