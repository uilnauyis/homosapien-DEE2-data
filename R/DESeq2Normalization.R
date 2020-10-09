#' DESeq2 normalization for the selected DEE2 data
#' 'In some experiment, there might be gene-dependent dependencies which vary 
#' across samples. For instance, GC-content bias or length bias might vary across 
#' samples coming from different labs or processed at different times...'
#'
#' @param dee2Data SummarizedExperiment
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom BiocGenerics estimateSizeFactors counts
#' @references http://www.sthda.com/english/wiki/rna-sequencing-data-analysis-counting-normalization-and-differential-expression#normalization-using-deseq2-size-factors
#' @references https://uclouvain-cbio.github.io/BSS2019/rnaseq_gene_summerschool_belgium_2019.html

DESeq2Normalization <- function(sExpr, counts.cutoff = 10) {
  sExpr <- .filterData(sExpr, counts.cutoff)

  # Create dESeq.ds from the summarizedExperiment object
  dESeq.ds <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(sExpr, "counts"),
    colData = SummarizedExperiment::colData(sExpr), 
    rowData = SummarizedExperiment::rowData(sExpr),
    design = ~1)

  # DESeq2 Default normalization method
  dESeq.dsDefault <- BiocGenerics::estimateSizeFactors(dESeq.ds)
  norm.counts <- BiocGenerics::counts(dESeq.dsDefault, normalized = TRUE)
  log.norm.counts <- log2(norm.counts + 1) 
  normalizedSExpr <- sExpr
  SummarizedExperiment::assay(normalizedSExpr, 'counts') <- log.norm.counts
  
  ## Return the normalized SummarizedExperiment object
  normalizedSExpr
}