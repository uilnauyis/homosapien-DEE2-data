library("SummarizedExperiment")
library("readODS")
library("scatterD3")
library(plotly)
library (DESeq2)

DESeq2Normalization <- function(sExpt, iqr.cutoff=0.1, plot=c("2d", "3d")) {
  
  ###########################################################################
  ## PART 0. Check function arguments
  ###########################################################################
  fun.main <- deparse(match.call()[[1]])
    .stopIfNotExpressionSetOrSummarizedExperiment(inputObj, 'inputObj', fun.main)
    .stopIfNotNumeric0to1(iqr.cutoff, 'min.diff.cutoff', fun.main)
    sExpt <- .tryMakeSummarizedExperimentFromExpressionSet(inputObj)
  
  ###########################################################################
  ## PART I. Filter samples according to phenoData
  ###########################################################################
  ## o phenoData table contains which samples should be used in the analysis
  ##    o the samples which should be used in the analysis will have an
  ##     assigned category sExpt@phenoData$category, as "standard" or "test"
  ##    o non-assigned samples with NA values will be ignored
  ## o NOTE
  ##   assigning NA values to samples is an easy way to eliminate samples
  ##   from the analysis, without having to remove them from all input tables
  ##   (eg removing from sExpt, pdata, calls)
  phenoData <- colData(sExpt)
  ## filter out samples with missing category and/or general cell type
  phenoData.sel <- .filterPheno(phenoData, fun.main, "na")
  
  ############################################################################
  ## PART II. Filter genes in dataset
  ############################################################################
  ## A. For the cosine scoring, keep only the samples in the selected category
  ## (either 'standard' or 'test')

  if (sum(phenoData.sel) == 0)
    stop(paste("No samples in selected category found, exiting function",
               fun.main))

  sExpt <- sExpt[, phenoData.sel]

  # remove genes without any counts
  sExpt <- sExpt[rowSums(assay(sExpt, "counts")) > 0 , ]

  # Create DESeq.ds from the summarizedExperiment object
  DESeq.ds <- DESeqDataSetFromMatrix(countData = (assay(sExpt, "counts") + 1),
    colData = colData(sExpt),
    rowData = rowData(sExpt),
    design = ~cell_type)

  # DESeq2 Default normalization method
  DESeq.dsDefault <- estimateSizeFactors(DESeq.ds)
  counts.sf_normalized <- counts(DESeq.dsDefault, normalized = TRUE)
  log.norm.counts <- log2(counts.sf_normalized + 1)

  log.norm.counts
}



selectedDEE2SExpr <- readRDS("./rdsFiles/selectedDEE2SExpr.rds")
DESeq2Normalization()
  
  
  pcaPlot(log.norm.counts, pdata, plot)
  pcaPlot(vst.norm.counts, pdata, plot)
