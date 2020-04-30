library("SummarizedExperiment")
library("readODS")
library("scatterD3")
library(plotly)

pcaAnalysis <- function(inputObj, iqr.cutoff=0.1, category.selected = c("standard", "test"),
                        plot=c("2d", "3d")) {
  
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
  pdata <- colData(sExpt)
  ## filter out samples with missing category and/or general cell type
  pdata.sel <- .filterPheno(pdata, fun.main, "na")
  
  ############################################################################
  ## PART II. Filter genes in dataset
  ############################################################################
  ## A. For the cosine scoring, keep only the samples in the selected category
  ## (either 'standard' or 'test')
  sel <- pdata.sel$category %in% category.selected
  if (sum(sel) == 0)
    stop(paste("No samples in selected category found, exiting function",
               fun.main))
  pdata.sel <- pdata.sel[sel, ]
  ynorm <- assay(sExpt[, sel], "counts")
  
  ## B. Filter out not variable probes,
  ##    o variance in terms of IQR of median expressions
  ##    o want to get rid of probes that are below a given IQR threshold
  
  ## Calculate IQR of the samples and filter by the given IQR threshold
  medIQR <- apply(ynorm, 1, IQR)
  selGenes <- medIQR >= quantile(medIQR, probs=1 - iqr.cutoff)
  if (sum(selGenes) == 0) {
    ## this should almost never happen
    stop(paste("No gene passed the IQR-based filtering, exiting function",
               fun.main))
  }
  
  ## We keep for score calculation the gene-filtered dataset with standards
  ## and test samples
  ynormIQR <- assay(sExpt[selGenes, rownames(pdata.sel)], "counts")
  
  ############################################################################
  ## PART III. PCA Analysis of the filtered data
  ############################################################################
  ##  1. Euclidean distance between categories
  ##  2. Cosine similarity between categories
  ## -DONT USE-3. cosine "score" between categories
  ##      o 0 < cosine score <1 where range is set to 1-min(cosine.similarity)
  ## What is this good for?
  ##   a. Visualization and exploratory data analysis
  ##   b. Generates values to be used for cell scoring
  
  ## Remove any rows that have zero variance
  sel <- apply(ynormIQR, 1, function(x) sd(x) != 0)
  ynormIQR <- ynormIQR[sel,]
  
  pca <- prcomp(t(ynormIQR), scale=TRUE)
  
  ## Get proportion of variance explained by PC1 and PC2
  pca.sum <- summary(pca)
  pc1 <- pca.sum$importance[2,1] * 100
  pc2 <- pca.sum$importance[2,2] * 100
  pc3 <- pca.sum$importance[2,3] * 100
  
  print(pc1)
  print(pc2)
  print(pc3)
    
  ## Extract coordinates from pca object
  pca.comp <- as.data.frame(pca$x[, c(1, 2, 3)])
  pca.comp <- cbind(pca.comp,
    cell_type = pdata.sel[, "cell_type"])
  
  tooltips <- paste(" <strong>", rownames(pca.comp),"</strong><br />", pca.comp$cell_type)

  if (plot == "2d") {
    ## 2D plot
    scatterD3(data = pca.comp, x=PC1, y=PC2, tooltip_text = tooltips,
            col_var=cell_type)
  } else {
    ## 3D plot
    fig <- plot_ly(pca.comp, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cell_type)
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = paste('PC1')),
                                     yaxis = list(title = ('PC2')),
                                     zaxis = list(title = ('PC3'))))
  
    fig
  }
}



selectedDEE2SExpr <- readRDS("./rdsFiles/selectedDEE2SExpr.rds")
pcaAnalysis(selectedDEE2SExpr, 0.1, "standard", "3d")
