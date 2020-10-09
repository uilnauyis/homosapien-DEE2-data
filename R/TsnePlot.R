#' tsne analysis of the count data
#'
#' @param sExpr  a SummarizedExperiment object that contains count data along 
#' with gene information and sample metadata.
#' @param plot the option of the plot of the PCA analysis result, which should
#' be either '2d' or '3d'
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom plotly layout plot_ly add_markers
#' @importFrom dplyr '%>%'
#' @importFrom matrixStats rowVars
#' @importFrom Rtsne Rtsne

TsnePlot <- function(sExpr, perplexity, steps, dim) {
  counts <- data.matrix(assay(sExpr, 'counts'))
  pdata <- colData(sExpr)
  
  # remove the gene of which the variance of counts accross samples equals to 0  
  counts <- counts[rowVars(counts) > 0, ]
  
  # Sets seed for reproducibility
  set.seed(41) 
  
  # Get the result of Tsne analysis
  tsne <- Rtsne(t(counts), check_duplicates = FALSE, perplexity = perplexity,
                dim = dim, steps = steps)
  tsne$Y <- as.data.frame(tsne$Y)
  rownames(tsne$Y) <- rownames(pdata)
  tsne.comp <- merge(tsne$Y, pdata, by=0, all=TRUE)
  
  tooltips <- paste(" <strong>", rownames(tsne.comp),"</strong><br />", 
                    tsne.comp$cell_type)
  
  ## Text on hover
  text.hover <- ~paste('Cell type:', cell_type, 
                       '<br>SRR Accession:', SRR_accession, 
                       '<br>Category:', category)
  
  symbol.factor <- ~category
  symbols = c('circle','x','o')
  
  
  ## 2D plot
  fig <- plot_ly(data = tsne.comp, x = ~V1, y = ~V2, color = ~cell_type, 
                  text = text.hover, symbol = symbol.factor, 
                  symbols = symbols)
  fig <- fig %>% layout(title="Tsne plot (factor1, factor2)",
                          xaxis = list(title = 'Factor 1'),
                          yaxis = list(title = 'Factor 2'))
  fig
}