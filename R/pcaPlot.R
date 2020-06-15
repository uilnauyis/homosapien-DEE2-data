#' PCA analysis of normalized count data and plot the result in 2D / 3D scatter
#' graphs
#'
#' @param sExpr  a SummarizedExperiment object that contains filtered and normalized
#'    count data along with gene information and sample metadata.
#' @param plot the type of plot of the PCA analysis result
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom stats prcomp
#' @importFrom scatterD3 scatterD3
#' @importFrom plotly layout plot_ly add_markers
#' @importFrom dplyr '%>%'
#' @importFrom matrixStats rowVars

pcaPlot <- function(sExpr, plot = c("3d", "2d")) {
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
  
  normalizedCount <- assay(sExpr, 'counts')
  pdata <- colData(sExpr)
  
  # remove the gene of which the variance of counts accross samples equals to 0  
  normalizedCount <- normalizedCount[rowVars(normalizedCount) > 0, ]

  ## Remove any rows that have zero variance
  pca <- prcomp(t(normalizedCount), scale=TRUE)
  
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
  pca.comp <- merge(pca.comp, pdata, by=0, all=TRUE)
  
  tooltips <- paste(" <strong>", rownames(pca.comp),"</strong><br />", 
                    pca.comp$cell_type)
  
  if (plot == "2d") {
    ## 2D plot
    scatterD3(data = pca.comp, x=PC1, y=PC2, tooltip_text = tooltips,
              col_var=cell_type)
    
    scatterD3(data = pca.comp, x=PC1, y=PC3, tooltip_text = tooltips,
              col_var=cell_type)
    
    scatterD3(data = pca.comp, x=PC2, y=PC3, tooltip_text = tooltips,
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