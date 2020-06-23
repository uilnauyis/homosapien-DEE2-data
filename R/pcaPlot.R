#' PCA analysis of normalized count data and plot the result in 2D / 3D scatter
#' graphs
#'
#' @param sExpr  a SummarizedExperiment object that contains filtered and normalized
#'    count data along with gene information and sample metadata.
#' @param plot the type of plot of the PCA analysis result
#' @export
#' @importFrom SummarizedExperiment colData assay
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
  
  normalizedCount <- SummarizedExperiment::assay(sExpr, 'counts')
  pdata <- SummarizedExperiment::colData(sExpr)
  
  # remove the gene of which the variance of counts accross samples equals to 0  
  normalizedCount <- normalizedCount[matrixStats::rowVars(normalizedCount) > 0, ]

  pca <- prcomp(t(normalizedCount), scale=TRUE)
  
  ## Get proportion of variance explained by PC1 and PC2
  pca.sum <- summary(pca)
  pc1ExplainedVar <- pca.sum$importance[2,1] * 100
  pc2ExplainedVar <- pca.sum$importance[2,2] * 100
  pc3ExplainedVar <- pca.sum$importance[2,3] * 100
  
  ## Extract coordinates from pca object
  pca.comp <- as.data.frame(pca$x[, c(1, 2, 3)])
  pca.comp <- merge(pca.comp, pdata, by=0, all=TRUE)
  
  tooltips <- paste(" <strong>", rownames(pca.comp),"</strong><br />", 
                    pca.comp$cell_type)
  
  if (plot == "2d") {

    ## 2D plot
    fig1 <- plotly::plot_ly(data =  pca.comp, x = ~PC1, y = ~PC2, color = ~cell_type)
    fig1 <- fig1 %>% plotly::layout(title="PCA plot (PC1, PC2)",
                            xaxis = list(title = paste('PC1(', pc1ExplainedVar, '%)', sep = '')),
                            yaxis = list(title = paste('PC2(', pc2ExplainedVar, '%)', sep = '')))

    fig2 <- plotly::plot_ly(data =  pca.comp, x = ~PC1, y = ~PC3, color = ~cell_type)
    fig2 <- fig2 %>% plotly::layout(title="PCA plot (PC1, PC3)",
                            xaxis = list(title = paste('PC1(', pc1ExplainedVar, '%)', sep = '')),
                            yaxis = list(title = paste('PC3(', pc3ExplainedVar, '%)', sep = '')))

    fig3 <- plotly::plot_ly(data =  pca.comp, x = ~PC2, y = ~PC3, color = ~cell_type)
    fig3 <- fig3 %>% plotly::layout(title="PCA plot (PC2, PC3)",
                            xaxis = list(title = paste('PC2(', pc2ExplainedVar, '%)', sep = '')),
                            yaxis = list(title = paste('PC3(', pc3ExplainedVar, '%)', sep = '')))

    return(list(fig1, fig2, fig3))
  } else {
    ## 3D plot
    fig <- plotly::plot_ly(pca.comp, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cell_type)
    fig <- fig %>% plotly::add_markers()
    fig <- fig %>% plotly::layout(scene = list(xaxis = list(title = paste('PC1(', pc1ExplainedVar, '%)', sep = '')),
                                       yaxis = list(title = paste('PC2(', pc2ExplainedVar, '%)', sep = '')),
                                       zaxis = list(title = paste('PC3(', pc3ExplainedVar, '%)', sep = ''))))
    
    fig
  }
}