#' PCA analysis of the count data and plot the result in 2D / 3D scatter graphs
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

PcaPlot <- function(sExpr, plot = '3d') {
  
  # Check plot option variable, which should be either '3d' or '2d' 
  stopifnot(plot == '2d' || plot == '3d')
  
  counts <- data.matrix(assay(sExpr, 'counts'))
  pdata <- colData(sExpr)
  
  # remove the gene of which the variance of counts accross samples equals to 0  
  counts <- counts[rowVars(counts) > 0, ]

  pca <- prcomp(t(counts), scale=TRUE)
  
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
  
  ## Text on hover
  text.hover <- ~paste('Cell type:', cell_type, 
         '<br>SRR Accession:', SRR_accession, 
         '<br>Category:', category)
  
  symbol.factor <- ~category
  symbols = c('circle','x','o')
  
  if (plot == "2d") {

    ## 2D plot
    fig1 <- plot_ly(data = pca.comp, x = ~PC1, y = ~PC2, color = ~cell_type, 
                    text = text.hover, symbol = symbol.factor, 
                    symbols = symbols)
    fig1 <- fig1 %>% layout(title="PCA plot (PC1, PC2)",
                            xaxis = list(title = paste('PC1(', pc1ExplainedVar, '%)', 
                                                       sep = '')),
                            yaxis = list(title = paste('PC2(', pc2ExplainedVar, '%)', 
                                                       sep = '')))

    fig2 <- plot_ly(data =  pca.comp, x = ~PC1, y = ~PC3, color = ~cell_type, 
                    text = text.hover, symbol = symbol.factor, 
                    symbols = symbols)
    fig2 <- fig2 %>% layout(title="PCA plot (PC1, PC3)",
                            xaxis = list(title = paste('PC1(', pc1ExplainedVar, '%)', 
                                                       sep = '')),
                            yaxis = list(title = paste('PC3(', pc3ExplainedVar, '%)', 
                                                       sep = '')))

    fig3 <- plot_ly(data =  pca.comp, x = ~PC2, y = ~PC3, color = ~cell_type, 
                    text = text.hover, symbol = symbol.factor, 
                    symbols = symbols)
    fig3 <- fig3 %>% layout(title="PCA plot (PC2, PC3)",
                            xaxis = list(title = paste('PC2(', pc2ExplainedVar, '%)', 
                                                       sep = '')),
                            yaxis = list(title = paste('PC3(', pc3ExplainedVar, '%)', 
                                                       sep = '')))

    return(list(fig1, fig2, fig3))
  } else {
    ## 3D plot
    fig <- plot_ly(pca.comp, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cell_type, 
                   text = text.hover, symbol = symbol.factor, 
                   symbols = symbols)
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = paste('PC1(', pc1ExplainedVar, '%)', sep = '')),
                                  yaxis = list(title = paste('PC2(', pc2ExplainedVar, '%)', sep = '')),
                                  zaxis = list(title = paste('PC3(', pc3ExplainedVar, '%)', sep = ''))))
    
    fig
  }
}