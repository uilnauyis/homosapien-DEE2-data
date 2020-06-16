#' Histogram Plots of Normalized Counts
#'
#' @param sExpr  SummarizedExperiment object that is created via DEE2 R interface 
#' @param plot the type of plot of the PCA analysis result
#' @export
#' 
NormalizedCountsHistogramPlot <- function(dee2Data, breaks = 100) {
    counts <- assay(dee2Data$normalizedSExpr, 'counts')

    sapply(colnames(counts),
        function(colName) {
            hist(counts[, colName],  
                main = colName,
                breaks = breaks)
    })
}