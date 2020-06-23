#' Histogram Plots of Normalized Counts
#'
#' @param sExpr  SummarizedExperiment object that is created via DEE2 R interface
#' @export
#' 
NormalizedCountsHistogramPlot <- function(sExpr) {
    counts <- assay(sExpr, 'counts')
    plot(density(counts[, 1]))
    sapply(colnames(counts[, 2:ncol(counts)]),
        function(colName) {
            lines(density(counts[, colName]))
    })
}