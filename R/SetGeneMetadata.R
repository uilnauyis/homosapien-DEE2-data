#' SetGeneMetadata
#'
#' @param sExpr SummarizedExperiment
#' @param species character
#' @export 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay assays
#' @importFrom biomaRt useMart listDatasets useDataset getBM
SetGeneMetadata <- function(sExpr, species) {
  datasetName <- paste(species, '_gene_ensembl', sep='')
  
  ensembl <- useMart("ensembl")
  datasets <- listDatasets(ensembl)
  
  # the dataset that we retrieve gene metadata from is expected to exists as one
  # of the biomaRt datasets
  stopifnot(datasetName %in% datasets$dataset)
  
  ensembl <- useDataset(datasetName, mart=ensembl)
  
  geneMetaData <- getBM(
    attributes=c('ensembl_gene_id',
                 'hgnc_symbol',
                 'external_gene_name', 
                 'entrezgene_id'), 
        filters = 'ensembl_gene_id', 
        values = rownames(sExpr), 
        mart = ensembl)
  
  rowData(sExpr) <- geneMetaData
  
  sExpr
}