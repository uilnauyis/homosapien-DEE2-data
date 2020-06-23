################################################################################
## FUNCTION: .filterData
################################################################################
#' .filterData function
#' 
#' @param sExpr raw sExpr data object 
#' @param counts.cutoff the threshold for filtering genes. For a specific gene,
#'    the total cound 
#' @return sExpr 
.filterData <- function(sExpr, counts.cutoff = 10) {
  # remove genes smaller than the cutoff
  sExpr <- sExpr[rowSums(assay(sExpr)) > counts.cutoff, ]

  sExpr
}

################################################################################
## FUNCTION: .downloadAbundanceTsv
################################################################################
#' .downloadAbundanceTsv function
#'
#' @param species
#' @param srrAccessions
#' @param outDir 
#' @return abundancefiles 
.downloadAbundanceTsv <- function(species, srrAccessions, outDir = NULL) {
    print('Downloading abundance files from DEE2 web API for each SRR accessions,
        This step could take a very long time if the number of SRR accessions is
        huge.')
    abundancefiles <- lapply(srrAccessions, function(srrAccession) {
        abundanceUrl <- file.path('http://dee2.io/data', 
            species, srrAccession, paste(srrAccession, '.ke.tsv.gz', sep = ""))
        abundanceName <- NULL
        if(is.null(outDir)){
            outDir <- tempfile()
            abundanceName <- paste(outDir, '.abundance.tsv.gz', sep = "")
        } else {
            abundanceName <- file.path(outDir, 
                paste(srrAccession, '.abundance.tsv.gz', sep = ""))
        }
        download.file(abundanceUrl, destfile=abundanceName)
        abundanceDf <- read.table(gzfile(abundanceName)) 
        abundanceDf <- abundanceDf[2:nrow(abundanceDf), ]
        colnames(abundanceDf) <- c('target_id', 'length', 'eff_length', 
            'est_counts', 'tpm')
        write.table(abundanceDf, 
            file=gzfile(abundanceName), 
            row.names=FALSE, 
            quote=FALSE, 
            sep='\t')
        abundanceName
    })
    abundancefiles <- file.path(abundancefiles)
    names(abundancefiles) <- srrAccessions
    abundancefiles
}
