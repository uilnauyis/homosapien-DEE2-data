test_that("When 'bulkDownloading' set to FALSE, 'DownloadDEE2Data' function 
          could get the same data for all given SRR accessions as when downloaded 
          for each accessions seperately using DEE2 legacy download", {
            SRRvec = c("SRR1783836", "SRR1783837", "SRR1783838", "SRR1999221", 
                       "SRR2153338", "SRR2153409", "SRR2153289")
            
            ## Disable bulk download option
            se <- DownloadDEE2Data('hsapiens', SRRvec, bulkDownloading = FALSE)
            
            ## Download the DEE2 data of sample 'SRR2153289' using DEE2 legacy
            ## download.
            listLegacy_SRR1783836 <- getDEE2(species = "hsapiens", 
                                             SRRvec = c("SRR1783836"), 
                                             legacy = TRUE)
            
            listLegacy_SRR1783837 <- getDEE2(species = "hsapiens", 
                                             SRRvec = c("SRR1783837"), 
                                             legacy = TRUE)
            
            listLegacy_SRR2153289 <- getDEE2(species = "hsapiens", 
                                             SRRvec = c("SRR2153289"), 
                                             legacy = TRUE)
            
            
            ## Check if counts equal for the first given SRR accession
            expect_equal(sum(listLegacy_SRR1783836$GeneCounts[, c('SRR1783836')]),
                         sum(assay(se)[, c('SRR1783836')]))
            
            ## Check if counts equal for the second given SRR accession
            expect_equal(sum(listLegacy_SRR1783837$GeneCounts[, c('SRR1783837')]),
                         sum(assay(se)[, c('SRR1783837')]))
            
            ## Check if counts equal for the last given SRR accession
            expect_equal(sum(listLegacy_SRR2153289$GeneCounts[, c('SRR2153289')]),
                         sum(assay(se)[, c('SRR2153289')]))
          })