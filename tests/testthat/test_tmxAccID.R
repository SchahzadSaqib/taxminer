library(taxminer)
library(testthat)
library(stringr)

### checks of the function retrieval accuracy
testthat::test_that("TaxIDs are successfully retrieved", {
    outfile <- "Fungi_noEnv.txt"
    if (file.exists(outfile)) {
    #Delete file if it exists
    file.remove(outfile)
  }
  txm_accIDs(text_query = "Fungi [organism] AND vagina",
             db2src = "nucleotide",
             out_name = outfile)
  results <- read.csv(outfile, header = FALSE)
  expect_true(nrow(results) > 0)
})

#### error handling checks
testthat::test_that("error handling when wrong database name is given", {
  
  expect_error(txm_accIDs(text_query = "Fungi [organism] AND vagina",
                          db2src = "NotADatabase",
                          out_name = outfile))
})

