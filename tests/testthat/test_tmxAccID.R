context("tests for the function txm_accIDs")
library(taxminer)
library(testthat)
library(stringr)

### checks of the function retrieval accuracy
test_that("a correct number of TaxID is retrieved from NCBI nr database", {
    outfile <- "Fungi_noEnv.txt"
    if (file.exists(outfile)) {
    #Delete file if it exists
    file.remove(outfile)
  }
  txm_accIDs(text_query = "Fungi [organism] AND vagina",
             db2src = "nucleotide",
             out_name = outfile)
  results <- read.csv(outfile, header = FALSE)
  
  expect_equal(nrow(results), 339)
})


test_that("a correct number of TaxID is retrieved from NCBI nr database", {
  outfile <- "Fungi_noEnv.txt"
  if (file.exists(outfile)) {
    #Delete file if it exists
    file.remove(outfile)
  }
  txm_accIDs(text_query = "this shouldn't retrieve anything",
             db2src = "nucleotide",
             out_name = outfile)
  results <- read.csv(outfile, header = FALSE)
  
  expect_equal(nrow(results), 0)
})


#### error handling checks
test_that("error handling when wrong database name is given", {
  
  expect_error(txm_accIDs(text_query = "Fungi [organism] AND vagina",
                          db2src = "NotADatabase",
                          out_name = outfile), 
               "Database not found", 
               ignore.case = TRUE)
})

test_that("error handling when no file name is given", {
  expect_error(txm_accIDs(text_query = "Fungi [organism] AND vagina",
                          db2src = "nucleotide"), 
               "Output file name missing",
               ignore.case = TRUE)
})

