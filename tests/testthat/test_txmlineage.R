context("tests for the function txm_lineage")
library(taxminer)
library(testthat)
library(stringr)

### checks of the function retrieval accuracy


test_that("a correct number of lineages are retrived from NCBI", {
  TaxId <- c(555079, 177437, 702113, 999552, 1547445, 110662, 585423, 349106, 351745, 
             186497, 314225, 314260, 1198449, 1248727, 1501269, 167879, 1225337, 375451, 388467, 216432) 
  my_testData <- data.frame(TaxId)
  
  results <- txm_lineage(my_testData, bindtoAcc = F)
  
  expect_equal(nrow(results), 20)
  expect_equal(ncol(results), 7)
})

test_that("a correct number of lineages are retrived from NCBI", {
  TaxId <- read.csv("marineOrgs-02-18.log") 
  my_testData <- data.frame(TaxId)
  
  results <- txm_lineage(my_testData, bindtoAcc = F)
  
  expect_equal(nrow(results), 4548)
  expect_equal(ncol(results), 7)
})


#### error handling checks

test_that("error handling when no taxID is given", {
  my_testData <-data.frame(TaxId=integer())
  
  expect_error(txm_lineage(my_testData, bindtoAcc = F), 
               "No TaxID provided", 
               ignore.case = TRUE)
})



