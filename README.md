
<!-- README.md is generated from README.Rmd. Please edit that file -->

# taxminer

<!-- badges: start -->
<!-- badges: end -->

Taxonomic annotations - BLAST alignment and text-mining based filtration
in R

## Installation

You can install the released version of taxminer from
[Github](https://github.com/) with:

``` r
devtools::install_github("taxminer")
```

## Example

``` r
library(taxminer)
dir.create("demo")
## extracting accession numbers
taxminer::txm_accIDs(text_query = "Fungi [organism] AND vagina",
                    db2src = "nucleotide",
                    out_name = "demo/Fungi_noEnv.seq")
#> Done! - Found 339 IDs
```

``` r
get_accIds <- readr::read_delim("demo/Fungi_noEnv.seq", delim = "\n", col_names = "AccIDs", col_types = readr::cols(AccIDs = readr::col_character()))
get_accIds
#> # A tibble: 339 x 1
#>    AccIDs    
#>    <chr>     
#>  1 LC601978.1
#>  2 MW559947.1
#>  3 CP048242.1
#>  4 CP048241.1
#>  5 CP048240.1
#>  6 CP048239.1
#>  7 CP048238.1
#>  8 CP048237.1
#>  9 CP048236.1
#> 10 CP048235.1
#> # ... with 329 more rows
```
