
<!-- README.md is generated from README.Rmd. Please edit that file -->

# taxminer

<!-- badges: start -->
<!-- badges: end -->

The goal of taxminer is to …

## Installation

You can install the released version of taxminer from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("taxminer")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(taxminer)
## extracting accession numbers
taxminer::txm_accIDs(text_query = "Fungi [organism] AND vagina NOT (environmental samples [organism] OR project)",
                    db2src = "nucleotide",
                    out_name = "demo/Fungi_noEnv.seq")
#> Done! - Found 310 IDs
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
test <- readr::read_delim("demo/Fungi_noEnv.seq", delim = "\n", col_names = "AccIDs", col_types = readr::cols(AccIDs = readr::col_character()))
test
#> # A tibble: 620 x 1
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
#> # ... with 610 more rows
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
