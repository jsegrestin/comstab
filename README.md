comstab
================
[![R-CMD-check](https://github.com/jsegrestin/comstab/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jsegrestin/comstab/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/608223971.svg)](https://zenodo.org/doi/10.5281/zenodo.10687005)

`comstab` is an R package that contains basic functions to apply the
unified framework for partitioning the drivers of stability of
ecological communities (Segrestin <i>et al.</i> 2024).

## Installation

You can install the stable version from CRAN with:

``` r
install.packages("comstab")
```

Alternatively, you can install the development version with:

``` r
require(devtools)
devtools::install_github("jsegrestin/comstab")
```

## Dependencies

Apart from base, `comstab` depends on `stats`, `graphics`, and `Ternary`
