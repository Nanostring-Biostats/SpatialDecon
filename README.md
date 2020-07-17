---
output: github_document
---


## Overview

The SpatialDecon library implements the SpatialDecon algorithm for mixed cell deconvolution in spatial gene expression datasets. (This algorithm also works in bulk expression profiling data.)
Details can be found in the SpatialDecon manuscript: ____.



## Guide to functions

#### Data preparation functions:

* "download_profile_matrix" Downloads any one of ~30 cell profile matrices compiled for deconvolution of diverse tissue types. 
* "derive_GeoMx_background_at_normalized_scale" Estimates the background expected from each data point in a GeoMx dataset. Accurate background estimation is key for SpatialDecon's accuracy. 

#### Deconvolution functions:

* "spatialdecon" runs the basic SpatialDecon algorithm. 
* it has numerous advanced options specified in the help file and demonstrated in the vignette.

#### Plotting functions:

* "TIL_barplot" is a convenient way to draw barplots of cell type abundance/ proportion
* "florets" is for plotting cells in space. For each data point specified in xy space, it draws a circular barplot showing the localized abundance of cell types. 

#### Post-deconvolution analyses:

* "reverseDecon" Models genes ~ decon-derived cell scores. It produces fitted values, residuals, and various metrics of how much genes depend on cell mixing. In studies where cell mixing is a dominant source of variance, these residuals aid interpretation. 


## Installation

``` r
install.packages("SpatialDecon")
```

## Getting started

``` r
library(SpatialDecon)
```


See the package's vignette for an example of its use. 
