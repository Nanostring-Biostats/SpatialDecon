
## Overview

The SpatialDecon library implements the SpatialDecon algorithm for mixed cell deconvolution in spatial gene expression datasets. (This algorithm also works in bulk expression profiling data.)

Details can be found in the SpatialDecon manuscript: Danaher & Kim (2020), "Advances in mixed cell deconvolution enable quantification of cell types in spatially-resolved gene expression data.""


## Guide to functions

#### Data preparation functions:

* "download_profile_matrix" Downloads any one of ~75 cell profile matrices compiled for deconvolution of diverse tissue types. 
* "create_profile_matrix" Creates custom profile matrix from single cell data: count matrix and cell type annotations.  
* "derive_GeoMx_background_at_normalized_scale" Estimates the background expected from each data point in a GeoMx dataset. Accurate background estimation is key for SpatialDecon's accuracy. 

#### SpatialDecon: the core deconvolution function:

spatialdecon runs the SpatialDecon algorithm for estimating mixed cell type abundance in the regions
 of spatially-resolved gene expression studies. 
 It is also appropriate to apply to bulk gene expression data.

Its minimal required input is:

* A normalized data matrix
* A matrix of expected background counts at each element of the normalized data matrix 
* A matrix of expected cell type expression profiles


spatialdecon has numerous advanced options specified in the help file and demonstrated in the vignette.
These include:

* Merge closely-related cell types
* Estimate cell abundance on the scale of absolute cell counts
* Use regions of pure tumor cells (or any other cell type missing from the cell profile matrix) to infer a missing cell type's profile and estimate it alongside the cells with known profiles.


#### Plotting functions:

* "TIL_barplot" is a convenient way to draw barplots of cell type abundance/ proportion
* "florets" is for plotting cells in space. For each data point specified in xy space, it draws a circular barplot showing the localized abundance of cell types. 

#### Post-deconvolution analyses:

* "reverseDecon" Models genes ~ decon-derived cell scores. It produces fitted values, residuals, and various metrics of how much genes depend on cell mixing. In studies where cell mixing is a dominant source of variance, these residuals aid interpretation. 


## Installation

``` r
devtools::install_github("Nanostring-Biostats/SpatialDecon",
                         ref = "master", 
                         build_vignettes = FALSE)
```

## Getting started

``` r
library(SpatialDecon)
```


See the package's vignette for an example of its use. 


## Other

* Depends: R >= 4.0, Mac, Unix or Windows
* Typical installation time: seconds
* Expected run time for vignette: 1 minute
* For reproducible code of all analyses in the SpatialDecon manuscript, see https://github.com/Nanostring-Biostats/SpatialDecon-manuscript-analyses
* license: GPL-3

## Computational benchmarking

The below memory usage and runtimes were gathered from applying spatialdecon to increasing numbers of GeoMx AOIs, using a 544-gene x 18 cell-type cell profile matrix:

| n  | memory (MB) | runtime (ms)  |
|---|---|---|
| 10  | 324  | 500  |
| 50  |  1598 | 2190  |
| 100  | 3130  | 4280  |
| 200  | 6266  | 8610  |
| 500  | 15244  | 21610  |
| 1000  | 30451  | 42580  |
| 5000  | 61522  | 83720  |
