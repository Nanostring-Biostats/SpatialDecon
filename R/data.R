#' Default colors for the cell types in the safeTME matrix
#'
#' A named vector of colors, giving colors for the cell types of the safeTME
#'  matrix.
#'
#' @format A named vector
"cellcols"

#' Small example GeoMx data 
#'
#' A miniature GeoMx dataset used by the spatialdecon examples. 
#'
#' @format A list with the following elements:
#'  \itemize{
#'  \item normalized: normalized data matrix
#'  \item raw: raw data matrix
#'  \item annot: AOI annotation data frame
#'  }
"mini_geomx_dataset"

#' Large example GeoMx data 
#'
#' A GeoMx dataset with dense AOIs gridded over a NSCLC tumor. Each AOI is split into tumor and microenvironment segments.
#'
#' @format GeoMxSet Object
"nsclc"



#' Mapping from granularly-defined cell populations to broaded cell populations
#'
#' Mapping from granularly-defined cell populations to broaded cell populations,
#'  for use by the convertCellTypes function.
#'
#' @format A list. Each element of the list contains the granular cell types
#'  that roll up
#'  to a single coarse cell type.
"safeTME.matches"


#' SafeTME matrix
#'
#' A matrix of expression profiles of 906 genes over 18 cell types.
#'
#' @format A matrix with 906 genes (rows) and 18 cell types (columns)
"safeTME"


#' Genes' biological variability in immune deconvolution from TCGA.
#'
#' Genes' biological SDs, as estimated from immune deconvolution from TCGA.
#' Used to weight genes in spatialdecon.
#'
#' @format A named vector giving SDs of 1179 genes.
"mean.resid.sd"

#' Mini human colon single cell dataset
#'
#' Random 250 cells and most informative genes (CV > 10) between cell types from
#' Kinchen, J. et al. Structural Remodeling of the Human Colonic 
#' Mesenchyme in Inflammatory Bowel Disease. Cell 175, 372-386.e17 (2018).
#'
#' @format A list with the following elements:
#'  \itemize{
#'  \item mtx: sparse count matrix
#'  \item annots: cell type annotation data frame
#'  }
"mini_singleCell_dataset"

