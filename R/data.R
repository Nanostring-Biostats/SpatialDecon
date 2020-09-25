#' Default colors for the cell types in the safeTME matrix
#'
#' A named vector of colors, giving colors for the cell types of the safeTME
#'  matrix.
#'
#' @format A named vector
"cellcols"


#' Default colors for the cell types in the safeTME matrix
#'
#' A named vector of colors, giving colors for the cell types of the safeTME
#'  matrix.
#'
#' @format A named vector
"mini_geomx_dataset"


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
