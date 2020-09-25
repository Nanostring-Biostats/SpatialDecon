# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
# Copyright (C) 2020, NanoString Technologies, Inc.
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.
# 530 Fairview Avenue N
# Seattle, WA 98109
# Tel: (888) 358-6266
# pdanaher@nanostring.com


#' SpatialDecon: A package for computating the notorious bar statistic.
#'
#' The SpatialDecon package estimates mixed cell type abundance in the regions
#' of spatially-resolved gene
#' expression studies, using the method of Danaher & Kim (2020), "Advances in
#'  mixed cell deconvolution enable
#' quantification of cell types in spatially-resolved gene expression data."
#' It is also appropriate to apply to bulk gene expression data.
#'
#' @section functions:
#' Functions to help set up deconvolution:
#' \itemize{
#'  \item derive_GeoMx_background Estimates the background levels from GeoMx
#'  experiments
#'  \item collapseCellTypes reformats deconvolution results to merge
#'  closely-related cell types
#'  \item download_profile_matrix Downloads a cell profile matrix.
#'  \item safeTME: a data object, a matrix of immune cell profiles for use in
#'   tumor-immune deconvolution.
#' }
#' Deconvolution functions:
#' \itemize{
#'  \item spatialdecon runs the core deconvolution function
#'  \item reverseDecon runs a transposed/reverse deconvolution problem, fitting
#'  the data as a function of cell abundance estimates.
#'   Used to measure genes' dependency on cell mixing and to calculate gene
#'    residuals from cell mixing.
#' }
#' Plotting functions:
#' \itemize{
#'  \item florets Plot cell abundance on a specified x-y space, with each point
#'   a cockscomb plot showing the cell abundances of that region/sample.
#'  \item TIL_barplot Plot abundances of tumor infiltrating lymphocytes (TILs)
#'   estimated from the safeTME cell profile matrix
#' }
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # run decon with bells and whistles:
#' res <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME,
#'   cellmerges = safeTME.matches,
#'   cell_counts = mini_geomx_dataset$annot$nuclei,
#'   is_pure_tumor = mini_geomx_dataset$annot$AOI.name == "Tumor"
#' )
#' @docType package
#' @name SpatialDecon-package
NULL
