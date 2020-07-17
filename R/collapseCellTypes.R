#SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression data
#Copyright (C) 2020, NanoString Technologies, Inc.
#    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.
#Contact us:
#NanoString Technologies, Inc.
#530 Fairview Avenue N
#Seattle, WA 98109
#Tel: (888) 358-6266
#pdanaher@nanostring.com


#' Collapse related cell types within a deconvolution result
#' 
#' Given the input of an SpatialDecon result output and a list of which cell types to combine,
#'  returns a reshaped deconvolution result object with the specified cell types merged.
#' @param fit The object (a list) returned by the SpatialDecon algorithm
#' @param matches A named vector specifying the cell types to be combined. 
#'  Each element of the vector is a target/ rolled-up cell type; element names give the original 
#'  cell types. 
#'