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



#' Barplot of abundance estimates
#' 
#' Draw barplot of the "betas" from a decon fit
#' 
#' @param mat Matrix of cell proportions or abundances, in the same dimensions output by spatialdecon
#'  (cells in rows, observations in columns). User is free to re-order columns/observations in 
#'  whatever order is best for display.
#' @param draw_legend Logical. If TRUE, the function draws a legend in a new plot frame.
#' @param main Title for barplot
#' @param ... Arguments passed to barplot()
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg = derive_GeoMx_background(
#'    norm = mini_geomx_dataset$normalized,
#'    probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'    negnames = "NegProbe")
#' # run basic decon:
#' res0 = spatialdecon(norm = mini_geomx_dataset$normalized, 
#'                     bg = mini_geomx_dataset$bg,
#'                     X = safeTME)
#' # run barplot:
#' TIL_barplot(mat = res0$beta)
#' # run barplot and draw a color legend
#' TIL_barplot(mat = res0$beta, draw_legend = TRUE)
#' @export
TIL_barplot = function(mat, draw_legend = FALSE, main = "", ...) {
  
  usecells = intersect(rownames(mat), names(SpatialDecon::cellcols))
  
  # draw barplot:
  graphics::barplot(mat[usecells, ], cex.lab = 1.5,
                    col = SpatialDecon::cellcols[rownames(mat)], border = NA,
                    las = 2, main = main, ...)
  
  # draw a legend:
  if (draw_legend) {
    graphics::frame()
    graphics::legend("center", 
                     fill = rev(SpatialDecon::cellcols[usecells]), 
                     legend = rev(usecells))
  }
  
}
