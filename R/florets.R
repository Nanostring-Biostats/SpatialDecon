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

#' Draw coxcomb plots as points in a graphics window
#'
#' Draws a scatterplot where each point is a circular barplot, intended to show
#' decon results
#'
#' @param x Vector of x coordinates
#' @param y Vector of y coordinates
#' @param b matrix or cell abundances, with columns aligned with the elements
#' of x and y
#' @param col vector of colors, aligned to the rows of b.
#' @param legendwindow Logical. If TRUE, the function draws a color legend in a
#'  new window
#' @param rescale.by.sqrt Logical, for whether to rescale b by its square root
#' to make value proportional to
#'  shape area, not shape length.
#' @param border Color of pie segment border, defauls to NA/none
#' @param add Logical. If TRUE, the function draws florets atop an existing
#' graphics device (TRUE) or call a new device (FALSE).
#' @param cex Floret size. Florets are scaled relative to the range of x and y;
#' this further scales up or down.
#' @param bty bty argument passed to plot()
#' @param xaxt xaxt argument passed to plot()
#' @param yaxt yaxt argument passed to plot()
#' @param xlab xlab, defaults to ""
#' @param ylab ylab, defaults to ""
#' @param ... additional arguments passed to plot()
#' @return Draws a coxcomb plot, returns no data.
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
#' # draw florets:
#' florets(
#'   x = mini_geomx_dataset$annot$x,
#'   y = mini_geomx_dataset$annot$y,
#'   b = res0$beta, cex = 2
#' )
#' @export
florets <- function(x, y, b, col = NULL, legendwindow = FALSE,
                    rescale.by.sqrt = TRUE, border = NA, add = FALSE, cex = 1,
                    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                    ...) {
    # rescale b by sqrt so magnitude is proportional to coxcomb area, not length
    if (rescale.by.sqrt) {
        b <- sqrt(b)
    }
    # make b a matrix:
    if (is.vector(b)) {
        b2 <- matrix(b, nrow = length(b))
        rownames(b2) <- names(b)
        b <- b2
        rm(b2)
    }
    # choose colors if not given:
    if ((length(col) == 0) &
        all(is.element(rownames(b), names(SpatialDecon::cellcols)))) {
        col <- SpatialDecon::cellcols[rownames(b)]
    }
    if ((length(col) == 0) &
        !all(is.element(rownames(b), names(SpatialDecon::cellcols)))) {
        manycols <- c(
            "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
            "#B3DE69", "#FCCDE5", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
            "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1B9E77", "#D95F02",
            "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
            sample(grDevices::colors(), 99)
        )
        col <- manycols[seq_len(nrow(b))]
        names(col) <- rownames(b)
    }
    # convert colors to matrix of the same dimension as b:
    if (length(col) == 1) {
        col <- matrix(col, nrow = nrow(b), ncol = ncol(b))
    }
    if (is.vector(col)) {
        col <- matrix(col, nrow = nrow(b), ncol = ncol(b))
    }

    # get radians:
    angles <- seq(0, 2 * pi, length.out = nrow(b) + 1)

    # scale b based on the range of x and y:
    maxrange <- max(diff(range(x, na.rm = TRUE)), diff(range(y, na.rm = TRUE)))
    b <- b * maxrange / mean(b, na.rm = TRUE) * 0.007 * cex

    # draw plot:
    if (!add) {
        graphics::plot(x, y,
            col = 0, bty = bty, xaxt = xaxt, yaxt = yaxt,
            xlab = xlab, ylab = ylab, ...
        )
    }

    # draw florets:
    if (nrow(b) > 1) {
        for (i in seq_len(length(x))) {
            for (j in seq_len(nrow(b))) {
                tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
                xt <- b[j, i] * cos(tempangles)
                yt <- b[j, i] * sin(tempangles)
                graphics::polygon(x[i] + c(0, xt), y[i] + c(0, yt),
                    col = col[j, i],
                    border = border, lwd = 0.5
                )
            }
        }
    }

    # if just one point, draw a full circle:
    if (nrow(b) == 1) {
        for (i in seq_len(length(x))) {
            for (j in seq_len(nrow(b))) {
                tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
                xt <- b[j, i] * cos(tempangles)
                yt <- b[j, i] * sin(tempangles)
                graphics::polygon(x[i] + xt, y[i] + yt,
                    col = col[j],
                    border = border, lwd = 0.5
                )
            }
        }
    }

    # draw a legend:
    if (legendwindow) {
        graphics::plot(0, 0,
            col = 0, xlim = c(-1, 1), ylim = c(-1, 1), xaxt = "n",
            yaxt = "n", xlab = "", ylab = "", ...
        )
        for (j in seq_len(length(angles))) {
            graphics::lines(c(0, 0.75 * cos(angles[j])), c(0, 0.75 * sin(angles[j])),
                col = col[j], lwd = 2
            )
            graphics::text(0.85 * cos(angles[j]), 0.85 * sin(angles[j]),
                rownames(b)[j],
                cex = 1.4
            )
        }
    }
}
