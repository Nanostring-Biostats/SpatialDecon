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


#' Identify outlier genes in a decon result
#'
#' Analyses a decon result's residuals to identify poorly-fit genes and flag
#' them for removal.
#'  Rule: flag anything with
#'
#' @param Y p-length expression vector or p * N expression matrix - the actual
#' (linear-scale) data
#' @param yhat Expectation of Y given the decon fit
#' @param resids Log2-scale residuals of Y vs. yhat
#' (log2(pmax(Y, 0.5)) - log2(yhat))
#' @param wts Matrix of data point weights, aligned to dimensions of resids
#' @param resid_thresh A scalar, sets a threshold on how extreme individual
#' data points' values
#'  can be (in log2 units) before getting flagged as outliers and set to NA.
#' @return a vector of names of poorly-fit genes
flagOutliers <- function(Y, yhat, resids, wts, resid_thresh = 3) {

    # get weighted resids:
    if (length(wts) == 0) {
        wres <- resids
    }
    if (length(wts) > 0) {
        wres <- resids * wts
    }

    # flag bad genes:
    outlier_genes <- c() # <-- this line makes it so no outlier genes are filtered
    # flag bad data points: (not doing anything for now)
    outlier_data_points <- abs(resids) > resid_thresh
    return(outlier_data_points)
}
