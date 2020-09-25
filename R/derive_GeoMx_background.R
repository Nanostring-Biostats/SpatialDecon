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


#' Derive background at the scale of the normalized data for GeoMx data
#'
#' Estimates per-datapoint background levels from a GeoMx experiment.
#' In studies with two or more probe pools, different probes will have different
#' background levels. This function provides a convenient way to account for
#' this phenomenon.
#'
#' @param norm Matrix of normalized data, genes in rows and segments in columns.
#'  Must include negprobes, and must have rownames.
#' @param probepool Vector of probe pool names for each gene, aligned to the
#' rows of "norm".
#' @param negnames Names of all negProbes in the dataset. Must be at least one
#'  neg.name within each probe pool.
#' @return A matrix of expected background values, in the same scale and
#'  dimensions as the "norm" argument.
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' @export
derive_GeoMx_background <- function(norm, probepool, negnames) {

    # check data input:
    if (nrow(norm) != length(probepool)) {
        stop("nrow(norm) != length(probepool)")
    }

    # initialize:
    bg <- norm * 0


    # fill in expected background at scale of normalized data:
    for (pool in unique(probepool)) {

        # get the pool's negProbes:
        tempnegs <- intersect(negnames, rownames(norm)[probepool == pool])
        if (length(tempnegs) == 0) {
            stop(paste0(pool, " probe pool didn't have any negprobes specified"))
        }
        tempnegfactor <- colMeans(norm[tempnegs, , drop = FALSE])

        # fill in the corresponding elements of bg:
        bg[probepool == pool, ] <-
            sweep(bg[probepool == pool, ], 2, tempnegfactor, "+")
    }
    return(bg)
}
