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



#' Compute gene weights
#'
#' Compute gene weights from pre-defined biological noise estimates and from
#'  raw count/ error-model-defined technical noise estimates
#' @param norm Matrix of data used in decon. Used to define the gene list and
#' the
#' shape of the output "wts" matrix.
#' @param raw Matrix of raw data. If provided, used to define technical noise
#' @param error.model Which error model to use. Defaults to "dsp"
#' @param weight.by.TIL.resid.sd If TRUE, then genes are weighted in part based
#' on their
#'  biological variability as estimated by their residual SD from decon
#'   performed on TCGA.
#' @return A matrix of weights, in the same dimension as norm
deriveWeights <- function(norm, raw = NULL, error.model = "dsp",
                          weight.by.TIL.resid.sd = FALSE) {

    # get tech SDs if raw data provided:
    if (length(raw) == 0) {
        sds.tech <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
    }
    if (length(raw) > 0) {
        sds.tech <- runErrorModel(
            counts = raw,
            platform = "dsp"
        )
    }

    # if the mean.resid.sd vector (which defines genes' biological SD) is in
    # the environment, get biological noise:
    if (!weight.by.TIL.resid.sd) {
        sds.bio <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
    }
    if (weight.by.TIL.resid.sd) {
        sds.bio <- matrix(NA, nrow(raw), ncol(raw), dimnames = dimnames(raw))
        for (gene in intersect(
            names(SpatialDecon::mean.resid.sd),
            rownames(sds.bio)
        )) {
            sds.bio[gene, ] <- SpatialDecon::mean.resid.sd[gene]
        }
        sds.bio <- replace(sds.bio, is.na(sds.bio), mean(sds.bio, na.rm = TRUE))
    }

    # define total SD, and invert to get weights
    sds.tot <- sqrt(sds.tech^2 + sds.bio^2)
    wts <- 1 / sds.tech
    return(wts)
}
