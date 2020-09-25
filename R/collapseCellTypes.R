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

#' Collapse related cell types within a deconvolution result
#'
#' Given the input of an SpatialDecon result output and a list of which cell
#' types to combine,
#'  returns a reshaped deconvolution result object with the specified cell
#'  types merged.
#' @param fit The object (a list) returned by the SpatialDecon algorithm
#' @param matching A list object holding the mapping from beta's cell names to
#' official cell names.
#'  See str(safeTME.matches) for an example.
#' @return A reshaped deconvolution result object
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
#' res1 <- collapseCellTypes(
#'     fit = res0,
#'     matching = SpatialDecon::safeTME.matches
#' )
#' @export
collapseCellTypes <- function(fit, matching) {

    # results object to hold the collapsed results:
    out <- fit

    # format matching list as a matrix to take a linear combination of beta:
    startingcellnames <- unlist(matching)
    A <- matrix(0, length(matching), nrow(fit$beta),
        dimnames = list(names(matching), rownames(fit$beta))
    )
    for (name in names(matching)) {
        cellnames <- matching[[name]]
        A[name, cellnames] <- 1
    }

    # apply A transformation to beta:
    for (name in c("beta", "prop_of_all", "prop_of_nontumor")) {
        if (is.element(name, names(fit))) {
            out[[name]] <- A[, startingcellnames] %*% fit[[name]][startingcellnames, ]
        }
    }

    # if Sigma provided, get vcov of beta2:
    if (is.element("sigmas", names(out))) {
        sigma <- fit$sigmas
        if (length(dim(sigma)) == 2) {
            out$sigmas <- A[, startingcellnames] %*%
                sigma[startingcellnames, startingcellnames, ] %*%
                t(A[, startingcellnames])
        }
        if (length(dim(sigma)) == 3) {
            out$sigmas <- array(NA,
                dim = c(nrow(A), nrow(A), dim(sigma)[3]),
                dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
            )
            for (i in seq_len(dim(sigma)[3])) {
                out$sigmas[, , i] <- A %*% sigma[, , i] %*% t(A)
            }
        }
    }

    # re-calculate p, se, t:
    if (is.element("beta", names(out)) & is.element("sigmas", names(out))) {
        # compute p-values
        tempbeta <- out$beta
        tempse <- tempp <- tempt <- tempbeta * NA
        for (i in seq_len(ncol(tempse))) {
            tempse[, i] <- suppressWarnings(sqrt(diag(out$sigmas[, , i])))
        }
        out$se <- tempse
        out$t <- (tempbeta / tempse)
        if (is.element("X", names(out))) {
            out$p <- 2 * (1 - stats::pt(out$t, df = nrow(fit$X) - ncol(fit$X) - 1))
        }
    }

    return(out)
}
