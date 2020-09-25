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

#' Reverse deconvolution
#'
#' Performs "reverse deconvolution", modelling each gene expression's ~
#' cell scores.
#' Returns a matrix of "fitted" expression values, a matrix of residuals,
#'  a matrix of
#' reverse decon coefficients for genes * cells.
#'
#' @param norm Matrix of normalized data, with genes in rows and observations
#' in columns
#' @param beta Matrix of cell abundance estimates, with cells in rows and
#' observations in columns.
#'  Columns are aligned to "norm".
#' @param epsilon All y and yhat values are thresholded up to this point when
#' performing decon.
#'  Essentially says, "ignore variability in counts below this threshold."
#' @return A list:
#' \itemize{
#' \item coefs, a matrix of coefficients for genes * cells, where element i,j is
#'  interpreted as
#' "every unit increase in cell score j is expected to increase expression of
#' gene i by _".
#' \item yhat, a matrix of fitted values, in the same dimension as norm
#' \item resids, a matrix of log2-scale residuals from the reverse decon fit,
#'  in the same
#'  dimension as norm
#' \item cors, a vector giving each gene's correlation between fitted and
#' observed expression
#' \item resid.sd, a vector of each gene's residual SD, a metric of how much
#' variability genes
#'  have independend of cell mixing.
#' }
#' @import logNormReg
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
#' # run reverse decon:
#' rdecon <- reverseDecon(
#'   norm = mini_geomx_dataset$norm,
#'   beta = res0$beta
#' )
#' @export
reverseDecon <- function(norm, beta, epsilon = NULL) {

    # remove cell types with no SD:
    beta <- beta[apply(beta, 1, stats::sd) > 0, ]
    # remove cell types that get removed by lm()
    # (presumably removed due to linear dependence)
    lm1 <- stats::lm(norm[1, ] ~ t(beta))
    beta <- beta[!is.na(lm1$coef[setdiff(names(lm1$coef), "(Intercept)")]), ,
        drop = FALSE
    ]

    # run reverse decon for all genes:
    rd <- function(y) {
        fit <- suppressWarnings(
            logNormReg::lognlm(y ~ t(beta),
                lik = FALSE,
                method = "L-BFGS-B",
                lower = rep(0, ncol(beta) + 1),
                upper = rep(Inf, ncol(beta) + 1),
                opt = "optim",
                control = list(maxit = 1000)
            )
        )
        return(fit$coefficients)
    }
    coefs <- t(apply(norm, 1, rd))
    colnames(coefs)[-1] <- rownames(beta)

    # get yhat
    yhat <- norm * NA
    for (ind in seq_len(ncol(yhat))) {
        yhat[, ind] <- coefs %*% c(1, beta[, ind])
    }

    # auto-select a reasonable epsilon if not provided
    if (length(epsilon) == 0) {
        epsilon <- stats::quantile(norm[norm > 0], 0.01)
    }

    # get resids:
    resids <- log2(pmax(norm, epsilon)) - log2(pmax(yhat, epsilon))

    # get summary stats:
    cors <- suppressWarnings(diag(stats::cor(t(norm), t(yhat))))
    resid.sd <- apply(resids, 1, stats::sd)

    out <- list(
        coefs = coefs, yhat = yhat, resids = resids, cors = cors,
        resid.sd = resid.sd
    )
    return(out)
}
