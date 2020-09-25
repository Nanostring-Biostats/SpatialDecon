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



#' Deconvolution using logNormReg package to run linear mean model and log error
#'  model
#'
#' Calls lognlm() to optimize the model.
#'
#' @param Y p-length expression vector or p * N expression matrix - the actual
#'  (linear-scale) data
#' @param X p * K Training matrix
#' @param bg scalar or matrix of expected background counts per data point.
#' @param weights The same as the weights argument used by lm
#' @param epsilon optional,  a very small non-zero number to use as a lower
#' threshold to make fits well-behaved
#' @param maxit Maximum number of iterations. Default 1000.
#' @return a list: beta (estimate), sigmas (covariance matrix of estimate,
#' derived by inverting the hessian from lognlm)
#' @import logNormReg
deconLNR <- function(Y, X, bg = 0, weights = NULL, epsilon = NULL,
                     maxit = 1000) {
    if (length(weights) == 0) {
        weights <- replace(Y, TRUE, 1)
    }
    if (is.vector(Y)) {
        Y <- matrix(Y, nrow = length(Y))
    }
    if (length(bg) == 1) {
        bg <- matrix(bg, nrow(Y), ncol(Y))
    }
    # choose "epsilon": a very small non-zero number to make fits well-behaved
    if (length(epsilon) == 0) {
        epsilon <- min(replace(Y, (Y == 0) & !is.na(Y), NA), na.rm = TRUE)
    }

    # matrix-like data for apply:
    mat <- rbind(Y, bg, weights)
    # fn to apply:
    fn <- function(zz) {
        # break into y, b, w:
        y <- zz[seq_len((length(zz)) / 3)]
        b <- zz[seq_len((length(zz)) / 3) + (length(zz) / 3)]
        wts <- zz[seq_len((length(zz)) / 3) + (length(zz) / 3) * 2]

        # remove NA data:
        use <- !is.na(y)
        y <- y[use]
        b <- b[use]
        Xtemp <- X[use, , drop = FALSE]
        wts <- wts[use]

        init <- rep(mean(y) / (mean(X) * ncol(X)), ncol(X))
        names(init) <- colnames(X)

        # run lognlm:
        fit <- lognlm(pmax(y, epsilon) ~ b + Xtemp - 1,
            lik = FALSE,
            weights = wts,
            start = c(1, init),
            method = "L-BFGS-B",
            lower = c(1, rep(0, ncol(Xtemp))),
            upper = c(1, rep(Inf, ncol(Xtemp))),
            opt = "optim",
            control = list(maxit = maxit)
        )
        fnout <- list(
            beta = fit$coefficients[-1],
            sigma = solve(fit$hessian)[-1, -1]
        )
        return(fnout)
    }
    # apply across all observations:
    fnlist <- apply(mat, 2, fn)
    # extract beta and sigmas:
    getbeta <- function(zz) {
        return(zz$beta)
    }
    getsigma <- function(zz) {
        return(zz$sigma)
    }

    beta <- vapply(X = fnlist, FUN = getbeta, FUN.VALUE = numeric(ncol(X)))
    rownames(beta) <- colnames(X)
    sigmas <- array(vapply(
        X = fnlist,
        FUN = getsigma,
        FUN.VALUE = numeric(ncol(X)^2)
    ),
    dim = c(ncol(X), ncol(X), ncol(Y)),
    dimnames = list(colnames(X), colnames(X), colnames(Y))
    )

    out <- list(beta = pmax(beta, 0), sigmas = sigmas)
    return(out)
}
