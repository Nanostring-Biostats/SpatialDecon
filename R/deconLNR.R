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
#' @keywords internal
#' @noRd
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
        fit <- lognlm3(pmax(y, epsilon) ~ b + Xtemp - 1,
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



# lognlm function from logNormReg v3.0
lognlm3 = function (formula, data, subset, weights, na.action, y = TRUE, 
                    start, model = TRUE, lik = TRUE, opt = c("nlminb", "optim"), 
                    contrasts = NULL, ...) 
{
    opt <- match.arg(opt)
    call <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    intercMt <- attr(mt, "intercept")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else stop("error in the design matrix")
    p <- ncol(X)
    if (!missing(start)) {
        if (lik && (length(start) != (p + 1))) 
            stop("if 'lik=TRUE', length(start) should have ncol(X)+1 values ")
        if (!lik && (length(start) != (p))) 
            stop("if 'lik=FALSE',  length(start) should have ncol(X) values  ")
    }
    attrContr <- attr(X, "contrasts")
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                          length(offset), NROW(Y)), domain = NA)
    }
    n <- nrow(X)
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (missing(start)) {
        b <- drop(solve(crossprod(X), crossprod(X, Y)))
        if (colnames(X)[1] == "(Intercept)") 
            b[1] <- max(b[1], median(Y))
        s <- log(mad(Y))
    }
    else {
        b <- start[1:p]
        s <- start[p + 1]
    }
    par0 <- if (lik) 
        c(b, s)
    else b
    obj <- lognlm.fit(X = X, y = Y, par = par0, lik = lik, opt = opt, 
                      offset = offset, ...)
    names(obj$coefficients) <- colnames(X)
    names(obj$s2) <- ""
    obj$call <- call
    if (y) 
        obj$y <- Y
    class(obj) <- "lognlm"
    obj$na.action <- attr(mf, "na.action")
    obj$offset <- offset
    obj$contrasts <- attr(X, "contrasts")
    obj$xlevels <- .getXlevels(mt, mf)
    obj$terms <- mt
    obj$opt <- opt
    obj$lik <- lik
    if (model) 
        obj$model <- mf
    obj
}