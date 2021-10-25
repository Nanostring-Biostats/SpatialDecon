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



# lognlm function from logNormReg v0.3
lognlm3 = function (formula, data, subset, weights, na.action, y = TRUE, 
                    start, model = TRUE, lik = TRUE, opt = c("nlminb", "optim"), 
                    contrasts = NULL, ...) {
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
    obj <- lognlm.fit3(X = X, y = Y, par = par0, lik = lik, opt = opt, 
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

# lognlm.fit function from logNormReg v0.3
lognlm.fit3 <-
    function(X, y, par, lik=TRUE, opt=c("nlminb","optim"), offset=NULL, ...){ #method=c("BFGS", "Nelder-Mead")
        #fitter function to estimate multiple linear regression with logNormal errors
        #X,y: design matrix and response
        #par: optional, starting values for 
        opt<-match.arg(opt)
        if(min(y)<=0) stop("Data can include only positive values")
        n<-length(y)
        p<-ncol(X)
        if(!missing(par)) {
            if(lik && length(par)!=(ncol(X)+1) ) stop("if 'par' is provided, it has to be length(par)= ncol(X)+1")
            if(!lik && length(par)!=ncol(X) ) stop("if 'par' is provided, it has to be length(par)= ncol(X)")
        }
        #if (!is.null(offset)) y <- y - offset
        if (is.null(offset)) offset<-0
        #y<-y-offset #Non funziona perche' poi le y possono essere negative..
        
        if(lik){
            
            mioNR <- function(start, X, y, h=.5, tol=1e-5,offs=0){
                par0<-start
                eps<-10
                it<-0
                while(eps>tol){
                    par1<- par0 - h* drop(solve(hess(par0,X,y,offs), lgrad(par0,X,y,offs)))
                    eps<- max(abs((par1-par0)/par0))
                    par0<-par1
                    it<-it+1
                    cat(it, par0, "\n")
                }
                par1
            }
            
            llik<-function(par,X,y,offs=0){
                b <- par[-length(par)]
                s <- par[length(par)]
                mu<- pmax(X %*% b,.0001)+offs
                mz<- log(mu)-s^2/2
                #r<- -sum(dlnorm(y,mz,s,log=TRUE))
                li<- -.5*log(s^2) -(log(y)-mz)^2/(2*s^2) #tolto "-log(y)" e poi dovrebbe essere "log(2*pi*s^2)".. non serve..
                r <- -sum(li)
                r
            }
            lgrad<-function(par,X,y,offs=0){
                b <- par[-length(par)]
                s <- par[length(par)]
                mu<- pmax(X %*% b, 0.0001)+offs
                mz<- log(mu)-s^2/2
                deriv.b<- t(X) %*% ((log(y)-mz)/mu)/s^2
                #questa  e' rispetto a s2:
                #deriv.s<- -n/(2*s^2) - (s^2*sum((log(y)-mz))-sum((log(y)-mz)^2))/(2*s^4) 
                #rispetto a s (Attenzione su qualche dataset usando la deriv s (e non s2) la stima di s viene negativa!!!)
                deriv.s<- -n/s - (s^2*sum((log(y)-mz))-sum((log(y)-mz)^2))/(s^3) 
                #--
                r<-c(-deriv.b, -deriv.s)
                #r<-grad(llik, c(b,s), X=X,y=y)
                r
            }
            hess<-function(par, X, y, offs=0){
                #hessiana della logLik 
                n<-length(y)
                b <- par[-length(par)]
                s <- par[length(par)]
                mu<- pmax(drop(X %*% b), 0.0001) + offs
                mz<- log(mu)-s^2/2
                w <- drop((mz - log(y) - 1)/mu^2)
                #i w possono essere negativi..
                #r <- crossprod(X*sqrt(w))/(s^2)
                H.bb<-t(X) %*%diag(w) %*%X/(s^2)
                H.bs<- - (2/s^3) * drop(t(X) %*% ((log(y) - log(mu))/mu))
                H.ss<-n/(s^2)- (3/s^4)*sum((log(y)-log(mu))^2)-n/4
                r<- -cbind(rbind(H.bb, H.bs), c(H.bs, H.ss))
                r
            }
            
            if(missing(par)) {
                #b<-solve(crossprod(X),crossprod(X,log(y)))
                b<-solve(crossprod(X),crossprod(X,y))
                s<- log(sqrt(sum((y-drop(X%*%b))^2)/n))
                #s<- sqrt(sum((log(y)-drop(X%*%b))^2)/n)
            } else {
                b<-par[1:p]
                s<-par[p+1]
            }
            opz<-list(...)
            #browser()
            if(opt=="optim"){
                opz$par<-c(b,s)
                opz$fn<- llik
                opz$gr<-lgrad
                opz$offs<- quote(offset)
                opz$X<-quote(X)
                opz$y<-quote(y)
                opz$hessian<-TRUE
                if(is.null(opz$method)) opz$method<-"BFGS"
                o<-do.call(optim, opz )
                ll<-o$value
                hh<-o$hessian
            } else {
                #browser()
                #      o<-mioNR(c(b,s),X,y)
                opz$start<- c(b, s)
                opz$objective<- quote(llik)
                opz$gradient<- quote(lgrad)
                opz$hessian<- quote(hess)
                opz$offs<- quote(offset)
                opz$X<-quote(X)
                opz$y<-quote(y)
                o<-do.call(nlminb, opz)
                ll<-o$objective
                hh<-hess(o$par,X,y)
            }
            gradd<-lgrad(o$par,X=X,y=y,offs=offset)
            est.b<-o$par[1:p]
            est.s2<-o$par[(p+1)]^2
            
        } else {
            ###### Min distance
            dmin<-function(par,X,y,offs=0){
                b <- par
                mu<- pmax(drop(X %*% b),.0001)+offs
                r<- sum((log(y)-log(mu))^2)
                r
            }
            lgrad<-function(par,X,y,offs=0){
                b <- par                      
                mu<- pmax(drop(X %*% b),.0001) +offs
                deriv.b<- -2*t(X)%*%((log(y)-log(mu))/mu)
                r<-deriv.b
                r
            }
            if(missing(par)) {
                #b<-solve(crossprod(X),crossprod(X,log(y)))
                b<-solve(crossprod(X),crossprod(X, y))
            } else {
                b<-par[1:p]
            }
            
            opz<-list(...)
            opz$par<-b
            opz$fn<- dmin
            opz$gr<-lgrad
            opz$offs<- quote(offset)
            opz$X<-X
            opz$y<-y
            opz$hessian<-TRUE
            if(is.null(opz$method)) opz$method<-"BFGS"
            o<-do.call(optim, opz )
            gradd<-lgrad(o$par, X=X, y=y)
            est.b<-o$par[1:p]
            est.s2<- sum((log(y) - log(drop(X %*% est.b)+offset))^2)/n
            ll<- -o$value
            hh<- o$hessian
        }
        
        if(o$convergence!=0) warning("Unsuccessful convergence") 
        
        est<-c(est.b, est.s2)
        fits<-drop(X%*%est.b)+offset
        r<-list(coefficients=est.b, loglik= -ll, s2=est.s2, fitted.values=fits, residuals=y-fits, 
                grad=drop(gradd) , hessian=hh, convergence=o$convergence)
        return(r)
    }

