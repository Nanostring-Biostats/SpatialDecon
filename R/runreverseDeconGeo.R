#' Run reversedecon on a NanostringGeomxSet object
#'
#' A wrapper for applying reversedecon to the Spatial data element in a NanostringGeomxSet object.
#'
#' @param object A NanostringGeomxSet object.
#' @param norm_elt normalized data element in assayData.
#' @param beta Matrix of cell abundance estimates, with cells in columns and observations in rows.
#'  Columns are aligned to "norm".
#' @param epsilon All y and yhat values are thresholded up to this point when performing decon.
#'  Essentially says, "ignore variability in counts below this threshold."
#' @return A list:
#' \itemize{
#'  \item coefs, a matrix of coefficients for genes * cells, where element i,j is interpreted as
#'               "every unit increase in cell score j is expected to increase expression of gene i by _".
#'  \item yhat, a matrix of fitted values, in the same dimension as norm
#'  \item resids, a matrix of log2-scale residuals from the reverse decon fit, in the same
#'                dimension as norm
#'  \item cors, a vector giving each gene's correlation between fitted and observed expression
#'  \item resid.sd, a vector of each gene's residual SD, a metric of how much variability genes
#'                  have independend of cell mixing.
#' }
#' @import Biobase assayData
#' @examples
#' data(target_demoData)
#' # run basic decon:
#' res0 <- runSpatialdeconGeomx(object = target_demoData,
#'                       norm_elt = "neg_norm",
#'                       raw_elt = "exprs")
#'
#' # run reverse decon:
#' runReverseDeconGeomx(object = target_demoData,
#' norm_elt = "neg_norm",
#' beta = res0$beta)
#' )
#'
#' @export
#'


setGeneric("runReverseDeconGeomx", signature = "object",
           function(object, ...) standardGeneric("runReverseDeconGeomx"))

setMethod("runReverseDeconGeomx",  "NanoStringGeoMxSet",
          function(object, norm_elt, beta, epsilon = NULL){

            # prep components:
            norm <- as.matrix(Biobase::assayDataElement(object, elt = norm_elt))

            # run runReverseDeconGeomx:
            result <- runReverseDeconGeomx(object = norm,
                                      beta = t(beta),
                                      epsilon = epsilon)


            # append results to the object
            # add coefs
            Biobase::fData(object)[["coefs"]] <- matrix(NA, nrow = nrow(object), ncol = ncol(result$coefs),
                                                       dimnames = list(featureNames(object), colnames(result$coefs)))
            Biobase::fData(object)[["coefs"]][rownames(result$coefs), ] <- result$coefs


            # add yhat matrix
            yhat_mat <- matrix(NA, nrow = nrow(object), ncol = ncol(object),
                               dimnames = dimnames(object))
            yhat_mat[rownames(result$yhat), colnames(result$yhat)] <- result$yhat
            Biobase::assayDataElement(object, "yhat") <- yhat_mat


            # add resids matrix
            resids_mat <- matrix(NA, nrow = nrow(object), ncol = ncol(object),
                                 dimnames = dimnames(object))
            resids_mat[rownames(result$resids), colnames(result$resids)] <- result$resids
            Biobase::assayDataElement(object, "resids") <- resids_mat


            # add cors
            Biobase::fData(object)[["cors"]] <- NA
            Biobase::fData(object)[["cors"]][match(names(result$cors), featureNames(object), nomatch = 0)] <- result$cors


            # add resid.sd
            Biobase::fData(object)[["resid.sd"]] <- NA
            Biobase::fData(object)[["resid.sd"]][match(names(result$resid.sd), featureNames(object), nomatch = 0)] <- result$resid.sd

            return(object)

            })


setMethod("runReverseDeconGeomx",  "matrix",
          function(object, beta, epsilon = NULL){

          # run reverseDecon
          res <- reverseDecon(norm = object,
                                beta = beta,
                                epsilon = epsilon)

          return(res)

          })
