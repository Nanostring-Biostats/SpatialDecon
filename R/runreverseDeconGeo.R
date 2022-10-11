#' Run Reversedecon
#' 
#' Runs reversedecon from an S4 object
#' @param object An S4 object such as a GeoMxSet object
#' @param ... Arguments passed to reversedecon
#' 
#'

setGeneric("runReverseDecon", signature = "object",
           function(object, ...) standardGeneric("runReverseDecon"))

#' Run reversedecon on a NanostringGeomxSet object
#'
#' A wrapper for applying reversedecon to a NanostringGeomxSet object.
#'
#'
#' @param norm_elt normalized data element in assayData.
#' @param beta Matrix of cell abundance estimates, with cells in columns and observations in rows.
#'  Columns are aligned to "norm".
#' @param epsilon All y and yhat values are thresholded up to this point when performing decon.
#'  Essentially says, "ignore variability in counts below this threshold."
#' @return a valid GeoMx S4 object including the following items:
#' \itemize{
#'  \item in fData
#'  \itemize{
#'  \item coefs, a matrix of coefficients for genes * cells, where element i,j is interpreted as
#'               "every unit increase in cell score j is expected to increase expression of gene i by _".
#'  \item cors, a vector giving each gene's correlation between fitted and observed expression
#'  \item resid.sd, a vector of each gene's residual SD, a metric of how much variability genes
#'                  have independend of cell mixing.
#'  }
#'  \item in assayData
#'  \itemize{
#'    \item yhat, a matrix of fitted values, in the same dimension as norm
#'    \item resids, a matrix of log2-scale residuals from the reverse decon fit, in the same
#'                dimension as norm
#'  }
#'  
#' }
#' @importFrom Biobase assayData
#' @examples
#' library(GeomxTools)
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' demoData <- shiftCountsOne(demoData)
#' target_demoData <- aggregateCounts(demoData)
#' 
#' target_demoData <- normalize(target_demoData, "quant")
#'                 
#' # run basic decon:
#' res0 <- runspatialdecon(object = target_demoData,
#'                         norm_elt = "exprs_norm",
#'                         raw_elt = "exprs")
#'
#' # run reverse decon:
#' target_demoData <- runReverseDecon(object = target_demoData,
#'                                    norm_elt = "exprs_norm",
#'                                    beta = pData(res0)$beta)
#'
#' @export
#' @rdname runReverseDecon
setMethod("runReverseDecon",  "NanoStringGeoMxSet",
          function(object, norm_elt = NULL, beta, epsilon = NULL){
              
              if(is.null(norm_elt)){
                  stop("norm_elt must be set")
              }
              
              if (!is.element(norm_elt, names(object@assayData))) {
                  stop(paste(norm_elt, "is not an element in assaysData slot"))
              }
              
              # prep components:
              norm <- as.matrix(Biobase::assayDataElement(object, elt = norm_elt))
              
              # run runReverseDeconGeomx:
              result <- reverseDecon(norm = norm,
                                     beta = t(beta),
                                     epsilon = epsilon)
              
              
              # append results to the object
              # add coefs
              Biobase::fData(object)[["coefs"]] <- matrix(NA, nrow = nrow(object), ncol = ncol(result$coefs),
                                                          dimnames = list(Biobase::featureNames(object), colnames(result$coefs)))
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
              Biobase::fData(object)[["cors"]][match(names(result$cors), Biobase::featureNames(object), nomatch = 0)] <- result$cors
              
              
              # add resid.sd
              Biobase::fData(object)[["resid.sd"]] <- NA
              Biobase::fData(object)[["resid.sd"]][match(names(result$resid.sd), Biobase::featureNames(object), nomatch = 0)] <- result$resid.sd
              
              return(object)
              
          })

