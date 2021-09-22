#' Run mergeTumorIntoX on a NanostringGeomxSet object
#'
#' A wrapper for applying mergeTumorIntoX to the Spatial data element in a NanostringGeomxSet object.
#'
#' @param object A NanostringGeomxSet object.
#' @param bg matrix of expected background, on the scale of norm.
#' @param pure_tumor_ids Vector identifying columns of norm that are pure tumor.
#'                       Can be indices, logicals or column names.
#' @param X The training matrix
#' @param K the number of clusters to fit
#' @param norm_elt norm data element in assayData
#' @return an updated X matrix with new columns, "tumor.1", "tumor.2", ...
#' @importFrom Biobase assayData
#' @examples
#' data(target_demoData)
#' data(safeTME)
#' tumor.ids <- as.logical(sample(x = c("TRUE","FALSE"), size = 88, replace = TRUE))
#' safeTME.with.tumor <- runMergeTumorIntoX(object = target_demoData,
#'                                          X = safeTME,
#'                                          K = 3,
#'                                          pure_tumor_ids = tumor.ids,
#'                                          norm_elt = "neg_norm")
#'
#' @export
#'


setGeneric("runMergeTumorIntoX", signature = "object",
           function(object, ...) standardGeneric("runMergeTumorIntoX"))

setMethod("runMergeTumorIntoX", "NanoStringGeoMxSet",
          function(object, X, K = 10, pure_tumor_ids, norm_elt){

            #norm
            norm <- Biobase::assayDataElement(object, elt = norm_elt)

            # estimate background
            bg <- derive_GeoMx_background(norm = norm,
                                          # access the probe pool information from the feature metadata
                                          probepool = Biobase::fData(object)$Module,
                                          # access the names of the negative control probes
                                          negnames = Biobase::fData(negativeControlSubset(object))$TargetName)

            # run runMergeTumorIntoX:
            result <- runMergeTumorIntoX(object = norm,
                                      bg = bg,
                                      pure_tumor_ids = pure_tumor_ids,
                                      X = X,
                                      K = K)


            # append result to the object
            # add tumor_res matrix
            Biobase::fData(object)[["tumor_res"]] <- matrix(NA, nrow = nrow(object), ncol = ncol(result),
                                                    dimnames = list(featureNames(object), colnames(result)))
            Biobase::fData(object)[["tumor_res"]][rownames(result), ] <- result

            return(object)

          })

setMethod("runMergeTumorIntoX", "matrix",
          function(object, bg, pure_tumor_ids, X, K = 10){

            # run mergeTumorIntoX
            res <- mergeTumorIntoX(norm = object,
                                   bg = bg,
                                   pure_tumor_ids = pure_tumor_ids,
                                   X = X,
                                   K = K)

            return(res)

          })
