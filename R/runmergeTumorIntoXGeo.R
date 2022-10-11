#' Run MergeTumorIntoX
#' 
#' Runs mergeTumorIntoX from an S4 object
#' @param object An S4 object such as a GeoMxSet object
#' @param ... Arguments passed to mergeTumorIntoX
#' @return updated X matrix with new columns, "tumor.1", "tumor.2", ...
#'

setGeneric("runMergeTumorIntoX", signature = "object",
           function(object, ...) standardGeneric("runMergeTumorIntoX"))

#' Run mergeTumorIntoX on a NanostringGeomxSet object
#'
#' A wrapper for applying mergeTumorIntoX to a NanostringGeomxSet object.
#'
#' @param pure_tumor_ids Vector identifying columns of norm that are pure tumor.
#'                       Can be indices, logicals or column names.
#' @param X The training matrix
#' @param K the number of clusters to fit
#' @param norm_elt norm data element in assayData
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
#' data(safeTME)
#' tumor.ids <- as.logical(sample(x = c("TRUE","FALSE"), size = 88, replace = TRUE))
#' safeTME.with.tumor <- runMergeTumorIntoX(object = target_demoData,
#'                                          X = safeTME,
#'                                          K = 3,
#'                                          pure_tumor_ids = tumor.ids,
#'                                          norm_elt = "exprs_norm")
#'
#' @importFrom methods is
#' @export
#' @rdname runMergeTumorIntoX

setMethod("runMergeTumorIntoX", "NanoStringGeoMxSet",
          function(object, X, K = 10, pure_tumor_ids = NULL, norm_elt = NULL){
              
              if(is.null(norm_elt)){
                  stop("norm_elt must be set")
              }
              
              if(is.null(pure_tumor_ids)){
                  stop("pure_tumor_ids must be set")
              }
              
              if (!is.element(norm_elt, names(object@assayData))) {
                  stop(paste(norm_elt, "is not an element in assaysData slot"))
              }
              
              if (!any(pure_tumor_ids %in% Biobase::sampleNames(object)) & 
                  !is(pure_tumor_ids,"logical") & !is(pure_tumor_ids, "integer")) {
                  stop(paste("No pure_tumor_ids match sample names in GeoMxSet object"))
              }
              
              #norm
              norm <- Biobase::assayDataElement(object, elt = norm_elt)
              
              # estimate background
              bg <- derive_GeoMx_background(norm = norm,
                                            # access the probe pool information from the feature metadata
                                            probepool = Biobase::fData(object)$Module,
                                            # access the names of the negative control probes
                                            negnames = Biobase::fData(object)$TargetName[Biobase::fData(object)$Negative])
              
              # run runMergeTumorIntoX:
              result <- mergeTumorIntoX(norm = norm,
                                        bg = bg,
                                        pure_tumor_ids = pure_tumor_ids,
                                        X = X,
                                        K = K)
              
              # # append result to the object
              # # add tumor_res matrix
              # Biobase::fData(object)[["tumor_res"]] <- matrix(NA, nrow = nrow(object), ncol = ncol(result),
              #                                         dimnames = list(Biobase::featureNames(object), colnames(result)))
              # Biobase::fData(object)[["tumor_res"]][rownames(result), ] <- result
              
              return(result)
              
          })

