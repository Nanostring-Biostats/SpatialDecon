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


#' Download a cell profile matrix
#'
#' Download a cell profile matrix from the online library
#'
#' @param species species of profile matrix 
#' @param age_group age_group of profile matrix, if fetal mouse please add the 
#' developmental stage separated with /, i.e. Fetal/E14.5 
#' @param matrixname name of profile matrix 
#' @return A cell profile matrix, suggested cell groups, and paper metadata
#' @details Valid matrices can be found on the github site 
#' \url{https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/NewProfileMatrices}
#' @examples
#' download_profile_matrix(species = "Human", age_group = "Adult", matrixname = "Colon_HCA")
#' head(profile_matrix)
#' print(cellGroups)
#' print(metadata)
#' @importFrom utils read.csv 
#' @importFrom repmis source_data
#' @export
download_profile_matrix <- function(species, age_group, matrixname) {

    # check formatting:
    if (length(species) > 1) {
        stop("specify just one species")
    }
    if (length(age_group) > 1) {
        stop("specify just one age")
    }
    if (length(matrixname) > 1) {
        stop("specify just one matrixname")
    }
    
    valid_species <- c("Human", "Mouse")
    
    if (!species %in% valid_species) {
        stop(paste0("Species input is invalid; must be \"", paste(valid_species, 
                                                                  collapse = "\" or \""), "\" (case sensitive)"))
    }
    
    if(species == "Human"){
        valid_ages <- c("Adult", "COVID-Infected", "Fetal")
    }else{
        valid_ages <- c("Adult", "Fetal/E14.5", "Fetal/E9.5-13.5", "Neonatal")
    }
    
    if (!age_group %in% valid_ages) {
        stop(paste0("Age input is invalid; must be \"", paste(valid_ages, collapse = "\" or \""), "\" (case sensitive)"))
    }
    
    metadata <- read.csv(paste0("https://raw.github.com/Nanostring-Biostats/CellProfileLibrary/NewProfileMatrices/", species, "/",
                                    species, "_datasets_metadata.csv"), header = T, sep = ",")
    
    librarynames <- paste0(metadata$Tissue, "_", metadata$Profile.Matrix)

    
    if (!is.element(matrixname, librarynames)) {
        stop(paste0(matrixname, " is not an expected cell profile matrix name. Did you mean \"", 
                    paste(librarynames[agrep(matrixname, librarynames)], collapse = "\" or \""), "\"?"))
    }
    
    matrixname <- paste(species, age_group, matrixname, sep = "/")

    suppressMessages(source_data(paste0("https://raw.github.com/Nanostring-Biostats/CellProfileLibrary/NewProfileMatrices/", 
                                        matrixname, ".RData?raw=True"), 
                cache = FALSE, rdata = TRUE, envir = globalenv()))
    
    assign("profile_matrix", as.matrix(profile_matrix), envir = globalenv())
}

