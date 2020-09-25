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
#' @param matrixname A name
#' @return A cell profile matrix
#' @details Valid values for the matrixname argument include:
#' \itemize{
#' \item Airway_Epithelium
#' \item Atlas_Adult_Retina_10x
#' \item Census_Adult_Immune_10x
#' \item Census_Newborn_Blood_10x
#' \item Diff_Fetal_Neuron_SS2
#' \item FetalMaternal_Adult_Blood_10x
#' \item FetalMaternal_Adult_Blood_SS2
#' \item FetalMaternal_Adult_Decidua_10x
#' \item FetalMaternal_Adult_Decidua_SS2
#' \item FetalMaternal_Fetal_Placenta_10x
#' \item Human_brain
#' \item Human_Cell_Landscape
#' \item IBD_Adult_Colon_10x
#' \item Landscape_Adult_Liver_10x
#' \item Lung_plus_neutrophils
#' \item Mouse_Brain
#' \item Profiling_Adult_BoneMarrow_10x
#' \item Reprogram_Embryo_Dendritic_10x
#' \item Sensitivity_Adult_Esophagus_10x
#' \item Sensitivity_Adult_Lung_10x
#' \item Sensitivity_Adult_Spleen_10x
#' \item Somatic_Adult_Pancreas_SS2
#' \item SpatioTemporal_Adult_Kidney_10x
#' \item SpatioTemporal_Fetal_Kidney_10x
#' \item Tcell_Adult_Blood_10x
#' \item Tcell_Adult_BoneMarrow_10x
#' \item Tcell_Adult_Lung_10x
#' \item Tcell_Adult_LymphNode_10x
#' }
#' @examples
#' X <- download_profile_matrix(matrixname = "Human_brain")
#' head(X)
#' @export
download_profile_matrix <- function(matrixname) {

    # check formatting:
    if (length(matrixname) > 1) {
        stop("specify just one matrixname")
    }

    librarynames <- c(
        "Airway_Epithelium", "Atlas_Adult_Retina_10x", "Census_Adult_Immune_10x",
        "Census_Newborn_Blood_10x", "Diff_Fetal_Neuron_SS2",
        "FetalMaternal_Adult_Blood_10x", "FetalMaternal_Adult_Blood_SS2",
        "FetalMaternal_Adult_Decidua_10x", "FetalMaternal_Adult_Decidua_SS2",
        "FetalMaternal_Fetal_Placenta_10x", "Human_brain", "Human_Cell_Landscape",
        "IBD_Adult_Colon_10x", "Landscape_Adult_Liver_10x",
        "Lung_plus_neutrophils", "Mouse_Brain", "Profiling_Adult_BoneMarrow_10x",
        "Reprogram_Embryo_Dendritic_10x", "Sensitivity_Adult_Esophagus_10x",
        "Sensitivity_Adult_Lung_10x", "Sensitivity_Adult_Spleen_10x",
        "Somatic_Adult_Pancreas_SS2", "SpatioTemporal_Adult_Kidney_10x",
        "SpatioTemporal_Fetal_Kidney_10x", "Tcell_Adult_Blood_10x",
        "Tcell_Adult_BoneMarrow_10x", "Tcell_Adult_Lung_10x",
        "Tcell_Adult_LymphNode_10x"
    )
    if (!is.element(matrixname, librarynames)) {
        warning(paste0(matrixname, " is not an expected cell profile matrix name."))
    }

    X <- as.matrix(utils::read.csv(paste0(
        "https://raw.githubusercontent.com/patrickjdanaher/cell-profile-library/master/profile_matrices/",
        matrixname, ".csv"
    ), row.names = 1))

    X <- X[rowSums(X) > 0, ]
    return(X)
}
