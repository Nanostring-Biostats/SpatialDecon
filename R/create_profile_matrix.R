# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
# Copyright (C) 2020, NanoString Technologies, Inc.
#        This program is free software: you can redistribute it and/or modify it
#        under the terms of the GNU General Public License as published by the Free
#        Software Foundation, either version 3 of the License, or (at your option)
#        any later version.
#        This program is distributed in the hope that it will be useful, but WITHOUT
#        ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#        FITNESS FOR A PARTICULAR PURPOSE.    See the GNU General Public License for
#        more details.
#        You should have received a copy of the GNU General Public License along
#        with this program.    If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.
# 530 Fairview Avenue N
# Seattle, WA 98109
# Tel: (888) 358-6266
# pdanaher@nanostring.com


#' Create Custom Cell Profile Matrix
#'
#' Create custom cell profile matrix using single cell data. The average gene expression for each cell type is returned.
#'
#' @param mtx gene x cell count matrix
#' @param cellAnnots cell annotations with cell type and cell name as columns
#' @param cellTypeCol column containing cell type
#' @param cellNameCol column containing cell ID/name
#' @param matrixName name of final profile matrix
#' @param outDir path to desired output directory, set to NULL if matrix should not be written
#' @param geneList gene list to filter profile matrix to 
#' @param normalize Should data be normalized? (TRUE/FALSE) if TRUE data will be normalize using total gene count
#' @param scalingFactor what should all values be multiplied by for final matrix, set to 1 if no scaling is wanted
#' @param minCellNum minimum number of cells of one type needed to create profile, exclusive 
#' @param minGenes minimum number of genes expressed in a cell, exclusive
#' @param discardCellTypes should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.
#' @return A custom cell profile matrix genes (rows) by cell types (columns), matrix gets written to disk and outDir
#' @examples
#' cellNames <- paste0("Cell", seq_len(1500))
#' geneNames <- paste0("Gene", seq_len(1500))
#' mtx <- matrix(data=sample(size = length(cellNames)*length(geneNames),
#'                           replace = TRUE,
#'                           x = c(0,seq_len(100)), 
#'                           prob = c(0.6784, rep(0.0075, 15), rep(0.005, 25),
#'                                    rep(0.002, 25), rep(0.001, 35))), 
#'                           ncol = length(cellNames), nrow = length(geneNames), 
#'                           dimnames = list(geneNames, cellNames))
#' cellAnnots <- as.data.frame(cbind(CellID=cellNames, 
#'                                   cellType=sample(size = length(cellNames), 
#'                                                   replace = TRUE,
#'                                                   x = c("A", "B", "C", "D"),
#'                                                   prob = c(0.1, 0.4, 0.3, 0.2))))
#' table(cellAnnots$cellType)
#' profile_matrix <- create_profile_matrix(mtx = mtx,    
#'                                         cellAnnots = cellAnnots, 
#'                                         cellTypeCol = "cellType",
#'                                         cellNameCol = "CellID",
#'                                         minGenes = 10,
#'                                         scalingFactor = 1)
#' head(profile_matrix)
#' @importFrom methods is
#' @export

create_profile_matrix <- function(mtx, cellAnnots, cellTypeCol, cellNameCol,
                                  matrixName = "Custom", outDir = "./", 
                                  geneList = NULL, normalize = FALSE,
                                  scalingFactor = 5, minCellNum = 15,
                                  minGenes = 100, discardCellTypes = FALSE) {
    
    # checking user input values
    if(is.null(mtx)){
        stop("count matrix is necessary")
    }
    if(!is.null(outDir)){
        if (!dir.exists(outDir) ) {
            stop("Output directory is not valid")
        }
    }
    if (is.null(cellAnnots)){
        stop("Cell Annotations are needed")
    }
    
    if(is.null(cellTypeCol)){
        stop("cellTypeCol must not be NULL")
    }else if (!cellTypeCol %in% colnames(cellAnnots)){
        stop("cellTypeCol not in cellAnnots")
    }
    
    if(is.null(cellNameCol)){
        stop("cellNameCol must not be NULL")
    }else if (!cellNameCol %in% colnames(cellAnnots)){
        stop("cellNameCol not in cellAnnots")
    }
    
    if(!is(normalize, "logical")){
        warning("normalize not a boolean, continuing with assumption that data should not be normalized")
        normalize <- TRUE
    }
    if(!is(discardCellTypes, "logical")){
        warning("discardCellTypes not a boolean, continuing with default of discarding cell types")
        discardCellTypes <- TRUE
    }
    
    if(!is(scalingFactor, "numeric")){
        warning("scalingFactor not a numeric, continuing with default value of 5")
        scalingFactor <- 5
    }
    if(!is(minCellNum, "numeric")){
        warning("minCellNum not a numeric, continuing with default value of 15")
        minCellNum <- 15
    }
    if(!is(minGenes, "numeric")){
        warning("minGenes not a numeric, continuing with default value of 100")
        minGenes <- 100
    }
    
    #make a sparse matrix
    mtx <- Matrix::Matrix(as.matrix(mtx), sparse = TRUE) 
    
    cellTypes <- NULL
    
    #read in cell type annotation file
    #get cell types 
    cellTypes <- cellAnnots[[cellTypeCol]]
    #assign cell name to type
    names(cellTypes) <- cellAnnots[[cellNameCol]]
    
    if(is.null(cellTypes)){
        stop("cellAnnots and/or cellTypeCol arguments are incorrectly formatted. 
                 cellAnnots should be a data frame, and cellTypeCol should give the 
                 name of the column holding each cell's cell type.")
    }
    
    rm(cellAnnots)
    
    # mtx <- as.data.frame(mtx)
    
    if(!any(names(cellTypes) %in% colnames(mtx)) & 
       any(names(cellTypes) %in% rownames(mtx))){
        print("Transposing Matrix")
        mtx<- t(mtx)
    }
    
    if(!any(names(cellTypes) %in% colnames(mtx))){
        stop(paste("cellNameCol names does not match count matrix column names", 
                   "matrix cell names:", colnames(mtx)[1], "annots cell names:", names(cellTypes)[1]))
    }else if(!all(names(cellTypes) %in% colnames(mtx))){
        missing <- length(which(!names(cellTypes) %in% colnames(mtx)))
        warning(paste("not all cellNameCol names are in count matrix;", missing, "cells are missing"))
    }
    
    if(discardCellTypes == TRUE){
        #remove cells with no cell type assignment
        w2rm <- which(is.na(cellTypes) | tolower(cellTypes) %in% c("unspecified", "unknown", "not available")) 
        w2rm <- unique(c(w2rm, grep(pattern = "doublet|dividing|low q|filtered|mitotic", x = tolower(cellTypes))))
        if(length(w2rm) > 0){
            cellTypes <- cellTypes[-w2rm]
        }
    }
    
    #normalize data if necessary 
    if(normalize == TRUE){
        print("Normalizing Matrix")
        med <- median(Matrix::colSums(mtx))
        cols <- colnames(mtx)
        rows <- rownames(mtx)
        mtx <- Matrix::Matrix(sweep(mtx, 2, Matrix::colSums(mtx), "/") * med, 
                              sparse = TRUE) 
        
        colnames(mtx) <- cols
        rownames(mtx) <- rows
        
        rm(cols,rows)
    }
    
    atlas <- NULL
    
    #get all unique cell types
    CTs <- unique(cellTypes)
    
    print("Creating Atlas")
    
    #change to apply() if bioconductor requires it
    for(i in CTs){
        #print log of progress
        print(paste(which(CTs == i), "/", length(CTs), ":", i))
        
        #get cell names for this cell type
        cellsType <- names(cellTypes)[which(cellTypes == i)]
        #confirm cells are in matrix
        cellsType <- cellsType[which(cellsType %in% colnames(mtx))]
        
        if(length(cellsType) > minCellNum){
            
            if(length(cellsType) > 1){
                #remove cells with low gene expression
                cellsType <- cellsType[which(Matrix::colSums(mtx[,cellsType] > 0) > minGenes)]
            }else{
                cellsType <- cellsType[which(sum(mtx[,cellsType] > 0) > minGenes)]
            }
            
            if(length(cellsType) > minCellNum){
                #get average expression if there are enough cells for cell type
                if(length(cellsType) > minCellNum & length(cellsType) != 1){
                    atlas <- as.data.frame(cbind(atlas, Matrix::rowMeans(mtx[,cellsType], na.rm = TRUE)))
                    colnames(atlas)[ncol(atlas)] <- i
                }else{
                    atlas <- as.data.frame(cbind(atlas, mtx[,cellsType]))
                    colnames(atlas)[ncol(atlas)] <- i
                }
            }else{
                warning(paste("\n", i, "was dropped from matrix because it didn't have enough viable cells based on current filtering thresholds. 
                                        If this cell type is necessary consider changing minCellNum or minGenes\n"))
            }
        }else{
            warning(paste("\n", i, "was dropped from matrix because it didn't have enough viable cells based on current filtering thresholds. 
                                        If this cell type is necessary consider changing minCellNum or minGenes\n"))
        }
    }
    
    rm(cellTypes)
    
    numCellTypesExpr <- 1
    
    #subset to genes expressed in at least a user defined number of cell type(s)
    if(ncol(atlas) == 1){
        w2kp <- which(Matrix::rowSums(atlas > 0) >= numCellTypesExpr)
        cols <- colnames(atlas)
        rows <- rownames(atlas)[w2kp]
        
        atlas <- as.matrix(atlas[w2kp,])
        
        colnames(atlas) <- cols
        rownames(atlas) <- rows
        
        rm(cols, rows, w2kp)
    }else{
        atlas <- atlas[which(Matrix::rowSums(atlas > 0) >= numCellTypesExpr),]
    }
    
    #scale data
    atlas <- atlas * scalingFactor
    
    if(!is.null(geneList)){
        if(any(geneList %in% rownames(atlas))){
            #subset to genes in panel
            atlas <- atlas[rownames(atlas) %in% geneList,]
        }else{
            warning("geneList genes do not match genes in matrix, no filtering done")
        }
    }
    
    #ensure cell types don't contain ","
    colnames(atlas) <- gsub(pattern = ",", replacement = "-", x = colnames(atlas))
    
    if(!is.null(outDir)){
        #write profile matrix
        write.table(atlas, file = paste0(outDir, "/", matrixName, "_profileMatrix.csv"), 
                    row.names = TRUE, col.names = NA, quote = FALSE, sep = ",")
    }
    
    return(as.matrix(atlas))
}

