#!/usr/bin/env R

# Author: Sean Maden

#' meanratiosParam-class
#'
#' class definition for meanratiosParam, which uses 
#' DeconvoBuddies::get_mean_ratio2().
#' 
#' @param assayName Name of expression matrix in SingleCellExperiment assays 
#' (e.g. "counts").
#' @param singleCellExperiment Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @param cellTypeVariable Name of cell type variable in SingleCellExperiment 
#' coldata.
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- randomSingleCellExperiment()
#' newParam <- meanratiosParam(singleCellExperiment=singleCellExperimentExample, cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#' 
#' @details Main constructor for class \linkS4class{meanratiosParam}.
#' @rdname meanratiosParam-class
#' @seealso \linkS4class{typemarkersParam}
#' @aliases 
#' MeanratiosParam-class, MeanRatiosParam-class
#' 
#' @references 
#' Louise Huuki-Meyers. DeconvoBuddies. (2023). R package. GitHub URL https://github.com/LieberInstitute/DeconvoBuddies.
#' 
setClass("meanratiosParam", contains="typemarkersParam", 
         slots = c(assayName="character", 
                   singleCellExperiment="SingleCellExperiment",
                   cellTypeVariable="character"))

#' Make new object of class meanratiosParam
#'
#' Main constructor for class \linkS4class{meanratiosParam}.
#'
#' @param assayName Name of expression matrix in SingleCellExperiment assays 
#' (e.g. "counts").
#' @param singleCellExperiment Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @param cellTypeVariable Name of cell type variable in SingleCellExperiment 
#' coldata.
#' @param markersPerType Number of top markers to get per cell type.
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#'
#' @details Main class for mapping arguments to the mean ratios method 
#' implemented as \code{DeconvoBuddies::get_mean_ratio2()}.
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- randomSingleCellExperiment()
#' newParam <- meanratiosParam(singleCellExperiment=singleCellExperimentExample, cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#' markers
#' 
#' @references 
#' Louise Huuki-Meyers. DeconvoBuddies. (2023). R package. GitHub URL https://github.com/LieberInstitute/DeconvoBuddies.
#' 
#' @returns Object of class \linkS4class{meanratiosParam}
#' @seealso \linkS4class{typemarkersParam}
#' @export
meanratiosParam <- function(singleCellExperiment, assayName = "counts", 
                            cellTypeVariable = "cellType", 
                            markersPerType = 20, 
                            returnInfo = FALSE) {
  new("meanratiosParam", 
      singleCellExperiment=singleCellExperiment, assayName=assayName, 
      cellTypeVariable=cellTypeVariable, markersPerType=markersPerType, 
      returnInfo=returnInfo)
}

#' Cell type markers method for meanratiosParam
#'
#' Defines the typemarkers method for \linkS4class{meanratiosParam}.
#'
#' @importFrom DeconvoBuddies get_mean_ratio2
#' @importFrom dplyr %>%
#'
#' @param object An object of class \linkS4class{meanratiosParam}.
#' 
#' @returns Either a vector of gene markers, or a list containing such a vector 
#' with the original method outputs.
#' 
#' @details Takes an object of class \linkS4class{meanratiosParam} as input, 
#' returning either a vector of cell type gene markers, or (if 
#' \code{return.info == TRUE}) a list containing such a vector along with 
#' original function outputs.
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- randomSingleCellExperiment()
#' newParam <- meanratiosParam(singleCellExperiment=singleCellExperimentExample, cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#' markers
#'
#' @references 
#' Louise Huuki-Meyers. DeconvoBuddies. (2023). R package. GitHub URL https://github.com/LieberInstitute/DeconvoBuddies.
#' 
#' @export
setMethod("typemarkers", signature(object = "meanratiosParam"), function(object){
  singleCellExperiment <- object[["singleCellExperiment"]]
  cellTypeVariable <- object[["cellTypeVariable"]]
  assayName <- object[["assayName"]]
  markersPerType <- object[["markersPerType"]]
  # get marker results
  markerTable <- 
    DeconvoBuddies::get_mean_ratio2(
      sce=singleCellExperiment, cellType_col=cellTypeVariable, 
      assay_name=assayName, add_symbol=FALSE) |> as.data.frame()
  # filter top markers
  uniqueCellTypesVector <- singleCellExperiment[[cellTypeVariable]] |> 
    as.character() |> unique()
  
  topMarkersList <- lapply(uniqueCellTypesVector, function(uniqueCellType){
    markerTable %>% 
      dplyr::filter(cellType.target == uniqueCellType) %>% 
      dplyr::arrange(rank_ratio) %>% 
      dplyr::top_n(n = markersPerType)
  })
  topMarkersTable <- do.call(rbind, topMarkersList)
  topMarkersVector <- topMarkersTable$gene
  ## parse returnInfo
  returnList <- topMarkersVector %>% unique()
  if(object[["returnInfo"]]){
    returnList <- list(
      markers=topMarkersVector, result.info=topMarkerTable, metadata=object)}
  return(returnList)
})

#' Show generic behavior for object of class \linkS4class{meanratiosParam}
#' @param object An object of class \linkS4class{meanratiosParam} (see 
#' \code{?meanratiosParam}).
#' @details Method for behavior of show generic when called for object of class 
#' \linkS4class{meanratiosParam}
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- randomSingleCellExperiment()
#' newParam <- meanratiosParam(singleCellExperiment=singleCellExperimentExample, 
#' cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#' markers
#' 
#' @returns Shows object summaries.
#' 
#' @references 
#' Louise Huuki-Meyers. DeconvoBuddies. (2023). R package. GitHub URL https://github.com/LieberInstitute/DeconvoBuddies.
#' 
#' @export
setMethod("show", signature(object="meanratiosParam"), function(object){
  show(object)
})