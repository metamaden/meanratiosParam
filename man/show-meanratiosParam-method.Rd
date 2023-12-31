% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanratiosParam-class.R
\name{show,meanratiosParam-method}
\alias{show,meanratiosParam-method}
\title{Show generic behavior for object of class \linkS4class{meanratiosParam}}
\usage{
\S4method{show}{meanratiosParam}(object)
}
\arguments{
\item{object}{An object of class \linkS4class{meanratiosParam} (see 
\code{?meanratiosParam}).}
}
\value{
Shows object summaries.
}
\description{
Show generic behavior for object of class \linkS4class{meanratiosParam}
}
\details{
Method for behavior of show generic when called for object of class 
\linkS4class{meanratiosParam}
}
\examples{
exampleList <- getDeconvolutionExampleData()
singleCellExperimentExample <- randomSingleCellExperiment()
newParam <- meanratiosParam(singleCellExperiment=singleCellExperimentExample, 
cellTypeVariable="celltype", markersPerType=5)
markers <- typemarkers(newParam)
markers

}
\references{
Louise Huuki-Meyers. DeconvoBuddies. (2023). R package. GitHub URL https://github.com/LieberInstitute/DeconvoBuddies.
}
