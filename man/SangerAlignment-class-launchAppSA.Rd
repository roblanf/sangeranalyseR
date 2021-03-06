% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodSangerAlignment.R
\docType{methods}
\name{SangerAlignment-class-launchAppSA}
\alias{SangerAlignment-class-launchAppSA}
\alias{launchAppSA,SangerAlignment-method}
\title{launchAppSA}
\usage{
\S4method{launchAppSA}{SangerAlignment}(object, outputDir = NULL)
}
\arguments{
\item{object}{A SangerAlignment S4 instance.}

\item{outputDir}{The output directory of the saved new SangerContig S4 instance.}
}
\value{
A \code{shiny.appobj} object.
}
\description{
A SangerAlignment method which launches Shiny app for SangerAlignment instance.
}
\examples{
data("sangerAlignmentData")
RShinySA <- launchAppSA(sangerAlignmentData)
}
