### ============================================================================
### Self-defined constructor for SangerMergeReads
### ============================================================================
#' @description
#'
#' @param forwardReadFileName .
#' @param reverseReadFileName .
#' @param cutoffQualityScore
#' @param slidingWindowSize
#'
#' @return SangerMergeReads
#' @export
#' @author Kuan-Hao Chao
#' @example
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' samplesRegExp <- "ACHL"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory     = inputFilesParentDir,
#'                                  readsRegularExp     = samplesRegExp,
#'                                  cutoffQualityScore  = 50L,
#'                                  slidingWindowSize   = 8L)
SangerConsensusRead <- function(parentDirectory = character(0),
                                readsRegularExp = character(0),
                                cutoffQualityScore   = 20L,
                                slidingWindowSize    = 5L) {
    newConsensusReads <- new("SangerConsensusRead",
                             parentDirectory    = parentDirectory,
                             readsRegularExp    = readsRegularExp,
                             cutoffQualityScore = cutoffQualityScore,
                             slidingWindowSize  = slidingWindowSize)
    return(newConsensusReads)
}



### ============================================================================
### Self-defined constructor for SangerSingleRead
### ============================================================================
#' @description
#'
#' @param readFeature .
#' @param readFileName .
#' @param cutoffQualityScore .
#' @param slidingWindowSize .
#'
#' @return SangerSingleRead
#' @export
#' @author Kuan-Hao Chao
#' @example
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
#' A_chloroticaSingleRead <-
#'        SangerSingleRead(readFeature         = "ForwardRead",
#'                         readFileName        = A_chloroticaFdReadFN,
#'                         cutoffQualityScore  = 60L,
#'                         slidingWindowSize   = 8L)
SangerSingleRead <- function(readFeature = character(0),
                             readFileName = character(0),
                             cutoffQualityScore   = 20L,
                             slidingWindowSize    = 5L) {
    newSingleReads <- new("SangerSingleRead",
                          readFeature         = readFeature,
                          readFileName        = readFileName,
                          cutoffQualityScore  = cutoffQualityScore,
                          slidingWindowSize   = slidingWindowSize)
    return(newSingleReads)
}
