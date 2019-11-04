### ============================================================================
### Self-defined constructor for SangerMergeReads
### ============================================================================
#' @description
#'
#' @param forwardReadFileName .
#' @param consenesusReadName .
#' @param suffixForwardRegExp .
#' @param suffixReverseRegExp .
#' @param cutoffQualityScore .
#' @param slidingWindowSize .
#' @param refAminoAcidSeq .
#' @param minReadsNum .
#' @param minReadLength .
#' @param minFractionCall .
#' @param maxFractionLost .
#' @param geneticCode .
#' @param acceptStopCodons .
#' @param readingFrame .
#' @param processorsNum .
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
SangerConsensusRead <- function(parentDirectory        = character(0),
                                consenesusReadName     = character(0),
                                suffixForwardRegExp    = character(0),
                                suffixReverseRegExp    = character(0),
                                cutoffQualityScore     = 20,
                                slidingWindowSize      = 5,
                                refAminoAcidSeq        = "",
                                minReadsNum            = 2,
                                minReadLength          = 20,
                                minFractionCall        = 0.5,
                                maxFractionLost        = 0.5,
                                geneticCode            = GENETIC_CODE,
                                acceptStopCodons       = TRUE,
                                readingFrame           = 1,
                                processorsNum          = 1) {
    newConsensusReads <- new("SangerConsensusRead",
                             parentDirectory        = parentDirectory,
                             consenesusReadName     = consenesusReadName,
                             suffixForwardRegExp    = suffixForwardRegExp,
                             suffixReverseRegExp    = suffixReverseRegExp,
                             cutoffQualityScore     = cutoffQualityScore,
                             slidingWindowSize      = slidingWindowSize,
                             refAminoAcidSeq        = refAminoAcidSeq,
                             minReadsNum            = minReadsNum,
                             minReadLength          = minReadLength,
                             minFractionCall        = minFractionCall,
                             maxFractionLost        = maxFractionLost,
                             geneticCode            = geneticCode,
                             acceptStopCodons       = acceptStopCodons,
                             readingFrame           = readingFrame,
                             processorsNum          = processorsNum)
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
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
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
