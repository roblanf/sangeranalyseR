### ============================================================================
### Constructor for SangerMergeReads
### ============================================================================
#' @example
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
#' A_chloroticaRvReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R.ab1")
#' A_chloroticaRead <- new("SangerMergeReads",
#'                         forwardReadFileName = A_chloroticaFdReadFN,
#'                         reverseReadFileName = A_chloroticaRvReadFN,
#'                         cutoffQualityScore  = 50L,
#'                         slidingWindowSize   = 8L)
SangerMergeReads <- function(forwardReadSangerseq = character(0),
                             reverseReadSangerseq = character(0),
                             cutoffQualityScore   = 20L,
                             slidingWindowSize    = 5L) {
    A_chloroticaRead <- new("SangerMergeReads",
                            forwardReadFileName = forwardReadSangerseq,
                            reverseReadFileName = A_chloroticaRvReadFN,
                            cutoffQualityScore  = 50L,
                            slidingWindowSize   = 8L)
}
