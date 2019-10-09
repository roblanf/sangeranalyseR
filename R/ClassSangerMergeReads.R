#' @title SangerMergeReads
#'
#' @description  An S4 class for storing reads of sanger sequencing.
#'
#' @slot forwardReadSangerseq .
#' @slot reverseReadSangerseq .
#'
#' @name SangerMergeReads-class
#'
#' @rdname SangerMergeReads-class
#'
#' @exportClass SangerMergeReads
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerSingleRead.R
#' @examples
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
setClass("SangerMergeReads",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerMergeReads'
         ### -------------------------------------------------------------------
         representation(
             forwardReadSangerseq    = "sangerSingleRead",
             reverseReadSangerseq    = "sangerSingleRead"
         ),
)

### ============================================================================
### Overwrite initialize for 'SangerMergeReads' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerMergeReads",
          function(.Object, ...,
                   forwardReadFileName  = forwardReadFileName,
                   reverseReadFileName  = reverseReadFileName,
                   forwardReadSangerseq = new("sangerSingleRead"),
                   reverseReadSangerseq = new("sangerSingleRead"),
                   cutoffQualityScore   = 20L,
                   slidingWindowSize    = 5L) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()
    if (!file.exists(forwardReadFileName)) {
        msg <- paste("\n'", forwardReadFileName, "'",
                     " foward read file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (!file.exists(reverseReadFileName)) {
        msg <- paste("\n'", reverseReadFileName, "'",
                     " reverse read file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }


    ### ------------------------------------------------------------------------
    ### Prechecking success. Start to create forward and reverse reads.
    ### ------------------------------------------------------------------------
    if (length(errors) == 0) {
        forwardReadSangerseq = new("sangerSingleRead",
                                   readFeature         = "ForwardRead",
                                   readFileName        = forwardReadFileName,
                                   cutoffQualityScore  = cutoffQualityScore,
                                   slidingWindowSize   = slidingWindowSize)
        reverseReadSangerseq = new("sangerSingleRead",
                                   readFeature         = "ReverseRead",
                                   readFileName        = reverseReadFileName,
                                   cutoffQualityScore  = cutoffQualityScore,
                                   slidingWindowSize   = slidingWindowSize)
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   forwardReadSangerseq = forwardReadSangerseq,
                   reverseReadSangerseq = reverseReadSangerseq)
})
