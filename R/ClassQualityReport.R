#' @title QualityReport
#'
#' @description  An S4 class for quality report for a SangerSingleRead S4 object
#'
#' @slot readFeature .
#' @slot qualityPhredScores .
#' @slot qualityBaseScore .
#' @slot trimmedStartPos .
#' @slot trimmedFinishPos .
#' @slot cutoffQualityScore .
#' @slot slidingWindowSize .
#'
#' @name QualityReport-class
#'
#' @rdname QualityReport-class
#'
#' @exportClass QualityReport
#' @author Kuan-Hao Chao
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
#' A_chloroticaRead <-
#'        SangerSingleRead(readFeature         = "ForwardRead",
#'                         readFileName        = A_chloroticaFdReadFN,
#'                         cutoffQualityScore  = 60L,
#'                         slidingWindowSize   = 8L)
#' "@@"(A_chloroticaRead, QualityReport)
setClass("QualityReport",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             readFeature             = "character",
             qualityPhredScores      = "numeric",
             qualityBaseScore        = "numeric",
             rawSeqLength            = "numeric",
             trimmedSeqLength        = "numeric",
             trimmedStartPos        = "numeric",
             trimmedFinishPos       = "numeric",
             rawSecondaryPeakNum     = "numeric",
             trimmedSecondaryPeakNum = "numeric",
             rawMeanQualityScore     = "numeric",
             trimmedMeanQualityScore = "numeric",
             rawMinQualityScore      = "numeric",
             trimmedMinQualityScore  = "numeric",
             cutoffQualityScore      = "numeric",
             slidingWindowSize       = "numeric"
         ),
)

### ============================================================================
### Overwrite initialize for QualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "QualityReport",
          function(.Object, ...,
                   readFeature         = character(0),
                   qualityPhredScores  = numeric(0),
                   cutoffQualityScore  = 20,
                   slidingWindowSize   = 5) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              if (identical(readFeature, character(0))) {
                  msg <- paste("\nYou must assign value to 'readFeature'\n")
                  errors <- c(errors, msg)
              }
              if (length(errors) == 0) {
                  ### ----------------------------------------------------------
                  ### Prechecking success.
                  ### ----------------------------------------------------------

                  ### ----------------------------------------------------------
                  ### Quality Trimming (Using slideing window VERSION 1)
                  ### ----------------------------------------------------------
                  # calculate base score
                  # Calculate probability error per base (through column)
                  #     ==> Q = -10log10(P)
                  qualityBaseScore <- 10** (qualityPhredScores / (-10.0))
                  trimmingPos <- inside_calculate_trimming(qualityBaseScore,
                                                           cutoffQualityScore,
                                                           slidingWindowSize)
                  trimmedStartPos <- trimmingPos[1]
                  trimmedFinishPos <- trimmingPos[2]
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature         = readFeature,
                             qualityPhredScores  = qualityPhredScores,
                             qualityBaseScore    = qualityBaseScore,
                             trimmedStartPos    = trimmedStartPos,
                             trimmedFinishPos   = trimmedFinishPos,
                             cutoffQualityScore  = cutoffQualityScore,
                             slidingWindowSize   = slidingWindowSize)
          })
