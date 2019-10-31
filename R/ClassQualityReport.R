#' @title QualityReport
#'
#' @description  An S4 class for quality report for a SangerSingleRead S4 object
#'
#' @slot readFeature .
#' @slot qualityPhredScores .
#' @slot qualityBaseScore .
#' @slot trimmingStartPos .
#' @slot trimmingFinishPos .
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
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
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
             trimmingStartPos        = "numeric",
             trimmingFinishPos       = "numeric",
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
                   qualityPhredScores = qualityPhredScores,
                   qualityBaseScore    = 0,
                   trimmingStartPos    = 0,
                   trimmingFinishPos   = 0,
                   cutoffQualityScore  = 20,
                   slidingWindowSize   = 5) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
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
                  trimmingStartPos <- trimmingPos[1]
                  trimmingFinishPos <- trimmingPos[2]
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature         = readFeature,
                             qualityPhredScores = qualityPhredScores,
                             qualityBaseScore    = qualityBaseScore,
                             trimmingStartPos    = trimmingStartPos,
                             trimmingFinishPos   = trimmingFinishPos,
                             cutoffQualityScore  = cutoffQualityScore,
                             slidingWindowSize   = slidingWindowSize)
          })
