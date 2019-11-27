#' @title QualityReport
#'
#' @description  An S4 class for quality report for a SangerSingleRead S4 object
#'
#' @slot readFeature .
#' @slot qualityPhredScores .
#' @slot qualityBaseScore .
#' @slot rawSeqLength .
#' @slot trimmedSeqLength .
#' @slot trimmedStartPos .
#' @slot trimmedFinishPos .
#' @slot rawMeanQualityScore .
#' @slot trimmedMeanQualityScore .
#' @slot rawMinQualityScore .
#' @slot trimmedMinQualityScore .
#' @slot remainingRatio .
#' @slot TrimmingMethod .
#' @slot M1TrimmingCutoff .
#' @slot M2CutoffQualityScore .
#' @slot M2SlidingWindowSize .
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
#'        SangerSingleRead(readFeature           = "ForwardRead",
#'                         readFileName          = A_chloroticaFdReadFN,
#'                         TrimmingMethod        = "M2",
#'                         M1TrimmingCutoff      = NULL,
#'                         M2CutoffQualityScore  = 40,
#'                         M2SlidingWindowSize   = 10)
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
             trimmedStartPos         = "numeric",
             trimmedFinishPos        = "numeric",
             # rawSecondaryPeakNum     = "numeric",
             # trimmedSecondaryPeakNum = "numeric",
             rawMeanQualityScore     = "numeric",
             trimmedMeanQualityScore = "numeric",
             rawMinQualityScore      = "numeric",
             trimmedMinQualityScore  = "numeric",
             remainingRatio         = "numeric",
             TrimmingMethod          = "character",
             M1TrimmingCutoff        = "numericORNULL",
             M2CutoffQualityScore    = "numericORNULL",
             M2SlidingWindowSize     = "numericORNULL"
         ),
)

setClassUnion("numericORNULL", c("numeric", "NULL"))

### ============================================================================
### Overwrite initialize for QualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "QualityReport",
          function(.Object, ...,
                   readFeature           = character(0),
                   qualityPhredScores    = numeric(0),
                   TrimmingMethod        = "M1",
                   M1TrimmingCutoff      = 0.0001,
                   M2CutoffQualityScore  = NULL,
                   M2SlidingWindowSize   = NULL) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              if (identical(readFeature, character(0))) {
                  msg <- paste("\nYou must assign value to 'readFeature'\n")
                  errors <- c(errors, msg)
              }

              if (length(qualityPhredScores) == 0) {
                  msg <- paste("\n'qualityPhredScores' size cannot be zero.\n")
                  errors <- c(errors, msg)
              }

              ### --------------------------------------------------------------
              ### Input parameter prechecking for TrimmingMethod.
              ### --------------------------------------------------------------
              errors <- checkTrimParam(TrimmingMethod,
                                       M1TrimmingCutoff,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize,
                                       errors)

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

                  if (TrimmingMethod == "M1") {

                  } else if (TrimmingMethod == "M2") {
                      trimmingPos <-
                          inside_calculate_trimming(qualityPhredScores,
                                                    qualityBaseScore,
                                                    M2CutoffQualityScore,
                                                    M2SlidingWindowSize)
                  }
                  rawSeqLength <- trimmingPos[1]
                  rawMeanQualityScore <- trimmingPos[2]
                  rawMinQualityScore <- trimmingPos[3]
                  trimmedStartPos <- trimmingPos[4]
                  trimmedFinishPos <- trimmingPos[5]
                  trimmedSeqLength <- trimmingPos[6]
                  trimmedMeanQualityScore <- trimmingPos[7]
                  trimmedMinQualityScore <- trimmingPos[8]
                  remainingRatio <- trimmingPos[9]
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature             = readFeature,
                             qualityPhredScores      = qualityPhredScores,
                             qualityBaseScore        = qualityBaseScore,
                             rawSeqLength            = rawSeqLength,
                             trimmedSeqLength        = trimmedSeqLength,
                             trimmedStartPos         = trimmedStartPos,
                             trimmedFinishPos        = trimmedFinishPos,
                             rawMeanQualityScore     = rawMeanQualityScore,
                             trimmedMeanQualityScore = trimmedMeanQualityScore,
                             rawMinQualityScore      = rawMinQualityScore,
                             trimmedMinQualityScore  = trimmedMinQualityScore,
                             remainingRatio         = remainingRatio,
                             TrimmingMethod          = TrimmingMethod,
                             M1TrimmingCutoff        = M1TrimmingCutoff,
                             M2CutoffQualityScore    = M2CutoffQualityScore,
                             M2SlidingWindowSize     = M2SlidingWindowSize)
          })
