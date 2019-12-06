#' @title QualityReport
#'
#' @description  An S4 class for quality report for a SangerSingleRead S4 object
#'
#' @slot qualityPhredScores .
#' @slot qualityBaseScores .
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
#'                              "Allolobophora_chlorotica",
#'                              "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
#' A_chloroticaRead <-
#'        SangerSingleRead(readFeature           = "Forward Read",
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
             TrimmingMethod          = "character",
             M1TrimmingCutoff        = "numericORNULL",
             M2CutoffQualityScore    = "numericORNULL",
             M2SlidingWindowSize     = "numericORNULL",
             qualityPhredScores      = "numeric",
             qualityBaseScores       = "numeric",
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
             remainingRatio         = "numeric"
         ),
)

### ============================================================================
### Overwrite initialize for QualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "QualityReport",
          function(.Object, ...,
                   qualityPhredScores = numeric(0),
                   TrimmingMethod        = "M1",
                   M1TrimmingCutoff      = 0.0001,
                   M2CutoffQualityScore  = NULL,
                   M2SlidingWindowSize   = NULL) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              errors <- checkQualityPhredScores (qualityPhredScores, errors)

              ##### ------------------------------------------------------------
              ##### Input parameter prechecking for TrimmingMethod.
              ##### ------------------------------------------------------------
              errors <- checkTrimParam(TrimmingMethod,
                                       M1TrimmingCutoff,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize,
                                       errors)
              if (length(errors) == 0) {
                  ### ----------------------------------------------------------
                  ### Prechecking success.
                  ### ----------------------------------------------------------

                  # calculate base score
                  # Calculate probability error per base (through column)
                  #     ==> Q = -10log10(P)
                  qualityBaseScores <- 10** (qualityPhredScores / (-10.0))

                  ### ----------------------------------------------------------
                  ### Use 'qualityPhredScores' & 'qualityBaseScores' for
                  ###     trimming ==> Not the raw one !!!
                  ### ----------------------------------------------------------
                  if (TrimmingMethod == "M1") {
                      ### ------------------------------------------------------
                      ### Quality Trimming (Using Logarithmic Scale
                      ###                   Sliding Window Trimming METHOD 1)
                      ### ------------------------------------------------------
                      trimmingPos <-
                          M1inside_calculate_trimming(qualityPhredScores,
                                                      qualityBaseScores,
                                                      M1TrimmingCutoff)
                  } else if (TrimmingMethod == "M2") {
                      ### ------------------------------------------------------
                      ### Quality Trimming (Using Logarithmic
                      ###                   Scale Trimming  METHOD 2)
                      ### ------------------------------------------------------
                      trimmingPos <-
                          M2inside_calculate_trimming(qualityPhredScores,
                                                      qualityBaseScores,
                                                      M2CutoffQualityScore,
                                                      M2SlidingWindowSize)
                  }
                  rawSeqLength <- trimmingPos[["rawSeqLength"]]
                  rawMeanQualityScore <- trimmingPos[["rawMeanQualityScore"]]
                  rawMinQualityScore <- trimmingPos[["rawMinQualityScore"]]
                  trimmedStartPos <- trimmingPos[["trimmedStartPos"]]
                  trimmedFinishPos <- trimmingPos[["trimmedFinishPos"]]
                  trimmedSeqLength <- trimmingPos[["trimmedSeqLength"]]
                  trimmedMeanQualityScore <-
                      trimmingPos[["trimmedMeanQualityScore"]]
                  trimmedMinQualityScore <-
                      trimmingPos[["trimmedMinQualityScore"]]
                  remainingRatio <- trimmingPos[["remainingRatio"]]
                  qualityScoresID = "Before Basecall"
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             qualityPhredScores      = qualityPhredScores,
                             qualityBaseScores       = qualityBaseScores,
                             rawSeqLength            = rawSeqLength,
                             trimmedSeqLength        = trimmedSeqLength,
                             trimmedStartPos         = trimmedStartPos,
                             trimmedFinishPos        = trimmedFinishPos,
                             rawMeanQualityScore     = rawMeanQualityScore,
                             trimmedMeanQualityScore = trimmedMeanQualityScore,
                             rawMinQualityScore      = rawMinQualityScore,
                             trimmedMinQualityScore  = trimmedMinQualityScore,
                             remainingRatio          = remainingRatio,
                             TrimmingMethod          = TrimmingMethod,
                             M1TrimmingCutoff        = M1TrimmingCutoff,
                             M2CutoffQualityScore    = M2CutoffQualityScore,
                             M2SlidingWindowSize     = M2SlidingWindowSize)
          })
