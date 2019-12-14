#' @title QualityReport
#'
#' @description  An S4 class storing quality related inputs and results in a SangerRead S4 object.
#'
#' @slot qualityPhredScoresRaw The Phred quality scores of each base pairs from abif class PCON.2 in sangerseqR package before base calling.
#' @slot qualityPhredScores The Phred quality scores of each base pairs after base calling.
#' @slot qualityBaseScores The probability of incorrect base call of each base pairs. They are calculated from \code{qualityPhredScores}.
#' @slot rawSeqLength The number of nucleotides of raw primary DNA sequence.
#' @slot trimmedSeqLength The number of nucleotides of trimeed  primary DNA sequence.
#' @slot trimmedStartPos The base pair index of trimming start point from 5' end of the sequence.
#' @slot trimmedFinishPos The base pair index of trimming finish point from 3' end of the sequence.
#' @slot rawMeanQualityScore The mean quality score of the primary sequence after base calling. In other words, it is the mean of \code{qualityPhredScores}.
#' @slot trimmedMeanQualityScore The mean quality score of the trimmed primary sequence after base calling.
#' @slot rawMinQualityScore The minimum quality score of the primary sequence after base calling.
#' @slot trimmedMinQualityScore The minimum quality score of the trimmed primary sequence after base calling.
#' @slot remainingRatio The remaining sequence length ratio after trimming.
#' @slot TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @slot M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @slot M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @slot M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
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
#'        SangerRead(readFeature           = "Forward Read",
#'                         readFileName          = A_chloroticaFdReadFN,
#'                         TrimmingMethod        = "M2",
#'                         M1TrimmingCutoff      = NULL,
#'                         M2CutoffQualityScore  = 20,
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
             qualityPhredScoresRaw   = "numeric",
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
                   qualityPhredScoresRaw = numeric(0),
                   qualityPhredScores    = numeric(0),
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
                             qualityPhredScoresRaw   = qualityPhredScoresRaw,
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
