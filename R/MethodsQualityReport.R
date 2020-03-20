preQualityBasePlot <- function(object) {
    trimmedStartPos = object@trimmedStartPos
    trimmedFinishPos = object@trimmedFinishPos
    qualityPhredScores = object@qualityPhredScores
    readLen = length(qualityPhredScores)
    qualityPlotDf<- data.frame(1:length(qualityPhredScores),
                               qualityPhredScores)
    colnames(qualityPlotDf) <- c("Index", "Score")
    x <- list(
        title = "Base Pair Index"
    )
    y <- list(
        title = "Phred Quality Score"
    )
    p <- QualityBasePlotly(trimmedStartPos, trimmedFinishPos,
                           readLen, qualityPlotDf, x,  y)
    return(p)
}

### ============================================================================
### Plotting quality of each base pair for "QualityReport" S4 object
### ============================================================================
#' A QualityReport method which creates quality base interactive plot.
#'
#' @title qualityBasePlot
#' @name QualityReport-class-qualityBasePlot
#' @aliases qualityBasePlot,QualityReport-method
#'
#' @param object A QualityReport S4 instance.
#'
#' @examples
#' data("qualityReportData")
#' qualityBasePlot(qualityReportData)
setMethod("qualityBasePlot",  "QualityReport", function(object){
    plotting <- preQualityBasePlot(object)
    plotting
})

## =============================================================================
## Updating quality parameters for QualityReport object.
## =============================================================================
#' A QualityReport method which updates quality base interactive plot.
#'
#' @title updateQualityParam
#' @name QualityReport-class-updateQualityParam
#' @aliases updateQualityParam,QualityReport-method
#'
#' @param object A QualityReport S4 instance.
#' @param TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#'
#' @examples
#' data("qualityReportData")
#' updateQualityParam(qualityReportData,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 30,
#'                    M2SlidingWindowSize    = 15)
setMethod("updateQualityParam",  "QualityReport",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL) {
              ##### ------------------------------------------------------------
              ##### Input parameter prechecking for TrimmingMethod.
              ##### ------------------------------------------------------------
              errors <- character()
              errors <- checkTrimParam(TrimmingMethod,
                                       M1TrimmingCutoff,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize,
                                       errors)
              if (length(errors) == 0) {
                  qualityBaseScores <- object@qualityBaseScores
                  qualityPhredScores <- object@qualityPhredScores
                  if (TrimmingMethod == "M1") {
                      trimmingPos <-
                          M1inside_calculate_trimming(qualityPhredScores,
                                                      qualityBaseScores,
                                                      M1TrimmingCutoff)
                  } else if (TrimmingMethod == "M2") {
                      trimmingPos <-
                          M2inside_calculate_trimming(qualityPhredScores,
                                                      M2CutoffQualityScore,
                                                      M2SlidingWindowSize)
                  }
                  ### ----------------------------------------------------------
                  ### Updating QualityReport quality parameters
                  ### ----------------------------------------------------------
                  object@TrimmingMethod <- TrimmingMethod
                  object@M1TrimmingCutoff <- M1TrimmingCutoff
                  object@M2CutoffQualityScore <- M2CutoffQualityScore
                  object@M2SlidingWindowSize <- M2SlidingWindowSize

                  object@rawSeqLength <-
                      trimmingPos[["rawSeqLength"]]
                  object@rawMeanQualityScore <-
                      trimmingPos[["rawMeanQualityScore"]]
                  object@rawMinQualityScore <-
                      trimmingPos[["rawMinQualityScore"]]
                  object@trimmedStartPos <-
                      trimmingPos[["trimmedStartPos"]]
                  object@trimmedFinishPos <-
                      trimmingPos[["trimmedFinishPos"]]
                  object@trimmedSeqLength <-
                      trimmingPos[["trimmedSeqLength"]]
                  object@trimmedMeanQualityScore <-
                      trimmingPos[["trimmedMeanQualityScore"]]
                  object@trimmedMinQualityScore <-
                      trimmingPos[["trimmedMinQualityScore"]]
                  object@remainingRatio <-
                      trimmingPos[["remainingRatio"]]
                  return(object)
              } else {
                  stop(errors)
              }
})
