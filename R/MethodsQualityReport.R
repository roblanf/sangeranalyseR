### ============================================================================
### Preplotting quality for each base for "QualityReport" S4 object
### ============================================================================
setMethod("preQualityBasePlot",  "QualityReport", function(object, readFeature){
    trimmedStartPos = object@trimmedStartPos
    trimmedFinishPos = object@trimmedFinishPos
    qualityPhredScores = object@qualityPhredScores
    readLen = length(qualityPhredScores)
    qualityPlotDf<- data.frame(1:length(qualityPhredScores),
                               qualityPhredScores)
    colnames(qualityPlotDf) <- c("Index", "Score")
    x <- list(
        title = "Base Pair Index"
        # titlefont = f
    )
    y <- list(
        title = "Phred Quality Score"
        # titlefont = f
    )
    p <- QualityBasePlotly(trimmedStartPos, trimmedFinishPos,
                           readLen, qualityPlotDf, x,  y)
})

### ============================================================================
### Plotting quality for each base for "QualityReport" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaRead.RDdata")
#' qualityBasePlot(A_chloroticaSingleRead@QualityReport)
setMethod("qualityBasePlot",  "QualityReport", function(object){
    plotting <- preQualityBasePlot(object)
    plotting
})

## =============================================================================
## Updating quality parameters for QualityReport object.
## =============================================================================
#' @example
#' load("data/A_chloroticaRead.RDdata")
#' QualityReport <- A_chloroticaSingleRead@QualityReport
#' trimmingRatioPlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' "@@"(QualityReport, TrimmingMethod)
#' "@@"(QualityReport, M1TrimmingCutoff)
#' "@@"(QualityReport, M2CutoffQualityScore)
#' "@@"(QualityReport, M2SlidingWindowSize)
#'
#' QualityReport <- updateQualityParam(QualityReport, "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' "@@"(QualityReport, TrimmingMethod)
#' "@@"(QualityReport, M1TrimmingCutoff)
#' "@@"(QualityReport, M2CutoffQualityScore)
#' "@@"(QualityReport, M2SlidingWindowSize)
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
