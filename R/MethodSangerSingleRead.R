### ============================================================================
### Plotting trimmed and remaining ratio for "SangerSingleRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead)
setMethod("trimmingRatioPlot",  "SangerSingleRead", function(object){
    QualityReportObject = object@QualityReport
    plotting <- preTrimmingRatioPlot(QualityReportObject)
    plotting
})



### ============================================================================
### Plotting quality for each base for "SangerSingleRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' qualityBasePlot(A_chloroticaSingleRead)
setMethod("qualityBasePlot",  "SangerSingleRead", function(object){
    QualityReportObject = object@QualityReport
    plotting <- preQualityBasePlot(QualityReportObject)
    plotting
})

## =============================================================================
## Updating quality parameters for SangerSingleRead object.
## =============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead)
#' qualityBasePlot(A_chloroticaSingleRead)
#' A_chloroticaSingleRead@QualityReport@TrimmingMethod
#' A_chloroticaSingleRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaSingleRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@M2SlidingWindowSize
#'
#' A_chloroticaSingleRead <- updateQualityParam(A_chloroticaSingleRead,
#'                                              "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(A_chloroticaSingleRead)
#' qualityBasePlot(A_chloroticaSingleRead)
#' A_chloroticaSingleRead@QualityReport@TrimmingMethod
#' A_chloroticaSingleRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaSingleRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@M2SlidingWindowSize
setMethod("updateQualityParam",  "SangerSingleRead",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
              ### --------------------------------------------------------------
              ### Updating SangerSingleRead quality parameters
              ### --------------------------------------------------------------
              object@QualityReport <- updateQualityParam(object@QualityReport,
                                                         TrimmingMethod,
                                                         M1TrimmingCutoff,
                                                         M2CutoffQualityScore,
                                                         M2SlidingWindowSize)
              return(object)
          })

setMethod("MakeBaseCalls", "SangerSingleRead",
          function(obj, signalRatioCutoff=.33) {
              traceMatrix <- obj@traceMatrix
              peakPosMatrixRaw <- obj@peakPosMatrixRaw
              qualityPhredScoresRaw <- obj@QualityReport@qualityPhredScoresRaw
              qualityBaseScoresRaw <- obj@QualityReport@qualityBaseScoresRaw

              MBCResult <- MakeBaseCallsInside (traceMatrix, peakPosMatrixRaw,
                                                qualityPhredScoresRaw,
                                                qualityBaseScoresRaw,
                                                signalRatioCutoff=
                                                    signalRatioCutoff)






              obj@QualityReport@qualityPhredScores <-
                  MBCResult[["qualityPhredScores"]]

              obj@peakPosMatrix <- MBCResult[["peakPosMatrix"]]
              obj@peakAmpMatrix <- MBCResult[["peakAmpMatrix"]]
              obj@primarySeqID <- MBCResult[["primarySeqID"]]
              obj@primarySeq <- MBCResult[["primarySeq"]]
              obj@secondarySeqID <- MBCResult[["secondarySeqID"]]
              obj@secondarySeq <- MBCResult[["secondarySeq"]]
              return(obj)
          })
