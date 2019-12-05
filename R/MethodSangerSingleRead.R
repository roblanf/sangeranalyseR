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
              peakPosMatrix <- obj@peakPosMatrix
              QualityReport <- obj@QualityReport
              MBCResult <- MakeBaseCallsInside (traceMatrix, peakPosMatrix,
                                                QualityReport,
                                                signalRatioCutoff=.33)
              obj@QualityReport <- MBCResult[["QualityReport"]]
              obj@peakPosMatrixBC <- MBCResult[["peakPosMatrixBC"]]
              obj@peakAmpMatrixBC <- MBCResult[["peakAmpMatrixBC"]]
              obj@primarySeqID <- MBCResult[["primarySeqID"]]
              obj@primarySeqBC <- MBCResult[["primarySeqBC"]]
              obj@secondarySeqID <- MBCResult[["secondarySeqID"]]
              obj@secondarySeqBC <- MBCResult[["secondarySeqBC"]]
              return(obj)
          })
