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
#' A_chloroticaSingleRead@QualityReport@cutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@slidingWindowSize
#'
#' A_chloroticaSingleRead <- updateQualityParam(A_chloroticaSingleRead, 20L, 5L)
#'
#' trimmingRatioPlot(A_chloroticaSingleRead)
#' qualityBasePlot(A_chloroticaSingleRead)
#' A_chloroticaSingleRead@QualityReport@cutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@slidingWindowSize
setMethod("updateQualityParam",  "SangerSingleRead",
          function(object,
                   cutoffQualityScore = 20L,
                   slidingWindowSize  = 5L){
              ### --------------------------------------------------------------
              ### Updating SangerSingleRead quality parameters
              ### --------------------------------------------------------------
              object@QualityReport <- updateQualityParam(object@QualityReport,
                                                         cutoffQualityScore,
                                                         slidingWindowSize)
              return(object)
          })
