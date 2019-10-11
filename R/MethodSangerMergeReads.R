### ============================================================================
### Plotting trimmed and remaining ratio for "SangerMergeReads" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaMergeReads.RDdata")
#' trimmingRatioPlot(A_chloroticaMergeReads)
setMethod("trimmingRatioPlot",  "SangerMergeReads", function(object){
    ### ------------------------------------------------------------------------
    ### Trimmed ratio plot for forward read
    ### ------------------------------------------------------------------------
    fdQualityReportObject = object@forwardReadSangerseq@QualityReport
    fdPlotting <- preTrimmingRatioPlot(fdQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Trimmed ratio plot for reverse read
    ### ------------------------------------------------------------------------
    rvQualityReportObject = object@reverseReadSangerseq@QualityReport
    rvPlotting <- preTrimmingRatioPlot(rvQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Grid plotting for forward and reverse reads
    ### ------------------------------------------------------------------------
    grid.arrange(fdPlotting, rvPlotting, ncol=2)
})



### ============================================================================
### Plotting quality for each base for "SangerMergeReads" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaMergeReads.RDdata")
#' qualityBasePlot(A_chloroticaMergeReads)
setMethod("qualityBasePlot",  "SangerMergeReads", function(object){
    ### ------------------------------------------------------------------------
    ### Quality base plot for forward read
    ### ------------------------------------------------------------------------
    fdQualityReportObject = object@forwardReadSangerseq@QualityReport
    fdPlotting <- preQualityBasePlot(fdQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Grid plotting for forward and reverse reads
    ### ------------------------------------------------------------------------
    rvQualityReportObject = object@reverseReadSangerseq@QualityReport
    rvPlotting <- preQualityBasePlot(rvQualityReportObject)

    grid.arrange(fdPlotting, rvPlotting, ncol=2)
})


## ============================================================================
## Updating quality parameters for SangerMergeReads object.
## ============================================================================
#' @example
#' load("data/A_chloroticaMergeReads.RDdata")
#' trimmingRatioPlot(A_chloroticaMergeReads)
#' qualityBasePlot(A_chloroticaMergeReads)
#' A_chloroticaMergeReads@forwardReadSangerseq@QualityReport@cutoffQualityScore
#' A_chloroticaMergeReads@forwardReadSangerseq@QualityReport@slidingWindowSize
#'
#' A_chloroticaMergeReads <- updateQualityParam(A_chloroticaMergeReads, 20L, 5L)
#'
#' trimmingRatioPlot(A_chloroticaMergeReads)
#' qualityBasePlot(A_chloroticaMergeReads)
#' A_chloroticaMergeReads@forwardReadSangerseq@QualityReport@cutoffQualityScore
#' A_chloroticaMergeReads@forwardReadSangerseq@QualityReport@slidingWindowSize
setMethod("updateQualityParam",  "SangerMergeReads",
          function(object,
                   cutoffQualityScore = 20L,
                   slidingWindowSize  = 5L){
    ### ------------------------------------------------------------------------
    ### Updating forward read quality parameters
    ### ------------------------------------------------------------------------
    # qualityBaseScoreFD <-
    #     object@forwardReadSangerseq@QualityReport@qualityBaseScore
    object@forwardReadSangerseq <-
        updateQualityParam(object@forwardReadSangerseq,
                           cutoffQualityScore,
                           slidingWindowSize)

    ### ------------------------------------------------------------------------
    ### Updating reverse read quality parameters
    ### ------------------------------------------------------------------------
    object@reverseReadSangerseq <-
        updateQualityParam(object@reverseReadSangerseq,
                           cutoffQualityScore,
                           slidingWindowSize)
    return(object)
})
