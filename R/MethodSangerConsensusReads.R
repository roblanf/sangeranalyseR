### ============================================================================
### Plotting trimmed and remaining ratio for "SangerConsensusRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaConsensusRead.RDdata")
#' trimmingRatioPlot(A_chloroticaConsensusRead)
setMethod("trimmingRatioPlot",  "SangerConsensusRead", function(object){
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
### Plotting quality for each base for "SangerConsensusRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaConsensusRead.RDdata")
#' qualityBasePlot(A_chloroticaConsensusRead)
setMethod("qualityBasePlot",  "SangerConsensusRead", function(object){
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
## Updating quality parameters for SangerConsensusRead object.
## ============================================================================
#' @example
#' load("data/A_chloroticaConsensusRead.RDdata")
#' trimmingRatioPlot(A_chloroticaConsensusRead)
#' qualityBasePlot(A_chloroticaConsensusRead)
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@TrimmingMethod
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M1TrimmingCutoff
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2CutoffQualityScore
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2SlidingWindowSize
#'
#' A_chloroticaConsensusRead <-
#'                 updateQualityParam(A_chloroticaConsensusRead,
#'                                    "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(A_chloroticaConsensusRead)
#' qualityBasePlot(A_chloroticaConsensusRead)
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@TrimmingMethod
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M1TrimmingCutoff
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2CutoffQualityScore
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2SlidingWindowSize
setMethod("updateQualityParam",  "SangerConsensusRead",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
    ### ------------------------------------------------------------------------
    ### Updating forward read quality parameters
    ### ------------------------------------------------------------------------
    # qualityBaseScoreFD <-
    #     object@forwardReadSangerseq@QualityReport@qualityBaseScore
    object@forwardReadSangerseq <-
        updateQualityParam(object@forwardReadSangerseq,
                           TrimmingMethod,
                           M1TrimmingCutoff,
                           M2CutoffQualityScore,
                           M2SlidingWindowSize)

    ### ------------------------------------------------------------------------
    ### Updating reverse read quality parameters
    ### ------------------------------------------------------------------------
    object@reverseReadSangerseq <-
        updateQualityParam(object@reverseReadSangerseq,
                           TrimmingMethod,
                           M1TrimmingCutoff,
                           M2CutoffQualityScore,
                           M2SlidingWindowSize)
    return(object)
})
