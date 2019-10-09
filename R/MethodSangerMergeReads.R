### ============================================================================
### Plotting trimmed and remaining ratio for "SangerMergeReads" S4 object
### ============================================================================
setMethod("trimmingRatioPlot",  "SangerMergeReads", function(object){
    ### ------------------------------------------------------------------------
    ### Trimmed ratio plot for forward read
    ### ------------------------------------------------------------------------
    fdQualityReportObject = object@forwardReadSangerseq@qualityReport
    fdPlotting <- preTrimmingRatioPlot(fdQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Trimmed ratio plot for reverse read
    ### ------------------------------------------------------------------------
    rvQualityReportObject = object@reverseReadSangerseq@qualityReport
    rvPlotting <- preTrimmingRatioPlot(rvQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Grid plotting for forward and reverse reads
    ### ------------------------------------------------------------------------
    grid.arrange(fdPlotting, rvPlotting, ncol=2)
})



### ============================================================================
### Plotting quality for each base for "SangerMergeReads" S4 object
### ============================================================================
setMethod("qualityBasePlot",  "SangerMergeReads", function(object){
    ### ------------------------------------------------------------------------
    ### Quality base plot for forward read
    ### ------------------------------------------------------------------------
    fdQualityReportObject = object@forwardReadSangerseq@qualityReport
    fdPlotting <- preQualityBasePlot(fdQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Grid plotting for forward and reverse reads
    ### ------------------------------------------------------------------------
    rvQualityReportObject = object@reverseReadSangerseq@qualityReport
    rvPlotting <- preQualityBasePlot(rvQualityReportObject)

    grid.arrange(fdPlotting, rvPlotting, ncol=2)
})


## ============================================================================
##
## ============================================================================
setMethod("updateQualityParam",  "SangerMergeReads",
          function(object,
                   cutoffQualityScore = 20,
                   slidingWindowSize  = 5){
    ### ------------------------------------------------------------------------
    ### Updating forward read quality parameters
    ### ------------------------------------------------------------------------
    qualityBaseScoreFD <-
        object@forwardReadSangerseq@qualityReport@qualityBaseScore
    trimmingPosFD <- inside_calculate_trimming(qualityBaseScoreFD,
                                               cutoffQualityScore,
                                               slidingWindowSize)

    object@forwardReadSangerseq@qualityReport@cutoffQualityScore <-
        cutoffQualityScore
    object@forwardReadSangerseq@qualityReport@slidingWindowSize <-
        slidingWindowSize
    object@forwardReadSangerseq@qualityReport@trimmingStartPos <-
        trimmingPosFD[1]
    object@forwardReadSangerseq@qualityReport@trimmingFinishPos <-
        trimmingPosFD[2]

    ### ------------------------------------------------------------------------
    ### Updating reverse read quality parameters
    ### ------------------------------------------------------------------------
    qualityBaseScoreRV <-
        object@reverseReadSangerseq@qualityReport@qualityBaseScore
    trimmingPosRV <- inside_calculate_trimming(qualityBaseScoreRV,
                                               cutoffQualityScore,
                                               slidingWindowSize)
    trimmingStartPosRV <- trimmingPosRV[1]
    trimmingFinishPosRV <- trimmingPosRV[2]

    object@reverseReadSangerseq@qualityReport@cutoffQualityScore <- cutoffQualityScore
    object@reverseReadSangerseq@qualityReport@slidingWindowSize <- slidingWindowSize
    object@reverseReadSangerseq@qualityReport@trimmingStartPos <-
        trimmingPosRV[1]
    object@reverseReadSangerseq@qualityReport@trimmingFinishPos <-
        trimmingPosRV[2]
    return(object)
})
