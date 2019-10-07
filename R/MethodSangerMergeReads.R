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
