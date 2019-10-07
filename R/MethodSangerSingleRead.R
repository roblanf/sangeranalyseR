### ============================================================================
### Plotting trimmed and remaining ratio for "sangerSingleRead" S4 object
### ============================================================================
setMethod("trimmingRatioPlot",  "sangerSingleRead", function(object){
    qualityReportObject = object@qualityReport
    plotting <- preTrimmingRatioPlot(qualityReportObject)
    plotting
})



### ============================================================================
### Plotting quality for each base for "sangerSingleRead" S4 object
### ============================================================================
setMethod("qualityBasePlot",  "sangerSingleRead", function(object){
    qualityReportObject = object@qualityReport
    plotting <- preQualityBasePlot(qualityReportObject)
    plotting
})
