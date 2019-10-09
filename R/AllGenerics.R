### ============================================================================
### Defined in QualityReport
### ============================================================================
#' @export
setGeneric("preTrimmingRatioPlot", function(object) {
    standardGeneric("preTrimmingRatioPlot")
})

### ============================================================================
### Defined in QualityReport
### ============================================================================
#' @export
setGeneric("preQualityBasePlot", function(object) {
    standardGeneric("preQualityBasePlot")
})


### ============================================================================
### Defined in QualityReport, sangerSingleRead, sangerMergeReads
### ============================================================================
#' @export
setGeneric("trimmingRatioPlot", function(object) {
    standardGeneric("trimmingRatioPlot")
})

### ============================================================================
### Defined in QualityReport, sangerSingleRead, sangerMergeReads
### ============================================================================
#' @export
setGeneric("qualityBasePlot", function(object) {
    standardGeneric("qualityBasePlot")
})

### ============================================================================
### Defined in QualityReport, sangerSingleRead, sangerMergeReads
### ============================================================================
#' @export
setGeneric("updateQualityParam", function(object,
                                          cutoffQualityScore = 20L,
                                          slidingWindowSize  = 5L) {
    standardGeneric("updateQualityParam")
})


