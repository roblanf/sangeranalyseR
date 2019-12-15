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
#' @export
setGeneric("preQualityBasePlot", function(object,
                                          readFeature = "Forward Read") {
    standardGeneric("preQualityBasePlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerContig
### ============================================================================
#' @export
setGeneric("trimmingRatioPlot", function(object) {
    standardGeneric("trimmingRatioPlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerContig
### ============================================================================
#' @export
setGeneric("qualityBasePlot", function(object) {
    standardGeneric("qualityBasePlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerContig
### ============================================================================
#' @export
setGeneric("updateQualityParam", function(object,
                                          TrimmingMethod        = "M1",
                                          M1TrimmingCutoff      = 0.0001,
                                          M2CutoffQualityScore  = NULL,
                                          M2SlidingWindowSize   = NULL) {
    standardGeneric("updateQualityParam")
})

### ============================================================================
### Defined in SangerRead
### ============================================================================
#' @export
setGeneric("MakeBaseCalls", function(obj, signalRatioCutoff = 0.33) {
    standardGeneric("MakeBaseCalls")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' @export
setGeneric("writeFastaSA", function(obj,outputDir = NULL,
                                    compress  = FALSE,
                                    compression_level = NA,
                                    selection = "all") {
    standardGeneric("writeFastaSA")
})
#' @export
setGeneric("writeFastaSC", function(obj,outputDir = NULL,
                                    compress  = FALSE,
                                    compression_level = NA,
                                    selection = "all") {
    standardGeneric("writeFastaSC")
})
#' @export
setGeneric("writeFastaSR", function(obj,outputDir = NULL,
                                  compress  = FALSE,
                                  compression_level = NA) {
    standardGeneric("writeFastaSR")
})

### ============================================================================
### Defined in SangerContig, SangerAlignment
### ============================================================================
#' @export
setGeneric("launchAppSC", function(obj, outputDir = NULL) {
    standardGeneric("launchAppSC")
})

#' @export
setGeneric("launchAppSA", function(obj, outputDir = NULL) {
    standardGeneric("launchAppSA")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' @export
setGeneric("generateReportSR", function(obj, outputDir = NULL, ...) {
    standardGeneric("generateReportSR")
})
#' @export
setGeneric("generateReportSC", function(obj, outputDir = NULL,
                                        includeSangerRead = TRUE, ...) {
    standardGeneric("generateReportSC")
})
#' @export
setGeneric("generateReportSA", function(obj, outputDir = NULL,
                                        includeSangerContig = TRUE,
                                        includeSangerRead = TRUE) {
    standardGeneric("generateReportSA")
})


#' @export
setClassUnion("numericORNULL", c("numeric", "NULL"))
