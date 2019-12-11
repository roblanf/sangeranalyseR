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
### Defined in QualityReport, SangerRead, SangerMergeReads
### ============================================================================
#' @export
setGeneric("trimmingRatioPlot", function(object) {
    standardGeneric("trimmingRatioPlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerMergeReads
### ============================================================================
#' @export
setGeneric("qualityBasePlot", function(object) {
    standardGeneric("qualityBasePlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerMergeReads
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
### Create report
### ============================================================================
#' @export
setGeneric("createReport", function(obj) {
    standardGeneric("createReport")
})

### ============================================================================
### Write FASTA
### ============================================================================
#' @export
setGeneric("writeFASTA", function(obj,outputDir = tempdir(),
                                  compress  = FALSE,
                                  compression_level = NA) {
    standardGeneric("writeFASTA")
})
#' @export
setGeneric("writeFASTA", function(obj,outputDir = tempdir(),
                                  compress  = FALSE,
                                  compression_level = NA,
                                  selection = "all") {
    standardGeneric("writeFASTA")
})


### ============================================================================
### Run Shiny app
### ============================================================================
#' @export
setGeneric("launchAppSangerContig", function(obj, outputDir = NULL) {
    standardGeneric("launchAppSangerContig")
})

#' @export
setGeneric("launchAppSangerAlignment", function(obj, outputDir = NULL) {
    standardGeneric("launchAppSangerAlignment")
})

#' @export
setGeneric("generateReport", function(obj, outputDir = NULL) {
    standardGeneric("generateReport")
})

#' @export
setClassUnion("numericORNULL", c("numeric", "NULL"))

#' #' @export
#' setClassUnion("DNAStringORNULL", c("DNAString", "NULL"))
#'
#' #' @export
#' setClassUnion("matrixORNULL", c("matrix", "NULL"))
