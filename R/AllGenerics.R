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

### ============================================================================
### Write FASTA
### ============================================================================
#' @export
setGeneric("writeFASTA", function(obj,outputDir = tempdir(),
                                  compress  = FALSE,
                                  compression_level = NA,
                                  selection = "all") {
    standardGeneric("writeFASTA")
})



#' @export
setClassUnion("numericORNULL", c("numeric", "NULL"))

#' #' @export
#' setClassUnion("DNAStringORNULL", c("DNAString", "NULL"))
#'
#' #' @export
#' setClassUnion("matrixORNULL", c("matrix", "NULL"))
