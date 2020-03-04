### ============================================================================
### Defined in QualityReport
### ============================================================================
#' Method preTrimmingRatioPlot
#' @name preTrimmingRatioPlot
#' @rdname preTrimmingRatioPlot-methods
#' @exportMethod preTrimmingRatioPlot
setGeneric("preTrimmingRatioPlot", function(object) {
    standardGeneric("preTrimmingRatioPlot")
})

### ============================================================================
### Defined in QualityReport
### ============================================================================
#' Method preQualityBasePlot
#' @name preQualityBasePlot
#' @rdname preQualityBasePlot-methods
#' @exportMethod preQualityBasePlot
setGeneric("preQualityBasePlot", function(object) {
    standardGeneric("preQualityBasePlot")
})

#' Method preQualityBasePlot
#' @name preQualityBasePlot
#' @rdname preQualityBasePlot-methods
#' @exportMethod preQualityBasePlot
setGeneric("preQualityBasePlot", function(object,
                                          readFeature = "Forward Read") {
    standardGeneric("preQualityBasePlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerContig
### ============================================================================
#' Method trimmingRatioPlot
#' @name trimmingRatioPlot
#' @rdname trimmingRatioPlot-methods
#' @exportMethod trimmingRatioPlot
setGeneric("trimmingRatioPlot", function(object) {
    standardGeneric("trimmingRatioPlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerContig
### ============================================================================
#' Method qualityBasePlot
#' @name qualityBasePlot
#' @rdname qualityBasePlot-methods
#' @exportMethod qualityBasePlot
setGeneric("qualityBasePlot", function(object) {
    standardGeneric("qualityBasePlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerContig
### ============================================================================
#' Method updateQualityParam
#' @name updateQualityParam
#' @rdname updateQualityParam-methods
#' @exportMethod updateQualityParam
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
#' Method MakeBaseCalls
#' @name MakeBaseCalls
#' @rdname MakeBaseCalls-methods
#' @exportMethod MakeBaseCalls
setGeneric("MakeBaseCalls", function(obj, signalRatioCutoff = 0.33) {
    standardGeneric("MakeBaseCalls")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' Method writeFastaSA
#' @name writeFastaSA
#' @rdname writeFastaSA-methods
#' @exportMethod writeFastaSA
setGeneric("writeFastaSA", function(obj,outputDir = NULL,
                                    compress  = FALSE,
                                    compression_level = NA,
                                    selection = "all") {
    standardGeneric("writeFastaSA")
})
#' Method writeFastaSC
#' @name writeFastaSC
#' @rdname writeFastaSC-methods
#' @exportMethod writeFastaSC
setGeneric("writeFastaSC", function(obj,outputDir = NULL,
                                    compress  = FALSE,
                                    compression_level = NA,
                                    selection = "all") {
    standardGeneric("writeFastaSC")
})
#' Method writeFastaSR
#' @name writeFastaSR
#' @rdname writeFastaSR-methods
#' @exportMethod writeFastaSR
setGeneric("writeFastaSR", function(obj,outputDir = NULL,
                                  compress  = FALSE,
                                  compression_level = NA) {
    standardGeneric("writeFastaSR")
})

### ============================================================================
### Defined in SangerContig, SangerAlignment
### ============================================================================
#' Method launchAppSC
#' @name launchAppSC
#' @rdname launchAppSC-methods
#' @exportMethod launchAppSC
setGeneric("launchAppSC", function(obj, outputDir = NULL) {
    standardGeneric("launchAppSC")
})

#' Method launchAppSA
#' @name launchAppSA
#' @rdname launchAppSA-methods
#' @exportMethod launchAppSA
setGeneric("launchAppSA", function(obj, outputDir = NULL) {
    standardGeneric("launchAppSA")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' Method generateReportSR
#' @name generateReportSR
#' @rdname generateReportSR-methods
#' @exportMethod generateReportSR
setGeneric("generateReportSR", function(obj, outputDir = NULL, ...) {
    standardGeneric("generateReportSR")
})
#' Method generateReportSC
#' @name generateReportSC
#' @rdname generateReportSC-methods
#' @exportMethod generateReportSC
setGeneric("generateReportSC", function(obj, outputDir = NULL,
                                        includeSangerRead = TRUE, ...) {
    standardGeneric("generateReportSC")
})
#' Method generateReportSA
#' @name generateReportSA
#' @rdname generateReportSA-methods
#' @exportMethod generateReportSA
setGeneric("generateReportSA", function(obj, outputDir = NULL,
                                        includeSangerContig = TRUE,
                                        includeSangerRead = TRUE) {
    standardGeneric("generateReportSA")
})

#' ClassUnion numericORNULL
#' @name numericORNULL
#' @rdname numericORNULL-ClassUnion
#' @exportClass numericORNULL
setClassUnion("numericORNULL", c("numeric", "NULL"))

#' ClassUnion abifORNULL
#' @name abifORNULL
#' @rdname abifORNULL-ClassUnion
#' @exportClass abifORNULL
setClassUnion("abifORNULL", c("abif", "NULL"))

#' ClassUnion characterORNULL
#' @name characterORNULL
#' @rdname characterORNULL-ClassUnion
#' @exportClass characterORNULL
setClassUnion("characterORNULL", c("character", "NULL"))

