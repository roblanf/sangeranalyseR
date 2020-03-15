#' @description ClassUnion numericORNULL
#' @title S4 Class Union slot types
#' @name numericORNULL
#' @exportClass numericORNULL
#' @rdname slot-type
setClassUnion("numericORNULL", c("numeric", "NULL"))

#' @description ClassUnion abifORNULL
#' @title S4 Class Union slot types
#' @name abifORNULL
#' @exportClass abifORNULL
#' @rdname slot-type
setClassUnion("abifORNULL", c("abif", "NULL"))

#' @description  ClassUnion characterORNULL
#' @title S4 Class Union slot types
#' @name characterORNULL
#' @exportClass characterORNULL
#' @rdname slot-type
setClassUnion("characterORNULL", c("character", "NULL"))

### ============================================================================
### Defined in QualityReport
### ============================================================================
#' Method preTrimmingRatioPlot
#' @name preTrimmingRatioPlot
#' @rdname preTrimmingRatioPlot-methods
#'
#' @param object object
#'
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
#'
#' @param object object
#'
#' @exportMethod preQualityBasePlot
setGeneric("preQualityBasePlot", function(object) {
    standardGeneric("preQualityBasePlot")
})

### ============================================================================
### Defined in QualityReport, SangerRead, SangerContig
### ============================================================================
#' Method trimmingRatioPlot
#' @name trimmingRatioPlot
#' @rdname trimmingRatioPlot-methods
#'
#' @param object object
#'
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
#'
#' @param object object
#'
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
#'
#' @param object object
#' @param TrimmingMethod TrimmingMethod
#' @param M1TrimmingCutoff M1TrimmingCutoff
#' @param M2CutoffQualityScore M2CutoffQualityScore
#' @param M2SlidingWindowSize M2SlidingWindowSize
#'
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
#'
#' @param object object
#' @param signalRatioCutoff signalRatioCutoff
#'
#' @exportMethod MakeBaseCalls
setGeneric("MakeBaseCalls", function(object, signalRatioCutoff = 0.33) {
    standardGeneric("MakeBaseCalls")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' Method writeFastaSA
#' @name writeFastaSA
#' @rdname writeFastaSA-methods
#'
#' @param object object
#' @param outputDir outputDir
#' @param compress compress
#' @param compression_level compression_level
#' @param selection selection
#'
#' @exportMethod writeFastaSA
setGeneric("writeFastaSA", function(object,outputDir = NULL,
                                    compress  = FALSE,
                                    compression_level = NA,
                                    selection = "all") {
    standardGeneric("writeFastaSA")
})
#' Method writeFastaSC
#' @name writeFastaSC
#' @rdname writeFastaSC-methods
#'
#' @param object object
#' @param outputDir outputDir
#' @param compress compress
#' @param compression_level compression_level
#' @param selection selection
#'
#' @exportMethod writeFastaSC
setGeneric("writeFastaSC", function(object,outputDir = NULL,
                                    compress  = FALSE,
                                    compression_level = NA,
                                    selection = "all") {
    standardGeneric("writeFastaSC")
})
#' Method writeFastaSR
#' @name writeFastaSR
#' @rdname writeFastaSR-methods
#'
#' @param object object
#' @param outputDir outputDir
#' @param compress compress
#' @param compression_level compression_level
#'
#' @exportMethod writeFastaSR
setGeneric("writeFastaSR", function(object,outputDir = NULL,
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
#'
#' @param object object
#' @param outputDir outputDir
#'
#' @exportMethod launchAppSC
setGeneric("launchAppSC", function(object, outputDir = NULL) {
    standardGeneric("launchAppSC")
})

#' Method launchAppSA
#' @name launchAppSA
#' @rdname launchAppSA-methods
#'
#' @param object object
#' @param outputDir outputDir
#'
#' @exportMethod launchAppSA
setGeneric("launchAppSA", function(object, outputDir = NULL) {
    standardGeneric("launchAppSA")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' Method generateReportSR
#' @name generateReportSR
#' @rdname generateReportSR-methods
#'
#' @param object object
#' @param outputDir outputDir
#' @param ... ...
#'
#' @exportMethod generateReportSR
setGeneric("generateReportSR", function(object, outputDir = NULL, ...) {
    standardGeneric("generateReportSR")
})
#' Method generateReportSC
#' @name generateReportSC
#' @rdname generateReportSC-methods
#'
#' @param object object
#' @param outputDir outputDir
#' @param includeSangerRead includeSangerRead
#' @param ... ...
#'
#' @exportMethod generateReportSC
setGeneric("generateReportSC", function(object, outputDir = NULL,
                                        includeSangerRead = TRUE, ...) {
    standardGeneric("generateReportSC")
})
#' Method generateReportSA
#' @name generateReportSA
#' @rdname generateReportSA-methods
#'
#' @param object object
#' @param outputDir outputDir
#' @param includeSangerContig includeSangerContig
#' @param includeSangerRead includeSangerRead
#'
#' @exportMethod generateReportSA
setGeneric("generateReportSA", function(object, outputDir = NULL,
                                        includeSangerContig = TRUE,
                                        includeSangerRead = TRUE) {
    standardGeneric("generateReportSA")
})
