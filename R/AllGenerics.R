setClassUnion("numericORNULL", c("numeric", "NULL"))

setClassUnion("abifORNULL", c("abif", "NULL"))

setClassUnion("characterORNULL", c("character", "NULL"))

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
#' @examples
#' data(qualityReport)
#' data(sangerReadF)
#' qualityBasePlot(qualityReport)
#' qualityBasePlot(sangerReadF)
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
#' @examples
#' data(qualityReport)
#' data(sangerReadF)
#' data(sangerContig)
#' data(sangerAlignment)
#' updateQualityParam(qualityReport,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
#' updateQualityParam(sangerReadF,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
#' updateQualityParam(sangerContig,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
#' updateQualityParam(sangerAlignment,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
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
#' @examples
#' data(sangerReadF)
#' MakeBaseCalls(sangerReadF, signalRatioCutoff = 0.22)
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
#' @examples
#' data(sangerAlignment)
#' writeFastaSA(sangerAlignment)
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
#' @examples
#' data(sangerContig)
#' writeFastaSC(sangerContig)
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
#' @examples
#' data(sangerReadF)
#' writeFastaSR(sangerReadF)
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
#' @examples
#' data(sangerContig)
#' \dontrun{
#' launchAppSC(sangerContig)}
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
#' @examples
#' data(sangerAlignment)
#' \dontrun{
#' launchAppSA(sangerAlignment)}
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
#' @examples
#' data(sangerReadF)
#' \dontrun{
#' generateReportSR(sangerReadF)}
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
#' @examples
#' data(sangerContig)
#' \dontrun{
#' generateReportSC(sangerContig)}
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
#' @examples
#' data(sangerAlignment)
#' \dontrun{
#' generateReportSA(sangerAlignment)}
setGeneric("generateReportSA", function(object, outputDir = NULL,
                                        includeSangerContig = TRUE,
                                        includeSangerRead = TRUE) {
    standardGeneric("generateReportSA")
})
