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
#' @param object A QualityReport or SangerRead S4 instance
#'
#' @return A quality plot.
#'
#' @exportMethod qualityBasePlot
#' @examples
#' data(qualityReportData)
#' data(sangerReadFData)
#' qualityBasePlot(qualityReportData)
#' qualityBasePlot(sangerReadFData)
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
#' @param object A QualityReport, SangerRead, SangerContig, or SangerAlignment S4 instance.
#' @param TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#' @param ... Further updateQualityParam-related parameters.
#'
#' @return A \code{QualityReport}, \code{SangerRead}, \code{SangerContig}, or \code{SangerAlignment} instance.
#'
#' @exportMethod updateQualityParam
#' @examples
#' data(qualityReportData)
#' data(sangerReadFData)
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' updateQualityParam(qualityReportData,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
#' updateQualityParam(sangerReadFData,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
#' updateQualityParam(sangerContigData,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
#' updateQualityParam(sangerAlignmentData,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)}
setGeneric("updateQualityParam", function(object,
                                          TrimmingMethod        = "M1",
                                          M1TrimmingCutoff      = 0.0001,
                                          M2CutoffQualityScore  = NULL,
                                          M2SlidingWindowSize   = NULL,
                                          ...) {
    standardGeneric("updateQualityParam")
})

### ============================================================================
### Defined in SangerRead
### ============================================================================
#' Method MakeBaseCalls
#' @name MakeBaseCalls
#' @rdname MakeBaseCalls-methods
#'
#' @param object A SangerRead S4 instance.
#' @param signalRatioCutoff The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is \code{0.33}.
#'
#' @return A \code{SangerRead} instance.
#'
#' @exportMethod MakeBaseCalls
#' @examples
#' data(sangerReadFData)
#' MakeBaseCalls(sangerReadFData, signalRatioCutoff = 0.22)
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
#' @param object A SangerAlignment S4 instance.
#' @param outputDir The output directory of generated FASTA files.
#' @param compress Like for the \code{save} function in base R, must be \code{TRUE} or \code{FALSE} (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip". This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param compression_level This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param selection This value can be \code{all}, \code{contigs_alignment}, \code{contigs_unalignment} or \code{all_reads}. It generates reads and contigs FASTA files.
#'
#' @return The output directory of FASTA files.
#'
#' @exportMethod writeFastaSA
#' @examples
#' data(sangerAlignmentData)
#' writeFastaSA(sangerAlignmentData)
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
#' @param object A SangerContig S4 instance.
#' @param outputDir The output directory of generated FASTA files.
#' @param compress Like for the \code{save} function in base R, must be \code{TRUE} or \code{FALSE} (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip". This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param compression_level This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param selection This value can be \code{all}, \code{reads_alignment}, \code{reads_unalignment} or \code{contig}. It generates reads and the contig FASTA files.
#'
#' @return The output directory of FASTA files.
#'
#' @exportMethod writeFastaSC
#' @examples
#' data(sangerContigData)
#' writeFastaSC(sangerContigData)
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
#' @param object A SangerRead S4 instance.
#' @param outputDir The output directory of the generated FASTA file.
#' @param compress Like for the \code{save} function in base R, must be \code{TRUE} or \code{FALSE} (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip". This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param compression_level This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#'
#' @return The output absolute path to the FASTA file.
#'
#' @exportMethod writeFastaSR
#' @examples
#' data(sangerReadFData)
#' writeFastaSR(sangerReadFData)
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
#' @param object A SangerContig S4 instance.
#' @param outputDir The output directory of the saved new SangerContig S4 instance.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#'     
#' @return A \code{shiny.appobj} object.
#'
#' @exportMethod launchAppSC
#' @examples
#' data(sangerContigData)
#' \dontrun{
#' launchAppSC(sangerContigData)}
setGeneric("launchAppSC", function(object, outputDir = NULL, colors = "default") {
    standardGeneric("launchAppSC")
})

#' Method launchAppSA
#' @name launchAppSA
#' @rdname launchAppSA-methods
#'
#' @param object A SangerAlignment S4 instance.
#' @param outputDir The output directory of the saved new SangerAlignment S4 instance.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#'
#' @return A \code{shiny.appobj} object.
#'
#' @exportMethod launchAppSA
#' @examples
#' data(sangerAlignmentData)
#' \dontrun{
#' launchAppSA(sangerAlignmentData)}
setGeneric("launchAppSA", function(object, outputDir = NULL, colors = "default") {
    standardGeneric("launchAppSA")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' Method readTable
#' @name readTable
#' @rdname readTable-methods
#'
#' @param object A SangerRead, SangerContig, or SangerAlignment S4 instance.
#' @param indentation The indentation for different level printing
#' @param ... Further generateReportSR-related parameters.
#'
#' @return None.
#'
#' @exportMethod readTable
#' @examples
#' data(sangerReadFData)
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' readTable(sangerReadFData)
#' readTable(sangerContigData)
#' readTable(sangerAlignmentData)
#' }
setGeneric("readTable", function(object, indentation = 0, ...) {
    standardGeneric("readTable")
})

### ============================================================================
### Defined in SangerRead, SangerContig, SangerAlignment
### ============================================================================
#' Method generateReportSR
#' @name generateReportSR
#' @rdname generateReportSR-methods
#'
#' @param object A SangerRead S4 instance.
#' @param outputDir The output directory of the generated HTML report.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#' @param ... Further generateReportSR-related parameters.
#'
#' @return The output absolute path to the SangerRead's HTML file.
#'
#' @exportMethod generateReportSR
#' @examples
#' data(sangerReadFData)
#' \dontrun{
#' generateReportSR(sangerReadFData)}
setGeneric("generateReportSR", function(object, outputDir = NULL, colors="default", ...) {
    standardGeneric("generateReportSR")
})

#' Method generateReportSC
#' @name generateReportSC
#' @rdname generateReportSC-methods
#'
#' @param object A SangerContig S4 instance.
#' @param outputDir The output directory of the generated HTML report.
#' @param includeSangerRead The parameter that decides whether to include SangerRead level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#' @param ... Further generateReportSC-related parameters.
#'
#' @return The output absolute path to the SangerContig's HTML file.
#'
#' @exportMethod generateReportSC
#' @examples
#' data(sangerContigData)
#' \dontrun{
#' generateReportSC(sangerContigData)}
setGeneric("generateReportSC", function(object, outputDir = NULL,
                                        includeSangerRead = TRUE, colors="default", ...) {
    standardGeneric("generateReportSC")
})

#' Method generateReportSA
#' @name generateReportSA
#' @rdname generateReportSA-methods
#'
#' @param object A SangerAlignment S4 instance.
#' @param outputDir The output directory of the generated HTML report.
#' @param includeSangerContig The parameter that decides whether to include SangerContig level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param includeSangerRead The parameter that decides whether to include SangerRead level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#' @param ... Further generateReportSA-related parameters.
#'
#' @return The output absolute path to the SangerAlignment's HTML file.
#'
#' @exportMethod generateReportSA
#' @examples
#' data(sangerAlignmentData)
#' \dontrun{
#' generateReportSA(sangerAlignmentData)}
setGeneric("generateReportSA", function(object, outputDir = NULL,
                                        includeSangerContig = TRUE,
                                        includeSangerRead = TRUE, colors="default", ...) {
    standardGeneric("generateReportSA")
})
