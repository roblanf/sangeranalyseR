### ============================================================================
### Plotting quality for each base for "SangerRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaRead.RDdata")
#' qualityBasePlot(A_chloroticaRead)
setMethod("qualityBasePlot",  "SangerRead", function(object){
    QualityReportObject = object@QualityReport
    plotting <- preQualityBasePlot(QualityReportObject)
    plotting
})

## =============================================================================
## Updating quality parameters for SangerRead object.
## =============================================================================
#' @example
#' load("data/A_chloroticaRead.RDdata")
#' trimmingRatioPlot(A_chloroticaRead)
#' qualityBasePlot(A_chloroticaRead)
#' A_chloroticaRead@QualityReport@TrimmingMethod
#' A_chloroticaRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaRead@QualityReport@M2SlidingWindowSize
#'
#' A_chloroticaRead <- updateQualityParam(A_chloroticaRead,
#'                                              "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(A_chloroticaRead)
#' qualityBasePlot(A_chloroticaRead)
#' A_chloroticaRead@QualityReport@TrimmingMethod
#' A_chloroticaRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaRead@QualityReport@M2SlidingWindowSize
setMethod("updateQualityParam",  "SangerRead",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
              ### --------------------------------------------------------------
              ### Updating SangerRead quality parameters
              ###   Trimming parameters is checked in 'QualityReport' method
              ### --------------------------------------------------------------
              errors <- character()
              errors <- checkTrimParam(TrimmingMethod,
                                       M1TrimmingCutoff,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize,
                                       errors)
              if (length(errors) == 0) {
                  object@QualityReport <-
                      updateQualityParam(object@QualityReport,
                                         TrimmingMethod,
                                         M1TrimmingCutoff,
                                         M2CutoffQualityScore,
                                         M2SlidingWindowSize)
                  return(object)
              } else {
                  stop(errors)
              }
          })

setMethod("MakeBaseCalls", "SangerRead",
          function(obj, signalRatioCutoff) {
              errors <- character(0)
              errors <- checkSignalRatioCutoff(signalRatioCutoff, errors)
              if (length(errors) == 0) {
                  traceMatrix <- obj@traceMatrix
                  peakPosMatrixRaw <- obj@peakPosMatrixRaw
                  qualityPhredScoresRaw <- obj@QualityReport@qualityPhredScoresRaw
                  readFeature <- obj@readFeature
                  MBCResult <-
                      MakeBaseCallsInside (traceMatrix, peakPosMatrixRaw,
                                           qualityPhredScoresRaw,
                                           signalRatioCutoff, readFeature)
                  obj@peakPosMatrix <- MBCResult[["peakPosMatrix"]]
                  obj@peakAmpMatrix <- MBCResult[["peakAmpMatrix"]]
                  obj@primarySeq <- MBCResult[["primarySeq"]]
                  obj@secondarySeq <- MBCResult[["secondarySeq"]]

                  AASeqResult <-
                      calculateAASeq (obj@primarySeq, obj@geneticCode)
                  obj@primaryAASeqS1 <- AASeqResult[["primaryAASeqS1"]]
                  obj@primaryAASeqS2 <- AASeqResult[["primaryAASeqS2"]]
                  obj@primaryAASeqS3 <- AASeqResult[["primaryAASeqS3"]]
                  obj@ChromatogramParam@signalRatioCutoff <- signalRatioCutoff
                  return(obj)
              } else {
                  stop(errors)
              }
})

setMethod("writeFastaSR", "SangerRead", function(obj, outputDir, compress,
                                                 compression_level) {
    message("Start writing '", obj@readFileName, "' to FASTA format ...")
    fastaFilename <- gsub(file_ext(basename(obj@readFileName)), "fa",
                          basename(obj@readFileName))
    outputFilename <- file.path(outputDir, fastaFilename)
    writeTarget <- DNAStringSet(obj@primarySeq)
    names(writeTarget) <- basename(obj@readFileName)
    writeXStringSet(writeTarget,
                    filepath = outputFilename,
                    compress = compress,
                    compression_level = compression_level)
    message("\n >> '", outputFilename, "' is written")
    return(outputFilename)
})
setMethod("generateReportSR", "SangerRead",
          function(obj, outputDir,
                   navigationContigFN = NULL, navigationAlignmentFN = NULL) {
    ### ------------------------------------------------------------------------
    ### Make sure the input directory is not NULL
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    readName <- sub('\\.ab1$', '', basename(obj@readFileName))
    outputDirSR <- file.path(outputDir, readName)
    ### ------------------------------------------------------------------------
    ### Make sure the directory is exist (SangerRead level)
    ### ------------------------------------------------------------------------
    if (!dir.exists(outputDirSR)) {
        suppressWarnings(dir.create(outputDirSR, recursive = TRUE))
    }
    rootDir <- system.file(package = "sangeranalyseR")
    originRmd <- file.path(rootDir, "rmd", "SangerRead_Report.Rmd")
    outputHtml <- file.path(outputDirSR, "SangerRead_Report.html")
    res <- render(input = originRmd,
                  output_dir = outputDirSR,
                  params = list(SangerRead = obj,
                                outputDir = outputDirSR,
                                navigationContigFN = navigationContigFN,
                                navigationAlignmentFN = navigationAlignmentFN))
    return(outputHtml)
})
