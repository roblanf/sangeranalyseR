### ============================================================================
### Plotting trimmed and remaining ratio for "SangerSingleRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead)
setMethod("trimmingRatioPlot",  "SangerSingleRead", function(object){
    QualityReportObject = object@QualityReport
    plotting <- preTrimmingRatioPlot(QualityReportObject)
    plotting
})



### ============================================================================
### Plotting quality for each base for "SangerSingleRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' qualityBasePlot(A_chloroticaSingleRead)
setMethod("qualityBasePlot",  "SangerSingleRead", function(object){
    QualityReportObject = object@QualityReport
    plotting <- preQualityBasePlot(QualityReportObject)
    plotting
})

## =============================================================================
## Updating quality parameters for SangerSingleRead object.
## =============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead)
#' qualityBasePlot(A_chloroticaSingleRead)
#' A_chloroticaSingleRead@QualityReport@TrimmingMethod
#' A_chloroticaSingleRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaSingleRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@M2SlidingWindowSize
#'
#' A_chloroticaSingleRead <- updateQualityParam(A_chloroticaSingleRead,
#'                                              "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(A_chloroticaSingleRead)
#' qualityBasePlot(A_chloroticaSingleRead)
#' A_chloroticaSingleRead@QualityReport@TrimmingMethod
#' A_chloroticaSingleRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaSingleRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@M2SlidingWindowSize
setMethod("updateQualityParam",  "SangerSingleRead",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
              ### --------------------------------------------------------------
              ### Updating SangerSingleRead quality parameters
              ### --------------------------------------------------------------
              object@QualityReport <- updateQualityParam(object@QualityReport,
                                                         TrimmingMethod,
                                                         M1TrimmingCutoff,
                                                         M2CutoffQualityScore,
                                                         M2SlidingWindowSize)
              return(object)
          })

setMethod("MakeBaseCalls", "SangerSingleRead",
          function(obj, signalRatioCutoff=.33) {
              traceMatrix <- obj@traceMatrix
              peakPosMatrixRaw <- obj@peakPosMatrixRaw
              qualityPhredScoresRaw <- obj@QualityReport@qualityPhredScoresRaw
              readFeature <- obj@readFeature
              MBCResult <- MakeBaseCallsInside (traceMatrix, peakPosMatrixRaw,
                                                qualityPhredScoresRaw,
                                                signalRatioCutoff, readFeature)
              # Length won't change, so we don't need to update !
              # obj@QualityReport@qualityPhredScores <-
              #     MBCResult[["qualityPhredScores"]]

              obj@peakPosMatrix <- MBCResult[["peakPosMatrix"]]
              obj@peakAmpMatrix <- MBCResult[["peakAmpMatrix"]]
              obj@primarySeq <- MBCResult[["primarySeq"]]
              obj@secondarySeq <- MBCResult[["secondarySeq"]]
              return(obj)
          })


setMethod("createReport", "SangerSingleRead",
          function(obj) {
})

setMethod("writeFASTA", "SangerSingleRead", function(obj, outputDir, compress,
                                                     compression_level) {
    message("Start writing '", obj@readFileName, "' to FASTA format ...")
    fastaFilename <- gsub(file_ext(basename(obj@readFileName)), "fa",
                          basename(obj@readFileName))
    outputFilename = file.path(outputDir, fastaFilename)
    writeTarget <- DNAStringSet(obj@primarySeq)
    names(writeTarget) <- basename(obj@readFileName)
    writeXStringSet(writeTarget,
                    filepath = outputFilename,
                    compress = compress,
                    compression_level = compression_level)
    message("\n >> '", outputFilename, "' is written")
})
