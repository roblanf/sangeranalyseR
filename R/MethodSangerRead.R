### ============================================================================
### Plotting quality for each base for "SangerRead" S4 object
### ============================================================================
#' @title qualityBasePlot
#' @name SangerRead-class-qualityBasePlot
#' @rdname SangerRead-class-qualityBasePlot
#' @examples
#' load("data/sangerRead.RData")
#' qualityBasePlot(sangerRead)
setMethod("qualityBasePlot",  "SangerRead", function(object){
    if (object@inputSource == "ABIF") {
        QualityReportObject = object@QualityReport
        plotting <- preQualityBasePlot(QualityReportObject)
        plotting
    } else if (object@inputSource == "FASTA") {
        message("SangerRead with 'FASTA' inputSource ",
                "cannot plot quality plots")
    }
})

## =============================================================================
## Updating quality parameters for SangerRead object.
## =============================================================================
#' @title updateQualityParam
#' @name SangerRead-class-updateQualityParam
#' @rdname SangerRead-class-updateQualityParam
#' @examples
#' load("data/sangerRead.RData")
#' updateQualityParam(sangerRead,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
setMethod("updateQualityParam",  "SangerRead",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
    if (object@inputSource == "ABIF") {
        ### --------------------------------------------------------------------
        ### Updating SangerRead quality parameters
        ###   Trimming parameters is checked in 'QualityReport' method
        ### --------------------------------------------------------------------
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
    } else if (object@inputSource == "FASTA") {
        message("SangerRead with 'FASTA' inputSource ",
                "cannot update quality parameters")
    }
})

## =============================================================================
## Base calling for primarySeq in SangerRead
## =============================================================================
#' @title MakeBaseCalls
#' @name SangerRead-class-MakeBaseCalls
#' @rdname SangerRead-class-MakeBaseCalls
#' @examples
#' load("data/sangerRead.RData")
#' MakeBaseCalls(sangerRead, signalRatioCutoff = 0.22)
setMethod("MakeBaseCalls", "SangerRead", function(obj, signalRatioCutoff) {
    if (obj@inputSource == "ABIF") {
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
    } else if (obj@inputSource == "FASTA") {
        message("SangerRead with 'FASTA' inputSource cannot do base calling")
    }
})

## =============================================================================
## Writing primary sequence into FASTA format
## =============================================================================
#' @title writeFastaSR
#' @name SangerRead-class-writeFastaSR
#' @rdname SangerRead-class-writeFastaSR
#' @examples
#' load("data/sangerRead.RData")
#' writeFastaSR(sangerRead, "/Users/chaokuan-hao/Desktop/sangeranalyseR_fasta/SangerRead")
setMethod("writeFastaSR", "SangerRead", function(obj, outputDir, compress,
                                                 compression_level) {
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
        message(">>> outputDir : ", outputDir)
    }
    message("Start writing '", obj@readFileName, "' to FASTA format ...")
    fastaFilename <- gsub(file_ext(basename(obj@readFileName)), "fa",
                          basename(obj@readFileName))
    outputFilename <- file.path(outputDir, fastaFilename)
    ### ------------------------------------------------------------------------
    ### Add trimming in ABIF file format
    ### ------------------------------------------------------------------------
    if(obj@inputSource == "ABIF") {
        trimmedStartPos <- obj@QualityReport@trimmedStartPos
        trimmedFinishPos <- obj@QualityReport@trimmedFinishPos
        targetSeq <- DNAString(substr(as.character(obj@primarySeq),
                                      trimmedStartPos+1, trimmedFinishPos))
    } else if (obj@inputSource == "FASTA") {
        targetSeq <- obj@primarySeq
    }
    writeTarget <- DNAStringSet(targetSeq)
    names(writeTarget) <- basename(obj@readFileName)
    writeXStringSet(writeTarget,
                    filepath = outputFilename,
                    compress = compress,
                    compression_level = compression_level)
    message("\n >> '", outputFilename, "' is written")
    return(outputFilename)
})

## =============================================================================
## Generating report for SangerRead
## =============================================================================
#' @title generateReportSR
#' @name SangerRead-class-generateReportSR
#' @rdname SangerRead-class-generateReportSR
#' @examples
#' load("data/sangerRead.RData")
#' generateReportSR(sangerRead)
setMethod("generateReportSR", "SangerRead",
          function(obj, outputDir,
                   navigationContigFN = NULL, navigationAlignmentFN = NULL) {
    # Another Rmd for SangerRead with FASTA file source
    ### ------------------------------------------------------------------------
    ### Make sure the input directory is not NULL
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
        message(">>> outputDir : ", outputDir)
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
