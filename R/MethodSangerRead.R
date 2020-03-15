### ============================================================================
### Plotting quality for each base for "SangerRead" S4 object
### ============================================================================
#' A SangerRead method which creates quality base interactive plot.
#'
#' @title qualityBasePlot
#' @name SangerRead-class-qualityBasePlot
#' @aliases qualityBasePlot,SangerRead-method
#'
#' @docType methods
#' @examples
#' \dontrun{qualityBasePlot(sangerReadF)}
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
#' A SangerRead method which updates QualityReport parameter inside the SangerRead.
#'
#' @title updateQualityParam
#' @name SangerRead-class-updateQualityParam
#' @aliases updateQualityParam,SangerRead-method
#'
#' @docType methods
#' @examples
#' \dontrun{updateQualityParam(sangerReadF,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)}
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
            AASeqResult    <- calculateAASeq (object@primarySeq,
                                              object@QualityReport@trimmedStartPos,
                                              object@QualityReport@trimmedFinishPos,
                                              object@geneticCode)
            object@primaryAASeqS1 <- AASeqResult[["primaryAASeqS1"]]
            object@primaryAASeqS2 <- AASeqResult[["primaryAASeqS2"]]
            object@primaryAASeqS3 <- AASeqResult[["primaryAASeqS3"]]
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
#' A SangerRead method which does base calling on SangerRead instance
#'
#' @title MakeBaseCalls
#' @name SangerRead-class-MakeBaseCalls
#' @aliases MakeBaseCalls,SangerRead-method
#'
#' @docType methods
#'
#' @examples
#' \dontrun{MakeBaseCalls(sangerReadF, signalRatioCutoff = 0.22)}
setMethod("MakeBaseCalls", "SangerRead", function(object, signalRatioCutoff) {
    if (object@inputSource == "ABIF") {
        errors <- character(0)
        errors <- checkSignalRatioCutoff(signalRatioCutoff, errors)
        if (length(errors) == 0) {
            traceMatrix <- object@traceMatrix
            peakPosMatrixRaw <- object@peakPosMatrixRaw
            qualityPhredScoresRaw <- object@abifRawData@data$PCON.2
            readFeature <- object@readFeature
            MBCResult <-
                MakeBaseCallsInside (traceMatrix, peakPosMatrixRaw,
                                     qualityPhredScoresRaw,
                                     signalRatioCutoff, readFeature)
            object@peakPosMatrix <- MBCResult[["peakPosMatrix"]]
            object@peakAmpMatrix <- MBCResult[["peakAmpMatrix"]]
            object@primarySeq <- MBCResult[["primarySeq"]]
            object@secondarySeq <- MBCResult[["secondarySeq"]]

            AASeqResult <-
                calculateAASeq (object@primarySeq,
                                object@QualityReport@trimmedStartPos,
                                object@QualityReport@trimmedFinishPos,
                                object@geneticCode)
            object@primaryAASeqS1 <- AASeqResult[["primaryAASeqS1"]]
            object@primaryAASeqS2 <- AASeqResult[["primaryAASeqS2"]]
            object@primaryAASeqS3 <- AASeqResult[["primaryAASeqS3"]]
            object@ChromatogramParam@signalRatioCutoff <- signalRatioCutoff
            return(object)
        } else {
            stop(errors)
        }
    } else if (object@inputSource == "FASTA") {
        message("SangerRead with 'FASTA' inputSource cannot do base calling")
    }
})

## =============================================================================
## Writing primary sequence into FASTA format
## =============================================================================
#' A SangerRead method which writes the sequence into Fasta files.
#'
#' @title writeFastaSR
#' @name SangerRead-class-writeFastaSR
#' @aliases writeFastaSR,SangerRead-method
#'
#' @docType methods
#' @examples
#' \dontrun{writeFastaSR(sangerReadF, "/Users/chaokuan-hao/Desktop/sangeranalyseR_fasta/SangerRead")}
setMethod("writeFastaSR", "SangerRead", function(object, outputDir, compress,
                                                 compression_level) {
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
        message(">>> outputDir : ", outputDir)
    }
    message("Start writing '", object@readFileName, "' to FASTA format ...")
    fastaFilename <- gsub(file_ext(basename(object@readFileName)), "fa",
                          basename(object@readFileName))
    outputFilename <- file.path(outputDir, fastaFilename)
    ### ------------------------------------------------------------------------
    ### Add trimming in ABIF file format
    ### ------------------------------------------------------------------------
    if(object@inputSource == "ABIF") {
        trimmedStartPos <- object@QualityReport@trimmedStartPos
        trimmedFinishPos <- object@QualityReport@trimmedFinishPos
        targetSeq <- DNAString(substr(as.character(object@primarySeq),
                                      trimmedStartPos+1, trimmedFinishPos))
    } else if (object@inputSource == "FASTA") {
        targetSeq <- object@primarySeq
    }
    ### ------------------------------------------------------------------------
    ### When writing out FASTA, reverse read need to reverseComplement back ~
    ### ------------------------------------------------------------------------
    if (object@readFeature == "Reverse Read") {
        targetSeq <- reverseComplement(targetSeq)
    }
    writeTarget <- DNAStringSet(targetSeq)
    names(writeTarget) <- basename(object@readFileName)
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
#' A SangerRead method which generates final reports of the SangerRead instance.
#'
#' @title generateReportSR
#' @name SangerRead-class-generateReportSR
#' @aliases generateReportSR,SangerRead-method
#'
#' @docType methods
#' @examples
#' \dontrun{generateReportSR(sangerReadF)}
setMethod("generateReportSR", "SangerRead",
          function(object, outputDir,
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
    if (object@inputSource == "ABIF") {
        readName <- sub('\\.ab1$', '', basename(object@readFileName))
    } else if (object@inputSource == "FASTA") {
        readName <- sub('\\.fa$', '', basename(object@readFileName))
        readName <- sub('\\.fasta$', '', readName)
    }
    outputDirSR <- file.path(outputDir, readName)
    ### ------------------------------------------------------------------------
    ### Make sure the directory is exist (SangerRead level)
    ### ------------------------------------------------------------------------
    if (!dir.exists(outputDirSR)) {
        suppressWarnings(dir.create(outputDirSR, recursive = TRUE))
    }
    rootDir <- system.file(package = "sangeranalyseR")
    if (object@inputSource == "ABIF") {
        originRmd <- file.path(rootDir, "rmd", "SangerRead_Report_ab1.Rmd")
        outputHtml <- file.path(outputDirSR, "SangerRead_Report_ab1.html")
    } else if (object@inputSource == "FASTA") {
        originRmd <- file.path(rootDir, "rmd", "SangerRead_Report_fasta.Rmd")
        outputHtml <- file.path(outputDirSR, "SangerRead_Report_fasta.html")
    }
    res <- render(input = originRmd,
                  output_dir = outputDirSR,
                  params = list(SangerRead = object,
                                outputDir = outputDirSR,
                                navigationContigFN = navigationContigFN,
                                navigationAlignmentFN = navigationAlignmentFN))
    return(outputHtml)
})
