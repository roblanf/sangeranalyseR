### ============================================================================
### Plotting quality for each base for "SangerRead" S4 object
### ============================================================================
#' A SangerRead method which creates quality base interactive plot.
#'
#' @title qualityBasePlot
#' @name SangerRead-class-qualityBasePlot
#' @aliases qualityBasePlot,SangerRead-method
#'
#' @param object A SangerRead S4 instance.
#'
#' @return A quality plot.
#'
#' @examples
#' data("sangerReadFData")
#' \dontrun{
#' qualityBasePlot(sangerReadFData)}
setMethod("qualityBasePlot",  "SangerRead", function(object){
    if (object@inputSource == "ABIF") {
        plotting <- preQualityBasePlot(object@QualityReport)
        plotting
    } else if (object@inputSource == "FASTA") {
        log_info("SangerRead with 'FASTA' inputSource ",
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
#' @param object A SangerRead S4 instance.
#' @param TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#'
#' @return A SangerRead instance.
#'
#' @examples
#' data("sangerReadFData")
#' updateQualityParam(sangerReadFData,
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
            AASeqResult    <- calculateAASeq (object@primarySeq,
                                              object@QualityReport@trimmedStartPos,
                                              object@QualityReport@trimmedFinishPos,
                                              object@geneticCode)
            object@primaryAASeqS1 <- AASeqResult[["primaryAASeqS1"]]
            object@primaryAASeqS2 <- AASeqResult[["primaryAASeqS2"]]
            object@primaryAASeqS3 <- AASeqResult[["primaryAASeqS3"]]
            return(object)
        } else {
            log_error(errors)
        }
    } else if (object@inputSource == "FASTA") {
        log_info("SangerRead with 'FASTA' inputSource ",
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
#' @param object A SangerRead S4 instance.
#' @param signalRatioCutoff The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is \code{0.33}.
#'
#' @return A \code{SangerRead} instance.
#'
#' @examples
#' data("sangerReadFData")
#' newSangerReadFData <- MakeBaseCalls(sangerReadFData, signalRatioCutoff = 0.22)
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
            log_error(errors)
        }
    } else if (object@inputSource == "FASTA") {
        log_info("SangerRead with 'FASTA' inputSource cannot do base calling")
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
#' @param object A SangerRead S4 instance.
#' @param outputDir The output directory of the generated FASTA file.
#' @param compress Like for the \code{save} function in base R, must be \code{TRUE} or \code{FALSE} (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip". This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param compression_level This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#'
#' @return The output absolute path to the FASTA file.
#'
#' @examples
#' data("sangerReadFData")
#' writeFastaSR(sangerReadFData)
setMethod("writeFastaSR", "SangerRead", function(object, outputDir, compress,
                                                 compression_level) {
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
        log_info(">>> outputDir : ", outputDir)
    }
    log_info("Start writing '", object@readFileName, "' to FASTA format ...")
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
    log_info("\n >> '", outputFilename, "' is written")
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
#' @param object A SangerRead S4 instance.
#' @param outputDir The output directory of the generated HTML report.
#' @param navigationContigFN The internal parameter passed to HTML report. Users should not modify this parameter on their own.
#' @param navigationAlignmentFN The internal parameter passed to HTML report. Users should not modify this parameter on their own.
#'
#' @return The output absolute path to the SangerRead's HTML file.
#'
#' @examples
#' data("sangerReadFData")
#' \dontrun{
#' generateReportSR(sangerReadFData)}
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
        log_info(">>> outputDir : ", outputDir)
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

## =============================================================================
## Generating summary table for SangerRead instance
## =============================================================================
#' A SangerRead method which generates summary table for SangerRead instance
#'
#' @title readTable
#' @name SangerRead-class-readTable
#' @aliases readTable,readTable-method
#'
#' @param object A SangerRead S4 instance.
#'
#' @return None
#'
#' @examples
#' data(sangerReadFData)
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' readTable(sangerReadFData)
#' readTable(sangerContigData)
#' readTable(sangerAlignmentData)
#' }
setMethod("readTable", "SangerRead", function(object, indentation = 0) {
    space <- paste(rep(' ', indentation), collapse = "")
    inputSource <- object@inputSource
    readFeature <- object@readFeature
    readFileNameAbs <- object@readFileName
    readFileNameBase <- basename(object@readFileName)
    if (object@inputSource == "ABIF") {
        TrimmingMethod <- object@QualityReport@TrimmingMethod
        M1TrimmingCutoff <- object@QualityReport@M1TrimmingCutoff
        M2CutoffQualityScore <- object@QualityReport@M2CutoffQualityScore
        M2SlidingWindowSize <- object@QualityReport@M2SlidingWindowSize
        rawSeqLength <- object@QualityReport@rawSeqLength
        trimmedSeqLength <- object@QualityReport@trimmedSeqLength
        rawMeanQualityScore <- object@QualityReport@rawMeanQualityScore
        trimmedMeanQualityScore <- object@QualityReport@trimmedMeanQualityScore
        primaryDNA <- as.character(object@primarySeq)
        secondaryDNA <- as.character(object@secondarySeq)
        # passed secondary peak cutoff [yesn/no]
        signalRatioCutoff <- object@ChromatogramParam@signalRatioCutoff
        # read included in contig [yes/no]
        trimmedStartPos <- object@QualityReport@trimmedStartPos
        trimmedFinishPos <- object@QualityReport@trimmedFinishPos
        primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
        secondaryDNA <- substr(secondaryDNA, trimmedStartPos+1,trimmedFinishPos)
        if (TrimmingMethod == "M1") {
            cat(space, "SangerRead S4 instance\n",
                space, "                 Input source : ", inputSource, "\n",
                space, "                 Read feature : ", readFeature, "\n",
                space, "          Read fileName (abs) : ", readFileNameAbs, "\n",
                space, "         Read fileName (base) : ", readFileNameBase, "\n",
                space, "              Trimming method : ", TrimmingMethod, "\n",
                space, "              Trimming cutoff : ", M1TrimmingCutoff, "\n",
                space, "     Trimming start base pair : ", trimmedStartPos, "\n",
                space, "    Trimming finish base pair : ", trimmedFinishPos, "\n",
                space, "  Read length before trimming : ", rawSeqLength, "\n",
                space, "   Read length after trimming : ", trimmedSeqLength, "\n",
                space, " Read quality before trimming : ", rawMeanQualityScore, "\n",
                space, "  Read quality after trimming : ", trimmedMeanQualityScore,"\n",
                space, "        Secondary peak cutoff : ", signalRatioCutoff,"\n",
                space, "             Primary Sequence : ", primaryDNA, "\n",
                space, "           Secondary Sequence : ", secondaryDNA, "\n"
            )
        } else if (TrimmingMethod == "M2") {
            cat(space, "SangerRead S4 instance\n",
                space, "                 Input source : ", inputSource, "\n",
                space, "                 Read feature : ", readFeature, "\n",
                space, "          Read fileName (abs) : ", readFileNameAbs, "\n",
                space, "         Read fileName (base) : ", readFileNameBase, "\n",
                space, "              Trimming method : ", TrimmingMethod, "\n",
                space, "              Trimming cutoff : ", M2CutoffQualityScore, "\n",
                space, "      Trimming sliding window : ", M2SlidingWindowSize, "\n",
                space, "     Trimming start base pair : ", trimmedStartPos, "\n",
                space, "    Trimming finish base pair : ", trimmedFinishPos, "\n",
                space, "  Read length before trimming : ", rawSeqLength, "\n",
                space, "   Read length after trimming : ", trimmedSeqLength, "\n",
                space, " Read quality before trimming : ", rawMeanQualityScore, "\n",
                space, "  Read quality after trimming : ", trimmedMeanQualityScore,"\n",
                space, "        Secondary peak cutoff : ", signalRatioCutoff,"\n",
                space, "             Primary Sequence : ", primaryDNA, "\n",
                space, "           Secondary Sequence : ", secondaryDNA, "\n"
            )
        }
    } else if (object@inputSource == "FASTA") {
        fastaReadName <- object@fastaReadName
        seqLength <- length(object@primarySeq)
        primarySeq <- as.character(object@primarySeq)
        cat(space, "SangerRead S4 instance\n",
            space, "                 Input source : ", inputSource, "\n",
            space, "                 Read feature : ", readFeature, "\n",
            space, "          Read fileName (abs) : ", readFileNameAbs, "\n",
            space, "         Read fileName (base) : ", readFileNameBase, "\n",
            space, "              Fasta Read Name : ", fastaReadName, "\n",
            space, "                  Read length : ", seqLength, "\n",
            space, "             Primary Sequence : ", primarySeq, "\n"
        )
    }
})

# sangerContigData@inputSource
# sangerContigData@fastaFileName
# sangerContigData@namesConversionCSV
# sangerContigData@parentDirectory
# 
# # contig name,
# sangerContigData@contigName
# sangerContigData@suffixForwardRegExp
# sangerContigData@suffixReverseRegExp
# sangerContigData@minReadsNum
# # passed length cutoff [yes/no],
# sangerContigData@minReadLength
# sangerContigData@refAminoAcidSeq
# sangerContigData@minFractionCall
# sangerContigData@maxFractionLost
# sangerContigData@acceptStopCodons
# sangerContigData@readingFrame
# sangerContigData@contigSeq
# sangerContigData@alignment
# sangerContigData@differencesDF
# sangerContigData@distanceMatrix
# sangerContigData@dendrogram
# sangerContigData@indelsDF
# sangerContigData@stopCodonsDF
# sangerContigData@secondaryPeakDF
# 
# # number of indels,
# # passed indel cutoff [yes/no]
# sangerContigData@indelsDF
# 
# # secondary peak cutoff,
# # number of secondary peaks
# sangerContigData@secondaryPeakDF
# 
# # passed quality cutoff [yes/no],