#' @title SangerSingleRead
#'
#' @description  An S4 class extending sangerseq S4 class
#'
#' @slot readFeature .
#'
#' @name SangerSingleRead-class
#'
#' @rdname SangerSingleRead-class
#'
#' @exportClass SangerSingleRead
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "RBNII396-13[C_LepFolF,C_LepFolR]_F_1.ab1")
#' A_chloroticaSingleRead <- new("SangerSingleRead",
#'                               readFeature         = "ForwardRead",
#'                               readFileName        = A_chloroticaFdReadFN,
#'                               TrimmingMethod        = "M1",
#'                               M1TrimmingCutoff      = 0.0001,
#'                               M2CutoffQualityScore  = NULL,
#'                               M2SlidingWindowSize   = NULL)
setClass(
    "SangerSingleRead",
    ### -------------------------------------------------------------------
    ### Input type of each variable of 'SangerMergeReads'.
    ###     * Inherit from 'sangerseq' from sangerseqR.
    ### -------------------------------------------------------------------
    contains="sangerseq",
    slots=c(readFeature         = "character",
            readFileName        = "character",
            abifRawData         = "abif",
            QualityReport       = "QualityReport")
) -> SangerSingleRead


### ============================================================================
### Overwrite initialize for SangerSingleRead (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerSingleRead",
          function(.Object, ...,
                   readFeature          = character(0),
                   readFileName         = character(0),
                   primarySeqID         = primarySeqID,
                   primarySeq           = primarySeq,
                   secondarySeqID       = secondarySeqID,
                   secondarySeq         = secondarySeq,
                   traceMatrix          = traceMatrix,
                   peakPosMatrix        = peakPosMatrix,
                   peakAmpMatrix        = peakAmpMatrix,
                   abifRawData          = abifRawData,
                   QualityReport        = QualityReport,
                   TrimmingMethod       = "M1",
                   M1TrimmingCutoff     = 0.0001,
                   M2CutoffQualityScore = NULL,
                   M2SlidingWindowSize  = NULL) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              if (identical(readFeature, character(0))) {
                  msg <- paste("\nYou must assign value to 'readFeature'\n")
                  errors <- c(errors, msg)
              }
              if (!file.exists(readFileName)) {
                  cat ("readFileName", readFileName)
                  msg <- paste("\n'", readFileName, "'",
                               " foward read file does not exist.\n", sep = "")
                  errors <- c(errors, msg)
              }

              ### --------------------------------------------------------------
              ### Input parameter prechecking for TrimmingMethod.
              ### --------------------------------------------------------------
              errors <- checkTrimParam(TrimmingMethod,
                                       M1TrimmingCutoff,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize,
                                       errors)

              ### --------------------------------------------------------------
              ### Prechecking success. Start to create 'SangerSingleRead'
              ### --------------------------------------------------------------
              if (length(errors) == 0) {
                  message(readFeature, " read: Creating abif & sangerseq ...")
                  message("    Creating ", readFeature , " raw abif ...")
                  readRawAbif = read.abif(readFileName)
                  message("    Creating ", readFeature , " raw sangerseq ...")
                  readSangerseq = sangerseq(readRawAbif)

                  primarySeqID        = readSangerseq@primarySeqID
                  primarySeq          = readSangerseq@primarySeq
                  secondarySeqID      = readSangerseq@secondarySeqID
                  secondarySeq        = readSangerseq@secondarySeq
                  traceMatrix         = readSangerseq@traceMatrix
                  peakPosMatrix       = readSangerseq@peakPosMatrix
                  peakAmpMatrix       = readSangerseq@peakAmpMatrix
                  abifRawData         = readRawAbif
                  QualityReport <- new("QualityReport",
                                       readFeature = readFeature,
                                       qualityPhredScores =
                                           readRawAbif@data$PCON.2,
                                       TrimmingMethod = TrimmingMethod,
                                       M1TrimmingCutoff = M1TrimmingCutoff,
                                       M2CutoffQualityScore = M2CutoffQualityScore,
                                       M2SlidingWindowSize = M2SlidingWindowSize)
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature         = readFeature,
                             readFileName        = readFileName,
                             primarySeqID        = primarySeqID,
                             primarySeq          = primarySeq,
                             secondarySeqID      = secondarySeqID,
                             secondarySeq        = secondarySeq,
                             traceMatrix         = traceMatrix,
                             peakPosMatrix       = peakPosMatrix,
                             peakAmpMatrix       = peakAmpMatrix,
                             abifRawData         = abifRawData,
                             QualityReport       = QualityReport)
          })
