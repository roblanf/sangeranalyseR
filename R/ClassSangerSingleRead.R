#' @title sangerSingleRead
#'
#' @description  An S4 class extending sangerseq S4 class
#'
#' @slot readFeature .
#'
#' @name sangerSingleRead-class
#'
#' @rdname sangerSingleRead-class
#'
#' @exportClass sangerSingleRead
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
#' A_chloroticaRead <- new("sangerSingleRead",
#'                         readFeature         = "ForwardRead",
#'                         readFileName        = A_chloroticaFdReadFN,
#'                         cutoffQualityScore  = 60L,
#'                         slidingWindowSize   = 8L)
setClass(
    "sangerSingleRead",
    ### -------------------------------------------------------------------
    ### Input type of each variable of 'SangerMergeReads'.
    ###     * Inherit from 'sangerseq' from sangerseqR.
    ### -------------------------------------------------------------------
    contains="sangerseq",
    slots=c(readFeature         = "character",
            readFileName        = "character",
            abifRawData         = "abif",
            qualityReport       = "qualityReport")
) -> sangerSingleRead


### ============================================================================
### Overwrite initialize for sangerSingleRead (New constructor)
### ============================================================================
setMethod("initialize",
          "sangerSingleRead",
          function(.Object, ...,
                   readFeature         = character(0),
                   readFileName        = character(0),
                   primarySeqID        = primarySeqID,
                   primarySeq          = primarySeq,
                   secondarySeqID      = secondarySeqID,
                   secondarySeq        = secondarySeq,
                   traceMatrix         = traceMatrix,
                   peakPosMatrix       = peakPosMatrix,
                   peakAmpMatrix       = peakAmpMatrix,
                   abifRawData         = abifRawData,
                   qualityReport       = qualityReport,
                   cutoffQualityScore  = 20L,
                   slidingWindowSize   = 5L) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              if (identical(readFeature, character(0))) {
                  msg <- paste("\nYou must assign value to 'readFeature'\n")
                  errors <- c(errors, msg)
              }
              if (!file.exists(readFileName)) {
                  msg <- paste("\n'", readFileName, "'",
                               " foward read file does not exist.\n", sep = "")
                  errors <- c(errors, msg)
              }

              ### --------------------------------------------------------------
              ### Prechecking success. Start to create 'sangerSingleRead'
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
                  qualityReport <- new("qualityReport",
                                       readFeature = readFeature,
                                       qualityPhredScores =
                                           readRawAbif@data$PCON.2,
                                       cutoffQualityScore = cutoffQualityScore,
                                       slidingWindowSize = slidingWindowSize)
                  cutoffQualityScore <- cutoffQualityScore
                  slidingWindowSize <- slidingWindowSize
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
                             qualityReport       = qualityReport)
          })
