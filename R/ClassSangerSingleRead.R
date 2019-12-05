#' @title SangerSingleRead
#'
#' @description  An S4 class extending sangerseq S4 class
#'
#' @slot readFeature .
#' @slot readFeature .
#' @slot readFileName .
#' @slot abifRawData .
#' @slot QualityReport .
#' @slot ChromatogramParam .
#' @slot primaryAASeq .
#' @slot geneticCode .
#' @slot primarySeqBC .
#' @slot secondarySeqBC .
#' @slot peakPosMatrixBC .
#' @slot peakAmpMatrixBC .
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
#'                               readFeature           = "Forward Read",
#'                               readFileName          = A_chloroticaFdReadFN,
#'                               geneticCode           = GENETIC_CODE,
#'                               TrimmingMethod        = "M2",
#'                               M1TrimmingCutoff      = NULL,
#'                               M2CutoffQualityScore  = 40,
#'                               M2SlidingWindowSize   = 10,
#'                               baseNumPerRow         = 100,
#'                               heightPerRow          = 200,
#'                               signalRatioCutoff     = 0.33,
#'                               showTrimmed           = TRUE)
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
            QualityReport       = "QualityReport",
            ChromatogramParam   = "ChromatogramParam",
            primaryAASeq        = "AAString",
            geneticCode         = "character",
            primarySeqBC        = "DNAStringORNULL",
            secondarySeqBC      = "DNAStringORNULL",
            peakPosMatrixBC     = "matrixORNULL",
            peakAmpMatrixBC     = "matrixORNULL"
            )
) -> SangerSingleRead


### ============================================================================
### Overwrite initialize for SangerSingleRead (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerSingleRead",
          function(.Object, ...,
                   readFeature          = character(0),
                   readFileName         = character(0),
                   geneticCode          = GENETIC_CODE,
                   TrimmingMethod       = "M1",
                   M1TrimmingCutoff     = 0.0001,
                   M2CutoffQualityScore = NULL,
                   M2SlidingWindowSize  = NULL,
                   baseNumPerRow        = 100,
                   heightPerRow         = 200,
                   signalRatioCutoff    = 0.33,
                   showTrimmed          = TRUE) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              errors <- checkReadFeature (readFeature, errors)
              errors <- checkReadFileName (readFileName, errors)
              errors <- checkGeneticCode (geneticCode, errors)

              ##### ------------------------------------------------------------
              ##### Input parameter prechecking for TrimmingMethod.
              ##### ------------------------------------------------------------
              errors <- checkTrimParam(TrimmingMethod,
                                       M1TrimmingCutoff,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize,
                                       errors)

              ##### ------------------------------------------------------------
              ##### Input parameter prechecking for ChromatogramParam
              ##### ------------------------------------------------------------
              errors <- checkBaseNumPerRow (baseNumPerRow, errors)
              errors <- checkHeightPerRow (baseNumPerRow, errors)
              errors <- checkSignalRatioCutoff (signalRatioCutoff, errors)
              errors <- checkShowTrimmed (showTrimmed, errors)

              ### --------------------------------------------------------------
              ### Prechecking success. Start to create 'SangerSingleRead'
              ### --------------------------------------------------------------
              if (length(errors) == 0) {
                  message(readFeature, " read: Creating abif & sangerseq ...")
                  message("    * Creating ", readFeature , " raw abif ...")
                  readRawAbif = read.abif(readFileName)
                  message("    * Creating ", readFeature , " raw sangerseq ...")
                  readSangerseq = sangerseq(readRawAbif)
                  primarySeqID        = readSangerseq@primarySeqID
                  secondarySeqID      = readSangerseq@secondarySeqID
                  if (readFeature == "Forward Read") {
                      primarySeq = readSangerseq@primarySeq
                      secondarySeq = readSangerseq@secondarySeq
                  } else if (readFeature == "Reverse Read") {
                      primarySeq = reverseComplement(
                          readSangerseq@primarySeq)
                      secondarySeq = reverseComplement(
                          readSangerseq@secondarySeq)
                  }
                  primaryAASeq        = suppressWarnings(translate(primarySeq,
                                                  genetic.code = geneticCode,
                                                  no.init.codon=TRUE,
                                                  if.fuzzy.codon="solve"))
                  traceMatrix         = readSangerseq@traceMatrix
                  peakPosMatrix       = readSangerseq@peakPosMatrix
                  peakAmpMatrix       = readSangerseq@peakAmpMatrix
                  abifRawData         = readRawAbif

                  ### ----------------------------------------------------------
                  ### Definition of 'PCON.1' & 'PCON.2'
                  ##### PCON.1: char => Per-base quality values (edited)
                  ##### PCON.2: char => Per-base quality values
                  ### ----------------------------------------------------------
                  QualityReport <- new("QualityReport",
                                       qualityPhredScores =
                                           readRawAbif@data$PCON.2,
                                       TrimmingMethod = TrimmingMethod,
                                       M1TrimmingCutoff = M1TrimmingCutoff,
                                       M2CutoffQualityScore = M2CutoffQualityScore,
                                       M2SlidingWindowSize = M2SlidingWindowSize)
                  ChromatogramParam <- new("ChromatogramParam",
                                           baseNumPerRow     = baseNumPerRow,
                                           heightPerRow      = heightPerRow,
                                           signalRatioCutoff = signalRatioCutoff,
                                           showTrimmed       = showTrimmed)

                  ### ----------------------------------------------------------
                  ### Before running MakeBaseCall, the parameters below are NULL
                  ###     'MakeBaseCallsInside' will be called during S4 object
                  ###     creation
                  ### ----------------------------------------------------------
                  MBCResult <- MakeBaseCallsInside (traceMatrix, peakPosMatrix,
                                                    QualityReport,
                                                    signalRatioCutoff =
                                                        signalRatioCutoff)
                  QualityReport <- MBCResult[["QualityReport"]]
                  peakPosMatrixBC <- MBCResult[["peakPosMatrixBC"]]
                  peakAmpMatrixBC <- MBCResult[["peakAmpMatrixBC"]]
                  primarySeqID <- MBCResult[["primarySeqID"]]
                  primarySeqBC <- MBCResult[["primarySeqBC"]]
                  secondarySeqID <- MBCResult[["secondarySeqID"]]
                  secondarySeqBC <- MBCResult[["secondarySeqBC"]]
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature         = readFeature,
                             readFileName        = readFileName,
                             geneticCode         = geneticCode,
                             primarySeqID        = primarySeqID,
                             primarySeq          = primarySeq,
                             primarySeqBC        = primarySeqBC,
                             secondarySeqID      = secondarySeqID,
                             secondarySeq        = secondarySeq,
                             secondarySeqBC      = secondarySeqBC,
                             primaryAASeq        = primaryAASeq,
                             traceMatrix         = traceMatrix,
                             peakPosMatrix       = peakPosMatrix,
                             peakPosMatrixBC     = peakPosMatrixBC,
                             peakAmpMatrix       = peakAmpMatrix,
                             peakAmpMatrixBC     = peakAmpMatrixBC,
                             abifRawData         = abifRawData,
                             QualityReport       = QualityReport,
                             ChromatogramParam   = ChromatogramParam)
          })
