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
#' @slot primarySeqRaw .
#' @slot secondarySeqRaw .
#' @slot peakPosMatrixRaw .
#' @slot peakAmpMatrixRaw .
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
            primaryAASeqS1      = "AAString",
            primaryAASeqS2      = "AAString",
            primaryAASeqS3      = "AAString",
            geneticCode         = "character",
            primarySeqRaw       = "DNAString",
            secondarySeqRaw     = "DNAString",
            peakPosMatrixRaw    = "matrix",
            peakAmpMatrixRaw    = "matrix"
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
                  primarySeqID = readSangerseq@primarySeqID
                  secondarySeqID = readSangerseq@secondarySeqID

                  ### ----------------------------------------------------------
                  ### With non-raw & raw primarySeq / secondarySeq
                  ### ----------------------------------------------------------
                  primarySeqRaw = readSangerseq@primarySeq
                  primarySeq = readSangerseq@primarySeq
                  secondarySeqRaw = readSangerseq@secondarySeq
                  secondarySeq = readSangerseq@secondarySeq

                  ### ----------------------------------------------------------
                  ### After Here, if 'MakeBaseCall' is called, we need to
                  ###    update parameters !!
                  ### 'primarySeq', 'secondarySeq',
                  ### 'traceMatrix', 'peakPosMatrix', 'peakAmpMatrix',
                  ### 'QualityReport@', 'ChromatogramParam@'
                  ### ----------------------------------------------------------
                  if (readFeature == "Reverse Read") {
                      primarySeqRaw =
                          reverseComplement(readSangerseq@primarySeq)
                      primarySeq =
                          reverseComplement(readSangerseq@primarySeq)
                      secondarySeqRaw =
                          reverseComplement(readSangerseq@secondarySeq)
                      secondarySeq =
                          reverseComplement(readSangerseq@secondarySeq)
                  }
                  primaryAASeqS1 =
                      suppressWarnings(translate(primarySeq,
                                                 genetic.code = geneticCode,
                                                 no.init.codon=TRUE,
                                                 if.fuzzy.codon="solve"))
                  DNASeqshift1 <-DNAString(substr(as.character(primarySeq),
                                                  2, length(primarySeq)))

                  primaryAASeqS2 =
                      suppressWarnings(translate(DNASeqshift1,
                                                 genetic.code = geneticCode,
                                                 no.init.codon=TRUE,
                                                 if.fuzzy.codon="solve"))
                  DNASeqshift2 <-
                      DNAString(substr(as.character(primarySeq),
                                       3, length(primarySeq)))
                  primaryAASeqS3 <-
                      suppressWarnings(translate(DNASeqshift2,
                                                 genetic.code = geneticCode,
                                                 no.init.codon=TRUE,
                                                 if.fuzzy.codon="solve"))

                  traceMatrixRaw      = readSangerseq@traceMatrix
                  traceMatrix         = readSangerseq@traceMatrix
                  peakPosMatrixRaw    = readSangerseq@peakPosMatrix
                  peakPosMatrix       = readSangerseq@peakPosMatrix
                  peakAmpMatrixRaw    = readSangerseq@peakAmpMatrix
                  peakAmpMatrix       = readSangerseq@peakAmpMatrix

                  abifRawData         = readRawAbif

                  ### ----------------------------------------------------------
                  ### Definition of 'PCON.1' & 'PCON.2'
                  ##### PCON.1: char => Per-base quality values (edited)
                  ##### PCON.2: char => Per-base quality values
                  ### ----------------------------------------------------------
                  QualityReport <- new("QualityReport",
                                       qualityPhredScoresRaw =
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
                  MBCResult <-
                      MakeBaseCallsInside (traceMatrixRaw, peakPosMatrixRaw,
                                           QualityReport@qualityPhredScoresRaw,
                                           QualityReport@qualityBaseScoresRaw,
                                           signalRatioCutoff=signalRatioCutoff)

                  # New to re-trimm again !!!!
                  QualityReport@qualityScoresID <-
                      MBCResult[["qualityScoresID"]]
                  QualityReport@qualityPhredScores <-
                      MBCResult[["qualityPhredScores"]]
                  QualityReport@qualityBaseScores <-
                      MBCResult[["qualityBaseScores"]]

                  peakPosMatrix <- MBCResult[["peakPosMatrix"]]
                  peakAmpMatrix <- MBCResult[["peakAmpMatrix"]]
                  primarySeqID <- MBCResult[["primarySeqID"]]
                  primarySeq <- MBCResult[["primarySeq"]]
                  secondarySeqID <- MBCResult[["secondarySeqID"]]
                  secondarySeq <- MBCResult[["secondarySeq"]]
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature         = readFeature,
                             readFileName        = readFileName,
                             geneticCode         = geneticCode,
                             primarySeqID        = primarySeqID,
                             primarySeqRaw       = primarySeqRaw,
                             primarySeq          = primarySeq,
                             secondarySeqID      = secondarySeqID,
                             secondarySeqRaw     = secondarySeqRaw,
                             secondarySeq        = secondarySeq,
                             primaryAASeqS1      = primaryAASeqS1,
                             primaryAASeqS2      = primaryAASeqS2,
                             primaryAASeqS3      = primaryAASeqS3,
                             traceMatrix         = traceMatrix,
                             peakPosMatrix       = peakPosMatrix,
                             peakPosMatrixRaw    = peakPosMatrixRaw,
                             peakAmpMatrix       = peakAmpMatrix,
                             peakAmpMatrixRaw    = peakAmpMatrixRaw,
                             abifRawData         = abifRawData,
                             QualityReport       = QualityReport,
                             ChromatogramParam   = ChromatogramParam)
          })
