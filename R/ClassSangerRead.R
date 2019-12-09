#' @title SangerRead
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
#' @name SangerRead-class
#'
#' @rdname SangerRead-class
#'
#' @exportClass SangerRead
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "RBNII396-13[C_LepFolF,C_LepFolR]_F_1.ab1")
#' A_chloroticaRead <- new("SangerRead",
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
    "SangerRead",
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
) -> SangerRead


### ============================================================================
### Overwrite initialize for SangerRead (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerRead",
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
              ### Prechecking success. Start to create 'SangerRead'
              ### --------------------------------------------------------------
              if (length(errors) == 0) {
                  message(readFeature, " read: Creating abif & sangerseq ...")
                  message("    * Creating ", readFeature , " raw abif ...")
                  abifRawData = read.abif(readFileName)
                  message("    * Creating ", readFeature , " raw sangerseq ...")
                  readSangerseq = sangerseq(abifRawData)
                  primarySeqID = readSangerseq@primarySeqID
                  secondarySeqID = readSangerseq@secondarySeqID

                  ### ----------------------------------------------------------
                  ### With non-raw & raw primarySeq / secondarySeq
                  ### ----------------------------------------------------------
                  primarySeqRaw = readSangerseq@primarySeq
                  primarySeq = readSangerseq@primarySeq
                  secondarySeqRaw = readSangerseq@secondarySeq
                  secondarySeq = readSangerseq@secondarySeq

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

                  traceMatrix         = readSangerseq@traceMatrix
                  peakPosMatrixRaw    = readSangerseq@peakPosMatrix
                  peakPosMatrix       = readSangerseq@peakPosMatrix
                  peakAmpMatrixRaw    = readSangerseq@peakAmpMatrix
                  peakAmpMatrix       = readSangerseq@peakAmpMatrix

                  # qualityBaseScoresRaw <- 10** (qualityPhredScoresRaw / (-10.0))

                  ### ----------------------------------------------------------
                  ### After Here, if 'MakeBaseCall' is called, we need to
                  ###    update parameters !!
                  ### ===== Update Once ========================================
                  ### 'primarySeqID', 'secondarySeqID'
                  ### 'primarySeq', 'QualityReport@qualityScoresID',
                  ### 'QualityReport@qualityPhredScores',
                  ### 'QualityReport@qualityBaseScores'
                  ### ===== Update everytime ===================================
                  ### 'secondarySeq',
                  ### 'peakPosMatrix', 'peakAmpMatrix',
                  ### ----------------------------------------------------------

                  ### ----------------------------------------------------------
                  ### Definition of 'PCON.1' & 'PCON.2'
                  ##### PCON.1: char => Per-base quality values (edited)
                  ##### PCON.2: char => Per-base quality values
                  ### ----------------------------------------------------------


                  ### ----------------------------------------------------------
                  ### Running 'MakeBaseCall' !! Remember to update parameters !!
                  ### ----------------------------------------------------------
                  MBCResult <-
                      MakeBaseCallsInside (traceMatrix, peakPosMatrixRaw,
                                           abifRawData@data$PCON.2,
                                           signalRatioCutoff, readFeature)

                  ### ==========================================================
                  ### Update Once (Only during creation)
                  ###    Basecall primary seq length will be same !
                  ### ==========================================================
                  qualityPhredScores <- MBCResult[["qualityPhredScores"]]
                  ### ----------------------------------------------------------
                  ##### 'QualityReport' creation
                  ### ----------------------------------------------------------
                  QualityReport <-
                      new("QualityReport",
                          qualityPhredScoresRaw = abifRawData@data$PCON.2,
                          qualityPhredScores = qualityPhredScores,
                          TrimmingMethod = TrimmingMethod,
                          M1TrimmingCutoff = M1TrimmingCutoff,
                          M2CutoffQualityScore = M2CutoffQualityScore,
                          M2SlidingWindowSize = M2SlidingWindowSize)

                  ### ==========================================================
                  ### Update everytime (whenever 'signalRatioCutoff' is changed)
                  ### ==========================================================
                  primarySeq <- MBCResult[["primarySeq"]]
                  secondarySeq <- MBCResult[["secondarySeq"]]
                  peakPosMatrix <- MBCResult[["peakPosMatrix"]]
                  peakAmpMatrix <- MBCResult[["peakAmpMatrix"]]
                  ### ----------------------------------------------------------
                  ##### 'QualityReport' & 'ChromatogramParam' creation
                  ### ----------------------------------------------------------
                  ChromatogramParam <- new("ChromatogramParam",
                                           baseNumPerRow     = baseNumPerRow,
                                           heightPerRow      = heightPerRow,
                                           signalRatioCutoff = signalRatioCutoff,
                                           showTrimmed       = showTrimmed)

                  AASeqResult <- calculateAASeq (primarySeq, geneticCode)
                  primaryAASeqS1 <- AASeqResult[["primaryAASeqS1"]]
                  primaryAASeqS2 <- AASeqResult[["primaryAASeqS2"]]
                  primaryAASeqS3 <- AASeqResult[["primaryAASeqS3"]]
                  ### ==========================================================
                  ### ==========================================================
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
