#' @title SangerRead
#'
#' @description  An S4 class extending sangerseq S4 class which corresponds to a single ABIF file in Sanger sequencing.
#'
#' @slot readFeature The direction of the Sanger read. The value must be \code{"Forward Read"} or \code{"Reverse Read"}.
#' @slot readFileName The filename of the target ABIF file.
#' @slot abifRawData A S4 class containing all fields in the ABIF file. It is defined in sangerseqR package.
#' @slot QualityReport A S4 class containing quality trimming related inputs and trimming results.
#' @slot ChromatogramParam A S4 class containing chromatogram inputs.
#' @slot primaryAASeqS1 A polypeptide translated from primary DNA sequence starting from the first nucleic acid.
#' @slot primaryAASeqS2 A polypeptide translated from primary DNA sequence starting from the second nucleic acid.
#' @slot primaryAASeqS3 A polypeptide translated from primary DNA sequence starting from the third nucleic acid.
#' @slot geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @slot primarySeqRaw The raw primary sequence from sangerseq class in sangerseqR package before base calling.
#' @slot secondarySeqRaw The raw secondary sequence from sangerseq class in sangerseqR package before base calling.
#' @slot peakPosMatrixRaw The raw peak position matrix from sangerseq class in sangerseqR package before base calling.
#' @slot peakAmpMatrixRaw The raw peak amplitude matrix from sangerseq class in sangerseqR package before base calling.
#'
#' @name SangerRead-class
#'
#' @rdname SangerRead-class
#'
#' @exportClass SangerRead
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R
#' @import sangerseqR
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdFN <- file.path(inputFilesPath,
#'                               "Allolobophora_chlorotica",
#'                               "ACHLO",
#'                               "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
#' sangerRead <- new("SangerRead",
#'                    readFeature           = "Forward Read",
#'                    readFileName          = A_chloroticaFdFN,
#'                    geneticCode           = GENETIC_CODE,
#'                    TrimmingMethod        = "M1",
#'                    M1TrimmingCutoff      = 0.0001,
#'                    M2CutoffQualityScore  = NULL,
#'                    M2SlidingWindowSize   = NULL,
#'                    baseNumPerRow         = 100,
#'                    heightPerRow          = 200,
#'                    signalRatioCutoff     = 0.33,
#'                    showTrimmed           = TRUE)
setClass(
    "SangerRead",
    ### ------------------------------------------------------------------------
    ### Input type of each variable of 'SangerContig'.
    ###     * Inherit from 'sangerseq' from sangerseqR.
    ### ------------------------------------------------------------------------
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
)

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
                  traceMatrix      <- readSangerseq@traceMatrix
                  peakPosMatrixRaw <- readSangerseq@peakPosMatrix
                  peakPosMatrix    <- readSangerseq@peakPosMatrix
                  peakAmpMatrixRaw <- readSangerseq@peakAmpMatrix
                  peakAmpMatrix    <- readSangerseq@peakAmpMatrix

                  ### ----------------------------------------------------------
                  ### Definition of 'PCON.1' & 'PCON.2'
                  ##### PCON.1: char => Per-base quality values (edited)
                  ##### PCON.2: char => Per-base quality values
                  ### ----------------------------------------------------------
                  ### ----------------------------------------------------------
                  ### 1. Running 'MakeBaseCall'! Remember to update parameters!
                  ### ----------------------------------------------------------
                  MBCResult <-
                      MakeBaseCallsInside (traceMatrix, peakPosMatrixRaw,
                                           abifRawData@data$PCON.2,
                                           signalRatioCutoff, readFeature)

                  ### ==========================================================
                  ### 2. Update Once (Only during creation)
                  ###    Basecall primary seq length will be same !
                  ### ==========================================================
                  qualityPhredScores <- MBCResult[["qualityPhredScores"]]
                  ### ----------------------------------------------------------
                  ##### 'QualityReport' creation
                  ### ----------------------------------------------------------
                  QualityReport <-
                      new("QualityReport",
                          qualityPhredScoresRaw = abifRawData@data$PCON.2,
                          qualityPhredScores    = qualityPhredScores,
                          TrimmingMethod        = TrimmingMethod,
                          M1TrimmingCutoff      = M1TrimmingCutoff,
                          M2CutoffQualityScore  = M2CutoffQualityScore,
                          M2SlidingWindowSize   = M2SlidingWindowSize)

                  ### ==========================================================
                  ### 3. Update everytime (whenever 'signalRatioCutoff' changed)
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

                  AASeqResult    <- calculateAASeq (primarySeq, geneticCode)
                  primaryAASeqS1 <- AASeqResult[["primaryAASeqS1"]]
                  primaryAASeqS2 <- AASeqResult[["primaryAASeqS2"]]
                  primaryAASeqS3 <- AASeqResult[["primaryAASeqS3"]]
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
