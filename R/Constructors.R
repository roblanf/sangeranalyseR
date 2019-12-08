### ============================================================================
### Self-defined constructor for AlignedConsensusSet
### ============================================================================
#' @description
#'
#' @param parentDirectory .
#' @param consensusReadName .
#' @param suffixForwardRegExp .
#' @param suffixReverseRegExp .
#' @param TrimmingMethod .
#' @param M1TrimmingCutoff .
#' @param M2CutoffQualityScore .
#' @param M2SlidingWindowSize .
#' @param baseNumPerRow .
#' @param heightPerRow .
#' @param signalRatioCutoff .
#' @param showTrimmed .
#' @param refAminoAcidSeq .
#' @param minReadsNum .
#' @param minReadLength .
#' @param minFractionCall .
#' @param maxFractionLost .
#' @param geneticCode .
#' @param acceptStopCodons .
#' @param readingFrame .
#' @param processorsNum .
#'
#' @return AlignedConsensusSet
#' @export
#' @author Kuan-Hao Chao
#' @example
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' SangerAlignedConsensusSet <- SangerAlignedConsensusSet(
#'                                parentDirectory       = rawDataDir,
#'                                suffixForwardRegExp   = suffixForwardRegExp,
#'                                suffixReverseRegExp   = suffixReverseRegExp,
#'                                refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                                TrimmingMethod        = "M2",
#'                                M1TrimmingCutoff      = NULL,
#'                                M2CutoffQualityScore  = 40,
#'                                M2SlidingWindowSize   = 10,
#'                                baseNumPerRow         = 100,
#'                                heightPerRow          = 200,
#'                                signalRatioCutoff     = 0.33,
#'                                showTrimmed           = TRUE)
SangerAlignedConsensusSet <- function(parentDirectory        = character(0),
                                      consensusReadName      = character(0),
                                      suffixForwardRegExp    = character(0),
                                      suffixReverseRegExp    = character(0),
                                      TrimmingMethod         = "M1",
                                      M1TrimmingCutoff       = 0.0001,
                                      M2CutoffQualityScore   = NULL,
                                      M2SlidingWindowSize    = NULL,
                                      baseNumPerRow          = 100,
                                      heightPerRow           = 200,
                                      signalRatioCutoff      = 0.33,
                                      showTrimmed            = TRUE,
                                      refAminoAcidSeq        = "",
                                      minReadsNum            = 2,
                                      minReadLength          = 20,
                                      minFractionCall        = 0.5,
                                      maxFractionLost        = 0.5,
                                      geneticCode            = GENETIC_CODE,
                                      acceptStopCodons       = TRUE,
                                      readingFrame           = 1,
                                      processorsNum          = 1) {
    newAlignedConsensusSet <- new("SangerAlignedConsensusSet",
                                  parentDirectory        = parentDirectory,
                                  consensusReadName     = consensusReadName,
                                  suffixForwardRegExp    = suffixForwardRegExp,
                                  suffixReverseRegExp    = suffixReverseRegExp,
                                  TrimmingMethod         = TrimmingMethod,
                                  M1TrimmingCutoff       = M1TrimmingCutoff,
                                  M2CutoffQualityScore   = M2CutoffQualityScore,
                                  M2SlidingWindowSize    = M2SlidingWindowSize,
                                  baseNumPerRow          = baseNumPerRow,
                                  heightPerRow           = heightPerRow,
                                  signalRatioCutoff      = signalRatioCutoff,
                                  showTrimmed            = showTrimmed,
                                  refAminoAcidSeq        = refAminoAcidSeq,
                                  minReadsNum            = minReadsNum,
                                  minReadLength          = minReadLength,
                                  minFractionCall        = minFractionCall,
                                  maxFractionLost        = maxFractionLost,
                                  geneticCode            = geneticCode,
                                  acceptStopCodons       = acceptStopCodons,
                                  readingFrame           = readingFrame,
                                  processorsNum          = processorsNum)
    return(newAlignedConsensusSet)
}



### ============================================================================
### Self-defined constructor for SangerConsensusRead
### ============================================================================
#' @description
#'
#' @param parentDirectory .
#' @param consensusReadName .
#' @param suffixForwardRegExp .
#' @param suffixReverseRegExp .
#' @param TrimmingMethod .
#' @param M1TrimmingCutoff .
#' @param M2CutoffQualityScore .
#' @param M2SlidingWindowSize .
#' @param baseNumPerRow .
#' @param heightPerRow .
#' @param signalRatioCutoff .
#' @param showTrimmed .
#' @param refAminoAcidSeq .
#' @param minReadsNum .
#' @param minReadLength .
#' @param minFractionCall .
#' @param maxFractionLost .
#' @param geneticCode .
#' @param acceptStopCodons .
#' @param readingFrame .
#' @param processorsNum .
#'
#' @return SangerConsensusRead
#' @export
#' @author Kuan-Hao Chao
#' @example
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' consensusReadName <- "RBNII395-13[C_LepFolF,C_LepFolR]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' A_chloroticConsensusReads <- SangerConsensusRead(
#'                                parentDirectory       = inputFilesParentDir,
#'                                consensusReadName    = consensusReadName,
#'                                suffixForwardRegExp   = suffixForwardRegExp,
#'                                suffixReverseRegExp   = suffixReverseRegExp,
#'                                TrimmingMethod        = "M2",
#'                                M1TrimmingCutoff      = NULL,
#'                                M2CutoffQualityScore  = 40,
#'                                M2SlidingWindowSize   = 10,
#'                                baseNumPerRow         = 100,
#'                                heightPerRow          = 200,
#'                                signalRatioCutoff     = 0.33,
#'                                showTrimmed           = TRUE)
SangerConsensusRead <- function(parentDirectory        = character(0),
                                consensusReadName     = character(0),
                                suffixForwardRegExp    = character(0),
                                suffixReverseRegExp    = character(0),
                                TrimmingMethod         = "M1",
                                M1TrimmingCutoff       = 0.0001,
                                M2CutoffQualityScore   = NULL,
                                M2SlidingWindowSize    = NULL,
                                baseNumPerRow          = 100,
                                heightPerRow           = 200,
                                signalRatioCutoff      = 0.33,
                                showTrimmed            = TRUE,
                                refAminoAcidSeq        = "",
                                minReadsNum            = 2,
                                minReadLength          = 20,
                                minFractionCall        = 0.5,
                                maxFractionLost        = 0.5,
                                geneticCode            = GENETIC_CODE,
                                acceptStopCodons       = TRUE,
                                readingFrame           = 1,
                                processorsNum          = 1) {
    newConsensusReads <- new("SangerConsensusRead",
                             parentDirectory        = parentDirectory,
                             consensusReadName     = consensusReadName,
                             suffixForwardRegExp    = suffixForwardRegExp,
                             suffixReverseRegExp    = suffixReverseRegExp,
                             TrimmingMethod         = TrimmingMethod,
                             M1TrimmingCutoff       = M1TrimmingCutoff,
                             M2CutoffQualityScore   = M2CutoffQualityScore,
                             M2SlidingWindowSize    = M2SlidingWindowSize,
                             baseNumPerRow          = baseNumPerRow,
                             heightPerRow           = heightPerRow,
                             signalRatioCutoff      = signalRatioCutoff,
                             showTrimmed            = showTrimmed,
                             refAminoAcidSeq        = refAminoAcidSeq,
                             minReadsNum            = minReadsNum,
                             minReadLength          = minReadLength,
                             minFractionCall        = minFractionCall,
                             maxFractionLost        = maxFractionLost,
                             geneticCode            = geneticCode,
                             acceptStopCodons       = acceptStopCodons,
                             readingFrame           = readingFrame,
                             processorsNum          = processorsNum)
    return(newConsensusReads)
}



### ============================================================================
### Self-defined constructor for SangerSingleRead
### ============================================================================
#' @description
#'
#' @param readFeature .
#' @param readFileName .
#' @param TrimmingMethod .
#' @param M1TrimmingCutoff .
#' @param M2CutoffQualityScore .
#' @param M2SlidingWindowSize .
#' @param baseNumPerRow .
#' @param heightPerRow .
#' @param signalRatioCutoff .
#' @param showTrimmed .
#'
#' @return SangerSingleRead
#' @export
#' @author Kuan-Hao Chao
#' @example
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <-
#'             file.path(inputFilesPath,
#'                       "Allolobophora_chlorotica",
#'                       "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
#' A_chloroticaSingleRead <-
#'        SangerSingleRead(readFeature           = "Forward Read",
#'                         readFileName          = A_chloroticaFdReadFN,
#'                         geneticCode           = GENETIC_CODE,
#'                         TrimmingMethod        = "M2",
#'                         M1TrimmingCutoff      = NULL,
#'                         M2CutoffQualityScore  = 40,
#'                         M2SlidingWindowSize   = 10,
#'                         baseNumPerRow         = 100,
#'                         heightPerRow          = 200,
#'                         signalRatioCutoff     = 0.33,
#'                         showTrimmed           = TRUE)
SangerSingleRead <- function(readFeature           = character(0),
                             readFileName          = character(0),
                             geneticCode           = GENETIC_CODE,
                             TrimmingMethod        = "M2",
                             M1TrimmingCutoff      = NULL,
                             M2CutoffQualityScore  = 40,
                             M2SlidingWindowSize   = 10,
                             baseNumPerRow         = 100,
                             heightPerRow          = 200,
                             signalRatioCutoff     = 0.33,
                             showTrimmed           = TRUE) {
    newSingleReads <- new("SangerSingleRead",
                          readFeature          = readFeature,
                          readFileName         = readFileName,
                          geneticCode          = geneticCode,
                          TrimmingMethod       = TrimmingMethod,
                          M1TrimmingCutoff     = M1TrimmingCutoff,
                          M2CutoffQualityScore = M2CutoffQualityScore,
                          M2SlidingWindowSize  = M2SlidingWindowSize,
                          baseNumPerRow        = baseNumPerRow,
                          heightPerRow         = heightPerRow,
                          signalRatioCutoff    = signalRatioCutoff,
                          showTrimmed          = showTrimmed)
    return(newSingleReads)
}
