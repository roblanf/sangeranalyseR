#' @title SangerAlignedConsensusSet
#'
#' @description  An S4 class for storing multiple single consensus reads to
#'
#' @slot parentDirectory .
#' @slot suffixForwardRegExp.
#' @slot suffixReverseRegExp .
#' @slot minFractionCallSCSet .
#' @slot maxFractionLostSCSet .
#' @slot geneticCode .
#' @slot refAminoAcidSeq .
#' @slot consensusReadsList .
#' @slot consensusReadSCSet .
#' @slot alignmentSCSet .
#' @slot alignmentTreeSCSet .
#'
#' @name SangerAlignedConsensusSet-class
#'
#' @rdname SangerAlignedConsensusSet-class
#'
#' @exportClass SangerAlignedConsensusSet
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerSingleRead.R
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' SangerAlignedConsensusSet <- new("SangerAlignedConsensusSet",
#'                                parentDirectory       = inputFilesParentDir,
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
setClass("SangerAlignedConsensusSet",
         # Users need to name their ab1 files in a systematic way. Here is the
         # regulation:
         #  1. Naming: XXXXX_F[0-9]*.ab1 / XXXXX_R[0-9]*.ab1
         #        For same consensus reads, XXXXX must be same.
         #  2. Users can set
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerAlignedConsensusSet'
         ### -------------------------------------------------------------------
         representation(
             parentDirectory             = "character",
             suffixForwardRegExp         = "character",
             suffixReverseRegExp         = "character",
             minFractionCallSCSet        = "numeric",
             maxFractionLostSCSet        = "numeric",
             geneticCode                 = "character",
             refAminoAcidSeq             = "character",
             consensusReadsList          = "list",
             consensusReadSCSet          = "DNAString",
             alignmentSCSet              = "DNAStringSet",
             alignmentTreeSCSet          = "phylo"
             ),
)

### ============================================================================
### Overwrite initialize for 'SangerConsensusRead' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerAlignedConsensusSet",
          function(.Object, ...,
                   parentDirectory        = parentDirectory,
                   suffixForwardRegExp    = "_[F]_[0-9]*.ab1",
                   suffixReverseRegExp    = "_[R]_[0-9]*.ab1",
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
                   minFractionCallSCSet   = 0.5,
                   maxFractionLostSCSet   = 0.5,
                   processorsNum          = 1) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()

    ### ------------------------------------------------------------------------
    ### 'parentDirectory' prechecking
    ### ------------------------------------------------------------------------
    errors <- checkParentDirectory (parentDirectory, errors)

    ### ------------------------------------------------------------------------
    ### Input parameter prechecking for TrimmingMethod.
    ### ------------------------------------------------------------------------
    errors <- checkTrimParam(TrimmingMethod,
                             M1TrimmingCutoff,
                             M2CutoffQualityScore,
                             M2SlidingWindowSize,
                             errors)

    ##### ----------------------------------------------------------------------
    ##### Input parameter prechecking for ChromatogramParam
    ##### ----------------------------------------------------------------------
    errors <- checkBaseNumPerRow (baseNumPerRow, errors)
    errors <- checkHeightPerRow (baseNumPerRow, errors)
    errors <- checkSignalRatioCutoff (signalRatioCutoff, errors)
    errors <- checkShowTrimmed (showTrimmed, errors)

    ##### ----------------------------------------------------------------------
    ##### Input parameter prechecking for ConsensusRead parameter
    ##### ----------------------------------------------------------------------
    errors <- checkMinReadsNum(minReadsNum, errors)
    errors <- checkMinReadLength(minReadLength, errors)
    errors <- checkMinFractionCall(minFractionCall, errors)
    errors <- checkMaxFractionLost(maxFractionLost, errors)
    errors <- checkGeneticCode(geneticCode, errors)
    errors <- checkAcceptStopCodons(acceptStopCodons, errors)
    errors <- checkReadingFrame(readingFrame, errors)

    ##### ----------------------------------------------------------------------
    ##### Input parameter prechecking for processorsNum
    ##### ----------------------------------------------------------------------
    errors <- checkProcessorsNum(processorsNum, errors)

    if (length(errors) == 0) {
        ### --------------------------------------------------------------------
        ### Automatically finding consensus read name by forward&reverse suffix
        ### --------------------------------------------------------------------
        parentDirFiles <- list.files(parentDirectory)
        forwardSelectInputFiles <- parentDirFiles[grepl(suffixForwardRegExp,
                                                        parentDirFiles)]
        reverseSelectInputFiles <- parentDirFiles[grepl(suffixReverseRegExp,
                                                        parentDirFiles)]

        # Find possible consensus Name for forward and reverse reads
        forwardConsensusName <-
            unlist(str_split(forwardSelectInputFiles, suffixForwardRegExp,
                             n = Inf, simplify = FALSE))[c(TRUE, FALSE)]
        reverseConsensusName <-
            unlist(str_split(reverseSelectInputFiles, suffixReverseRegExp,
                             n = Inf, simplify = FALSE))[c(TRUE, FALSE)]

        consensusReadsName <- union(forwardConsensusName, reverseConsensusName)
        consensusReadsNumber <- length(consensusReadsName)

        # Create consensusReads for all list of consensusReadsNumber
        ### --------------------------------------------------------------------
        ### Creating each SangerConsensusRead (store as SangerConsensusReadList)
        ### --------------------------------------------------------------------
        SangerConsensusReadList <-
            sapply(consensusReadsName,
                   function(eachConsRead) {
                       SangerConsensusRead(
                           parentDirectory, eachConsRead,
                           suffixForwardRegExp, suffixReverseRegExp,
                           TrimmingMethod, M1TrimmingCutoff,
                           M2CutoffQualityScore, M2SlidingWindowSize,
                           baseNumPerRow, heightPerRow, signalRatioCutoff,
                           showTrimmed, refAminoAcidSeq, minReadsNum,
                           minReadLength, minFractionCall, maxFractionLost,
                           geneticCode, acceptStopCodons,
                           readingFrame, processorsNum)
                   })

        acResult <-
            alignConsensusReads(SangerConsensusReadList,
                                geneticCode, refAminoAcidSeq,
                                minFractionCallSCSet, maxFractionLostSCSet,
                                processorsNum)
        consensus <- acResult[["consensus"]]
        aln <- acResult[["aln"]]
        aln.tree <- acResult[["aln.tree"]]
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory           = parentDirectory,
                   suffixForwardRegExp       = suffixForwardRegExp,
                   suffixReverseRegExp       = suffixReverseRegExp,
                   consensusReadsList        = SangerConsensusReadList,
                   minFractionCallSCSet      = minFractionCallSCSet,
                   maxFractionLostSCSet      = maxFractionLostSCSet,
                   geneticCode               = geneticCode,
                   consensusReadSCSet        = consensus,
                   refAminoAcidSeq           = refAminoAcidSeq,
                   alignmentSCSet            = aln,
                   alignmentTreeSCSet        = aln.tree
                   )
})



