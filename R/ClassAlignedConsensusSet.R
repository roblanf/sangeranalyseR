#' @title SangerAlignedConsensusSet
#'
#' @description  An S4 class for storing multiple single consensus reads to
#'
#' @slot parentDirectory
#' @slot suffixRegExp
#' @slot consensusReadsList
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
#'                      parentDirectory       = inputFilesParentDir,
#'                      suffixForwardRegExp   = suffixForwardRegExp,
#'                      suffixReverseRegExp   = suffixReverseRegExp,
#'                      TrimmingMethod        = "M1",
#'                      M1TrimmingCutoff      = 0.0001,
#'                      M2CutoffQualityScore  = NULL,
#'                      M2SlidingWindowSize   = NULL)
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
             parentDirectory           = "character",
             suffixForwardRegExp       = "character",
             suffixReverseRegExp       = "character",
             consensusReadsList        = "list"
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
                   TrimmingMethod        = "M1",
                   M1TrimmingCutoff      = 0.0001,
                   M2CutoffQualityScore  = NULL,
                   M2SlidingWindowSize   = NULL,
                   refAminoAcidSeq        = "",
                   minReadsNum            = 2,
                   minReadLength          = 20,
                   minFractionCall        = 0.5,
                   maxFractionLost        = 0.5,
                   geneticCode            = GENETIC_CODE,
                   acceptStopCodons       = TRUE,
                   readingFrame           = 1,
                   processorsNum          = 1) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()


    ### --------------------------------------------------------------
    ### Input parameter prechecking for TrimmingMethod.
    ### --------------------------------------------------------------
    errors <- checkTrimParam(TrimmingMethod,
                             M1TrimmingCutoff,
                             M2CutoffQualityScore,
                             M2SlidingWindowSize,
                             errors)
    errors <- checkMinReadsNum(minReadsNum, errors)
    errors <- checkMinReadLength(minReadLength, errors)
    errors <- checkMinFractionCall(minFractionCall, errors)
    errors <- checkMaxFractionLost(maxFractionLost, errors)
    errors <- checkReadingFrame(readingFrame, errors)
    errors <- checkGeneticCode(geneticCode, errors)

    ### ------------------------------------------------------------------------
    ### 'parentDirectory' prechecking
    ### ------------------------------------------------------------------------
    if (!file.exists(parentDirectory)) {
        msg <- paste("\n'", parentDirectory, "'",
                     " parent directory does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }

    ### ------------------------------------------------------------------------
    ### 'forwardAllReads' & 'reverseAllReads' files prechecking
    ### ------------------------------------------------------------------------
    parentDirFiles <- list.files(parentDirectory)
    forwardSelectInputFiles <- parentDirFiles[grepl(suffixForwardRegExp,
                                                    parentDirFiles)]
    reverseSelectInputFiles <- parentDirFiles[grepl(suffixReverseRegExp,
                                                    parentDirFiles)]

    # Find possible consensus Name for forward and reverse reads
    forwardConsensusName <- unlist(str_split(forwardSelectInputFiles, suffixForwardRegExp, n = Inf, simplify = FALSE))[c(TRUE, FALSE)]
    reverseConsensusName <- unlist(str_split(reverseSelectInputFiles, suffixReverseRegExp, n = Inf, simplify = FALSE))[c(TRUE, FALSE)]

    consensusReadsName <- union(forwardConsensusName, reverseConsensusName)
    consensusReadsNumber <- length(consensusReadsName)

    # Create consensusReads for all list of consensusReadsNumber

    SangerConsensusReadList <- sapply(consensusReadsName, function(eachConsRead) {
        SangerConsensusRead(parentDirectory, eachConsRead,
                            suffixForwardRegExp, suffixReverseRegExp,
                            TrimmingMethod, M1TrimmingCutoff,
                            M2CutoffQualityScore, M2SlidingWindowSize,
                            refAminoAcidSeq, minReadsNum, minReadLength,
                            minFractionCall, maxFractionLost, geneticCode,
                            acceptStopCodons, readingFrame, processorsNum)
    })

    if (length(errors) == 0) {

    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory           = parentDirectory,
                   suffixForwardRegExp       = suffixForwardRegExp,
                   suffixReverseRegExp       = suffixReverseRegExp,
                   consensusReadsList        = SangerConsensusReadList)
})

