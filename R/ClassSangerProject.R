#' @title SangerProject
#'
#' @description  An S4 class for storing multiple single consensus reads to
#'
#' @slot parentDirectory
#' @slot suffixRegExp
#' @slot consensusReadsList
#'
#' @name SangerProject-class
#'
#' @rdname SangerProject-class
#'
#' @exportClass SangerProject
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerSingleRead.R
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' forwardRegExp <- "^ACHLO([0-9]*)-09\\[LCO1490_t1,HCO2198_t1\\]_F.ab1$"
#' reverseRegExp <- "^ACHLO([0-9]*)-09\\[LCO1490_t1,HCO2198_t1\\]_R.ab1$"
setClass("SangerProject",
         # Users need to name their ab1 files in a systematic way. Here is the
         # regulation:
         #  1. Naming: XXXXX_F[0-9]*.ab1 / XXXXX_R[0-9]*.ab1
         #        For same consensus reads, XXXXX must be same.
         #  2. Users can set
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerConsensusRead'
         ### -------------------------------------------------------------------
         representation(
             parentDirectory           = "character",
             suffixRegExp              = "character",
             consensusReadsList        = "list"
             ),
)

### ============================================================================
### Overwrite initialize for 'SangerConsensusRead' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerProject",
          function(.Object, ...,
                   parentDirectory        = parentDirectory,
                   suffixForwardRegExp    = "_[F][0-9]*.ab1",
                   suffixReverseRegExp    = "_[R][0-9]*.ab1",
                   cutoffQualityScore     = 20,
                   slidingWindowSize      = 5,
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

    # errors <- checkCutoffQualityScore(cutoffQualityScore, errors)
    # errors <- checkSlidingWindowSize(slidingWindowSize, errors)
    # errors <- checkMinReadsNum(minReadsNum, errors)
    # errors <- checkMinReadLength(minReadLength, errors)
    # errors <- checkMinFractionCall(minFractionCall, errors)
    # errors <- checkMaxFractionLost(maxFractionLost, errors)
    # errors <- checkReadingFrame(readingFrame, errors)
    # errors <- checkGeneticCode(geneticCode, errors)

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
    forwardSelectInputFiles <- parentDirFiles[grepl(forwardReadsRegularExp,
                                                    parentDirFiles)]
    reverseSelectInputFiles <- parentDirFiles[grepl(reverseReadsRegularExp,
                                                    parentDirFiles)]
    forwardAllReads <- lapply(parentDirectory, file.path,
                              forwardSelectInputFiles)
    reverseAllReads <- lapply(parentDirectory, file.path,
                              reverseSelectInputFiles)

    forwardNumber <- length(forwardAllReads[[1]])
    reverseNumber <- length(reverseAllReads[[1]])
    if (length(errors) == 0) {

    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory        = parentDirectory,
                   forwardReadsRegularExp = forwardReadsRegularExp,
                   reverseReadsRegularExp = reverseReadsRegularExp,
                   forwardReadsList       = forwardReadsList,
                   reverseReadsList       = reverseReadsList,
                   minReadsNum            = minReadsNum,
                   minReadLength          = minReadLength,
                   refAminoAcidSeq        = refAminoAcidSeq,
                   minFractionCall        = minFractionCall,
                   maxFractionLost        = maxFractionLost,
                   geneticCode            = geneticCode,
                   acceptStopCodons       = acceptStopCodons,
                   readingFrame           = readingFrame,
                   consensusRead          = consensusGapfree,
                   differencesDF          = diffsDf,
                   alignment              = aln2,
                   distanceMatrix         = dist,
                   dendrogram             = dend,
                   indelsDF               = indels,
                   stopCodonsDF           = stopsDf,
                   secondaryPeakDF        = spDf)
})

