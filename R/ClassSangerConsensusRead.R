#' @title SangerConsensusRead
#'
#' @description  An S4 class for storing multiple single reads to build up new
#'  consensus read
#'
#' @slot parentDirectory
#' @slot consenesusReadName .
#' @slot suffixForwardRegExp .
#' @slot suffixReverseRegExp .
#' @slot forwardReadsList
#' @slot reverseReadsList
#' @slot minReadsNum
#' @slot minReadLength
#' @slot refAminoAcidSeq
#' @slot minFractionCall
#' @slot maxFractionLost
#' @slot geneticCode
#' @slot acceptStopCodons
#' @slot readingFrame
#' @slot consensusRead
#' @slot alignment
#' @slot differencesDF
#' @slot distanceMatrix
#' @slot dendrogram
#' @slot indelsDF
#' @slot stopCodonsDF
#' @slot secondaryPeakDF
#'
#' @name SangerConsensusRead-class
#'
#' @rdname SangerConsensusRead-class
#'
#' @exportClass SangerConsensusRead
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerSingleRead.R
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' consenesusReadName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory       = inputFilesParentDir,
#'                                  consenesusReadName    = consenesusReadName,
#'                                  suffixForwardRegExp   = suffixForwardRegExp,
#'                                  suffixReverseRegExp   = suffixReverseRegExp,
#'                                  TrimmingMethod        = "M2",
#'                                  M1TrimmingCutoff      = NULL,
#'                                  M2CutoffQualityScore  = 40,
#'                                  M2SlidingWindowSize   = 10,
#'                                  baseNumPerRow         = 100,
#'                                  signalRatioCutoff     = 0.33,
#'                                  showTrimmed           = TRUE)
setClass("SangerConsensusRead",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerConsensusRead'
         ### -------------------------------------------------------------------
         representation(
             parentDirectory           = "character",
             consenesusReadName        = "character",
             suffixForwardRegExp       = "character",
             suffixReverseRegExp       = "character",
             forwardReadsList          = "list",
             reverseReadsList          = "list",
             minReadsNum               = "numeric",
             minReadLength             = "numeric",
             refAminoAcidSeq           = "character",
             minFractionCall           = "numeric",
             maxFractionLost           = "numeric",
             geneticCode               = "character",
             acceptStopCodons          = "logical",
             readingFrame              = "numeric",
             consensusRead             = "DNAString",
             alignment                 = "DNAStringSet",
             differencesDF             = "data.frame",
             distanceMatrix            = "matrix",
             dendrogram                = "list",
             indelsDF                  = "data.frame",
             stopCodonsDF              = "data.frame",
             secondaryPeakDF           = "data.frame"
             ),
)

### ============================================================================
### Overwrite initialize for 'SangerConsensusRead' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerConsensusRead",
          function(.Object, ...,
                   parentDirectory        = parentDirectory,
                   consenesusReadName     = consenesusReadName,
                   suffixForwardRegExp    = suffixForwardRegExp,
                   suffixReverseRegExp    = suffixReverseRegExp,
                   TrimmingMethod        = "M1",
                   M1TrimmingCutoff      = 0.0001,
                   M2CutoffQualityScore  = NULL,
                   M2SlidingWindowSize   = NULL,
                   baseNumPerRow         = 100,
                   signalRatioCutoff     = 0.33,
                   showTrimmed           = TRUE,
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
    errors <- checkParentDirectory (parentDirectory, errors)

    ### ------------------------------------------------------------------------
    ### 'forwardAllReads' & 'reverseAllReads' files prechecking
    ### ------------------------------------------------------------------------
    parentDirFiles <- list.files(parentDirectory)
    consensusSubGroupFiles <- parentDirFiles[grepl(consenesusReadName,
                                                    parentDirFiles, fixed=TRUE)]
    forwardSelectInputFiles <- consensusSubGroupFiles[grepl(suffixForwardRegExp,
                                                    consensusSubGroupFiles)]
    reverseSelectInputFiles <- consensusSubGroupFiles[grepl(suffixReverseRegExp,
                                                    consensusSubGroupFiles)]
    forwardAllReads <- lapply(parentDirectory, file.path,
                              forwardSelectInputFiles)
    reverseAllReads <- lapply(parentDirectory, file.path,
                              reverseSelectInputFiles)

    forwardNumber <- length(forwardAllReads[[1]])
    reverseNumber <- length(reverseAllReads[[1]])
    # sapply to check all forwardAllReads files are exist.
    forwardAllErrorMsg <- sapply(c(forwardAllReads[[1]]), function(filePath) {
        if (!file.exists(filePath)) {
            msg <- paste("\n'", filePath, "' forward read file does ",
                         "not exist.\n", sep = "")
            return(msg)
        }
        return()
    })
    reverseAllErrorMsg <- sapply(c(reverseAllReads[[1]]), function(filePath) {
        if (!file.exists(filePath)) {
            msg <- paste("\n'", filePath, "'",
                         " reverse read file does not exist.\n", sep = "")
            return(msg)
        }
        return()
    })
    errors <- c(errors, unlist(forwardAllErrorMsg), use.names = FALSE)
    errors <- c(errors, unlist(reverseAllErrorMsg), use.names = FALSE)

    ### ------------------------------------------------------------------------
    ### Prechecking success. Start to create multiple reads.
    ### ------------------------------------------------------------------------
    if (length(errors) == 0) {
        # sapply to create SangerSingleRead list.

        forwardReadsList <- sapply(forwardAllReads[[1]], SangerSingleRead,
                                   readFeature = "Forward Read",
                                   TrimmingMethod = TrimmingMethod,
                                   M1TrimmingCutoff = M1TrimmingCutoff,
                                   M2CutoffQualityScore = M2CutoffQualityScore,
                                   M2SlidingWindowSize = M2SlidingWindowSize,
                                   baseNumPerRow = baseNumPerRow,
                                   signalRatioCutoff = signalRatioCutoff,
                                   showTrimmed = showTrimmed)

        reverseReadsList <- sapply(reverseAllReads[[1]], SangerSingleRead,
                                   readFeature = "Reverse Read",
                                   TrimmingMethod = TrimmingMethod,
                                   M1TrimmingCutoff = M1TrimmingCutoff,
                                   M2CutoffQualityScore = M2CutoffQualityScore,
                                   M2SlidingWindowSize = M2SlidingWindowSize,
                                   baseNumPerRow = baseNumPerRow,
                                   signalRatioCutoff = signalRatioCutoff,
                                   showTrimmed = showTrimmed)

        CSResult<-
            calculateConsensusRead (forwardReadsList, reverseReadsList,
                                    refAminoAcidSeq, minFractionCall,
                                    maxFractionLost, geneticCode,
                                    acceptStopCodons, readingFrame)

        consensusGapfree <- CSResult$consensusGapfree
        diffsDf <- CSResult$diffsDf
        aln2 <- CSResult$aln2
        dist <- CSResult$dist
        dend <- CSResult$dend
        indels <- CSResult$indels
        stopsDf <- CSResult$stopsDf
        spDf <- CSResult$spDf
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory        = parentDirectory,
                   consenesusReadName     = consenesusReadName,
                   suffixForwardRegExp    = suffixForwardRegExp,
                   suffixReverseRegExp    = suffixReverseRegExp,
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

