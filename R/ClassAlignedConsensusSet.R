#' @title SangerAlignedConsensusSet
#'
#' @description  An S4 class for storing multiple single consensus reads to
#'
#' @slot parentDirectory .
#' @slot suffixForwardRegExp .
#' @slot suffixReverseRegExp .
#' @slot consensusReadsList .
#' @slot SCconsensusRead .
#' @slot SCalignment .
#' @slot SCdifferencesDF .
#' @slot SCdistanceMatrix .
#' @slot SCdendrogram .
#' @slot SCindelsDF .
#' @slot SCstopCodonsDF .
#' @slot SCsecondaryPeakDF .
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
#'                      TrimmingMethod        = "M2",
#'                      M1TrimmingCutoff      = NULL,
#'                      M2CutoffQualityScore  = 40,
#'                      M2SlidingWindowSize   = 10,
#'                      baseNumPerRow         = 100,
#'                      heightPerRow          = 200,
#'                      signalRatioCutoff     = 0.33,
#'                      showTrimmed           = TRUE)
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
             consensusReadsList          = "list",
             SCconsensusRead             = "DNAString",
             SCalignment                 = "DNAStringSet",
             SCdifferencesDF             = "data.frame",
             SCdistanceMatrix            = "matrix",
             SCdendrogram                = "list",
             SCindelsDF                  = "data.frame",
             SCstopCodonsDF              = "data.frame",
             SCsecondaryPeakDF           = "data.frame"
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

    SangerConsensusReadList <- sapply(consensusReadsName,
                                      function(eachConsRead) {
        SangerConsensusRead(parentDirectory, eachConsRead,
                            suffixForwardRegExp, suffixReverseRegExp,
                            TrimmingMethod, M1TrimmingCutoff,
                            M2CutoffQualityScore, M2SlidingWindowSize,
                            baseNumPerRow, heightPerRow, signalRatioCutoff,
                            showTrimmed, refAminoAcidSeq, minReadsNum,
                            minReadLength, minFractionCall, maxFractionLost,
                            geneticCode, acceptStopCodons,
                            readingFrame, processorsNum)
    })

    message("Filtering readsets with < ", minReadsNum, " reads...")

    # sapply(consensusReadsName, function(eachCSName) {
    #     eachCSName
    #
    #     grepl(suffixForwardRegExp,
    #           parentDirFiles)
    # })
    #
    # SangerConsensusReadList


    # calculateConsensusRead (forwardReadsList, reverseReadsList,
    #                         refAminoAcidSeq, minFractionCall,
    #                         maxFractionLost, geneticCode,
    #                         acceptStopCodons, readingFrame)

















    SangerConsensusReadDNAList <- sapply(SangerConsensusReadList, function(SangerConsensusRead) {
        as.character(SangerConsensusRead@consensusRead)
    })

    SangerConsensusReadDNASet <- DNAStringSet(SangerConsensusReadDNAList)

    ### --------------------------------------------------------------------
    ### DNAStringSet storing forward & reverse reads ! (Origin)
    ### --------------------------------------------------------------------

    if(length(SangerConsensusReadDNASet) < 2) {
        error <- paste("\n'Valid abif files should be more than 2.\n",
                       sep = "")
        stop(error)
    }
    processorsNum <- getProcessors(processorsNum)

    ### --------------------------------------------------------------------
    ### Amino acid reference sequence CorrectFrameshifts correction
    ### --------------------------------------------------------------------
    if (refAminoAcidSeq != "") {
        message("Correcting frameshifts in reads using amino acid",
                "reference sequence")
        # My test refAminoAcidSeq data
        # no_N_string <- str_replace_all(SangerConsensusReadDNASet[1], "N", "T")
        # example.dna <- DNAStringSet(c(`IGHV1-18*01`=no_N_string))
        # refAminoAcidSeq <- translate(example.dna)
        # Verbose should be FALSE, but I get error when calling it
        corrected =
            CorrectFrameshifts(myXStringSet = SangerConsensusReadDNASet,
                               myAAStringSet = AAStringSet(refAminoAcidSeq),
                               geneticCode = geneticCode,
                               type = 'both',
                               processors = processorsNum)
        SangerConsensusReadDNASet = corrected$sequences
        indels = getIndelDf(corrected$indels)
        stops = as.numeric(unlist(mclapply(SangerConsensusReadDNASet, countStopSodons,
                                           readingFrame, geneticCode,
                                           mc.cores = processorsNum)))
        stopsDf = data.frame("read" = names(SangerConsensusReadDNASet),
                             "stop.codons" = stops)
        SangerConsensusReadDNASetLen = unlist(lapply(SangerConsensusReadDNASet, function(x) length(x)))
        SangerConsensusReadDNASet = SangerConsensusReadDNASet[which(SangerConsensusReadDNASetLen>0)]
    } else {
        indels = data.frame()
        stopsDf = data.frame()
    }
    if(length(SangerConsensusReadDNASet) < 2) {
        error <- paste("\n'After running 'CorrectFrameshifts' function, ",
                       "forward and reverse reads should be more than 2.\n",
                       sep = "")
        stop(error)
    }

    ### --------------------------------------------------------------------
    ### Reads with stop codons elimination
    ### --------------------------------------------------------------------
    ### ----------------------------------------------------------------
    ### Remove reads with stop codons
    ### ----------------------------------------------------------------
    if (!acceptStopCodons) {
        print("Removing reads with stop codons")
        if(refAminoAcidSeq == ""){ # otherwise we already did it above
            stops =
                as.numeric(unlist(mclapply(SangerConsensusReadDNASet,
                                           countStopSodons,
                                           readingFrame, geneticCode,
                                           mc.cores = processorsNum)))
            stopsDf = data.frame("read" = names(SangerConsensusReadDNASet),
                                 "stopCodons" = stops)
        }
        old_length = length(SangerConsensusReadDNASet)
        SangerConsensusReadDNASet = SangerConsensusReadDNASet[which(stops==0)]
        # Modify
        message(old_length - length(SangerConsensusReadDNASet),
                "reads with stop codons removed")
    }

    if(length(SangerConsensusReadDNASet) < 2) {
        error <- paste("\n'After removing reads with stop codons, ",
                       "forward and reverse reads should be more than 2.\n",
                       sep = "")
        stop(error)
    }

    ### --------------------------------------------------------------------
    ### Start aligning reads
    ### --------------------------------------------------------------------
    if (refAminoAcidSeq != "") {
        aln = AlignTranslation(SangerConsensusReadDNASet, geneticCode = geneticCode,
                               processors = processorsNum, verbose = FALSE)
    } else {
        aln = AlignSeqs(SangerConsensusReadDNASet,
                        processors = processorsNum, verbose = FALSE)
    }
    names(aln) = paste(1:length(aln), "Read",
                       basename(names(aln)), sep="_")
    consensus = ConsensusSequence(aln,
                                  minInformation = minFractionCall,
                                  includeTerminalGaps = TRUE,
                                  ignoreNonBases = TRUE,
                                  threshold = maxFractionLost,
                                  noConsensusChar = "-",
                                  ambiguity = TRUE
    )[[1]]

    diffs = mclapply(aln, nPairwiseDiffs,
                     subject = consensus, mc.cores = processorsNum)
    diffs = do.call(rbind, diffs)
    diffsDf = data.frame("name" = names(aln),
                         "pairwise.diffs.to.consensus" = diffs[,1],
                         "unused.chars" = diffs[,2])
    rownames(diffsDf) = NULL

    # get a dendrogram
    dist = DistanceMatrix(aln, correction = "Jukes-Cantor",
                          penalizeGapLetterMatches = FALSE,
                          processors = processorsNum, verbose = FALSE)
    dend = IdClusters(dist, type = "both",
                      showPlot = FALSE,
                      processors = processorsNum, verbose = FALSE)

    # add consensus to alignment
    aln2 = c(aln, DNAStringSet(consensus))
    names(aln2)[length(aln2)] = "Consensus"
    # strip gaps from consensus (must be an easier way!!)
    consensusGapfree = RemoveGaps(DNAStringSet(consensus))[[1]]

    # count columns in the alignment with >1 coincident secondary peaks
    spDf = countCoincidentSp(aln, processors = processorsNum)
    if (is.null(spDf)) {
        spDf = data.frame()
    }


















    if (length(errors) == 0) {

    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory           = parentDirectory,
                   suffixForwardRegExp       = suffixForwardRegExp,
                   suffixReverseRegExp       = suffixReverseRegExp,
                   consensusReadsList        = SangerConsensusReadList,
                   SCconsensusRead           = consensusGapfree,
                   SCdifferencesDF           = diffsDf,
                   SCalignment               = aln2,
                   SCdistanceMatrix          = dist,
                   SCdendrogram              = dend,
                   SCindelsDF                = indels,
                   SCstopCodonsDF            = stopsDf,
                   SCsecondaryPeakDF         = spDf)
})

