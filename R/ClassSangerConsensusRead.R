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
#'                                  cutoffQualityScore    = 20,
#'                                  slidingWindowSize     = 8)
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

    errors <- checkCutoffQualityScore(cutoffQualityScore, errors)
    errors <- checkSlidingWindowSize(slidingWindowSize, errors)
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
                                   cutoffQualityScore, slidingWindowSize)
        reverseReadsList <- sapply(reverseAllReads[[1]], SangerSingleRead,
                                   readFeature = "Reverse Read",
                                   cutoffQualityScore, slidingWindowSize)

        ### --------------------------------------------------------------------
        ### forward & reverse character reads list string creation
        ### --------------------------------------------------------------------
        fRDNAStringSet <- sapply(forwardReadsList, function(forwardRead) {
            as.character(primarySeq(forwardRead))
        })
        rRDNAStringSet <- sapply(reverseReadsList, function(reverseRead) {
            as.character(reverseComplement(primarySeq(reverseRead)))
        })

        ### --------------------------------------------------------------------
        ### DNAStringSet storing forward & reverse reads ! (Origin)
        ### --------------------------------------------------------------------
        frReadSet <- DNAStringSet(c(unlist(fRDNAStringSet),
                                          unlist(rRDNAStringSet)))
        frReadFeatureList <- c(rep("Forward Reads", length(fRDNAStringSet)),
                               rep("Reverse Reads", length(rRDNAStringSet)))

        if(length(frReadSet) < 2) {
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
            # no_N_string <- str_replace_all(frReadSet[1], "N", "T")
            # example.dna <- DNAStringSet(c(`IGHV1-18*01`=no_N_string))
            # refAminoAcidSeq <- translate(example.dna)
            # Verbose should be FALSE, but I get error when calling it
            corrected =
                CorrectFrameshifts(myXStringSet = frReadSet,
                                   myAAStringSet = AAStringSet(refAminoAcidSeq),
                                   geneticCode = geneticCode,
                                   type = 'both',
                                   processors = processorsNum)
            frReadSet = corrected$sequences
            indels = getIndelDf(corrected$indels)
            stops = as.numeric(unlist(mclapply(frReadSet, countStopSodons,
                                               readingFrame, geneticCode,
                                               mc.cores = processorsNum)))
            stopsDf = data.frame("read" = names(frReadSet),
                                 "stop.codons" = stops)
            frReadSetLen = unlist(lapply(frReadSet, function(x) length(x)))
            frReadSet = frReadSet[which(frReadSetLen>0)]
        } else {
            indels = data.frame()
            stopsDf = data.frame()
        }

        if(length(frReadSet) < 2) {
            error <- paste("\n'After running 'CorrectFrameshifts' function, ",
                           "forward and reverse reads should be more than 2.\n",
                           sep = "")
            stop(error)
        }

        ### --------------------------------------------------------------------
        ### Reads with stop codons elimination
        ### --------------------------------------------------------------------
        if (!acceptStopCodons) {
            ### ----------------------------------------------------------------
            ### Remove reads with stop codons
            ### ----------------------------------------------------------------
            print("Removing reads with stop codons")
            if(refAminoAcidSeq == ""){ # otherwise we already did it above
                stops =
                    as.numeric(unlist(mclapply(frReadSet,
                                               countStopSodons,
                                               readingFrame, geneticCode,
                                               mc.cores = processorsNum)))
                stopsDf = data.frame("read" = names(frReadSet),
                                     "stopCodons" = stops)
            }
            old_length = length(frReadSet)
            frReadSet = frReadSet[which(stops==0)]
            # Modify
            message(old_length - length(frReadSet),
                    "reads with stop codons removed")
        }

        if(length(frReadSet) < 2) {
            error <- paste("\n'After removing reads with stop codons, ",
                           "forward and reverse reads should be more than 2.\n",
                           sep = "")
            stop(error)
        }

        ### --------------------------------------------------------------------
        ### Start aligning reads
        ### --------------------------------------------------------------------
        if (refAminoAcidSeq != "") {
            aln = AlignTranslation(frReadSet, geneticCode = geneticCode,
                                   processors = processorsNum, verbose = FALSE)
        } else {
            aln = AlignSeqs(frReadSet,
                            processors = processorsNum, verbose = FALSE)
        }
        names(aln) = paste(1:length(aln), "Read",
                           basename(names(aln)), sep="_")

        print("Calling consensus sequence")
        consensus = ConsensusSequence(aln,
                                      minInformation = minFractionCall,
                                      includeTerminalGaps = TRUE,
                                      ignoreNonBases = TRUE,
                                      threshold = maxFractionLost,
                                      noConsensusChar = "-",
                                      ambiguity = TRUE
        )[[1]]

        print("Calculating differences between reads and consensus")
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

