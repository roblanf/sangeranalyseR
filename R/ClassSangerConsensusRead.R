#' @title SangerConsensusRead
#'
#' @description  An S4 class for storing multiple single reads to build up new
#'  consensus read
#'
#' @slot parentDirectory .
#' @slot readsRegularExp .
#' @slot cutoffQualityScore .
#' @slot slidingWindowSize .
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
#' samplesRegExp <- "ACHL"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory = inputFilesParentDir,
#'                                  readsRegularExp = samplesRegExp,
#'                                  cutoffQualityScore  = 50L,
#'                                  slidingWindowSize   = 8L)
setClass("SangerConsensusRead",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerConsensusRead'
         ### -------------------------------------------------------------------
         representation(
             parentDirectory    = "character",
             readsRegularExp    = "character",
             SangerReadsList    = "list"
         ),
)

### ============================================================================
### Overwrite initialize for 'SangerConsensusRead' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerConsensusRead",
          function(.Object, ...,
                   parentDirectory      = parentDirectory,
                   readsRegularExp      = readsRegularExp,
                   cutoffQualityScore   = 20L,
                   slidingWindowSize    = 5L) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()
    if (!file.exists(parentDirectory)) {
        msg <- paste("\n'", parentDirectory, "'",
                     " parent directory does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }

    parentDirFiles <- list.files(parentDirectory)
    selectInputFiles <- parentDirFiles[grepl(readsRegularExp, parentDirFiles)]
    allReads <- lapply(parentDirectory, file.path, selectInputFiles)
    readsNumber <- length(allReads[[1]])
    # sapply to check all files are exist.
    allErrorMsg <- sapply(c(allReads[[1]]), function(filePath) {
        if (!file.exists(filePath)) {
            msg <- paste("\n'", filePath, "'",
                         " read file does not exist.\n", sep = "")
            return(msg)
        }
        return()
    })
    errors <- c(errors, unlist(allErrorMsg), use.names = FALSE)

    ### ------------------------------------------------------------------------
    ### Prechecking success. Start to create multiple reads.
    ### ------------------------------------------------------------------------
    if (length(errors) == 0) {
        # sapply to create SangerSingleRead list.
        SangerSingleReadList <- sapply(allReads[[1]], SangerSingleRead,
               readFeature = "Reads", cutoffQualityScore, slidingWindowSize)


        # # Try to correct frameshifts in the input sequences
        # if(!is.null(ref.aa.seq)) {
        #
        #     print("Correcting frameshifts in reads using amino acid reference sequence")
        #     corrected = CorrectFrameshifts(myXStringSet = readset, myAAStringSet = AAStringSet(ref.aa.seq), geneticCode = genetic.code, type = 'both', processors = processors, verbose = FALSE)
        #     readset = corrected$sequences
        #     indels = get.indel.df(corrected$indels)
        #     stops = as.numeric(unlist(mclapply(readset, count.stop.codons, reading.frame, genetic.code, mc.cores = processors)))
        #     stops.df = data.frame("read" = names(readset), "stop.codons" = stops)
        #     readset.lengths = unlist(lapply(readset, function(x) length(x)))
        #     readset = readset[which(readset.lengths>0)]
        # }else{
        #     indels = NULL
        #     stops.df = NULL
        # }


        # # Remove reads with stop codons
        # if(accept.stop.codons == FALSE){
        #     print("Removing reads with stop codons")
        #     if(is.null(ref.aa.seq)){ # otherwise we already did it above
        #         stops = as.numeric(unlist(mclapply(readset, count.stop.codons, reading.frame, genetic.code, mc.cores = processors)))
        #         stops.df = data.frame("read" = names(readset), "stop.codons" = stops)
        #     }
        #     old_length = length(readset)
        #     readset = readset[which(stops==0)]
        #     print(sprintf("%d reads with stop codons removed", old_length - length(readset)))
        # }


        # SangerSingleRead list number checking
        if (length(SangerSingleReadList) < 2) {
            stop("There are only ", length(SangerSingleReadList), " reads",
                 " in this directory. (requires at least two reads)")
        }

        print("Aligning reads")


        # if(!is.null(ref.aa.seq)){
        #     aln = AlignTranslation(readset, geneticCode = genetic.code, processors = processors, verbose = FALSE)
        # }else{
        #     aln = AlignSeqs(readset, processors = processors, verbose = FALSE)
        # }





        SangerSingleReadNum <- length(SangerSingleReadList)
        SangerDNAStringList <- sapply(1:SangerSingleReadNum, function(i)
            as.character(primarySeq(SangerSingleReadList[[i]])))
        SangerDNAStringSet <- DNAStringSet(SangerDNAStringList)

        aln = AlignSeqs(SangerDNAStringSet, verbose = FALSE)


        if(is.null(names(aln))){
            names(aln) = paste("read", 1:length(aln), sep="_")
        }

        # call consensus
        print("Calling consensus sequence")
        consensus = ConsensusSequence(aln,
                                      # minInformation = minInformation,
                                      includeTerminalGaps = TRUE,
                                      ignoreNonBases = TRUE,
                                      # threshold = threshold,
                                      noConsensusChar = "-",
                                      ambiguity = TRUE
        )[[1]]

        # print("Calculating differences between reads and consensus")
        # diffs = mclapply(aln, n.pairwise.diffs, subject = consensus, mc.cores = processors)
        # diffs = do.call(rbind, diffs)
        # diffs.df = data.frame("name" = names(aln), "pairwise.diffs.to.consensus" = diffs[,1], "unused.chars" = diffs[,2])
        # rownames(diffs.df) = NULL
        #
        # # get a dendrogram
        # dist = DistanceMatrix(aln, correction = "Jukes-Cantor", penalizeGapLetterMatches = FALSE, processors = processors, verbose = FALSE)
        # dend = IdClusters(dist, type = "dendrogram", processors = processors, verbose = FALSE)
        #
        # # add consensus to alignment
        # aln2 = c(aln, DNAStringSet(consensus))
        # names(aln2)[length(aln2)] = "consensus"
        # # strip gaps from consensus (must be an easier way!!)
        # consensus.gapfree = RemoveGaps(DNAStringSet(consensus))[[1]]
        #
        # # count columns in the alignment with >1 coincident secondary peaks
        # sp.df = count.coincident.sp(aln, processors = processors)
        #
        # merged.read = list("consensus" = consensus.gapfree,
        #                    "alignment" = aln2,
        #                    "differences" = diffs.df,
        #                    "distance.matrix" = dist,
        #                    "dendrogram" = dend,
        #                    "indels" = indels,
        #                    "stop.codons" = stops.df,
        #                    "secondary.peak.columns" = sp.df)
        #
        # class(merged.read) = "merged.read"




























    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory = parentDirectory,
                   readsRegularExp = readsRegularExp,
                   SangerReadsList = SangerSingleReadList)
})
