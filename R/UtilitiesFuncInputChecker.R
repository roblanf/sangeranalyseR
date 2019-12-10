checkReadFileName <- function(readFileName, errors) {
    if (!file.exists(readFileName)) {
        cat ("readFileName", readFileName)
        msg <- paste("\n'", readFileName, "'",
                     " foward read file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
}

### ============================================================================
### QualityReport related: 'readFeature', 'qualityPhredScores'
### ============================================================================
checkReadFeature <- function(readFeature, errors) {
    if (readFeature != "Forward Read" && readFeature != "Reverse Read") {
        msg <- paste("\n'readFeature' must be
                     'Forward Read' or 'Reverse Read'\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkQualityPhredScores <- function(qualityPhredScores, errors) {
    if (length(qualityPhredScores) == 0) {
        msg <- paste("\n'qualityPhredScores'
                               length cannot be zero.\n")
        errors <- c(errors, msg)
    }
    if (!all(qualityPhredScores%%1 == 0)) {
        msg <- paste("\nAll elements in 'qualityPhredScores' vector
                               must be integer.\n")
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### Quality trimming related: 'TrimmingMethod', 'M1TrimmingCutoff',
###                           'M2CutoffQualityScore', 'M2SlidingWindowSize'
### ============================================================================
checkTrimParam <- function(TrimmingMethod, M1TrimmingCutoff,
                           M2CutoffQualityScore, M2SlidingWindowSize, errors) {
    if (TrimmingMethod == "M1") {
        if (!is.numeric(M1TrimmingCutoff)) {
            msg<- paste("\n'M1TrimmingCutoff' must be numeric",
                        "(You choose M1).\n")
            errors <- c(errors, msg)
        } else {
            # Ristriction about M1TrimmingCutoff !
            # if (M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
            #     M2CutoffQualityScore%%1!=0) {
            #     msg <- paste("\n'Your input M2CutoffQualityScore is: ",
            #                  M2CutoffQualityScore, "' is invalid.",
            #                  "'M2CutoffQualityScore' should",
            #                  "be between 0 and 60.\n", sep = "")
            #     errors <- c(errors, msg)
            # }
        }
        if (!is.null(M2CutoffQualityScore)) {
            msg<- paste("\n'M2CutoffQualityScore' must be null",
                        "(You choose M1).\n")
            errors <- c(errors, msg)
        }
        if (!is.null(M2SlidingWindowSize)) {
            msg<- paste("\n'M2SlidingWindowSize' must be null",
                        "(You choose M1).\n")
            errors <- c(errors, msg)
        }
    } else if (TrimmingMethod == "M2") {
        if (!is.null(M1TrimmingCutoff)) {
            msg<- paste("\n'M1TrimmingCutoff' must be null",
                        "(You choose M2).\n")
            errors <- c(errors, msg)
        }
        if (!is.numeric(M2CutoffQualityScore)) {
            msg<- paste("\n'M2CutoffQualityScore' must be numeric",
                        "(You choose M2).\n")
            errors <- c(errors, msg)
        } else {
            if (M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
                M2CutoffQualityScore%%1!=0) {
                msg <- paste("\n'Your input M2CutoffQualityScore is: ",
                             M2CutoffQualityScore, "' is invalid.",
                             "'M2CutoffQualityScore' should",
                             "be between 0 and 60.\n", sep = "")
                errors <- c(errors, msg)
            }
        }
        if (!is.numeric(M2SlidingWindowSize)) {
            msg<- paste("\n'M2SlidingWindowSize' must be numeric",
                        "(You choose M2).\n")
            errors <- c(errors, msg)
        } else {
            if (M2SlidingWindowSize > 20 || M2SlidingWindowSize < 0 ||
                M2SlidingWindowSize%%1!=0) {
                msg <- paste("\n'Your input M2SlidingWindowSize is: ",
                             M2SlidingWindowSize, "' is invalid.",
                             "'M2SlidingWindowSize' should",
                             "be between 0 and 20.\n", sep = "")
                errors <- c(errors, msg)
            }
        }
    } else {
        msg <- paste("\n'TrimmingMethod' must be 'M1' or 'M2'.\n")
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### ConsensusRead related: 'minReadsNum', 'minReadLength', 'minFractionCall'
###                        'maxFractionLost' prechecking
### ============================================================================
checkMinReadsNum <- function(minReadsNum, errors) {
    if (minReadsNum%%1!=0) {
        msg <- paste("\n'minReadsNum' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMinReadLength <- function(minReadLength, errors) {
    if (minReadLength%%1!=0) {
        msg <- paste("\n'minReadLength' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}



checkMinFractionCall <- function(minFractionCall, errors) {
    if (minFractionCall > 1 || minFractionCall < 0) {
        msg <- paste("\n'minFractionCall' must be between 0 and 1.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMaxFractionLost <- function(maxFractionLost, errors) {
    if (maxFractionLost > 1 || maxFractionLost < 0) {
        msg <- paste("\n'maxFractionLost' must be between 0 and 1.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}


checkGeneticCode <- function(geneticCode, errors) {
    if(!("*" %in% geneticCode)) {
        msg <- paste("\n'geneticCode' does not specify any stop codons.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}


checkReadingFrame <- function(readingFrame, errors) {
    if(!readingFrame %in% c(1,2,3)) {
        msg <- paste("\n'readingFrame' must be 1, 2, or 3.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkAcceptStopCodons <- function(acceptStopCodons, errors) {
    if (!is.logical(acceptStopCodons)) {
        msg <- paste("\n'acceptStopCodons' must be 'TRUE' or 'FALSE'\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkProcessorsNum <- function(processorsNum, errors) {
    if (!(processorsNum %% 1 == 0)) {
        msg <- paste("\n'processorsNum' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### 'parentDirectory' prechecking
### ============================================================================
checkParentDirectory <- function(parentDirectory, errors) {
    if (!dir.exists(parentDirectory)) {
        msg <- paste("\n'", parentDirectory, "'",
                     " parent directory does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
}

### ============================================================================
### 'baseNumPerRow', 'signalRatioCutoff', 'showTrimmed' prechecking
### ============================================================================
checkBaseNumPerRow <- function(baseNumPerRow, errors) {
    if (baseNumPerRow%%1!=0) {
        msg <- paste("\n'baseNumPerRow' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (baseNumPerRow < 0 || baseNumPerRow > 200) {
        msg <- paste("\n'baseNumPerRow' must be between 0 and 200.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkHeightPerRow <- function(heightPerRow, errors) {
    if (heightPerRow%%1!=0) {
        msg <- paste("\n'heightPerRow' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (heightPerRow < 50 || heightPerRow > 600) {
        msg <- paste("\n'heightPerRow' must be between 0 and 200.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### MakeBaseCalls Utilities function
### ============================================================================
checkSignalRatioCutoff <- function(signalRatioCutoff, errors) {
    if (signalRatioCutoff < 0 || signalRatioCutoff > 1) {
        msg <- paste("\n'signalRatioCutoff' must be between 0 and 1.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkShowTrimmed <- function(showTrimmed, errors) {
    if (!is.logical(showTrimmed)) {
        msg <- paste("\n'showTrimmed' must be between TRUE and FALSE.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}



### ============================================================================
### Quality trimming related parameter
### ============================================================================
M1inside_calculate_trimming <- function(qualityPhredScores,
                                        qualityBaseScores,
                                        M1TrimmingCutoff) {
    rawSeqLength <- length(qualityBaseScores)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    start = FALSE
    trimmedStartPos = 0
    qualityBaseScoresCutOff = M1TrimmingCutoff - qualityBaseScores
    ### ------------------------------------------------------------------------
    ### calculate cummulative score
    ### if cumulative value < 0, set it to 0
    ### the BioPython implementation always trims the first base,
    ### this implementation does not.
    ### ------------------------------------------------------------------------
    score = qualityBaseScoresCutOff[1]
    if(score < 0){
        score = 0
    }else{
        trimmedStartPos = 1
        start = TRUE
    }
    cummul_score = c(score)
    ### ------------------------------------------------------------------------
    ### trimmedStartPos = value when cummulative score is first > 0
    ### ------------------------------------------------------------------------
    ### ------------------------------------------------------------------------
    ### trimmedFinishPos = index of highest cummulative score,
    ### marking the end of sequence segment with highest cummulative score
    ### ------------------------------------------------------------------------
    for(i in 2:length(qualityBaseScoresCutOff)){
        score = cummul_score[length(cummul_score)] + qualityBaseScoresCutOff[i]
        if (score <= 0) {
            cummul_score = c(cummul_score, 0)
        }else{
            cummul_score = c(cummul_score, score)
            if(start == FALSE){
                trimmedStartPos = i
                start = TRUE
            }
        }
        trimmedFinishPos = which.max(cummul_score)
    }
    ### ------------------------------------------------------------------------
    ### fix an edge case, where all scores are worse than the cutoff
    ### in this case you wouldn't want to keep any bases at all
    ### ------------------------------------------------------------------------
    if(sum(cummul_score)==0){trimmedFinishPos = 0}
    if (trimmedFinishPos - trimmedStartPos == 0) {
        trimmedStartPos = 1
        trimmedFinishPos = 2
    }
    trimmedSeqLength = trimmedFinishPos - trimmedStartPos
    trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
    trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
    trimmedMinQualityScore <- min(trimmedQualityPhredScore)
    remainingRatio = trimmedSeqLength / rawSeqLength

    return(c("rawSeqLength" = rawSeqLength,
             "rawMeanQualityScore" = rawMeanQualityScore,
             "rawMinQualityScore" = rawMinQualityScore,
             "trimmedStartPos" = trimmedStartPos,
             "trimmedFinishPos" = trimmedFinishPos,
             "trimmedSeqLength" = trimmedSeqLength,
             "trimmedMeanQualityScore" = trimmedMeanQualityScore,
             "trimmedMinQualityScore" = trimmedMinQualityScore,
             "remainingRatio" = remainingRatio))
}

M2inside_calculate_trimming <- function(qualityPhredScores, qualityBaseScores,
                                        M2CutoffQualityScore, M2SlidingWindowSize) {
    rawSeqLength <- length(qualityBaseScores)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    qualityPbCutoff <- 10** (M2CutoffQualityScore / (-10.0))
    remainingIndex <- c()
    if (M2SlidingWindowSize > 20 || M2SlidingWindowSize < 0 ||
        M2SlidingWindowSize%%1!=0 ||
        M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
        M2CutoffQualityScore%%1!=0) {
        trimmedStartPos = NULL
        trimmedFinishPos = NULL
    } else {
        for (i in 1:(rawSeqLength-M2SlidingWindowSize+1)) {
            meanSLidingWindow <-
                mean(qualityBaseScores[i:(i+M2SlidingWindowSize-1)])
            if (meanSLidingWindow < qualityPbCutoff) {
                remainingIndex <- c(remainingIndex, i)
                # or ==> i + floor(M2SlidingWindowSize/3)
            }
        }
        trimmedStartPos = remainingIndex[1]
        trimmedFinishPos = remainingIndex[length(remainingIndex)]
        if (is.null(trimmedStartPos) || is.null(trimmedFinishPos)) {
            trimmedStartPos <- 1
            trimmedFinishPos <- 2
        }
        trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
        trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
        trimmedMinQualityScore <- min(trimmedQualityPhredScore)
        trimmedSeqLength = trimmedFinishPos - trimmedStartPos
        remainingRatio = trimmedSeqLength / rawSeqLength
    }
    return(list("rawSeqLength" = rawSeqLength,
                "rawMeanQualityScore" = rawMeanQualityScore,
                "rawMinQualityScore" = rawMinQualityScore,
                "trimmedStartPos" = trimmedStartPos,
                "trimmedFinishPos" = trimmedFinishPos,
                "trimmedSeqLength" = trimmedSeqLength,
                "trimmedMeanQualityScore" = trimmedMeanQualityScore,
                "trimmedMinQualityScore" = trimmedMinQualityScore,
                "remainingRatio" = remainingRatio))
}

### ============================================================================
### Calculating SangerContig
### ============================================================================
calculateContigSeq <- function(forwardReadsList, forwardReadList,
                               refAminoAcidSeq, minFractionCall,
                               maxFractionLost, geneticCode,
                               acceptStopCodons, readingFrame,
                               processorsNum) {
    ### ------------------------------------------------------------------------
    ### forward & reverse character reads list string creation
    ### ------------------------------------------------------------------------
    fRDNAStringSet <- sapply(forwardReadsList, function(forwardRead) {
        trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
        trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
        primaryDNA <- as.character(forwardRead@primarySeq)
        substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
    })
    rRDNAStringSet <- sapply(forwardReadList, function(reverseRead) {
        trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
        trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
        primaryDNA <- as.character(reverseRead@primarySeq)
        substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
    })

    ### --------------------------------------------------------------------
    ### DNAStringSet storing forward & reverse reads ! (Origin)
    # ### --------------------------------------------------------------------
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
    ### ----------------------------------------------------------------
    ### Remove reads with stop codons
    ### ----------------------------------------------------------------
    if (!acceptStopCodons) {
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
    return(list("consensusGapfree" = consensusGapfree,
                "diffsDf"          = diffsDf,
                "aln2"             = aln2,
                "dist"             = dist,
                "dend"             = dend,
                "indels"           = indels,
                "stopsDf"          = stopsDf,
                "spDf"             = spDf))
}

### ============================================================================
### Aligning SangerContigs (SangerAlignment)
### ============================================================================
alignContigs <- function(SangerContigList, geneticCode, refAminoAcidSeq,
                         minFractionCallSA, maxFractionLostSA, processorsNum) {
    ### ------------------------------------------------------------------------
    ### Creating SangerContigList DNAStringSet
    ### ------------------------------------------------------------------------
    SangerContigDNAList <-
        sapply(SangerContigList, function(SangerContig) {
            as.character(SangerContig@contigSeq)
        })
    SangerContigDNASet <- DNAStringSet(SangerContigDNAList)
    ### ------------------------------------------------------------------------
    ### Aligning consensus reads
    ### ------------------------------------------------------------------------
    if(length(SangerContigDNASet) > 1) {
        message("Aligning consensus reads ... ")
        if(refAminoAcidSeq != ""){
            aln = AlignTranslation(SangerContigDNASet,
                                   geneticCode = geneticCode,
                                   processors = processorsNum,
                                   verbose = FALSE)
        }else{
            aln = AlignSeqs(SangerContigDNASet, processors = processorsNum,
                            verbose = FALSE)
        }
        # Making a rough NJ tree. Labels are rows in the summary df
        neat.labels = match(names(aln),
                            as.character(names(SangerContigDNASet)))
        aln2 = aln
        names(aln2) = neat.labels

        aln.bin = as.DNAbin(aln2)

        aln.dist = dist.dna(aln.bin, pairwise.deletion = TRUE)
        # Making a rough NJ tree. Labels are rows in the summary df
        #    (If tree cannot be created ==> NULL)
        aln.tree = NULL
        try({
            aln.tree = bionjs(aln.dist)
            aln.tree$tip.label <- names(aln)
            # deal with -ve branches
            # This is not necessarily accurate, but it is good enough to
            # judge seuqences using the tree
            aln.tree$edge.length[which(aln.tree$edge.length<0)] =
                abs(aln.tree$edge.length[which(aln.tree$edge.length<0)])
        },
        silent = TRUE
        )

        # Get consensus read and add to alignment result
        consensus = ConsensusSequence(aln,
                                      minInformation = minFractionCallSA,
                                      includeTerminalGaps = TRUE,
                                      ignoreNonBases = TRUE,
                                      threshold = maxFractionLostSA,
                                      noConsensusChar = "-",
                                      ambiguity = TRUE)[[1]]
    } else {
        aln = NULL
        aln.tree = NULL
    }
    return(list("consensus" = consensus,
                "aln"       = aln,
                "aln.tree"  = aln.tree))
}


### ============================================================================
### chromatogram row number counting
### ============================================================================
chromatogramRowNum <- function(width, rawLength, trimmedLength, showTrimmed) {
    if (showTrimmed) {
        numplots = ceiling(rawLength / width)
    } else {
        numplots = ceiling(trimmedLength / width)
    }
}
