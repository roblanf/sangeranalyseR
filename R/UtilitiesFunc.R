### ============================================================================
### Global helper functions
### ============================================================================
getProcessors <- function(processors) {
    sysinf <- Sys.info()
    if (!is.null(sysinf)){
        os <- sysinf["sysname"]
        if (os == "Darwin")
            os <- "macos"
    } else {
        ## mystery machine
        os <- .Platform$OS.type
        if (grepl("^darwin", R.version$os))
            os <- "macos"
        if (grepl("linux-gnu", R.version$os))
            os <- "linux"
    }
    os <- tolower(os)
    if (os != "linux" && os != "osx") {
        ### --------------------------------------------------------------------
        ### Operating system is Window
        ### --------------------------------------------------------------------
        return(1)
    } else {
        ### --------------------------------------------------------------------
        ### Operating system is macOS or Linux
        ### --------------------------------------------------------------------
        if(is.null(processors)){
            processors = detectCores(all.tests = FALSE, logical = FALSE)
        }
        return(processors)
    }
}
suppressPlotlyMessage <- function(p) {
    suppressMessages(plotly_build(p))
}

# <---------------------------------------------------------------------------->

### ============================================================================
### SangerAlignment related helper functions
### ============================================================================
### ----------------------------------------------------------------------------
### Aligning SangerContigs (SangerAlignment)
### ----------------------------------------------------------------------------
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
            message('Before building!!')
            aln = AlignSeqs(SangerContigDNASet, processors = processorsNum,
                            verbose = FALSE)
            message('After building!!')
        }
        # Making a rough NJ tree. Labels are rows in the summary df

        if (length(aln) > 2) {
            neat.labels = match(names(aln),
                                as.character(names(SangerContigDNASet)))
            aln2 = aln
            names(aln2) = neat.labels

            aln.bin = as.DNAbin(aln2)

            aln.dist = dist.dna(aln.bin, pairwise.deletion = TRUE)
            # Making a rough NJ tree. Labels are rows in the summary df
            #    (If tree cannot be created ==> NULL)
            aln.tree = read.tree(text="();")
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
        } else {
            aln.tree = read.tree(text="();")
        }

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

# <---------------------------------------------------------------------------->

### ============================================================================
### SangerContig related helper functions
### ============================================================================
getIndelDf <- function(indelList){
    r = lapply(indelList, indelRow)
    indelDf = data.frame(matrix(unlist(r), byrow = TRUE, nrow = length(indelList)))
    indelDf = cbind(names(indelList), indelDf)
    names(indelDf) = c('read', 'insertions', 'deletions', 'distance')
    return(indelDf)
}
indelRow <- function(row){
    nIns = length(row$insertions)
    nDel = length(row$deletions)
    dist = row$distance
    return(c(nIns, nDel, dist))
}
nPairwiseDiffs <- function(pattern, subject){
    # pairwise differences assuming pattern and subject are aligned
    comp = compareStrings(pattern, subject)
    qs = str_count(comp, '\\?')
    ps = str_count(comp, '\\+')
    return(c(qs, ps))
}
countCoincidentSp <- function(aln, processorsNum){
    # make a data frame of columns in the alignment that have
    # more than one secondary peak
    is = 1:aln@ranges@width[1]
    r = mclapply(is, oneAmbiguousColumn, aln=aln, mc.cores = processorsNum)
    r = Filter(Negate(is.null), r)

    if(length(r)>0){
        r = as.data.frame(matrix(unlist(r), nrow=length(r), byrow=TRUE))
        names(r) = c('column.number', 'ambiguities', 'column')
        return(r)
    }else{
        return(NULL)
    }
}
oneAmbiguousColumn <- function(i, aln){
    ambiguous = names(IUPAC_CODE_MAP)[5:length(names(IUPAC_CODE_MAP))]
    col = as.character(subseq(aln, i, i))
    str = paste(col, sep="", collapse="")
    amb = sum(col %in% ambiguous)
    if(amb>1){
        return(c(i, amb, str))
    }
}
### ----------------------------------------------------------------------------
### Calculating SangerContig
### ----------------------------------------------------------------------------
calculateContigSeq <- function(inputSource, forwardReadList, reverseReadList,
                               refAminoAcidSeq, minFractionCall,
                               maxFractionLost, geneticCode,
                               acceptStopCodons, readingFrame,
                               processorsNum) {
    ### ------------------------------------------------------------------------
    ### forward & reverse character reads list string creation
    ### ------------------------------------------------------------------------
    fRDNAStringSet <- sapply(forwardReadList, function(forwardRead) {
        primaryDNA <- as.character(forwardRead@primarySeq)
        if (inputSource == "ABIF") {
            trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
            primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
        }
        return(primaryDNA)
    })
    rRDNAStringSet <- sapply(reverseReadList, function(reverseRead) {
        primaryDNA <- as.character(reverseRead@primarySeq)
        if (inputSource == "ABIF") {
            trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
            primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
        }
        return(primaryDNA)
    })
    ### ------------------------------------------------------------------------
    ### DNAStringSet storing forward & reverse reads ! (Origin)
    ### ------------------------------------------------------------------------
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

    ### ------------------------------------------------------------------------
    ### Amino acid reference sequence CorrectFrameshifts correction
    ### ------------------------------------------------------------------------
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
    ### ------------------------------------------------------------------------
    ### Reads with stop codons elimination
    ### ------------------------------------------------------------------------
    ### ------------------------------------------------------------------------
    ### Remove reads with stop codons
    ### ------------------------------------------------------------------------
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

    ### ------------------------------------------------------------------------
    ### Start aligning reads
    ### ------------------------------------------------------------------------
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
    spDf = countCoincidentSp(aln, processorsNum = processorsNum)
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
### ----------------------------------------------------------------------------
### MakeBaseCalls related function
### ----------------------------------------------------------------------------
MakeBaseCallsInside <- function(traceMatrix, peakPosMatrixRaw,
                                qualityPhredScoresRaw,
                                signalRatioCutoff, readFeature) {
    message("          * Making basecall !!")
    #get peaks for each base
    Apeaks <- getpeaks(traceMatrix[,1])
    Cpeaks <- getpeaks(traceMatrix[,2])
    Gpeaks <- getpeaks(traceMatrix[,3])
    Tpeaks <- getpeaks(traceMatrix[,4])

    #get window around primary basecall peaks
    primarypeaks <- peakPosMatrixRaw[,1]
    diffs <- diff(c(0,primarypeaks))
    starts <- primarypeaks - 0.5*diffs
    stops <- c(primarypeaks[1:(length(primarypeaks)-1)] +
                   0.5*diffs[2:length(diffs)],
               primarypeaks[length(diffs)] + 0.5*diffs[length(diffs)]
    )
    #hack for last peak. Just uses distance preceding peak
    #as distance after peak

    #Now get max peak value for each channel in each peak window.
    #If no peak return 0
    primary <- NULL
    secondary <- NULL
    tempPosMatrix <- matrix(nrow=length(starts), ncol=4)
    tempAmpMatrix <- matrix(nrow=length(starts), ncol=4)
    indexBaseCall <- c()
    for(i in 1:length(starts)) {
        Apeak <- peakvalues(Apeaks, starts[i], stops[i])
        Cpeak <- peakvalues(Cpeaks, starts[i], stops[i])
        Gpeak <- peakvalues(Gpeaks, starts[i], stops[i])
        Tpeak <- peakvalues(Tpeaks, starts[i], stops[i])
        if(is.na(Apeak[2]) &
           is.na(Cpeak[2]) &
           is.na(Gpeak[2]) &
           is.na(Tpeak[2])) {
            next #rare case where no peak found
        }
        ### --------------------------------------------------------------------
        ### My modification here: Tracking BaseCall index
        ###     Add qualtiy score when making basecall
        ### --------------------------------------------------------------------
        indexBaseCall <- c(indexBaseCall, i)
        signals <- c(Apeak[1], Cpeak[1], Gpeak[1], Tpeak[1])

        tempAmpMatrix[i,] <- signals
        # print(tempAmpMatrix[i,])

        positions <- c(Apeak[2], Cpeak[2], Gpeak[2], Tpeak[2])

        tempPosMatrix[i,] <- positions
        # print(tempPosMatrix[i,])


        signalratios <- signals/max(signals, na.rm=TRUE)
        Bases <- c("A", "C", "G", "T")
        Bases[signalratios < signalRatioCutoff] <- NA
        #sort by decreasing signal strength
        Bases <- Bases[order(signals, decreasing=TRUE)]
        positions <- positions[order(signals, decreasing=TRUE)]
        if(length(Bases[!is.na(Bases)]) == 4
           | length(Bases[!is.na(Bases)]) == 0) {
            primary <- c(primary, "N")
            secondary <- c(secondary, "N")
        } else if (length(Bases[!is.na(Bases)]) > 1) {
            primary <- c(primary, Bases[1])
            Bases2 <- Bases[2:4]
            sortedBase2<- sort(Bases2[!is.na(Bases2)])
            # message("Collapse: ", paste(sortedBase2, collapse=""))
            dicValue <- paste(sortedBase2, collapse="")
            secondaryLetter <- names(IUPAC_CODE_MAP[IUPAC_CODE_MAP == dicValue])
            # secondaryLetter <-
            #     mergeIUPACLetters(paste(sortedBase2, collapse=""))
            # message("secondaryLetter: ", secondaryLetter)
            secondary <- c(secondary, secondaryLetter)
        }
        else {
            primary <- c(primary, Bases[1])
            secondary <- c(secondary, Bases[1])
        }
    }
    if (readFeature == "Forward Read") {
        qualityPhredScores <- qualityPhredScoresRaw[indexBaseCall]
        primarySeq <- DNAString(paste(primary, collapse=""))
        secondarySeq <- DNAString(paste(secondary, collapse=""))
    } else if (readFeature == "Reverse Read") {
        qualityPhredScores <- rev(qualityPhredScoresRaw[indexBaseCall])
        primarySeq <- reverseComplement(DNAString(paste(primary, collapse="")))
        secondarySeq <- reverseComplement(DNAString(paste(secondary, collapse="")))
    }
    peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    message("          * Updating slots in 'SangerRead' instance !!")
    return(list("qualityPhredScores" = qualityPhredScores,
                "peakPosMatrix" = peakPosMatrix,
                "peakAmpMatrix" = peakAmpMatrix,
                "primarySeq" = primarySeq,
                "secondarySeq" = secondarySeq))
}
getpeaks <- function(trace) {
    r <- rle(trace)
    indexes <- which(rep(diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2,
                         times = r$lengths))
    cbind(indexes, trace[indexes])
}
peakvalues <- function(x, pstart, pstop) {
    region <- x[x[,1] > pstart & x[,1] < pstop, ,drop=FALSE]
    if (length(region[,1]) == 0) return(c(0, NA))
    else return(c(max(region[,2], na.rm=TRUE), region[which.max(region[,2]),1]))
}

### ----------------------------------------------------------------------------
### MakeBasecall secondary peak finding helper function
### ----------------------------------------------------------------------------
IUPAC_CODE_MAP <- c(
    A="A",
    C="C",
    G="G",
    T="T",
    M="AC",
    R="AG",
    W="AT",
    S="CG",
    Y="CT",
    K="GT",
    V="ACG",
    H="ACT",
    D="AGT",
    B="CGT",
    N="ACGT"
)

### ----------------------------------------------------------------------------
### chromatogram related function
### ----------------------------------------------------------------------------
chromatogramRowNum <- function(width, rawLength, trimmedLength, showTrimmed) {
    if (showTrimmed) {
        numplots = ceiling(rawLength / width)
    } else {
        numplots = ceiling(trimmedLength / width)
    }
}

# <---------------------------------------------------------------------------->

### ============================================================================
### SangerRead related helper functions
### ============================================================================
SetCharStyleList <- function(AASeqDF, selectChar, colorCode) {
    stopIndex <- AASeqDF %>% `==` (selectChar) %>% which()
    stopExcelIndex <- int2col(stopIndex)
    stopExcelIndexName <- paste0(stopExcelIndex, "1")
    styleList <-
        as.list(rep(paste('background-color:', colorCode,
                          "; font-weight: bold;"), length(stopExcelIndex)))
    if (length(stopIndex) != 0) {
        names(styleList) <- stopExcelIndexName
    }
    return(styleList)
}
SetAllStyleList <- function(AASeqDF, colorCode) {
    Index <- strtoi(names(AASeqDF))
    ExcelIndex <- int2col(Index)
    ExcelIndexName <- paste0(ExcelIndex, "1")
    styleList <-
        as.list(rep(paste('background-color:', colorCode,
                          "; font-weight: bold;"), length(ExcelIndex)))
    names(styleList) <- ExcelIndexName
    return(styleList)
}
countStopSodons <- function(sequence,
                            readingFrame = 1, geneticCode = GENETIC_CODE){
    l = length(sequence) + 1 - readingFrame
    if(l < 3){
        sprintf("Cannot calculate stop codons on sequence of length %d",
                " in reading frame %d",length(sequence), readingFrame)
        # return(NULL)
        error <- paste("\nCannot calculate stop codons on sequence of length ",
                       length(sequence), " in reading frame ", readingFrame,
                       ".\n", sep = "")
        stop(error)
    }
    # this comes almost straight from the BioStrings manual
    tri = trinucleotideFrequency(sequence[readingFrame:length(sequence)],step=3)
    names(tri) <- geneticCode[names(tri)]
    freqs = sapply(split(tri, names(tri)), sum)
    stops = freqs["*"]
    return(as.numeric(stops))
}
calculateAASeq <- function(primarySeq, trimmedStartPos,
                           trimmedFinishPos, geneticCode) {
    DNASeqshift0 <- DNAString(substr(as.character(primarySeq),
                                     1+trimmedStartPos, trimmedFinishPos))
    primaryAASeqS1 <-
        suppressWarnings(Biostrings::translate(DNASeqshift0,
                                   genetic.code = geneticCode,
                                   no.init.codon=TRUE,
                                   if.fuzzy.codon="solve"))

    DNASeqshift1 <- DNAString(substr(as.character(primarySeq),
                                     2+trimmedStartPos, trimmedFinishPos))
    primaryAASeqS2 <-
        suppressWarnings(Biostrings::translate(DNASeqshift1,
                                   genetic.code = geneticCode,
                                   no.init.codon=TRUE,
                                   if.fuzzy.codon="solve"))

    DNASeqshift2 <- DNAString(substr(as.character(primarySeq),
                                     3+trimmedStartPos, trimmedFinishPos))
    primaryAASeqS3 <-
        suppressWarnings(Biostrings::translate(DNASeqshift2,
                                   genetic.code = geneticCode,
                                   no.init.codon=TRUE,
                                   if.fuzzy.codon="solve"))
    return(list("primaryAASeqS1" = primaryAASeqS1,
                "primaryAASeqS2" = primaryAASeqS2,
                "primaryAASeqS3" = primaryAASeqS3))
}
### ----------------------------------------------------------------------------
### Quality trimming related parameter
### ----------------------------------------------------------------------------
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
    trimmedQualityPhredScore <-
        qualityPhredScores[(trimmedStartPos+1):trimmedFinishPos]
    trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
    trimmedMinQualityScore <- min(trimmedQualityPhredScore)
    remainingRatio = trimmedSeqLength / rawSeqLength

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
M2inside_calculate_trimming <-function(qualityPhredScores,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize) {
    rawSeqLength <- length(qualityPhredScores)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    if (M2SlidingWindowSize > 40 || M2SlidingWindowSize < 0 ||
        M2SlidingWindowSize%%1!=0 ||
        M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
        M2CutoffQualityScore%%1!=0) {
        trimmedStartPos = NULL
        trimmedFinishPos = NULL
    } else {
        ### ------------------------------------------------------------------------
        ### Find the trimming start point
        ###  First window that the average score is bigger than threshold score.
        ###  (Whole window will be kept)
        ### ------------------------------------------------------------------------
        totalThresholdScore <- M2CutoffQualityScore * M2SlidingWindowSize
        trimmedStartPos = 1
        trimmedFinishPos = 2
        for (i in 1:(rawSeqLength-M2SlidingWindowSize+1)) {
            totalScore <-
                sum(qualityPhredScores[i:(i+M2SlidingWindowSize-1)])
            if (totalScore > totalThresholdScore) {
                trimmedStartPos = i
                break
            }
        }
        qualityPhredScoresRev <- rev(qualityPhredScores)
        for (i in 1:(rawSeqLength-M2SlidingWindowSize+1)) {
            totalScore <-
                sum(qualityPhredScoresRev[i:(i+M2SlidingWindowSize-1)])
            if (totalScore > totalThresholdScore) {
                trimmedFinishPos = i
                break
            }
        }
        trimmedFinishPos <- length(qualityPhredScoresRev) - trimmedFinishPos + 1
        # for (i in (trimmedStartPos+M2SlidingWindowSize-1):(rawSeqLength-M2SlidingWindowSize+1)) {
        #     totalScore <-
        #         sum(qualityPhredScores[i:(i+M2SlidingWindowSize-1)])
        #     if (totalScore < totalThresholdScore) {
        #         # Keep all base pairs in the previous window.
        #         trimmedFinishPos = i + M2SlidingWindowSize - 2
        #         break
        #     }
        # }
        if (trimmedStartPos == (rawSeqLength-M2SlidingWindowSize+1) ||
            trimmedStartPos == trimmedFinishPos) {
            trimmedStartPos = 1
            trimmedFinishPos = 2
        }
        trimmedSeqLength = trimmedFinishPos - trimmedStartPos
        trimmedQualityPhredScore <-
            qualityPhredScores[(trimmedStartPos+1):trimmedFinishPos]
        trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
        trimmedMinQualityScore <- min(trimmedQualityPhredScore)
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

### ----------------------------------------------------------------------------
### Quality score base pair plot functions
### ----------------------------------------------------------------------------
QualityBasePlotly <- function(trimmedStartPos, trimmedFinishPos,
                              readLen, qualityPlotDf, x,  y) {
    p <- suppressPlotlyMessage(
        plot_ly(data=qualityPlotDf,
                x=~Index) %>%
            add_markers(y=~Score,
                        text = ~paste("BP Index : ",
                                      Index,
                                      '<sup>th</sup><br>Phred Quality Score :',
                                      Score),
                        name = 'Quality Each BP') %>%
            add_trace(x=seq(trimmedStartPos,
                            trimmedFinishPos,
                            len=trimmedFinishPos-trimmedStartPos+1),
                      y=rep(70, trimmedFinishPos-trimmedStartPos+1),
                      mode="lines", hoverinfo="text",
                      text=paste("Trimmed Reads BP length:",
                                 trimmedFinishPos-trimmedStartPos+1,
                                 "BPs <br>",
                                 "Trimmed Reads BP ratio:",
                                 round((trimmedFinishPos - trimmedStartPos+1)/
                                           readLen * 100,
                                       digits=2),
                                 "%"),
                      line = list(width = 12),
                      name = 'Trimmed Read') %>%
            add_trace(x=seq(0,readLen,len=readLen),
                      y=rep(80, readLen), mode="lines", hoverinfo="text",
                      text=paste("Whole Reads BP length:",
                                 readLen,
                                 "BPs <br>",
                                 "Trimmed Reads BP ratio: 100 %"),
                      line = list(width = 12),
                      name = 'Whole Read') %>%
            layout(xaxis = x, yaxis = y,
                   shapes = list(vline(trimmedStartPos),
                                 vline(trimmedFinishPos)),
                   legend = list(orientation = 'h',
                                 xanchor = "center",
                                 x = 0.5, y = 1.1)) %>%
            add_annotations(
                text = "Trimming Strat <br> BP Index",
                x = trimmedStartPos + 40,
                y = 15,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Trimming End <br> BP Index",
                x = trimmedFinishPos - 40,
                y = 15,
                showarrow=FALSE
            ))
    return(p)
}

vline <- function(x = 0, color = "red") {
    list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = x,
        x1 = x,
        line = list(color = color)
    )
}

