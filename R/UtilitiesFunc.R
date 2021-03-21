### ============================================================================
### Global helper functions
### ============================================================================
getProcessors <- function(processors = NULL) {
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
        vapply(SangerContigList, function(SangerContig) {
            as.character(SangerContig@contigSeq)
        }, FUN.VALUE = character(1))
    SangerContigDNASet <- DNAStringSet(SangerContigDNAList)
    ### ------------------------------------------------------------------------
    ### Aligning consensus reads
    ### ------------------------------------------------------------------------
    if(length(SangerContigDNASet) > 1) {
        log_info("Aligning consensus reads ... ")
        if(refAminoAcidSeq != ""){
            aln = AlignTranslation(SangerContigDNASet,
                                   geneticCode = geneticCode,
                                   processors = processorsNum,
                                   verbose = FALSE)
        }else{
            log_info('Before building!!')
            aln = AlignSeqs(SangerContigDNASet, processors = processorsNum,
                            verbose = FALSE)
            log_info('After building!!')
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
            tryCatch({
                aln.tree = bionjs(aln.dist)
                aln.tree$tip.label <- names(aln)
                # deal with -ve branches
                # This is not necessarily accurate, but it is good enough to
                # judge seuqences using the tree
                aln.tree$edge.length[which(aln.tree$edge.length<0)] =
                    abs(aln.tree$edge.length[which(aln.tree$edge.length<0)])
            }, warning = function(warning_condition) {
                log_info("The number of contigs is less than 3 or quality of reads ",
                        "are too low. 'Contigs Tree' cannot be created.")
            }, error = function(error_condition) {
                log_info("The number of contigs is less than 3 or quality of reads ",
                        "are too low. 'Contigs Tree' cannot be created.")
                aln.tree = read.tree(text="();")
            })
        } else {
            log_info("The number of contigs is less than 3 or quality of reads ",
                    "are too low. 'Contigs Tree' cannot be created.")
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
        consensus = NULL
        aln = NULL
        phyloT <- as.phylo(rtree(n = 2))
        phyloT$edge <- matrix(c(0,0,0,0), 2, 2)
        phyloT$tip.label <- NULL
        phyloT$edge.length <- NULL
        phyloT$Nnode <- NULL
        aln.tree = phyloT
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
    is = seq_len(aln@ranges@width[1])
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
                               processorsNum = NULL) {
    ### ------------------------------------------------------------------------
    ### forward & reverse character reads list string creation
    ### ------------------------------------------------------------------------
    fRDNAStringSet <- lapply(forwardReadList, function(forwardRead) {
        primaryDNA <- as.character(forwardRead@primarySeq)
        if (inputSource == "ABIF") {
            trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
            primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
        }
        return(primaryDNA)
    })
    rRDNAStringSet <- lapply(reverseReadList, function(reverseRead) {
        DNALen <- length(reverseRead@primarySeq)
        primaryDNA <- as.character(reverseComplement(reverseRead@primarySeq))
        if (inputSource == "ABIF") {
            trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
            ## Trimming on the original reads (not reverseComplement).
            primaryDNA <- substr(primaryDNA, DNALen - trimmedFinishPos + 1,
                                 DNALen - trimmedStartPos)
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
    # Read number in each contig can be 1!
    # if(length(frReadSet) < 2) {
    #     error <- paste("\n'Valid abif files should be more than 2.\n",
    #                    sep = "")
    #     log_error(error)
    # }
    processorsNum <- getProcessors(processorsNum)

    ### ------------------------------------------------------------------------
    ### Amino acid reference sequence CorrectFrameshifts correction
    ### ------------------------------------------------------------------------
    if (refAminoAcidSeq != "") {
        log_info("Correcting frameshifts in reads using amino acid",
                "reference sequence")
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
        log_error(error)
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
        log_info(old_length - length(frReadSet),
                "reads with stop codons removed")
    }

    if(length(frReadSet) < 2) {
        error <- paste("\n'After removing reads with stop codons, ",
                       "forward and reverse reads should be more than 2.\n",
                       sep = "")
        log_error(error)
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
    names(aln) = paste(seq_len(length(aln)), "Read",
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
                                signalRatioCutoff, readFeature, printLevel) {
    if (printLevel == "SangerRead") {
        log_info("          * Making basecall !!")
    }
    #get peaks for each base
    Apeaks <- getpeaks(traceMatrix[,1])
    Cpeaks <- getpeaks(traceMatrix[,2])
    Gpeaks <- getpeaks(traceMatrix[,3])
    Tpeaks <- getpeaks(traceMatrix[,4])
    
    #get window around primary basecall peaks
    primarypeaks <- peakPosMatrixRaw[,1]
    diffs <- diff(c(0,primarypeaks))
    starts <- primarypeaks - 0.5*diffs
    stops <- c(primarypeaks[seq_len((length(primarypeaks)-1))] +
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
    for(i in seq_len(length(starts))) {
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
        positions <- c(Apeak[2], Cpeak[2], Gpeak[2], Tpeak[2])
        tempPosMatrix[i,] <- positions
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
            dicValue <- paste(sortedBase2, collapse="")
            secondaryLetter <- names(IUPAC_CODE_MAP[IUPAC_CODE_MAP == dicValue])
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
        qualityPhredScores <- qualityPhredScoresRaw[indexBaseCall]
        primarySeq <- DNAString(paste(primary, collapse=""))
        secondarySeq <- DNAString(paste(secondary, collapse=""))
    }
    peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    if (printLevel == "SangerRead") {
        log_info("          * Updating slots in 'SangerRead' instance !!")
    }
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
        log_error(error)
    }
    # this comes almost straight from the BioStrings manual
    tri = trinucleotideFrequency(sequence[readingFrame:length(sequence)],step=3)
    names(tri) <- geneticCode[names(tri)]
    freqs = lapply(split(tri, names(tri)), sum)
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
        for (i in seq_len((rawSeqLength-M2SlidingWindowSize+1))) {
            totalScore <-
                sum(qualityPhredScores[i:(i+M2SlidingWindowSize-1)])
            if (totalScore > totalThresholdScore) {
                trimmedStartPos = i
                break
            }
        }
        qualityPhredScoresRev <- rev(qualityPhredScores)
        for (i in seq_len((rawSeqLength-M2SlidingWindowSize+1))) {
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
                text = "Trimming Start <br> BP Index",
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

SangerReadInnerTrimming <- function(SangerReadInst, inputSource) {
    primaryDNA <- as.character(SangerReadInst@primarySeq)
    if (inputSource == "ABIF") {
        trimmedStartPos <- 
            SangerReadInst@QualityReport@trimmedStartPos
        trimmedFinishPos <- 
            SangerReadInst@QualityReport@trimmedFinishPos
        primaryDNA <- 
            substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
    }
    return(primaryDNA)
}

#' @export
chromatogram_overwrite <- function(obj, trim5=0, trim3=0, 
                                   showcalls=c("primary", "secondary", "both", "none"), 
                                   width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
                                   filename=NULL, showtrim=FALSE, showhets=TRUE, colors="default") {
    if (colors == "default") {
        A_color = "green"
        T_color = "blue"
        C_color = "black"
        G_color = "red"
        unknown_color = "purple"
    } else if (colors == "cb_friendly") {
        A_color = rgb(0, 0, 0, max = 255)
        T_color = rgb(199, 199, 199, max = 255)
        C_color = rgb(0, 114, 178, max = 255)
        G_color = rgb(213, 94, 0, max = 255)
        unknown_color = rgb(204, 121, 167, max = 255)
    } else {
        A_color = colors[1]
        T_color = colors[2]
        C_color = colors[3]
        G_color = colors[4]
        unknown_color = colors[5]
    }
    originalpar <- par(no.readonly=TRUE)
    showcalls <- showcalls[1]
    traces <- obj@traceMatrix
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    basecalls1 <- basecalls1[1:length(aveposition)] 
    basecalls2 <- basecalls2[1:length(aveposition)] 
    if(showtrim == FALSE) {
        if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
        else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
        if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
        else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
        aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
    }
    indexes <- 1:length(basecalls1)
    trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all 
    #false if not trimmed
    if (!is.null(trim3)) {
        traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, 
                                nrow(traces))), ]
    }
    if (!is.null(trim5)) {
        offset <- max(c(1, aveposition[1] - 10))
        traces <- traces[offset:nrow(traces),]
        aveposition <- aveposition - (offset-1)
    }
    maxsignal <- apply(traces, 1, max)
    ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))           
    p <- c(0, aveposition, nrow(traces))
    midp <- diff(p)/2
    starts <- aveposition - midp[1:(length(midp)-1)]
    starthets <- starts
    starthets[basecalls1 == basecalls2] <- NA
    ends <- aveposition + midp[2:(length(midp))]
    endhets <- ends
    endhets[basecalls1 == basecalls2] <- NA
    starttrims <- starts
    starttrims[!trimmed] <- NA
    endtrims <- ends
    endtrims[!trimmed] <- NA
    
    colortranslate <- c(A=A_color, C=C_color, G=G_color, T=T_color)
    colorvector1 <- unname(colortranslate[basecalls1])
    colorvector1[is.na(colorvector1)] <- unknown_color
    colorvector2 <- unname(colortranslate[basecalls2])
    colorvector2[is.na(colorvector2)] <- unknown_color
    
    valuesperbase <- nrow(traces)/length(basecalls1)
    tracewidth <- width*valuesperbase
    breaks <- seq(1,nrow(traces), by=tracewidth)
    numplots <- length(breaks)
    if(!is.null(filename)) {
        if (! dir.exists(dirname(filename))) {
            dir.create(dirname(filename))
        }
        pdf(filename, width=8.5, height=height*numplots) 
    }
    par(mar=c(2,2,2,1), mfrow=c(numplots, 1))
    basecallwarning1 = 0
    basecallwarning2 = 0
    j = 1
    
    for(i in breaks) {
        range <- aveposition >= i & aveposition < (i+tracewidth)
        starthet <- starthets[range] - tracewidth*(j-1)
        starthet[starthet < 0] <- 0
        endhet <- endhets[range] - tracewidth*(j-1)
        endhet[endhet > tracewidth] <- tracewidth
        lab1 <- basecalls1[range]
        lab2 <- basecalls2[range]
        pos <- aveposition[range] - tracewidth*(j-1)
        colors1 <- colorvector1[range]
        colors2 <- colorvector2[range]
        starttrim <- starttrims[range] - tracewidth*(j-1)
        endtrim <- endtrims[range] - tracewidth*(j-1)
        plotrange <- i:min(i+tracewidth, nrow(traces))
        plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n", 
             bty="n", xlab="", yaxt="n", , xlim=c(1,tracewidth))
        if (showhets==TRUE) {
            rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')
        }
        if (showtrim==TRUE) {
            rect(starttrim, 0, endtrim, ylims[2], col='red', border='transparent', 
                 density=15)
        }
        lines(traces[plotrange,1], col=A_color)
        lines(traces[plotrange,2], col=T_color)
        lines(traces[plotrange,3], col=C_color)
        lines(traces[plotrange,4], col=G_color)
        mtext(as.character(which(range)[1]), side=2, line=0, cex=cex.mtext)
        
        for(k in 1:length(lab1)) {
            if (showcalls=="primary" | showcalls=="both") {
                if (is.na(basecalls1[1]) & basecallwarning1==0) {
                    warning("Primary basecalls missing")
                    basecallwarning1 = 1
                } 
                else if (length(lab1) > 0) {   
                    axis(side=3, at=pos[k], labels=lab1[k], col.axis=colors1[k], 
                         family="mono", cex=cex.base, line=ifelse(showcalls=="both", 0, 
                                                                  -1), tick=FALSE)
                }
            }
            if (showcalls=="secondary" | showcalls=="both") {
                if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
                    warning("Secondary basecalls missing")
                    basecallwarning2 = 1
                } 
                else if (length(lab2) > 0) { 
                    axis(side=3, at=pos[k], labels=lab2[k], col.axis=colors2[k], 
                         family="mono", cex=cex.base, line=-1, tick=FALSE)
                }
            }
        }
        j = j + 1
    }
    if(!is.null(filename)) {
        dev.off()
        cat(paste("Chromatogram saved to", filename, 
                  "in the current working directory"))
    }
    else par(originalpar)
}
