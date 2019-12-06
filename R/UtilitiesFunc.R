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


getIndelDf <- function(indelList){
    r = lapply(indelList, indelRow)
    indelDf = data.frame(matrix(unlist(r), byrow = T, nrow = length(indelList)))
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
        r = as.data.frame(matrix(unlist(r), nrow=length(r), byrow=T))
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


countStopSodons <- function(sequence, readingFrame = 1, geneticCode = GENETIC_CODE){
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
    tri = trinucleotideFrequency(sequence[readingFrame:length(sequence)], step=3)
    names(tri) <- geneticCode[names(tri)]
    freqs = sapply(split(tri, names(tri)), sum)
    stops = freqs["*"]
    return(as.numeric(stops))
}

suppressPlotlyMessage <- function(p) {
    suppressMessages(plotly_build(p))
}











### ============================================================================
### MakeBaseCalls related function
### ============================================================================
MakeBaseCallsInside <- function(traceMatrix, peakPosMatrixRaw,
                                qualityPhredScoresRaw,
                                signalRatioCutoff) {
    message("     * Making basecall !!")
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
        ### ----------------------------------------------------------
        ### My modification here: Tracking BaseCall index
        ###     Add qualtiy score when making basecall
        ### ----------------------------------------------------------
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
        }
        else if(length(Bases[!is.na(Bases)]) > 1) {
            primary <- c(primary, Bases[1])
            Bases2 <- Bases[2:4]
            secondary <- c(secondary,
                           mergeIUPACLetters(paste(sort(Bases2[!is.na(Bases2)]),
                                                   collapse="")))
        }
        else {
            primary <- c(primary, Bases[1])
            secondary <- c(secondary, Bases[1])
        }
    }

    qualityPhredScores <- qualityPhredScoresRaw[indexBaseCall]
    peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    primarySeq <- DNAString(paste(primary, collapse=""))
    secondarySeq <- DNAString(paste(secondary, collapse=""))
    message("     * Updating slots in 'SangerSingleRead' instance !!")

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

getStopList <- function(AASeqDF, selectChar, colorCode) {
    stopIndex <- AASeqDF %>% `==`(selectChar) %>% which()
    stopExcelIndex <- int2col(stopIndex)
    stopExcelIndexName <- paste0(stopExcelIndex, "1")
    styleList <- as.list(rep(paste('background-color:', colorCode, "; font-weight: bold;"), length(stopExcelIndex)))
    names(styleList) <- stopExcelIndexName
    return(styleList)
}
