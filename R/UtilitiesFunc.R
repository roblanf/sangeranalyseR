### ============================================================================
### Quality related: 'cutoffQualityScore' & 'slidingWindowSize' prechecking
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



#' Count stop codons in a DNA sequence
#'
#' @param sequence a DNAString object
#' @param readingFrame a number from 1 to 3 denoting the reading frame of the sequence
#' @param geneticCode Named character vector in the same format as GENETIC_CODE (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#'
#' @export count.stop.codons
#'
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

