#' Merge sequence reads into a single consensus sequence
#' 
#' This function attempts to merge any number of unaligned reads into a single consensus sequence. It calculates a number of statistics that should help you decide whether a sensible consensus sequence exists for your data.\cr
#'
#' @param readset a DNAStringSet object with at least 2 reads
#' @param ref.aa.seq an amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand.
#' @param minInformation minimum fraction of the sequences required to call a consensus sequence at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75implying that 3/4 of all reads must be present in order to call a consensus.
#' @param threshold Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.05, implying that each consensus base can ignore at most 5 percent of the information at a given position. 
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#' @param genetic.code Named character vector in the same format as GENETIC_CODE (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @param accept.stop.codons TRUE/FALSE. TRUE (the defualt): keep all reads, regardless of whether they have stop codons; FALSE: reject reads with stop codons. If FALSE is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
#' @param reading.frame 1, 2, or 3. Only used if accept.stop.codons == FALSE. This specifies the reading frame that is used to determine stop codons. If you use a ref.aa.seq, then the frame should always be 1, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame. 
#' 
#' @return A list with the following components:
#'          \enumerate{
#'              \item {consensus}: the consensus sequence from the merged reads.\cr
#'              \item {alignment}: the alignment of all the reads, with the called consensus sequence. E.g. if the list was called 'merged.reads', you could use BrowseSeqs(merged.reads$alignment) to view the alignment.
#'              \item {differences}: a data frame of the number of pairwise differences between each read and the consensus sequence, as well as the number of bases in each input read that did not contribute to the consensus sequence. Can assist in detecting incorrect reads, or reads with a lot of errors.
#'              \item {distance.matrix}: a distance matrix of genetic distances (corrected with the JC model) between all of the input reads.
#'              \item {dendrogram}: a dendrogram depicting the distance.matrix. E.g. if the list was called 'merged.reads', you could use plot(merged.reads$dendrogram) to see the dendrogram.
#'              \item {indels}: if you specified a reference sequence via ref.aa.seq, then this will be a data frame describing the number of indels and deletions that were made to each of the input reads in order to correct frameshift mutations.
#'          }
#'
#' @keywords merge, reads, alignment, sequence
#'
#' @export merge.reads
#'

merge.reads <- function(readset, ref.aa.seq = NULL, minInformation = 0.75, threshold = 0.05, processors = NULL, genetic.code = GENETIC_CODE, accept.stop.codons = TRUE, reading.frame = 1){

    # check input options
    processors = get.processors(processors)
    

    # this sometimes happens when we automate things
    if(length(readset) < 2) {return(NULL)}

    # Try to correct frameshifts in the input sequences 
    if(!is.null(ref.aa.seq)) {
        print("Correcting frameshifts in reads using amino acid reference sequence")
        corrected = CorrectFrameshifts(myXStringSet = readset, myAAStringSet = AAStringSet(ref.aa.seq), geneticCode = genetic.code, type = 'both', processors = processors, verbose = FALSE)
        readset = corrected$sequences
        indels = get.indel.df(corrected$indels)        
        stops = as.numeric(unlist(mclapply(readset, count.stop.codons, reading.frame, genetic.code, mc.cores = processors)))
        stops.df = data.frame("read" = names(readset), "stop.codons" = stops)

    }else{
        indels = NULL
        stops.df = NULL
    }

    # Remove reads with stop codons
    if(accept.stop.codons == FALSE){
        print("Removing reads with stop codons")
        if(is.null(ref.aa.seq)){ # otherwise we already did it above
            stops = as.numeric(unlist(mclapply(readset, count.stop.codons, reading.frame, genetic.code, mc.cores = processors)))
            stops.df = data.frame("read" = names(readset), "stop.codons" = stops)
        }
        old_length = length(readset)
        readset = readset[which(stops==0)]
        print(sprintf("%d reads with stop codons removed", old_length - length(readset)))
    }

    if(length(readset) < 2) {return(NULL)}

    print("Aligning reads")
    if(!is.null(ref.aa.seq)){
        aln = AlignTranslation(readset, geneticCode = genetic.code, processors = processors, verbose = FALSE)
    }else{
        aln = AlignSeqs(readset, processors = processors, verbose = FALSE)
    }

    if(is.null(names(aln))){
        names(aln) = paste("read", 1:length(aln), sep="_")
    }

    # call consensus
    print("Calling consensus sequence")
    consensus = ConsensusSequence(aln,
                                  minInformation = minInformation,
                                  includeTerminalGaps = TRUE,
                                  ignoreNonBases = FALSE,
                                  threshold = threshold,
                                  noConsensusChar = "-")[[1]]

    print("Calculating differences between reads and consensus")
    diffs = mclapply(aln, n.pairwise.diffs, subject = consensus, mc.cores = processors)
    diffs = do.call(rbind, diffs)
    diffs.df = data.frame("name" = names(aln), "pairwise.diffs.to.consensus" = diffs[,1], "unused.chars" = diffs[,2])
    rownames(diffs.df) = NULL

    # get a dendrogram
    dist = DistanceMatrix(aln, correction = "Jukes-Cantor", penalizeGapLetterMatches = FALSE, processors = processors, verbose = FALSE)
    dend = IdClusters(dist, asDendrogram = TRUE, processors = processors, verbose = FALSE)

    # add consensus to alignment
    aln2 = c(aln, DNAStringSet(consensus))
    names(aln2)[length(aln2)] = "consensus"
    # strip gaps from consensus (must be an easier way!!)
    consensus.gapfree = DNAString(paste(del.gaps(consensus), collapse = ''))

    return(list("consensus" = consensus.gapfree, 
                "alignment" = aln2, 
                "differences" = diffs.df, 
                "distance.matrix" = dist,
                "dendrogram" = dend,
                "indels" = indels,
                "stop.codons" = stops.df)) 
}


get.indel.df <- function(indel.list){

    r = lapply(indel.list, indel.row)
    indel.df = data.frame(matrix(unlist(r), byrow = T, nrow = length(indel.list)))
    indel.df = cbind(names(indel.list), indel.df)
    names(indel.df) = c('read', 'insertions', 'deletions', 'distance')    
    return(indel.df)

}

indel.row <- function(row){

    n.ins = length(row$insertions)
    n.del = length(row$deletions)
    dist = row$distance
    return(c(n.ins, n.del, dist))
}

n.pairwise.diffs <- function(pattern, subject){
    # pairwise differences assuming pattern and subject are aligned
    comp = compareStrings(pattern, subject)
    qs = str_count(comp, '\\?')
    ps = str_count(comp, '\\+')
    return(c(qs, ps))
}

