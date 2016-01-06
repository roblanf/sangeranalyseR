#' Automatically make consensus sequences by grouping .ab1 files by name.
#' 
#' @param input.folder The parent folder of all of the reads contained in ab1 files you wish to analyse. Subfolders will be scanned recursively.
#' @param forward.suffix the suffix of the filenames for forward reads, i.e. reads that do not need to be reverse-complemented. Include the full suffix, e.g. "forward.ab1".
#' @param reverse.suffix the suffix of the filenames for reverse reads, i.e. reads that *do* need to be reverse-complemented. Include the full suffix, e.g. "reverse.ab1".
#' @param min.reads The minimum number of reads required to make a consensus sequence, must be 2 or more (default 2). 
#' @param trim TRUE/FALSE trim sequences based on quality while creating readgroup? If TRUE, the trim.mott function is applied to each sequence before inclusion in the readgroup. Note, trimming only works if the raw data are stored in ab1 files with appropriate information.
#' @param trim.cutoff value passed to trim.mott as quality cutoff for sequencing trimming, only used if 'trim' == TRUE
#' @param max.secondary.peaks reads with more secondary peaks than this will not be included in the readset. The default (NULL) is to include all reads regardless of secondary peaks 
#' @param secondary.peak.ratio Only applies if max.secondary.peaks is not NULL. The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are counted. Those below the ratio are not. 
#' @param min.length reads shorter than this will not be included in the readset. The default (1) means that all reads with length of 1 or more will be included.
#' @param ref.aa.seq an amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand.
#' @param minInformation minimum fraction of the sequences required to call a consensus sequence at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75implying that 3/4 of all reads must be present in order to call a consensus.
#' @param threshold Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.05, implying that each consensus base can ignore at most 5 percent of the information at a given position. 
#' @param genetic.code Named character vector in the same format as GENETIC_CODE (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @param accept.stop.codons TRUE/FALSE. TRUE (the defualt): keep all reads, regardless of whether they have stop codons; FALSE: reject reads with stop codons. If FALSE is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
#' @param reading.frame 1, 2, or 3. Only used if accept.stop.codons == FALSE. This specifies the reading frame that is used to determine stop codons. If you use a ref.aa.seq, then the frame should always be 1, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame. 
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' @export make.consensus.seqs

make.consensus.seqs <- function(input.folder, forward.suffix, reverse.suffix, min.reads = 2, trim = TRUE, trim.cutoff = 0.0001, min.length = 1, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33, ref.aa.seq = NULL, minInformation = 0.75, threshold = 0.05, genetic.code = GENETIC_CODE, accept.stop.codons = TRUE, reading.frame = 1,  processors = NULL){

    processors = get.processors(processors)

    rs = make.readsets(input.folder = input.folder, 
                       forward.suffix = forward.suffix, 
                       reverse.suffix = reverse.suffix, 
                       trim = trim, 
                       trim.cutoff = trim.cutoff, 
                       min.length = min.length, 
                       max.secondary.peaks = max.secondary.peaks, 
                       secondary.peak.ratio = secondary.peak.ratio,
                       processors = processors
                       )

    readsets = rs$readsets
    summaries = rs$summaries

    readset.lengths = unlist(lapply(readsets, function(x) length(x)))
    valid.readsets = readsets[which(readset.lengths >= min.reads)]
    valid.readset.lengths = unlist(lapply(valid.readsets, function(x) length(x)))


    if(median(valid.readset.lengths) > length(readsets)){
        # better to do readgroups sequentially, but parallelise each
        mc.cores = 1 
    }else{
        # better to do readgroups in parallel, but sequentially within each
        mc.cores = processors
        processors = 1
    }

    merged.readsets = mclapply(valid.readsets,
                               merge.reads,
                               ref.aa.seq = ref.aa.seq, 
                               minInformation = minInformation, 
                               threshold = threshold, 
                               processors = processors, 
                               genetic.code = genetic.code, 
                               accept.stop.codons = accept.stop.codons, 
                               reading.frame = reading.frame,
                               mc.cores = mc.cores
                               )

    # which reads made it to the consensus sequence
    used.reads = unlist(lapply(merged.readsets, function(x) as.character(x$differences$name)))
    summaries$read.included.in.consensus = summaries$file.path %in% used.reads

    # a column for successful consensus sequence
    success = names(merged.readsets)
    success.indices = which(summaries$readset.name %in% success)
    summaries$consensus.name[success.indices] = as.character(summaries$readset.name[success.indices])


}