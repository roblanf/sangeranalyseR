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
    read.summaries = rs$read.summaries

    readset.lengths = unlist(lapply(readsets, function(x) length(x)))
    valid.readsets = readsets[which(readset.lengths >= min.reads)]
    valid.readset.lengths = unlist(lapply(valid.readsets, function(x) length(x)))


    if(median(valid.readset.lengths) > length(readsets)){
        # better to do readgroups sequentially, but parallelise each
        mc.cores = 1
        c.processors = processors
    }else{
        # better to do readgroups in parallel, but sequentially within each
        mc.cores = processors
        c.processors = 1
    }

    consensi = mclapply(valid.readsets,
                               merge.reads,
                               ref.aa.seq = ref.aa.seq, 
                               minInformation = minInformation, 
                               threshold = threshold, 
                               processors = c.processors, 
                               genetic.code = genetic.code, 
                               accept.stop.codons = accept.stop.codons, 
                               reading.frame = reading.frame,
                               mc.cores = mc.cores
                               )


    consensus.summaries = lapply(merged.readsets, summarise.merged.read)
    consensus.summaries = do.call(rbind, consensus.summaries)
    consensus.summaries = cbind("consensus.name" = row.names(consensus.summaries), consensus.summaries)
    row.names(consensus.summaries) = NULL
    consensus.summaries = data.frame(consensus.summaries)

    # which reads made it to the consensus sequence
    used.reads = unlist(lapply(merged.readsets, function(x) as.character(x$differences$name)))
    read.summaries$read.included.in.consensus = read.summaries$file.path %in% used.reads

    # a column for successful consensus sequence
    success = names(consensi)
    success.indices = which(read.summaries$readset.name %in% success)
    read.summaries$consensus.name[success.indices] = as.character(read.summaries$readset.name[success.indices])

    # Now we add more data from the read summaries
    used.read.summaries = filter(read.summaries, read.included.in.consensus==TRUE)    
    rsm = melt(used.read.summaries, id.vars = c("consensus.name", "folder.name", "file.name", "readset.name", "file.path", "read.included.in.readset", "read.included.in.consensus"))
    meds = dcast(rsm, consensus.name ~ variable, median)
    maxs = dcast(rsm, consensus.name ~ variable, max)
    mins = dcast(rsm, consensus.name ~ variable, min)
    more.summaries = data.frame("consensus.name" = meds$consensus.name,
                                "raw.secondary.peaks.min" = mins$raw.secondary.peaks, 
                                "raw.secondary.peaks.max" = maxs$raw.secondary.peaks,
                                "raw.secondary.peaks.med" = meds$raw.secondary.peaks,
                                "trimmed.secondary.peaks.min" = mins$trimmed.secondary.peaks, 
                                "trimmed.secondary.peaks.max" = maxs$trimmed.secondary.peaks,
                                "trimmed.secondary.peaks.med" = meds$trimmed.secondary.peaks,
                                "raw.mean.quality.min" = mins$raw.mean.quality, 
                                "raw.mean.quality.max" = maxs$raw.mean.quality,
                                "raw.mean.quality.med" = meds$raw.mean.quality,
                                "trimmed.mean.quality.min" = mins$trimmed.mean.quality, 
                                "trimmed.mean.quality.max" = maxs$trimmed.mean.quality,
                                "trimmed.mean.quality.med" = meds$trimmed.mean.quality
                               )
    consensus.summaries = merge(consensus.summaries, more.summaries, by = "consensus.name", sort = FALSE)

    # make the set of consensus sequences
    consensus.seqs = lapply(consensi, function(x) x$consensus)
    consensus.set  = DNAStringSet(consensus.seqs)

    # align the consensus sequences
    if(!is.null(ref.aa.seq)){
        aln = AlignTranslation(consensus.set, geneticCode = genetic.code, processors = processors, verbose = FALSE)
    }else{
        aln = AlignSeqs(consensus.set, processors = processors, verbose = FALSE)
    }

    # make a rough NJ tree. Labels are rows in the summary df
    aln.bin = as.DNAbin(aln)
    aln.dist = dist.dna(aln.bin)
    aln.tree = bionj(aln.dist)
    neat.labels = match(aln.tree$tip.label, 
                        as.character(consensus.summaries$consensus.name)
                        )
    aln.tree$tip.label = neat.labels

    return(list("read.summaries" = read.summaries, 
                "merged.reads" = consensi, 
                "consensus.summaries" = consensus.summaries, 
                "consensus.sequences" = consensus.set, 
                "consensus.alignment" = aln,
                "consensus.tree" = aln.tree
                )
            )

}