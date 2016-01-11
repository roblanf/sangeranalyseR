#' Create a detailed summary of a single ABI sequencing file
#' 
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export summarise.abi.file

summarise.abi.file <- function(seq.abif, trim.cutoff = 0.0001, secondary.peak.ratio = 0.33, output.folder = NA, prefix = "seq", processors = NULL){
 
    seq.sanger = sangerseq(seq.abif)

    # first we get the secondary peaks
    # note that the secondary peaks correspond to the seq.sanger object AFTER we
    # have called makeBaseCalls. And that this means the trim locations and the 
    # secondary peak locations do not match, since makeBaseCalls usually calls
    # fewer bases than the standard ABI calls.
    secondary.peaks.data = secondary.peaks(seq.sanger, secondary.peak.ratio, output.folder, prefix, processors = processors)
    secondary.peaks = secondary.peaks.data$secondary.peaks
    seq.sanger = secondary.peaks.data$read

    # now we trim the sequence
    trims = trim.mott(seq.abif, cutoff = trim.cutoff)
    qual = seq.abif@data$PCON.2
    qual.trimmed = qual[trims$start:trims$finish]

    # now we fix up the trim locations to correspond to the sangerseq primary seq object
    if(trims$start==1 && trims$finish==nchar(as.character(seq.abif@data$PBAS.2))){
        # there's nothing to do, because we didn't trim anything
        trim.start = 1
        trim.finish = length(primarySeq(seq.sanger))

    }else{
        trims.fixed = fix.trims(trims, seq.sanger, seq.abif, processors)
        trim.start = trims.fixed$start
        trim.finish = trims.fixed$finish
    }

    # get trimmed and untrimmed version of raw data
    seq.trimmed = seq.sanger@primarySeq[trim.start:trim.finish]
    secondary.peaks.trimmed = subset(secondary.peaks, position >= trim.start & position <= trim.finish)

    read.summary = c("raw.length"                       = length(seq.sanger@primarySeq), 
                     "trimmed.length"                   = length(seq.trimmed),
                     "trim.start"                       = trim.start,
                     "trim.finish"                      = trim.finish,
                     "raw.secondary.peaks"              = nrow(secondary.peaks),
                     "trimmed.secondary.peaks"          = nrow(secondary.peaks.trimmed),
                     "raw.mean.quality"                 = mean(qual),
                     "trimmed.mean.quality"             = mean(qual.trimmed),
                     "raw.min.quality"                  = min(qual),
                     "trimmed.min.quality"              = min(qual.trimmed)                     
                     )

    return(list("summary" = read.summary, "read" = seq.sanger))

}


fix.trims <- function(trims, seq.sanger, seq.abif, processors){

    # transfer trim locations from one sequence (denoted in the trims list, and which 
    # correspond to the seq.abif object to another 
    # the primarySeq(seq.sanger) from the seq.sanger object

    # 1. First we trim the original sequence
    original.seq = seq.abif@data$PBAS.2

    original.trimmed = substring(original.seq, trims$start, trims$finish)

    # 2. Align the original and recalled sequences
    recalled = primarySeq(seq.sanger, string = TRUE)
    seqs = DNAStringSet(c(original.trimmed, recalled))
    pa = AlignSeqs(seqs, iterations = 0, refinements = 0, verbose = FALSE, processors = processors)

    # 3. Get the sequence out, and find the first and last gaps.
    aligned.trimmed = as.character(pa[[1]])
    not.gaps = str_locate_all(aligned.trimmed, pattern = "[^-]")[[1]][,1]

    start = min(not.gaps)
    finish = max(not.gaps)

    if(start < 1){start = 1}
    if(finish > nchar(recalled)){finish = nchar(recalled)}

    return(list("start" = start, "finish" = finish))
}

