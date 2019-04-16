#' Create a detailed summary of a single ABI sequencing file
#' 
#' @param seq.abif an abif.seq s4 object from the sangerseqR package
#' @param trim.cutoff the cutoff at which you consider a base to be bad. This works on a logarithmic scale, such that if you want to consider a score of 10 as bad, you set cutoff to 0.1; for 20 set it at 0.01; for 30 set it at 0.001; for 40 set it at 0.0001; and so on. Contiguous runs of bases below this quality will be removed from the start and end of the sequence. Given the high quality reads expected of most modern ABI sequencers, the defualt is 0.0001.
#' @param secondary.peak.ratio the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not. 
#' @param output.folder If output.folder is NA (the default) no files are written. If a valid folder is provided, two files are written to that folder: a .csv file of the secondary peaks (see description below) and a .pdf file of the chromatogram.
#' @param file.prefix If output.folder is specified, this is the prefix which will be appended to the .csv and the .pdf file. The default is "seq".
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' @return A numeric vector including:
#'          \enumerate{
#'              \item {raw.length}: the length of the untrimmed sequence, note that this is the sequence after conversion to a sangerseq object, and then the recalling the bases with MakeBaseCalls from the sangerseqR package\cr
#'              \item {trimmed.length}: the length of the trimmed sequence, after trimming using trim.mott from this package and the parameter supplied to this function \cr
#'              \item {trim.start}: the start position of the good sequence, see trim.mott for more details\cr
#'              \item {trim.finish}: the finish position of the good sequence, see trim.mott for more details\cr
#'              \item {raw.secondary.peaks}: the number of secondary peaks in the raw sequence, called with the secondary.peaks function from this package and the parameters supplied to this function \cr
#'              \item {trimmed.secondary.peaks}: the number of secondary peaks in the trimmed sequence, called with the secondary.peaks function from this package and the parameters supplied to this function \cr
#'              \item {raw.mean.quality}: the mean quality score of the raw sequence \cr
#'              \item {trimmed.mean.quality}: the mean quality score of the trimmed sequence \cr
#'              \item {raw.min.quality}: the minimum quality score of the raw sequence \cr
#'              \item {trimmed.min.quality}: the minimum quality score of the trimmed sequence \cr
#'          }  
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

    if(length(qual.trimmed)==0){qual.trimmed = c(NA)} # so we can summarise later

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

    #print(qual.trimmed)
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

    if(trims$start == 0 & trims$finish == 0){
        # no need to do anything fancy here...
        return(trims)
    }

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

