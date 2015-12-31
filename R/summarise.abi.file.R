#' Create a detailed summary of a single ABI sequencing file
#' 
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export summarise.abi.file

summarise.abi.file <- function(seq.abif, trim.cutoff = 0.05, trim.segment = 20, secondary.peak.cutoff = 0.33, write.secondary.peak.files = FALSE, processors = NULL){
 
    seq.sanger = sangerseq(seq.abif)

    # first we get the secondary peaks
    if(write.secondary.peak.files == TRUE){
        output.folder = dirname(inputfile)
        prefix = basename(inputfile)
        secondary.peaks = secondary.peaks(seq.sanger, secondary.peak.cutoff, output.folder, prefix, processors = processors)
    }else if(write.secondary.peak.files == FALSE){
        secondary.peaks = secondary.peaks(seq.sanger, secondary.peak.cutoff, processors = processors)
    }else{
        stop("Unknown option for write.secondary.peak.files. Should be TRUE or FALSE")
    }

    # now we trim the sequence
    trims = trim.mott(seq.abif, cutoff = trim.cutoff)
    trim.start = trims["trim.start"][[1]]
    trim.finish = trims["trim.finish"][[1]]

    # get trimmed and untrimmed version of raw data
    seq.trimmed = seq.sanger@primarySeq[trim.start:trim.finish]
    qual = seq.abif@data$PCON.2
    qual.trimmed = qual[trim.start:trim.finish]
    secondary.peaks.trimmed = subset(secondary.peaks, position > trim.start && position < trim.finish)

    # get quality at secondary peak sites
    raw.secondary.peak.qual     = mean(qual[secondary.peaks$position])
    trimmed.secondary.peak.qual = mean(qual[secondary.peaks.trimmed$position])

    read.summary = c("raw.length"                       = length(seq.sanger@primarySeq), 
                     "trimmed.length"                   = length(seq.trimmed),
                     "trim.start"                       = trim.start,
                     "trim.finish"                      = trim.finish,
                     "raw.secondary.peaks"              = nrow(secondary.peaks),
                     "trimmed.secondary.peaks"          = nrow(secondary.peaks.trimmed),
                     "raw.mean.quality"                 = mean(qual),
                     "trimmed.mean.quality"             = mean(qual.trimmed),
                     "raw.min.quality"                  = min(qual),
                     "trimmed.min.quality"              = min(qual.trimmed),
                     "raw.secondary.peak.quality"       = raw.secondary.peak.qual,
                     "trimmed.secondary.peak.quality"   = trimmed.secondary.peak.qual
                     )

    return(read.summary)

}