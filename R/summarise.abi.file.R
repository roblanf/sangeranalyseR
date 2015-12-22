#' Create a detailed summary of a single ABI sequencing file
#' 
#' This function recursively scans a folder that contains abi files either in itself or its subfolders, and looks for secondary peaks in each sequence. It produces an annotated plot of the secondary peaks above a user-defined cutoff, as well as a table of the mismatches.
#' 
#' It recursively scanes the input folder for all abi files. For each abi file, it creates and returns a data frame and a plot of mismatches, and also writes both to the same file as the .ab1 file (as a .csv and .pdf respectively).
#'
#' @param input.folder the folder in which to search recursively for *.ab1 files.
#' @param cutoff the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not. 
#' @param write.files TRUE/FALSE. FALSE (the default): do not write any files to disk. TRUE: write plots and .csv files of secondary peaks to disk, in this case, files are named with the same filename as the .ab1 files, but with .csv and .pdf extensions, and are saved in the same folder as the .ab1 file from which they were created.
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#' 
#' @return output a dataframe of secondary peaks
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export summarise.abi.file

summarise.abi.file <- function(inputfile, trim.cutoff = 0.0001, secondary.peak.cutoff = 0.33, write.secondary.peak.files = FALSE, processors = NULL){
 
    seq.abif = read.abif(inputfile)
    seq = sangerseq(seq.abif)

    # first we get the secondary peaks
    if(write.secondary.peak.files == TRUE){
        output.folder = dirname(inputfile)
        prefix = basename(inputfile)
        secondary.peaks = secondary.peaks(seq, secondary.peak.cutoff, output.folder, prefix, processors = processors)
    }else if(write.secondary.peak.files == FALSE){
        secondary.peaks = secondary.peaks(seq, secondary.peak.cutoff, processors = processors)
    }else{
        stop("Unknown option for write.secondary.peak.files. Should be TRUE or FALSE")
    }

    # now we trim the sequence
    trims = trim.mott(seq.abif, cutoff = trim.cutoff)
    trim.start = trims["trim_start"][[1]]
    trim.finish = trims["trim_finish"][[1]]


    # get trimmed and untrimmed version of raw data
    seq.trimmed = seq@primarySeq[trim.start:trim.finish]
    qual = seq.abif@data$PCON.2
    qual.trimmed = qual[trim.start:trim.finish]
    secondary.peaks.trimmed = subset(secondary.peaks, position > trim.start && position < trim.end )

    return(c("raw length"               = length(seq@primarySeq), 
             "trimmed length"           = length(seq.trimmed), 
             "raw secondary peaks"      = nrow(secondary.peaks),
             "trimmed secondary peaks"  = nrow(secondary.peaks.trimmed),
             "raw mean quality"         = mean(qual),
             "trimmed mean quality"     = mean(qual.trimmed),
             "raw min quality"          = min(qual),
             "trimmed min quality"      = min(qual.trimmed)
             ))

}