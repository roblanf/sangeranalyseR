#' Check for secondary peaks in a sangerseq object
#' 
#' This function finds and reports secondary peaks in a sangerseq object. It returns a table of secondary peaks, and optionally saves an annotated chromatogram and a csv file of the peak locations.
#' 
#' @param s a sangerseq s4 object from the sangerseqR package
#' @param ratio the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not. 
#' @param output.folder If output.folder is NA (the default) no files are written. If a valid folder is provided, two files are written to that folder: a .csv file of the secondary peaks (see description below) and a .pdf file of the chromatogram.
#' @param file.prefix If output.folder is specified, this is the prefix which will be appended to the .csv and the .pdf file. The default is "seq".
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#' 
#' @return A list with two elements:
#'          \enumerate{
#'              \item {secondary.peaks}: a data frame with one row per secondary peak above the ratio, and three columns: "position" is the position of the secondary peak relative to the primary sequence; "primary.basecall" is the primary base call; "secondary.basecall" is the secondary basecall. \cr
#'              \item {read}: the input sangerseq s4 object after having the makeBaseCalls() function from sangerseqR applied to it. This re-calls the primary and secondary bases in the sequence, and resets a lot of the internal data.
#'          }
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export secondary.peaks
#'

secondary.peaks <- function(s, ratio = 0.33, output.folder = NA, file.prefix = "seq", processors = NULL){
  
    # make secondary basecalls, and align them to the original sequence
    basecalls = makeBaseCalls(s, ratio = ratio)

    primary = primarySeq(basecalls, string = TRUE)
    secondary = secondarySeq(basecalls, string = TRUE)

    # perhaps we don't need to align...
    #seqs = DNAStringSet(c(primary, secondary))

    # these seuqences should be VERY similar...
    #pa = AlignSeqs(seqs, iterations = 0, refinements = 0, verbose = FALSE, processors = processors)


    # NB: it would seem to make more sense to use mismatchTable here, 
    # but I recoded it this way because mismatchTable had a bug.
    comp = compareStrings(primary, secondary)
    diffs = str_locate_all(pattern ='\\?',comp)[[1]][,1]
    primary.vector = strsplit(primary, split="")[[1]]
    secondary.vector = strsplit(secondary, split="")[[1]]

    primary.basecall    = primary.vector[diffs]
    secondary.basecall  = secondary.vector[diffs]

    r = data.frame("position" = diffs, "primary.basecall" = primary.basecall, "secondary.basecall" = secondary.basecall)

    if(!is.na(output.folder)){
        if(dir.exists(output.folder)){
            chromname = paste(file.prefix, "_", "chromatogram.pdf", sep='')
            tablename = paste(file.prefix, "_", "secondary_peaks.csv", sep='')
            chrom = chromatogram(basecalls, height = 2, showcalls = 'both', filename = file.path(output.folder, chromname))
            write.csv(r, file = file.path(output.folder, tablename))
        }else{
            warning(sprintf("Couldn't find directory '%s', no files saved", output.folder))
        }
    }

    return(list("secondary.peaks" = r, "read" = basecalls))
  
}
