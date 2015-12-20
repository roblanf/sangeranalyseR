#' Check for secondary peaks in a sangerseq object
#' 
#' This function scans finds and reports secondary peaks in a sangerseq object. It returns a table of secondary peaks, and optionally saves an annotated chromatogram and a csv file of the peak locations.
#' 
#' @param seq a sangerseq s4 object from the sangerseqR package
#' @param cutoff the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not. 
#' @param output.folder If output.folder is NA (the default) no files are written. If a valid folder is provided, two files are written to that folder: a .csv file of the secondary peaks (see description below) and a .pdf file of the chromatogram.
#' @param file.prefix If output.folder is specified, this is the prefix which will be appended to the .csv and the .pdf file. The default is "seq".
#' 
#' @return a data frame with one row per secondary peak above the cutoff, and three columns: "position" is the position of the secondary peak relative to the primary sequence; "primary.basecall" is the primary base call; "secondary.basecall" is the secondary basecall. 
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export secondary.peaks
#'

secondary.peaks <- function(seq, cutoff = 0.33, output.folder = NA, file.prefix = "seq"){
  
    # make secondary basecalls, and align them to the original sequence
    basecalls = makeBaseCalls(s, ratio = cutoff)
    pa = pairwiseAlignment(primarySeq(basecalls), secondarySeq(basecalls), type = "global-local")
    mismatches = mismatchTable(pa)

    r = data.frame("position" = mismatches$PatternStart, "primary.basecall" = mismatches$PatternSubstring, "secondary.basecall" = mismatches$SubjectSubstring)


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

    return(r)
  
}
