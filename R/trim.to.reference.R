#' Trim an alignment to the limits of a supplied reference sequence
#'
#' This function profile aligns the alignment to the
#' reference sequence. It then trims the alignment to the limits of the 
#' reference sequence, and returns the trimmed alignment.
#'
#' @param alignment a DNAStringSet object
#' @param reference a DNA sequence as a DNAString sequence
#'
#' @export trim.to.reference

trim.to.reference <- function(alignment, reference){

    if(class(alignment)!='DNAStringSet'){ stop("alignment must be a DNAStringSet object")}
    if(class(reference)!='DNAString'){ stop("reference must be a DNAString object")}

    ref = DNAStringSet(reference)
    names(ref) = 'ref'
    aln = AlignProfiles(alignment, ref)

    ref.aligned = as.character(aln['ref'])
    not.gaps = str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    ref.start = min(not.gaps)
    ref.finish = max(not.gaps)

    aln.trimmed = subseq(alignment, ref.start, ref.finish)

    return(aln.trimmed)    

}
