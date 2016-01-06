#' Make a readgroup object
#' 
#' @param fwd.fnames a list of full file paths to forward reads from ab1 files (i.e. those that do not need to be reverse-complemented). 
#' @param rev.fnames a list of full file paths to reverse reads from ab1 files (i.e. those that *do* need to be reverse-complemented). 
#' @param trim TRUE/FALSE trim sequences based on quality while creating readgroup? If TRUE, the trim.mott function is applied to each sequence before inclusion in the readgroup. Note, trimming only works if the raw data are stored in ab1 files with appropriate information.
#' @param trim.cutoff value passed to trim.mott as quality cutoff for sequencing trimming, only used if 'trim' == TRUE
#' @param trim.segment value passed to trim.mott as minimum seuqence length for sequence trimming, only used if 'trim' == TRUE
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' A set of unaligned reads as a DNAstringset object, names are the input file paths.
#'
#' @export make.readset


make.readset <- function(fwd.fnames, rev.fnames, trim = TRUE, trim.cutoff = 0.05, trim.segment = 20, processors = NULL){

    processors = get.processors(processors)

    fwd.reads = mclapply(fwd.fnames, loadread, trim, trim.cutoff, trim.segment, revcomp = FALSE, processors = 1, mc.cores = processors)
    rev.reads = mclapply(rev.fnames, loadread, trim, trim.cutoff, trim.segment, revcomp = TRUE, processors = 1, mc.cores = processors)

    names(fwd.reads) = fwd.fnames
    names(rev.reads) = rev.fnames

    readset = DNAStringSet(c(fwd.reads, rev.reads))

    return(readset)

}

loadread <- function(fname, trim, trim.cutoff, trim.segment, revcomp, processors){

    read.abi = read.abif(fname)
    read.sanger = makeBaseCalls(sangerseq(read.abi))

    if(trim == TRUE){
        trims = trim.mott(read.abi, cutoff = trim.cutoff, segment = trim.segment)
        trims.fixed = fix.trims(trims, read.sanger, read.abi, processors)
        trim.start = trims.fixed$start
        trim.finish = trims.fixed$finish
    }else if(trim == FALSE){
        trim.start = 1
        trim.finish = length(read.sanger@primarySeq)
    }

    read.final = read.sanger@primarySeq[trim.start:trim.finish]

    if(revcomp == FALSE){
        return(read.final)
    }else if(revcomp == TRUE){
        return(reverseComplement(read.final))
    }


    return(read.final)

}

