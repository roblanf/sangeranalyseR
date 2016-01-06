#' Make a readgroup object
#' 
#' @param fwd.fnames a list of full file paths to forward reads from ab1 files (i.e. those that do not need to be reverse-complemented). 
#' @param rev.fnames a list of full file paths to reverse reads from ab1 files (i.e. those that *do* need to be reverse-complemented). 
#' @param trim TRUE/FALSE trim sequences based on quality while creating readgroup? If TRUE, the trim.mott function is applied to each sequence before inclusion in the readgroup. Note, trimming only works if the raw data are stored in ab1 files with appropriate information.
#' @param trim.cutoff value passed to trim.mott as quality cutoff for sequencing trimming, only used if 'trim' == TRUE
#' @param max.secondary.peaks reads with more secondary peaks than this will not be included in the readset. The default (NULL) is to include all reads regardless of secondary peaks 
#' @param secondary.peak.ratio Only applies if max.secondary.peaks is not NULL. The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are counted. Those below the ratio are not. 
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' A set of unaligned reads as a DNAstringset object, names are the input file paths.
#'
#' @export make.readset


make.readset <- function(fwd.fnames, rev.fnames, trim = TRUE, trim.cutoff = 0.05, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33, processors = NULL){

    processors = get.processors(processors)

    fwd.dat = mclapply(fwd.fnames, loadread, trim, trim.cutoff, revcomp = FALSE, max.secondary.peaks = max.secondary.peaks, secondary.peak.ratio = secondary.peak.ratio, processors = 1, mc.cores = processors)
    rev.dat = mclapply(rev.fnames, loadread, trim, trim.cutoff, revcomp = TRUE,  max.secondary.peaks = max.secondary.peaks, secondary.peak.ratio = secondary.peak.ratio, processors = 1, mc.cores = processors)


    fwd.reads = lapply(fwd.dat, function(x) x[["read"]])
    rev.reads = lapply(rev.dat, function(x) x[["read"]])
    names(fwd.reads) = fwd.fnames
    names(rev.reads) = rev.fnames

    # remove the NULL reads
    all.reads = c(fwd.reads, rev.reads)
    all.reads = Filter(Negate(is.null), all.reads)
    readset = DNAStringSet(all.reads)

    # build the summary data frame just as in summarise.abi.folder
    fwd.summaries = lapply(fwd.dat, function(x) x[["summary"]])
    rev.summaries = lapply(rev.dat, function(x) x[["summary"]])
    all.summaries = c(fwd.summaries, rev.summaries)
    all.summaries = do.call(rbind, all.summaries)
    abi.fnames = unlist(c(fwd.fnames, rev.fnames))
    folder.names = basename(dirname(abi.fnames))
    file.names = basename(abi.fnames)
    all.summaries = cbind.data.frame("file.path" = as.character(abi.fnames), "folder.name" = as.character(folder.names), "file.name" = file.names, all.summaries, stringsAsFactors = FALSE)


    return(list("readset" = readset, "summaries" = all.summaries))

}

loadread <- function(fname, trim, trim.cutoff, revcomp, max.secondary.peaks, secondary.peak.ratio, processors){

    read.abi = read.abif(fname)

    s = summarise.abi.file(read.abi, trim.cutoff, secondary.peak.ratio, write.secondary.peak.files = FALSE, processors = processors)

    summary = s$summary
    read.sanger = s$read

    if(trim == TRUE){
        trim.start = summary["trim.start"]
        trim.finish = summary["trim.finish"]
        sp = summary["trimmed.secondary.peaks"]

    }else if(trim == FALSE){
        trim.start = 1
        trim.finish = length(read.sanger@primarySeq)
        sp = summary["raw.secondary.peaks"]
    }

    # Return NULL if the read fails filters
    read = read.sanger@primarySeq[trim.start:trim.finish]

    if(!is.null(max.secondary.peaks)){
        if(sp > max.secondary.peaks){
            read = NULL
        }
    }

    if(!is.null(read)) {
        if(revcomp == TRUE){
            read = reverseComplement(read)
        }
    }
    return(list('read' = read, summary = summary))

}

