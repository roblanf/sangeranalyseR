#' Automatically load readsets by grouping .ab1 files by name.
#' 
#' Load all reads recursively, then group them by name after removing prefix/suffix. Return a list of lists of filenames of input reads.
#'
#' @param secondary.peak.ratio Only applies if max.secondary.peaks is not NULL. The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are counted. Those below the ratio are not. 
#'
#' @export make.readsets

make.readsets <- function(input.folder, forward.suffix, reverse.suffix, trim = TRUE, trim.cutoff = 0.0001, processors = NULL, min.length = 1, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33){

    processors = get.processors(processors)

    abi.files = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)

    # get a set of unique filenames after removing suffixes
    group.dataframe = get.group.dataframe(abi.files, forward.suffix, reverse.suffix)
    groups = unique(group.dataframe$group)

    # load full file paths for readgroups based on unique filenames 
    readset.fnames = mclapply(groups, 
                        get.readgroup.fnames,
                        abi.files = abi.files, 
                        forward.suffix = forward.suffix,
                        reverse.suffix = reverse.suffix,
                        mc.cores = processors)

    # how we parallelise depends on how many readgroups there are
    if(length(readset.fnames[[1]]) > length(readset.fnames)){
        # better to do readgroups sequentially, but parallelise each
        mc.cores = 1 
    }else{
        # better to do readgroups in parrelle, but sequentially within each
        mc.cores = processors
        processors = 1
    }

    rs = mclapply(readset.fnames, make.readset.from.list,
                          trim = trim,
                          trim.cutoff = trim.cutoff,
                          max.secondary.peaks = max.secondary.peaks,
                          secondary.peak.ratio = secondary.peak.ratio,
                          min.length = min.length,
                          processors = processors,                          
                          mc.cores = mc.cores
                          )

    names(rs) = groups

    readsets = lapply(rs, function(x) x[["readset"]])
    summaries = lapply(rs, function(x) x[["summaries"]])
    summaries = do.call(rbind, summaries)    
    rownames(summaries) = NULL

    summaries$group = group.dataframe$group

    # which reads did we end up using for the readsets?
    used.reads = unlist(lapply(readsets, function(x) names(x)))
    summaries$read.included = summaries$file.path %in% used.reads

    return(list("readsets" = readsets, "summaries" = summaries))

}


make.readset.from.list <- function(fnames, trim, trim.cutoff, max.secondary.peaks, secondary.peak.ratio, min.length, processors){
    # this just unpacks the fnames and sends them on, so I can use mclapply
    print(fnames)

    fwd.reads = fnames$forward.reads
    rev.reads = fnames$reverse.reads

    readset = make.readset(fwd.reads, rev.reads,
                           trim = trim,
                           trim.cutoff = trim.cutoff,
                           max.secondary.peaks = max.secondary.peaks,
                           secondary.peak.ratio = secondary.peak.ratio,
                           min.length = min.length,
                           processors = processors)

    return(readset)

}

get.group.dataframe <- function(fname.list, forward.suffix, reverse.suffix){

    files.cleaned = fname.list

    # try with and without the .ab1 after the suffixes, just in case
    files.cleaned = str_replace(files.cleaned, forward.suffix, replacement = "")
    files.cleaned = str_replace(files.cleaned, reverse.suffix, replacement = "")

    return(data.frame("file.path" = fname.list, "group" = files.cleaned))

}


load.sangerseqs <- function(filenames){
    seqs = lapply(filenames, readsangerseq)
    return(seqs)
}


filter.reads <- function(blah){
    #filter all the reads based on input criteria passed by user
}

get.readgroup <- function(readgroup.fnames, processors){

    # need to filter the reads as we go here...

    rg = make.readgroup(readgroup.fnames$forward.reads, 
                        readgroup.fnames$reverse.reads,
                        processors = processors)

    return(rg)
}

get.readgroup.fnames <- function(group, abi.files, forward.suffix, reverse.suffix){

    # we need to use teh base stri function to use literal strings
    # in case there are special characters in the filename
    indices = which(!is.na(stri_match_first_regex(abi.files, pattern = as.character(group), opts_regex = list("literal" = TRUE))))

    filenames = abi.files[indices]

    fwd.fnames = filenames[which(!is.na(stri_match_first_regex(filenames, pattern = forward.suffix, opts_regex = list("literal" = TRUE))))]
    rev.fnames = filenames[which(!is.na(stri_match_first_regex(filenames, pattern = reverse.suffix, opts_regex = list("literal" = TRUE))))]

    readgroup.fnames = list("forward.reads" = fwd.fnames, "reverse.reads" = rev.fnames)

    return(readgroup.fnames)

}