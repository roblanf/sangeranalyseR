#' Automatically load readsets by grouping .ab1 files by name.
#' 
#' Load all reads recursively, then group them by name after removing prefix/suffix. Return a list of lists of filenames of input reads.
#'
#' @export make.readsets


make.readsets <- function(input.folder, forward.suffix, reverse.suffix, processors = NULL, min.length = NULL, max.secondary.peaks = NULL){

    processors = get.processors(processors)

    abi.files = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)

    # get a set of unique filenames after removing suffixes
    group.dataframe = get.group.dataframe(abi.files, forward.suffix, reverse.suffix)
    groups = unique(group.dataframe$group)

    # load full file paths for readgroups based on unique filenames 
    # TODO: could be more efficient using split/apply
    readgroup.fnames = mclapply(groups, 
                        get.readgroup.fnames,
                        abi.files = abi.files, 
                        forward.suffix = forward.suffix,
                        reverse.suffix = reverse.suffix,
                        mc.cores = processors)

    # how we parallelise depends on how many readgroups there are
    if(length(readgroup.fnames[[1]]) > length(readgroup.fnames)){
        # better to do readgroups sequentially, but parallelise each
        mc.cores = 1 
    }else{
        # better to do readgroups in parrelle, but sequentially within each
        mc.cores = processors
        processors = 1
    }

    readgroups = mclapply(readgroup.fnames, get.readgroup, mc.cores = mc.cores, processors = processors)

    names(readgroups) = groups

    return(readgroups)

}

get.group.dataframe <- function(fname.list, forward.suffix, reverse.suffix){

    files.cleaned = fname.list
    files.cleaned = str_replace(files.cleaned, fwd.suffix, replacement = "")
    files.cleaned = str_replace(files.cleaned, rev.suffix, replacement = "")

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

    indices = which(!is.na(str_match(string = abi.files, pattern = group)))

    filenames = abi.files[indices]

    fwd.fnames = filenames[which(!is.na(str_match(filenames, forward.suffix)))]
    rev.fnames = filenames[which(!is.na(str_match(filenames, reverse.suffix)))]

    readgroup.fnames = list("forward.reads" = fwd.fnames, "reverse.reads" = rev.fnames)

    return(readgroup.fnames)

}