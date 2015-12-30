#' Group ABI reads by recursively reading them from an input folder
#' 
#' Load all reads recursively, then group them by name after removing prefix/suffix. Return a list of lists of filenames of input reads.
#'
#' @export load.readgroups


load.readgroups <- function(input.folder, forward.suffix, reverse.suffix, processors = NULL){

    if(is.null(processors)) { processors = detectCores(all.tests = FALSE, logical = FALSE)}

    abi.files = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)
    files.clean = abi.files
    files.clean = str_replace(files.clean, fwd.suffix, replacement = "")
    files.clean = str_replace(files.clean, rev.suffix, replacement = "")
    groups = unique(files.clean)
    readgroup.fnames = lapply(groups, 
                        get.readgroup.fnames,
                        abi.files = abi.files, 
                        forward.suffix = forward.suffix,
                        reverse.suffix = reverse.suffix)

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



load.sangerseqs <- function(filenames){
    seqs = lapply(filenames, readsangerseq)
    return(seqs)
}


get.readgroup <- function(readgroup.fnames, processors){

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