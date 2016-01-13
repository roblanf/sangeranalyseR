#' Automatically load readsets by grouping .ab1 files by name.
#' 
#' @param input.folder The parent folder of all of the reads contained in ab1 files you wish to analyse. Subfolders will be scanned recursively.
#' @param forward.suffix the suffix of the filenames for forward reads, i.e. reads that do not need to be reverse-complemented. Include the full suffix, e.g. "forward.ab1".
#' @param reverse.suffix the suffix of the filenames for reverse reads, i.e. reads that *do* need to be reverse-complemented. Include the full suffix, e.g. "reverse.ab1".
#' @param trim TRUE/FALSE trim sequences based on quality while creating readgroup? If TRUE, the trim.mott function is applied to each sequence before inclusion in the readgroup. Note, trimming only works if the raw data are stored in ab1 files with appropriate information.
#' @param trim.cutoff value passed to trim.mott as quality cutoff for sequencing trimming, only used if 'trim' == TRUE
#' @param max.secondary.peaks reads with more secondary peaks than this will not be included in the readset. The default (NULL) is to include all reads regardless of secondary peaks 
#' @param secondary.peak.ratio Only applies if max.secondary.peaks is not NULL. The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are counted. Those below the ratio are not. 
#' @param min.length reads shorter than this will not be included in the readset. The default (1) means that all reads with length of 1 or more will be included.
#' @param processors The number of processors to use, or NULL (the default) for all available processors
#'
#' @export make.readsets

make.readsets <- function(input.folder, forward.suffix, reverse.suffix, trim = TRUE, trim.cutoff = 0.0001, processors = NULL, min.length = 1, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33){

    processors = get.processors(processors)

    print("Looking for .ab1 files...")
    abi.files = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)

    # get a set of unique filenames after removing suffixes
    print(sprintf("Grouping %d files into sets by filename...", length(abi.files)))
    group.dataframe = get.group.dataframe(abi.files, forward.suffix, reverse.suffix)
    groups = unique(group.dataframe$group)

    print(sprintf("%d of the files matched either the forward or reverse suffixes", nrow(group.dataframe)))

    print(sprintf("Grouped these into %d sets of files", length(groups)))

    # load full file paths for readgroups based on unique filenames 
    print("Loading .ab1 files...")
    readset.fnames = mclapply(groups, 
                        get.readgroup.fnames,
                        group.dataframe = group.dataframe,
                        forward.suffix = forward.suffix,
                        reverse.suffix = reverse.suffix, 
                        mc.cores = processors)

    # how we parallelise depends on how many readgroups there are
    if(length(readset.fnames[[1]]) > length(readset.fnames)){
        # better to do readgroups sequentially, but parallelise each
        mc.cores = 1 
    }else{
        # better to do readgroups in parallel, but sequentially within each
        mc.cores = processors
        processors = 1
    }

    print("Constructing readsets...")
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

    print("Building read summaries...")
    readsets = lapply(rs, function(x) x[["readset"]])
    summaries = lapply(rs, function(x) x[["read.summaries"]])
    summaries = do.call(rbind, summaries)    
    rownames(summaries) = NULL

    summaries$readset.name = group.dataframe$group


    return(list("readsets" = readsets, "read.summaries" = summaries))

}


make.readset.from.list <- function(fnames, trim, trim.cutoff, max.secondary.peaks, secondary.peak.ratio, min.length, processors){
    # this just unpacks the fnames and sends them on, so I can use mclapply
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
    f.matches = str_match(files.cleaned, forward.suffix)
    f.indices = which(!is.na(f.matches))
    r.matches = str_match(files.cleaned, reverse.suffix)
    r.indices = which(!is.na(r.matches))

    #ONLY keep files that do match the suffixes:
    keep = c(f.indices, r.indices)
    files.cleaned = files.cleaned[keep]

    # try with and without the .ab1 after the suffixes, just in case
    files.cleaned = str_replace(files.cleaned, forward.suffix, replacement = "")
    files.cleaned = str_replace(files.cleaned, reverse.suffix, replacement = "")

    return(data.frame("file.path" = fname.list[keep], "group" = files.cleaned))

}


load.sangerseqs <- function(filenames){
    seqs = lapply(filenames, readsangerseq)
    return(seqs)
}


get.readgroup.fnames <- function(group, group.dataframe, forward.suffix, reverse.suffix){

    readgroup.fnames = as.character(group.dataframe$file.path[which(group.dataframe$group == group)])

    f.matches = str_match(readgroup.fnames, forward.suffix)
    f.indices = which(!is.na(f.matches))
    fwd.fnames = as.character(group.dataframe$file.path[f.indices])

    r.matches = str_match(readgroup.fnames, reverse.suffix)
    r.indices = which(!is.na(r.matches))
    rev.fnames = as.character(group.dataframe$file.path[r.indices])

    readgroup.fnames = list("forward.reads" = fwd.fnames, "reverse.reads" = rev.fnames)

    return(readgroup.fnames)

}