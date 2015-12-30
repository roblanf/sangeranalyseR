#' Make a readgroup object
#' 
#' A readgroup object contains a set of reads that all correspond to a single underlying sequence. The function takes two lists of filenames: those for forward reads and those for reverse reads. It then efficiently reads in all of the files, and returns a single object which contains two lists of forward and reverse sequences corresponding to the forward and reverse filenames respectively.
#'
#' @export make.readgroup


make.readgroup <- function(fwd.fnames, rev.fnames, name = NULL, processors = NULL){

    if(is.null(processors)) { processors = detectCores(all.tests = FALSE, logical = FALSE)}

    fwd.reads = mclapply(fwd.fnames, readsangerseq, mc.cores = processors)
    rev.reads = mclapply(rev.fnames, readsangerseq, mc.cores = processors)

    readgroup = list(
                    list("forward.reads" = fwd.reads,
                    "reverse.reads" = rev.reads))

    names(readgroup) = name

    class(readgroup) = "readgroup"

    return(readgroup)

}
