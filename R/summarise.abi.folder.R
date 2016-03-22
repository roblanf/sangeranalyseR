#' Create detailed summaries of all ABI sequencing reads in a folder (recursively searched)
#' 
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export summarise.abi.folder



summarise.abi.folder <- function(input.folder, trim.cutoff = 0.0001, secondary.peak.ratio = 0.33, write.secondary.peak.files = FALSE, processors = NULL){

    processors = get.processors(processors)

    print("Looking for .ab1 files...")
    abi.fnames = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)

    print(sprintf(("Found %d .ab1 files..."), length(abi.fnames)))


    print("Loading reads...")
    abi.seqs = mclapply(abi.fnames, read.abif, mc.cores = processors)

    print("Calculating read summaries...")
    # now make a data.frame of summaries of all the files
    summaries.dat = mclapply(abi.seqs, 
                         summarise.abi.file,
                         trim.cutoff = trim.cutoff,
                         secondary.peak.ratio = secondary.peak.ratio,
                         processors = 1,
                         mc.cores = processors  
                         )

    print("Cleaning up")
    summaries = mclapply(summaries.dat, function(x) x[["summary"]], mc.cores = processors)
    summaries = do.call(rbind, summaries)

    reads = mclapply(summaries.dat, function(x) x[["read"]], mc.cores = processors)

    folder.names = basename(dirname(abi.fnames))
    file.names = basename(abi.fnames)

    summaries = cbind.data.frame("file.path" = as.character(abi.fnames), "folder.name" = as.character(folder.names), "file.name" = file.names, summaries, stringsAsFactors = FALSE)

    return(list("summaries" = summaries, "reads" = reads))

}