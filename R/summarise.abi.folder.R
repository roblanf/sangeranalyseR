#' Create detailed summaries of all ABI sequencing reads in a folder (recursively searched)
#' 
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export summarise.abi.folder



summarise.abi.folder <- function(input.folder, trim.cutoff = 0.05, trim.segment = 20, secondary.peak.cutoff = 0.33, write.secondary.peak.files = FALSE, processors = NULL){

    processors = get.processors(processors)

    abi.fnames = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)

    abi.seqs = mclapply(abi.fnames, read.abif, mc.cores = processors)

    # now make a data.frame of summaries of all the files
    summaries = mclapply(abi.seqs, 
                         summarise.abi.file,
                         trim.cutoff = trim.cutoff,
                         trim.segment = trim.segment,
                         secondary.peak.cutoff = secondary.peak.cutoff,
                         write.secondary.peak.files = FALSE,
                         processors = 1,
                         mc.cores = processors  
                         )

    summaries = do.call(rbind, summaries)

    folder.names = basename(dirname(abi.fnames))
    file.names = basename(abi.fnames)

    summaries = cbind.data.frame("file.path" = as.character(abi.fnames), "folder.name" = as.character(folder.names), "file.name" = file.names, summaries, stringsAsFactors = FALSE)

    return(summaries)

}