#' Produce a one-row summary of a merged read
#' 
#' @param merged.read a merged.read object produced by the merge.reads() function
#'
#' @return A one-row summary of the merged read
#'
#' @export summarise.merged.read


summarise.merged.read <- function(merged.read){

    m = merged.read
    reads = m$alignment[1:(length(m$alignment)-1)]
    read.lens = unlist(lapply(reads, function(x) length(DNAString(paste(del.gaps(x), collapse = '')))))

    # the NAs allow us to take min/max/med and get NA back
    # TODO: reduce these to refer to only the reads that made it...
    insertions = m$indels$insertions
    if(is.null(insertions)){ insertions = NA }
    deletions = m$indels$deletions
    if(is.null(deletions)){ deletions = NA }
    stops = m$stop.codons$stop.codons
    if(is.null(stops)){ stops = NA }


    summary = list(
                "consensus.len"     = length(m$consensus),
                "n.reads"           = length(reads),
                "distance.min"      = min(m$distance.matrix),
                "distance.max"      = max(m$distance.matrix),
                "distance.med"      = median(m$distance.matrix),
                "readlen.min"       = min(read.lens),
                "readlen.max"       = max(read.lens),
                "readlen.med"       = median(read.lens),
                "insertions.min"    = min(insertions),
                "insertions.max"    = max(insertions),
                "insertions.med"    = median(insertions),
                "deletions.min"     = min(deletions),
                "deletions.max"     = max(deletions),
                "deletions.med"     = median(deletions),
                "stops.min"         = min(stops),
                "stops.max"         = max(stops),
                "stops.med"         = median(stops)
                )

    return(summary)

}
