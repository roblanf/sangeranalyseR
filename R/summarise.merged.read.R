#' Produce a one-row summary of a merged read
#' 
#' @param merged.read a merged.read object produced by the merge.reads() function
#'
#' @return A numeric vector including (note that insertions, deletions, and stop codons will only be called if the read was merged with a reference amino acid sequence):
#'          \enumerate{
#'              \item {consensus.len}: the length of the consensus sequence\cr
#'              \item {n.reads}: the number of reads used to build the consensus\cr
#'              \item {distance.min}: the minimum distance between two reads \cr
#'              \item {distance.max}: the maximum distance between two reads \cr
#'              \item {distance.med}: the median distance between pairs of reads \cr
#'              \item {readlen.min}: the length of the smallest read \cr
#'              \item {readlen.max}: the length of the longest read \cr
#'              \item {readlen.med}: the median of all read lengths \cr
#'              \item {insertions.min}: the smallest number of insertions inferred in a read \cr
#'              \item {insertions.max}: the larges number of insertions inferred in a read \cr
#'              \item {insertions.med}: the median number of insertions inferred in a read \cr
#'              \item {deletions.min}: the smallest number of deletions inferred in a read \cr
#'              \item {deletions.max}: the larges number of deletions inferred in a read \cr
#'              \item {deletions.med}: the median number of deletions inferred in a read \cr
#'              \item {stops.min}: the smallest number of stop codons inferred in a read \cr
#'              \item {stops.max}: the larges number of stop codons inferred in a read \cr
#'              \item {stops.med}: the median number of stop codons inferred in a read \cr
#'          }  
#'
#' @export summarise.merged.read


summarise.merged.read <- function(merged.read){

    if(class(merged.read != 'merged.read')){ stop("merged.read must be a merged.read object")}

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
