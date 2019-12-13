IUPAC_CODE_MAP <- c(
    A="A",
    C="C",
    G="G",
    T="T",
    M="AC",
    R="AG",
    W="AT",
    S="CG",
    Y="CT",
    K="GT",
    V="ACG",
    H="ACT",
    D="AGT",
    B="CGT",
    N="ACGT"
)

mergeIUPACLettersInside <- function(x)
{
    message("Inside 'mergeIUPACLettersInside': ")
    if (!is.character(x) || any(is.na(x)) || any(nchar(x) == 0))
        stop("'x' must be a vector of non-empty character strings")
    x <- CharacterList(strsplit(toupper(x), "", fixed=TRUE))
    yy <- unname(IUPAC_CODE_MAP[unlist(x, use.names=FALSE)])
    if (any(is.na(yy)))
        stop("some strings in 'x' contain non IUPAC letters")
    yy <- CharacterList(strsplit(yy, "", fixed=TRUE))
    y <- unstrsplit(sort(unique(regroupBySupergroupInside(yy, x))))
    names(IUPAC_CODE_MAP)[match(y, IUPAC_CODE_MAP)]
    message("IUPAC_CODE_MAP: ", IUPAC_CODE_MAP)
    message("y: ", y)
    return(y)
}


regroupBySupergroupInside <- function(x, supergroups)
{
    message("Inside 'regroupBySupergroupInside': ")
    supergroups <- PartitioningByEnd(supergroups)
    x_breakpoints <- end(PartitioningByEnd(x))
    ans_breakpoints <- x_breakpoints[end(supergroups)]
    nleading0s <- length(supergroups) - length(ans_breakpoints)
    if (nleading0s != 0L)
        ans_breakpoints <- c(rep.int(0L, nleading0s), ans_breakpoints)
    ans_partitioning <- PartitioningByEnd(ans_breakpoints,
                                          names=names(supergroups))
    if (is(x, "PartitioningByEnd"))
        return(ans_partitioning)
    message("unlist(x, use.names=FALSE): ", unlist(x, use.names=FALSE))
    message('ans_partitioning@end: ', ans_partitioning@end)
    reGroup <- relist(unlist(x, use.names=FALSE), skeleton = ans_partitioning@end)
    return(reGroup)
}
