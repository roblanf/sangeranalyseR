#' Trim an XStringSet object to the limits of a supplied reference sequence
#'
#' This function aligns all of the sequences in the XStringSet object along with
#' the reference sequence. It them trims the alignment to the limits of the 
#' reference sequence, and returns a trimmed alignment
#'
#' @param stringset an XStringSet object
#' @param reference an XString object of the same type as the XStringSet object
#'
#' @export trim.to.reference
