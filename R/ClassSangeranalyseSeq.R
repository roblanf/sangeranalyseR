#' @title SangeranalyseSeq
#'
#' @description  An S4 class extending sangerseq S4 class
#'
#' @slot qualityReport .
#'
#' @name SangeranalyseSeq-class
#'
#' @rdname SangeranalyseSeq-class
#'
#' @exportClass SangeranalyseSeq
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R
#' @examples
setClass(
    "SangeranalyseSeq",
    contains="sangerseq",
    slots=c(qualityReport="qualityReport")
) -> SangeranalyseSeq
