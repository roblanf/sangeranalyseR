#' @title sangerSingleRead
#'
#' @description  An S4 class extending sangerseq S4 class
#'
#' @slot qualityReport .
#'
#' @name sangerSingleRead-class
#'
#' @rdname sangerSingleRead-class
#'
#' @exportClass sangerSingleRead
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R
#' @examples
setClass(
    "sangerSingleRead",
    contains="sangerseq",
    slots=c(rawData        = "abif",
            qualityReport  = "qualityReport")
) -> sangerSingleRead
