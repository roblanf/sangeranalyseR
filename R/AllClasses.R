#' @title SangerReads
#'
#' @description  An S4 class for storing reads of sanger sequencing.
#'
#' @slot forward.read .
#' @slot reverse.read .
#' @slot consensus.read .
#' @slot quality.matrix .
#'
#' @name SangerReads-class
#'
#' @rdname SangerReads-class
#'
#' @exportClass SangerReads
#' @author Kuan-Hao Chao
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaForwardReadFN <- file.path(inputFilesPath,
#'                                        "Allolobophora_chlorotica",
#'                                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
#' A_chloroticaReverseReadFN <- file.path(inputFilesPath,
#'                                        "Allolobophora_chlorotica",
#'                                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R.ab1")
#' A_chloroticaRead <- new("SangerReads", forwardReadFileName = A_chloroticaForwardReadFN,
#'                                 reverseReadFileName = A_chloroticaReverseReadFN)
setClass("SangerReads",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             forwardReadFileName = "character",
             reverseReadFileName = "character",
             forwardRead         = "abif",
             reverseRead         = "abif"
         ),
)

# Constructor for SangerReads
setMethod("initialize", "SangerReads", function(.Object, ..., forwardReadFileName = forwardReadFileName, reverseReadFileName = reverseReadFileName, forwardRead=new("abif"), reverseRead=new("abif")) {
    ## do work of initialization
    errors <- character()
    if (!file.exists(forwardReadFileName)) {
        msg <- paste("\n'", forwardReadFileName, "'",
                     " foward read file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (!file.exists(reverseReadFileName)) {
        msg <- paste("\n'", reverseReadFileName, "'",
                     " reverse read file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) {
        message("Creating forwardRead ...\n")
        forwardRead = read.abif(forwardReadFileName)
        message("Creating reverseRead ...\n")
        reverseRead = read.abif(reverseReadFileName)
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   forwardReadFileName = forwardReadFileName,
                   reverseReadFileName = reverseReadFileName,
                   forwardRead=forwardRead,
                   reverseRead=reverseRead)
})
