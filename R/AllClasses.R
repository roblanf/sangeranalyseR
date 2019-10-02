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
             forwardReadFileName  = "character",
             reverseReadFileName  = "character",
             forwardReadAbif      = "abif",
             reverseReadAbif      = "abif",
             forwardReadSangerseq = "sangerseq",
             reverseReadSangerseq = "sangerseq"
         ),
)

### ============================================================================
### Overwrite initialize for SangerReads (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerReads",
          function(.Object, ...,
                   forwardReadFileName = forwardReadFileName,
                   reverseReadFileName = reverseReadFileName,
                   forwardReadAbif=new("abif"),
                   reverseReadAbif=new("abif"),
                   forwardReadSangerseq=new("sangerseq"),
                   reverseReadSangerseq=new("sangerseq")) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
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
        ### ------------------------------------------------------------------------
        ### Prechecking success. Start abif & sangerseq creation.
        ### ------------------------------------------------------------------------
        message("Forward read: Creating abif & sangerseq ...")
        message("    Creating forwardReadAbif ...")
        forwardReadAbif = read.abif(forwardReadFileName)
        message("    Creating forwardReadSangerseq ...\n")
        forwardReadSangerseq = sangerseq(forwardReadAbif)

        message("Reverse read: Creating abif & sangerseq ...")
        message("    Creating reverseReadAbif ...")
        reverseReadAbif = read.abif(reverseReadFileName)
        message("    Creating reverseReadSangerseq ...\n")
        reverseReadSangerseq = sangerseq(reverseReadAbif)
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   forwardReadFileName = forwardReadFileName,
                   reverseReadFileName = reverseReadFileName,
                   forwardReadAbif=forwardReadAbif,
                   reverseReadAbif=reverseReadAbif,
                   forwardReadSangerseq=forwardReadSangerseq,
                   reverseReadSangerseq=reverseReadSangerseq)
})
