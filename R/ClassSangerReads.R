#' @title SangerReads
#'
#' @description  An S4 class for storing reads of sanger sequencing.
#'
#' @slot forwardReadFileName .
#' @slot reverseReadFileName .
#' @slot forwardReadRawAbif .
#' @slot reverseReadRawAbif .
#' @slot forwardReadSangerseq .
#' @slot reverseReadSangerseq .
#'
#' @name SangerReads-class
#'
#' @rdname SangerReads-class
#'
#' @exportClass SangerReads
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangeranalyseSeq.R
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
#' A_chloroticaRvReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R.ab1")
#' A_chloroticaRead <- new("SangerReads",
#'                         forwardReadFileName = A_chloroticaFdReadFN,
#'                         reverseReadFileName = A_chloroticaRvReadFN)
setClass("SangerReads",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             forwardReadFileName     = "character",
             reverseReadFileName     = "character",
             forwardReadRawAbif      = "abif",
             reverseReadRawAbif      = "abif",
             forwardReadSangerseq    = "SangeranalyseSeq",
             reverseReadSangerseq    = "SangeranalyseSeq"
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
                   forwardReadRawAbif=new("abif"),
                   reverseReadRawAbif=new("abif"),
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
        ### --------------------------------------------------------------------
        ### Prechecking success. Start abif & sangerseq creation.
        ### --------------------------------------------------------------------
        message("Forward read: Creating abif & sangerseq ...")
        message("    Creating forwardReadRawAbif ...")
        forwardReadRawAbif = read.abif(forwardReadFileName)
        message("    Creating forwardReadSangerseq ...\n")
        forwardReadSangerseq = sangerseq(forwardReadRawAbif)
        forwardReadSangerseq <- as(forwardReadSangerseq,"SangeranalyseSeq")
        forwardQualityRep <- new("qualityReport",
                                 qualityScoreNumeric =
                                     forwardReadRawAbif@data$PCON.2)
        forwardReadSangerseq@qualityReport <- forwardQualityRep


        message("Reverse read: Creating abif & sangerseq ...")
        message("    Creating reverseReadRawAbif ...")
        reverseReadRawAbif = read.abif(reverseReadFileName)
        message("    Creating reverseReadSangerseq ...\n")
        reverseReadSangerseq = sangerseq(reverseReadRawAbif)
        reverseReadSangerseq <- as(reverseReadSangerseq,"SangeranalyseSeq")
        reverseQualityRep <- new("qualityReport",
                                 qualityScoreNumeric =
                                     reverseReadRawAbif@data$PCON.2)
        reverseReadSangerseq@qualityReport <- reverseQualityRep

    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   forwardReadFileName = forwardReadFileName,
                   reverseReadFileName = reverseReadFileName,
                   forwardReadRawAbif=forwardReadRawAbif,
                   reverseReadRawAbif=reverseReadRawAbif,
                   forwardReadSangerseq=forwardReadSangerseq,
                   reverseReadSangerseq=reverseReadSangerseq)
})
