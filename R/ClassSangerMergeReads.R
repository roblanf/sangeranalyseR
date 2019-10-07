#' @title SangerMergeReads
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
#' @name SangerMergeReads-class
#'
#' @rdname SangerMergeReads-class
#'
#' @exportClass SangerMergeReads
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerSingleRead.R
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
#' A_chloroticaRvReadFN <- file.path(inputFilesPath,
#'                                   "Allolobophora_chlorotica",
#'                                   "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R.ab1")
#' A_chloroticaRead <- new("SangerMergeReads",
#'                         forwardReadFileName = A_chloroticaFdReadFN,
#'                         reverseReadFileName = A_chloroticaRvReadFN)
setClass("SangerMergeReads",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             forwardReadFileName     = "character",
             reverseReadFileName     = "character",
             forwardReadSangerseq    = "sangerSingleRead",
             reverseReadSangerseq    = "sangerSingleRead"
         ),
)

### ============================================================================
### Overwrite initialize for SangerReads (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerMergeReads",
          function(.Object, ...,
                   forwardReadFileName = forwardReadFileName,
                   reverseReadFileName = reverseReadFileName,
                   forwardReadSangerseq=new("sangerSingleRead"),
                   reverseReadSangerseq=new("sangerSingleRead")) {
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
        forwardReadSangerseq <- as(forwardReadSangerseq,"sangerSingleRead")
        forwardQualityRep <- new("qualityReport",
                                 qualityScoreNumeric =
                                     forwardReadRawAbif@data$PCON.2)
        forwardReadSangerseq@rawData <- forwardReadRawAbif
        forwardReadSangerseq@qualityReport <- forwardQualityRep


        message("Reverse read: Creating abif & sangerseq ...")
        message("    Creating reverseReadRawAbif ...")
        reverseReadRawAbif = read.abif(reverseReadFileName)
        message("    Creating reverseReadSangerseq ...\n")
        reverseReadSangerseq = sangerseq(reverseReadRawAbif)
        reverseReadSangerseq <- as(reverseReadSangerseq,"sangerSingleRead")
        reverseQualityRep <- new("qualityReport",
                                 qualityScoreNumeric =
                                     reverseReadRawAbif@data$PCON.2)
        reverseReadSangerseq@rawData <- reverseReadRawAbif
        reverseReadSangerseq@qualityReport <- reverseQualityRep

    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   forwardReadFileName  = forwardReadFileName,
                   reverseReadFileName  = reverseReadFileName,
                   forwardReadSangerseq = forwardReadSangerseq,
                   reverseReadSangerseq = reverseReadSangerseq)
})
