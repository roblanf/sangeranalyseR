#' @title SangerConsensusRead
#'
#' @description  An S4 class for storing multiple single reads to build up new
#'  consensus read
#'
#' @slot parentDirectory .
#' @slot readsRegularExp .
#' @slot cutoffQualityScore .
#' @slot slidingWindowSize .
#'
#' @name SangerConsensusRead-class
#'
#' @rdname SangerConsensusRead-class
#'
#' @exportClass SangerConsensusRead
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerSingleRead.R
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' samplesRegExp <- "ACHL"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory = inputFilesParentDir,
#'                                  readsRegularExp = samplesRegExp,
#'                                  cutoffQualityScore  = 50L,
#'                                  slidingWindowSize   = 8L)
setClass("SangerConsensusRead",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerConsensusRead'
         ### -------------------------------------------------------------------
         representation(
             parentDirectory    = "character",
             readsRegularExp    = "character",
             SangerReadsList    = "list"
         ),
)

### ============================================================================
### Overwrite initialize for 'SangerConsensusRead' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerConsensusRead",
          function(.Object, ...,
                   parentDirectory      = parentDirectory,
                   readsRegularExp      = readsRegularExp,
                   cutoffQualityScore   = 20L,
                   slidingWindowSize    = 5L) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()
    if (!file.exists(parentDirectory)) {
        msg <- paste("\n'", parentDirectory, "'",
                     " parent directory does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }

    parentDirFiles <- list.files(parentDirectory)
    selectInputFiles <- parentDirFiles[grepl(readsRegularExp, parentDirFiles)]
    allReads <- lapply(parentDirectory, file.path, selectInputFiles)
    readsNumber <- length(allReads[[1]])
    # sapply to check all files are exist.
    allErrorMsg <- sapply(c(allReads[[1]]), function(filePath) {
        if (!file.exists(filePath)) {
            msg <- paste("\n'", filePath, "'",
                         " read file does not exist.\n", sep = "")
            return(msg)
        }
        return()
    })
    errors <- c(errors, unlist(allErrorMsg), use.names = FALSE)

    ### ------------------------------------------------------------------------
    ### Prechecking success. Start to create multiple reads.
    ### ------------------------------------------------------------------------
    if (length(errors) == 0) {
        # sapply to create SangerSingleRead list.
        SangerSingleReadList <- sapply(allReads[[1]], SangerSingleRead,
               readFeature = "Reads", cutoffQualityScore, slidingWindowSize)
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory = parentDirectory,
                   readsRegularExp = readsRegularExp,
                   SangerReadsList = SangerSingleReadList)
})
