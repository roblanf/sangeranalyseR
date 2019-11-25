#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' consenesusReadName <- "RBNII395-13[C_LepFolF,C_LepFolR]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory       = inputFilesParentDir,
#'                                  consenesusReadName    = consenesusReadName,
#'                                  suffixForwardRegExp   = suffixForwardRegExp,
#'                                  suffixReverseRegExp   = suffixReverseRegExp,
#'                                  cutoffQualityScore    = 20,
#'                                  slidingWindowSize     = 5)
#' RShiny <- launchAppConsensusRead(list(A_chloroticConsensusReads))
launchAppConsensusRead <- function(SangerConsensusRead, directory = NULL) {
    ### ------------------------------------------------------------------------
    ### Checking SangerConsensusRead input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    if (is.null(directory)) {
        directory <- tempdir()
        suppressWarnings(dir.create(directory))
    }
    if (dir.exists(directory)) {
        shinyOptions(SangerConsensusRead = SangerConsensusRead)
        shinyOptions(shinyDirectory = directory)
        newSangerConsensusRead <- shinyApp(consensusReadUI, consensusReadServer,
                                           options = SangerConsensusRead)
        return(newSangerConsensusRead)
    } else {
        stop("'", directory, "' is not valid. Please check again")
    }
}

#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' SangerAlignedConsensusSet <- new("SangerAlignedConsensusSet",
#'                      parentDirectory       = inputFilesParentDir,
#'                      suffixForwardRegExp   = suffixForwardRegExp,
#'                      suffixReverseRegExp   = suffixReverseRegExp,
#'                      cutoffQualityScore    = 20,
#'                      slidingWindowSize     = 8)
#' RShiny <- launchAppAlignedConsensusSet(list(SangerAlignedConsensusSet))
launchAppAlignedConsensusSet <- function(SangerAlignedConsensusSet, directory = NULL) {
    ### ------------------------------------------------------------------------
    ### Checking AlignedConsensusSet input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    if (is.null(directory)) {
        directory <- tempdir()
        suppressWarnings(dir.create(directory))
    }
    if (dir.exists(directory)) {
        shinyOptions(SangerAlignedConsensusSet = SangerAlignedConsensusSet)
        shinyOptions(shinyDirectory = directory)
        newSangerAlignedConsensusSet <- shinyApp(alignedConsensusSetUI, alignedConsensusSetServer,
                                                 options = SangerAlignedConsensusSet)
        return(newSangerAlignedConsensusSet)
    } else {
        stop("'", directory, "' is not valid. Please check again")
    }
    # shinyOptions(SangerAlignedConsensusSet = SangerAlignedConsensusSet)
    # newSangerConsensusRead <- shinyApp(alignedConsensusSetUI, alignedConsensusSetServer,
    #                                    options = SangerConsensusRead)
    # return(newSangerConsensusRead)
}

