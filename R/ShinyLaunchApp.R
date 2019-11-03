#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' consenesusReadName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory       = inputFilesParentDir,
#'                                  consenesusReadName    = consenesusReadName,
#'                                  suffixForwardRegExp   = suffixForwardRegExp,
#'                                  suffixReverseRegExp   = suffixReverseRegExp,
#'                                  cutoffQualityScore    = 20,
#'                                  slidingWindowSize     = 8)
#' RShiny <- launchAppConsensusRead(list(A_chloroticConsensusReads))
launchAppConsensusRead <- function(SangerConsensusRead) {
    ### ------------------------------------------------------------------------
    ### Checking SangerConsensusRead input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    shinyOptions(SangerConsensusReadSet = SangerConsensusRead)
    newSangerConsensusRead <- shinyApp(consensusUI, consensusServer, options = SangerConsensusRead)
    return(newSangerConsensusRead)
}

#' @export
launchAppSangerProject <- function() {
    shinyApp(ui, server)
}
