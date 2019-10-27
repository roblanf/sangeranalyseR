#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' samplesRegExp <- ".ab1"
#' B_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory = inputFilesParentDir,
#'                                  readsRegularExp = samplesRegExp,
#'                                  cutoffQualityScore  = 50L,
#'                                  slidingWindowSize   = 8L)
#' RShiny <- launchAppConsensusRead(list(B_chloroticConsensusReads))
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
