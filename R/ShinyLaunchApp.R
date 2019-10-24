#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' samplesRegExp <- "ACHL"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory = inputFilesParentDir,
#'                                  readsRegularExp = samplesRegExp,
#'                                  cutoffQualityScore  = 50L,
#'                                  slidingWindowSize   = 8L)
#' consensusReadsList <- list(A_chloroticConsensusReads)
#' launchAppConsensusRead(consensusReadsList)
launchAppConsensusRead <- function(SangerConsensusRead) {
    ### ------------------------------------------------------------------------
    ### Checking SangerConsensusRead input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    shinyOptions(SangerConsensusReadSet = SangerConsensusRead)
    shinyApp(consensusUI, consensusServer, options = SangerConsensusRead)
}

#' @export
launchAppSangerProject <- function() {
    shinyApp(ui, server)
}
