#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' forwardRegExp <- "^ACHLO([0-9]*)-09\\[LCO1490_t1,HCO2198_t1\\]_F.ab1$"
#' reverseRegExp <- "^ACHLO([0-9]*)-09\\[LCO1490_t1,HCO2198_t1\\]_R.ab1$"
#' A_chloroticConsensusReads <- new("SangerConsensusRead",
#'                                  parentDirectory       = inputFilesParentDir,
#'                                  forwardReadsRegularExp= forwardRegExp,
#'                                  reverseReadsRegularExp= reverseRegExp,
#'                                  cutoffQualityScore    = 50L,
#'                                  slidingWindowSize     = 8L)
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
