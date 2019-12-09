#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' consensusReadName <- "RBNII395-13[C_LepFolF,C_LepFolR]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' A_chloroticConsensusReads <- SangerConsensusRead(
#'                                  parentoutputDir       = inputFilesParentDir,
#'                                  consensusReadName    = consensusReadName,
#'                                  suffixForwardRegExp   = suffixForwardRegExp,
#'                                  suffixReverseRegExp   = suffixReverseRegExp,
#'                                  refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                                  TrimmingMethod        = "M1",
#'                                  M1TrimmingCutoff      = 0.0001,
#'                                  M2CutoffQualityScore  = NULL,
#'                                  M2SlidingWindowSize   = NULL,
#'                                  baseNumPerRow         = 80,
#'                                  heightPerRow          = 200,
#'                                  signalRatioCutoff     = 0.33,
#'                                  showTrimmed           = FALSE)
#' RShinyCS <- launchAppConsensusRead(list(A_chloroticConsensusReads))
launchAppConsensusRead <- function(SangerConsensusRead, outputDir = NULL) {
    ### ------------------------------------------------------------------------
    ### Checking SangerConsensusRead input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir))
    }
    if (dir.exists(outputDir)) {
        shinyOptions(SangerConsensusRead = SangerConsensusRead)
        shinyOptions(shinyoutputDir = outputDir)
        newSangerConsensusRead <- shinyApp(consensusReadUI, consensusReadServer,
                                           options = SangerConsensusRead)
        return(newSangerConsensusRead)
    } else {
        stop("'", outputDir, "' is not valid. Please check again")
    }
}

#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' SangerAlignedConsensusSet <- new("SangerAlignedConsensusSet",
#'                                  parentoutputDir       = rawDataDir,
#'                                  suffixForwardRegExp   = suffixForwardRegExp,
#'                                  suffixReverseRegExp   = suffixReverseRegExp,
#'                                  refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                                  TrimmingMethod        = "M1",
#'                                  M1TrimmingCutoff      = 0.001,
#'                                  M2CutoffQualityScore  = NULL,
#'                                  M2SlidingWindowSize   = NULL,
#'                                  baseNumPerRow         = 80,
#'                                  heightPerRow          = 200,
#'                                  signalRatioCutoff     = 0.23,
#'                                  showTrimmed           = FALSE)
#' RShinyCSSet <- launchAppAlignedConsensusSet(list(SangerAlignedConsensusSet))
launchAppAlignedConsensusSet <- function(SangerAlignedConsensusSet,
                                         outputDir = NULL) {
    ### ------------------------------------------------------------------------
    ### Checking AlignedConsensusSet input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir))
    }
    if (dir.exists(outputDir)) {
        shinyOptions(SangerAlignedConsensusSet = SangerAlignedConsensusSet)
        shinyOptions(shinyoutputDir = outputDir)
        newSangerAlignedConsensusSet <- shinyApp(alignedConsensusSetUI,
                                                 alignedConsensusSetServer,
                                                 options =
                                                     SangerAlignedConsensusSet)
        return(newSangerAlignedConsensusSet)
    } else {
        stop("'", outputDir, "' is not valid. Please check again")
    }
}

