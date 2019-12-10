#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' contigName <- "RBNII395-13[C_LepFolF,C_LepFolR]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' A_chloroticContig <- SangerContig(
#'                                  parentDirectory       = inputFilesParentDir,
#'                                  contigName            = contigName,
#'                                  suffixForwardRegExp   = suffixForwardRegExp,
#'                                  suffixReverseRegExp   = suffixReverseRegExp,
#'                                  TrimmingMethod        = "M2",
#'                                  M1TrimmingCutoff      = NULL,
#'                                  M2CutoffQualityScore  = 40,
#'                                  M2SlidingWindowSize   = 10,
#'                                  baseNumPerRow         = 100,
#'                                  heightPerRow          = 200,
#'                                  signalRatioCutoff     = 0.33,
#'                                  showTrimmed           = TRUE)
#' RShinyCS <- launchAppSangerContig(list(A_chloroticContig))
launchAppSangerContig <- function(SangerContig, outputDir = NULL) {
    ### ------------------------------------------------------------------------
    ### Checking SangerContig input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir))
    }
    if (dir.exists(outputDir)) {
        shinyOptions(sangerContig = SangerContig)
        shinyOptions(shinyDirectory = outputDir)
        newSangerContig <- shinyApp(SangerContigUI, SangerContigServer,
                                           options = SangerContig)
        return(newSangerContig)
    } else {
        stop("'", outputDir, "' is not valid. Please check again")
    }
}

#' @export
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' newAlignment <- SangerAlignment(
#'                     parentDirectory       = rawDataDir,
#'                     suffixForwardRegExp   = suffixForwardRegExp,
#'                     suffixReverseRegExp   = suffixReverseRegExp,
#'                     refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                     TrimmingMethod        = "M2",
#'                     M1TrimmingCutoff      = NULL,
#'                     M2CutoffQualityScore  = 40,
#'                     M2SlidingWindowSize   = 10,
#'                     baseNumPerRow         = 100,
#'                     heightPerRow          = 200,
#'                     signalRatioCutoff     = 0.33,
#'                     showTrimmed           = TRUE)
#' RShinyCSSet <- launchAppSangerAlignment(list(newAlignment))
launchAppSangerAlignment <- function(SangerAlignment,
                                         outputDir = NULL) {
    ### ------------------------------------------------------------------------
    ### Checking SangerAlignment input parameter is a list containing
    ### one S4 object.
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir))
    }
    if (dir.exists(outputDir)) {
        shinyOptions(sangerAlignment = SangerAlignment)
        shinyOptions(shinyDirectory = outputDir)
        newSangerAlignment <- shinyApp(SangerAlignmentUI, SangerAlignmentServer,
                                       options = SangerAlignment)
        return(newSangerAlignment)
    } else {
        stop("'", outputDir, "' is not valid. Please check again")
    }
}

