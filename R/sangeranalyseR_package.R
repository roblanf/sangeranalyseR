#' sangeranalyseR-package
#'
#' @name sangeranalyseR
#' @importFrom sangerseqR read.abif sangerseq primarySeq chromatogram
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @importFrom shiny shinyApp
#' @importFrom shinydashboard
#' @importFrom data.table tstrsplit
#' @importFrom shinyjs useShinyjs html
#' @importFrom plotly
#' @importFrom DECIPHER AlignSeqs AlignTranslation ConsensusSequence
#'               CorrectFrameshifts DistanceMatrix IdClusters RemoveGaps
#' @importFrom Biostrings DNAStringSet GENETIC_CODE trinucleotideFrequency
#'               reverseComplement translate
#' @importFrom DT dataTableOutput
#' @importFrom zeallot
#' @importFrom excelR
#' @importFrom shinycssloaders withSpinner
#' @importFrom ggdendro ggdendrogram
#' @importFrom ape BIONJ


#' @export
setOldClass("phylo")
