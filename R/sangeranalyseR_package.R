#' sangeranalyseR-package
#'
#' @name sangeranalyseR
#' @importFrom sangerseqR sangerseq read.abif primarySeq chromatogram
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @import shiny
#' @importFrom shinydashboard dashboardPage
#' @importFrom data.table tstrsplit
#' @importFrom shinyjs useShinyjs html
#' @importFrom plotly orca
#' @importFrom DECIPHER AlignSeqs AlignTranslation ConsensusSequence
#'               CorrectFrameshifts DistanceMatrix IdClusters RemoveGaps
#' @importFrom Biostrings DNAStringSet GENETIC_CODE trinucleotideFrequency
#'               reverseComplement translate writeXStringSet translate
#' @importFrom DT dataTableOutput
#' @importFrom zeallot
#' @importFrom excelR
#' @importFrom shinycssloaders withSpinner
#' @importFrom ggdendro ggdendrogram
#' @importFrom ape BIONJ
#' @importFrom shinyWidgets actionBttn
#' @importFrom openxlsx int2col
#' @importFrom tools file_ext
#' @importFrom rmarkdown render
#' @importFrom kableExtra
#' @importFrom seqinr read.fasta
