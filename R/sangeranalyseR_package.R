#' sangeranalyseR-package
#'
#' @name sangeranalyseR
#' @importFrom stringr str_extract
#' @importFrom sangerseqR sangerseq read.abif primarySeq chromatogram
#' @importFrom gridExtra grid.arrange
#' @import shiny
#' @importFrom shinydashboard dashboardPage
#' @importFrom data.table tstrsplit
#' @importFrom shinyjs useShinyjs html
#' @importFrom plotly plot_ly
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
