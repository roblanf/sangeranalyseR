#' sangeranalyseR-package
#'
#' @name sangeranalyseR
#' @importFrom stringr str_extract str_count
#' @importFrom sangerseqR sangerseq read.abif primarySeq chromatogram
#' @importFrom gridExtra grid.arrange
#' @importFrom shinydashboard renderMenu menuSubItem menuItem sidebarMenu
#'             updateTabItems box valueBox dashboardPage dashboardHeader
#'             dashboardSidebar sidebarMenuOutput dashboardBody
#' @importFrom shiny icon tags tagList isolate observeEvent getShinyOption
#'             reactiveValues renderUI h1 fluidRow column h4 uiOutput htmlOutput
#'             plotOutput sliderInput numericInput h3 checkboxInput stopApp
#'             showNotification includeHTML renderPlot textInput p
#'             removeNotification actionButton HTML textOutput verbatimTextOutput
#' @importFrom data.table tstrsplit
#' @importFrom shinyjs useShinyjs html
#' @importFrom plotly plot_ly add_markers add_trace layout add_annotations plotlyOutput renderPlotly plotly_build
#' @importFrom DECIPHER AlignSeqs AlignTranslation ConsensusSequence
#'               CorrectFrameshifts DistanceMatrix Treeline RemoveGaps BrowseSeqs
#' @importFrom Biostrings DNAString DNAStringSet AAStringSet GENETIC_CODE trinucleotideFrequency
#'               reverseComplement translate writeXStringSet translate subseq
#' @importFrom pwalign compareStrings
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom shinycssloaders withSpinner
#' @importFrom ggdendro ggdendrogram
#' @importFrom ape bionjs as.DNAbin dist.dna read.tree
#' @importFrom shinyWidgets actionBttn
#' @importFrom openxlsx int2col
#' @importFrom tools file_ext
#' @importFrom rmarkdown render
#' @importFrom excelR excelTable excelOutput renderExcel
#' @importFrom seqinr read.fasta
#' @importFrom parallel mclapply detectCores
#' @importFrom BiocParallel bplapply bpparam bpnworkers SerialParam MulticoreParam SnowParam
#' @importFrom Rcpp sourceCpp
#' @useDynLib sangeranalyseR, .registration = TRUE
#' @importFrom methods new is isVirtualClass setClass setClassUnion setGeneric
#'             setMethod setValidity slotNames validObject callNextMethod
#'             slot slot<-
#' @importFrom utils read.csv write.csv head tail capture.output data
#' @importFrom stats setNames IQR quantile aggregate
#' @importFrom stringr str_split
#' @importFrom Biostrings AAString
#' @importFrom plotly "%>%"
#' @importFrom grDevices colorRamp dev.off pdf rgb
#' @importFrom graphics axis lines mtext par rect
#' @importFrom ape as.phylo rtree
#' @importFrom shiny shinyApp shinyOptions
#' @importFrom S4Vectors isEmpty
#' @importFrom BiocGenerics width
#' @import logger
NULL

if (getRversion() >= "2.15.1") {
    # `<<-` assignments inside the Shiny servers walk past `getShinyOption(..)`
    # values that are seeded by launchAppSA / launchAppSC at startup. Declare
    # the names so R CMD check doesn't warn about "no visible binding".
    utils::globalVariables(c(
        "NEW_SANGER_ALIGNED_CONSENSUS_READ_SET",
        "NEW_SANGER_CONTIG"
    ))
}
