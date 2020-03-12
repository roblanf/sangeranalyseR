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
#'               CorrectFrameshifts DistanceMatrix IdClusters RemoveGaps BrowseSeqs
#' @importFrom Biostrings DNAString DNAStringSet AAStringSet GENETIC_CODE trinucleotideFrequency
#'               reverseComplement translate writeXStringSet translate compareStrings subseq
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom shinycssloaders withSpinner
#' @importFrom ggdendro ggdendrogram
#' @importFrom ape bionjs as.DNAbin dist.dna read.tree
#' @importFrom shinyWidgets actionBttn
#' @importFrom openxlsx int2col
#' @importFrom tools file_ext
#' @importFrom rmarkdown render
#' @importFrom kableExtra %>%
#' @importFrom excelR excelTable excelOutput renderExcel
#' @importFrom seqinr read.fasta
#' @importFrom parallel mclapply detectCores
NULL
