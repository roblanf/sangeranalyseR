### ============================================================================
### R shiny contigSeq server function
### ============================================================================
SangerContigServer <- function(input, output, session) {
    # Suppress Warning
    options(warn = -1)
    ### ------------------------------------------------------------------------
    ### SangerContig-related parameters initialization.
    ### ------------------------------------------------------------------------
    SangerContig <- getShinyOption("sangerContig")
    shinyDirectory <- getShinyOption("shinyDirectory")
    SangerContig <- SangerContig[[1]]
    trimmingMethod <- SangerContig@trimmingMethodSC
    if (trimmingMethod == "M1") {
        trimmingMethodName = "Method 1:
                               'Modified Mott Trimming'"
    } else if (trimmingMethod == "M2") {
        trimmingMethodName = "Method 2:
                               'Trimmomatics Sliding Window Trimming'"
    }

    ### ------------------------------------------------------------------------
    ### SangerRead-related parameters initialization.
    ### ------------------------------------------------------------------------
    forwardReadNum <- length(SangerContig@forwardReadList)
    reverseReadNum <- length(SangerContig@reverseReadList)
    SangerReadNum <- forwardReadNum + reverseReadNum
    # readFeature
    forwardReadFeature <- lapply(1:forwardReadNum, function(i)
        paste0(i, " ",
               SangerContig@forwardReadList[[i]]@readFeature))
    reverseReadFeature <- lapply(1:reverseReadNum, function(i)
        paste0(i, " ",
               SangerContig@reverseReadList[[i]]@readFeature))
    # readFileName (basename) (Fixed)
    forwardReadBFN <- lapply(1:forwardReadNum, function(i)
        basename(SangerContig@forwardReadList[[i]]@readFileName))
    reverseReadBFN <- lapply(1:reverseReadNum, function(i)
        basename(SangerContig@reverseReadList[[i]]@readFileName))
    SangerReadBFN <- c(forwardReadBFN, reverseReadBFN)

    ### ------------------------------------------------------------------------
    ### SangerContig reactiveValue
    ### ------------------------------------------------------------------------
    contigParam <-
        reactiveValues(
            contigSeq       = SangerContig@contigSeq,
            differencesDF   = SangerContig@differencesDF,
            alignment       = as.character(SangerContig@alignment),
            distanceMatrix  = SangerContig@distanceMatrix,
            dendrogram      = SangerContig@dendrogram,
            indelsDF        = SangerContig@indelsDF,
            stopCodonsDF    = SangerContig@stopCodonsDF,
            secondaryPeakDF = SangerContig@secondaryPeakDF)
    ### ------------------------------------------------------------------------
    ### SingleRead reactiveValue
    ### ------------------------------------------------------------------------
    sequenceParam <- reactiveValues(primarySeq = "",
                                    secondarySeq = "",
                                    primaryAASeqS1 = "",
                                    primaryAASeqS2 = "",
                                    primaryAASeqS3 = "")
    trimmedRV <- reactiveValues(rawSeqLength            = 0,
                                rawMeanQualityScore     = 0,
                                rawMinQualityScore      = 0,
                                trimmedStartPos         = 0,
                                trimmedFinishPos        = 0,
                                trimmedSeqLength        = 0,
                                trimmedMeanQualityScore = 0,
                                trimmedMinQualityScore  = 0,
                                remainingRatio          = 0)
    trimmedParam <- reactiveValues(M1TrimmingCutoff     = 0,
                                   M2CutoffQualityScore = 0,
                                   M2SlidingWindowSize  = 0)
    ChromatogramParam <- reactiveValues(baseNumPerRow     = 0,
                                        heightPerRow      = 0,
                                        signalRatioCutoff = 0,
                                        showTrimmed       = TRUE)

    ############################################################################
    ### output$ for all UI page
    ############################################################################
    ### ------------------------------------------------------------------------
    ### dynamic side menu bar
    ### ------------------------------------------------------------------------
    dynamicMenuSideBarSC(input, output, session,
                         forwardReadNum, reverseReadNum,
                         forwardReadFeature, reverseReadFeature)

    ############################################################################
    ### Main page switch
    ############################################################################
    output$SangerContig_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (input$sidebar_menu == "Sanger Contig Overview") {
            log_info(">>>>>>>> Inside '", input$sidebar_menu, "'")
            ### ----------------------------------------------------------------
            ### Dynamic page navigation: SangerContig overview page
            ### ----------------------------------------------------------------
            fluidRow(
                useShinyjs(),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Basic Information: "),
                                   style = "font-size: 26px;
                                   font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 2px hidden #A9A9A9;")),
                    fluidRow(
                        column(width = 12,
                               actionBttn("recalculateButton",
                                          "Re-calculate Contig",
                                          icon = icon("calculator"),
                                          style = "simple", color = "danger",
                                          block = TRUE, size = "lg")
                        ),
                        column(12,
                               tags$hr(
                                   style = ("border-top: 2px hidden #A9A9A9;")),
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Output Directory: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(shinyDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                "Raw ABI Parent Directory:"),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerContig@parentDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Contig Name: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerContig@contigName),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Trimming Method: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(trimmingMethodName),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Forward Suffix RegExp: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerContig@suffixForwardRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Forward Read Number: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(forwardReadNum),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Reverse Suffix RegExp: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerContig@suffixReverseRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Reverse Read Number: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(reverseReadNum),
                               )
                        ),
                    ),
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Contig Parameters",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(4,
                               uiOutput("SCMinReadsNum") ,
                        ),
                        column(4,
                               uiOutput("SCMinReadLength")  ,
                        ),
                        column(4,
                               uiOutput("SCMinFractionCall") ,
                        ),
                        column(4,
                               uiOutput("SCMaxFractionLost") ,
                        ),
                        column(4,
                               uiOutput("SCAcceptStopCodons") ,
                        ),
                        column(4,
                               uiOutput("SCReadingFrame") ,
                        ),
                    ),
                    box(title = tags$p("Genetic Code Data Frame",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 2,
                               tags$p("Tri-nucleotide:",
                                      style = "font-size: 15px;
                                      font-weight: bold;"),
                               tags$p("Amino Acid : ",
                                      style = "font-size: 15px;
                                      font-weight: bold;"),
                               tags$p("('*' : stop codon) ",
                                      style = "font-size: 12px;
                                       font-weight: italic;"),
                        ),
                        column(width = 10,
                               excelOutput("geneticCodeDF",
                                           width = "100%", height = "50"),
                               style = paste("height:100%; ",
                                             "overflow-y: hidden;",
                                             "overflow-x: scroll;")
                        ),
                    ),
                    uiOutput("SCrefAminoAcidSeq") ,
                ),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Contig Results: "),
                                   style = "font-size: 26px;
                                       font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Alignment",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               htmlOutput("contigAlignmentHTML"),
                        ),
                    ),
                    box(title = tags$p("Differences Data frame",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SCDifferencesDFUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Dendrogram",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               plotOutput("dendrogramPlot"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        ),
                        column(width = 12,
                               tags$hr(
                                   style = ("border-top: 4px hidden #A9A9A9;")),
                        ),
                        column(width = 12,
                               dataTableOutput("dendrogramDF"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Samples Distance",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               # plot()
                               uiOutput("SCDistanceMatrixPlotUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        ),
                        column(width = 12,
                               tags$hr(
                                   style = ("border-top: 4px hidden #A9A9A9;")),
                        ),
                        column(width = 12,
                               uiOutput("SCDistanceMatrixUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Indels Data frame",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SCIndelsDFUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Stop Codons Data frame",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SCStopCodonsDFUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    )
                ),
            )
        } else if (!is.na(strtoi(readIndex)) &&
                   (directionParam == "Forward" ||
                    directionParam == "Reverse")) {
            ### ----------------------------------------------------------------
            ### Dynamic page navigation: SangerRead page
            ### ----------------------------------------------------------------
            log_info(">>>>>>>> Inside '", input$sidebar_menu, "'")
            if (directionParam == "Forward") {
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerContig@
                                     forwardReadList[[readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerContig@
                                     forwardReadList[[readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerContig@
                                     forwardReadList[[readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerContig@
                                     forwardReadList[[readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerContig@
                                     forwardReadList[[readIndex]]@primaryAASeqS3)
                ChromatogramParam[["baseNumPerRow"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    ChromatogramParam@baseNumPerRow
                ChromatogramParam[["heightPerRow"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    ChromatogramParam@heightPerRow
                ChromatogramParam[["signalRatioCutoff"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    ChromatogramParam@signalRatioCutoff
                ChromatogramParam[["showTrimmed"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    ChromatogramParam@showTrimmed
                trimmedParam[["M1TrimmingCutoff"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@M1TrimmingCutoff
                trimmedParam[["M2CutoffQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@M2CutoffQualityScore
                trimmedParam[["M2SlidingWindowSize"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@M2SlidingWindowSize
                trimmedRV[["rawSeqLength"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerContig@forwardReadList[[readIndex]]@
                              QualityReport@remainingRatio * 100, 2)
                SSReadBFN <- basename(
                    SangerContig@forwardReadList[[readIndex]]@readFileName)
                SSReadAFN <- SangerContig@
                    forwardReadList[[readIndex]]@readFileName
            } else if (directionParam == "Reverse") {
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerContig@
                                     reverseReadList[[readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerContig@
                                     reverseReadList[[readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerContig@
                                     reverseReadList[[readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerContig@
                                     reverseReadList[[readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerContig@
                                     reverseReadList[[readIndex]]@primaryAASeqS3)
                ChromatogramParam[["baseNumPerRow"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    ChromatogramParam@baseNumPerRow
                ChromatogramParam[["heightPerRow"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    ChromatogramParam@heightPerRow
                ChromatogramParam[["signalRatioCutoff"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    ChromatogramParam@signalRatioCutoff
                ChromatogramParam[["showTrimmed"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    ChromatogramParam@showTrimmed
                trimmedParam[["M1TrimmingCutoff"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    M1TrimmingCutoff
                trimmedParam[["M2CutoffQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    M2CutoffQualityScore
                trimmedParam[["M2SlidingWindowSize"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    M2SlidingWindowSize
                trimmedRV[["rawSeqLength"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerContig@reverseReadList[[readIndex]]@
                              QualityReport@remainingRatio * 100, 2)
                SSReadBFN <- basename(
                    SangerContig@reverseReadList[[readIndex]]@readFileName)
                SSReadAFN <- SangerContig@
                    reverseReadList[[readIndex]]@readFileName
            }
            fluidRow(
                useShinyjs(),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Raw File: "),
                                   style = "font-size: 26px;
                                             font-weight: bold;"),
                    solidHeader = TRUE,
                    status = "success", width = 12,
                    h1(paste0(SSReadBFN)),
                    tags$h5(paste("( full path:", SSReadAFN, ")"),
                            style = "font-style:italic")),
                box(title =
                        tags$p(tagList(icon("dot-circle"),
                                       "DNA & Amino Acid Sequence
                                       (Before Trimming):"),
                                   style = "font-size: 26px;
                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    column(width = 12,
                           tags$p(tagList(icon("bars"),
                                          "Primary Sequence"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("primarySeqDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Secondary Sequence"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("secondSeqDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Quality Phred Score"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("qualityScoreDF",
                                       width = "100%", height = "50"),
                           style = paste("overflow-y: hidden;",
                                         "overflow-x: scroll;")
                    ),
                ),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Quality Report: "),
                                   style = "font-size: 26px;
                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p(tagList(icon("arrow-circle-right"),
                                               "Trimming Parameters Input"),
                                       style = "font-size: 24px;
                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        fluidRow(
                            column(width = 12,
                                   uiOutput("TrimmingMethodSelectionOutput"),
                            ),
                        ),
                        column(width = 12,
                               uiOutput("TrimmingMethodUI"),
                        ),
                        actionBttn("startTrimmingButton",
                                   "Apply Trimming Parameters",
                                   style = "simple", color = "danger",
                                   block = TRUE, size = "lg")
                    ),
                    box(title = tags$p(tagList(icon("arrow-circle-left"),
                                               "Trimmed Result Output"),
                                       style = "font-size: 24px;
                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        fluidRow(
                            box(title = tags$p("Before Trimming",
                                               style = "font-size: 21px;
                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       column(4,
                                              uiOutput("rawSeqLength") ,
                                       ),
                                       column(4,
                                              uiOutput("rawMeanQualityScore"),
                                       ),
                                       column(4,
                                              uiOutput("rawMinQualityScore"),
                                       ),
                                ),
                            ),
                        ),
                        fluidRow(
                            box(title = tags$p("After Trimming",
                                               style = "font-size: 21px;
                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       column(4,
                                              uiOutput("trimmedSeqLength"),
                                       ),
                                       column(4,
                                              uiOutput("trimmedMeanQualityScore"),
                                       ),
                                       column(4,
                                              uiOutput("trimmedMinQualityScore"),
                                       ),
                                ),

                                column(width = 12,
                                       column(4,
                                              uiOutput("trimmedStartPos") ,
                                       ),
                                       column(4,
                                              uiOutput("trimmedFinishPos") ,
                                       ),
                                       column(4,
                                              uiOutput("remainingRatio") ,
                                       )
                                ),
                            ),
                        ),
                        tags$hr(
                            style = ("border-top: 6px double #A9A9A9;")),
                        fluidRow(
                            box(title = tags$p("Quality Trimming Plot",
                                               style = "font-size: 21px;
                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                plotlyOutput("qualityQualityBasePlot") %>%
                                    withSpinner()),
                        ),
                    ),
                ),
                box(title =
                        tags$p(tagList(icon("dot-circle"),
                                       "DNA Sequence
                                       (After Trimming):"),
                               style = "font-size: 26px;
                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    column(width = 12,
                           tags$p(tagList(icon("bars"),
                                          "Primary Sequence"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("primarySeqTrimmedDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Secondary Sequence"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("secondSeqTrimmedDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Quality Phred Score"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("qualityScoreTrimmedDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "AA Sequence 1"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("PrimAASeqS1DF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "AA Sequence 2"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("PrimAASeqS2DF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "AA Sequence 3"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("PrimAASeqS3DF",
                                       width = "100%", height = "50"),
                           style = paste("overflow-y: hidden;",
                                         "overflow-x: scroll;")
                    )
                ),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Chromatogram: "),
                                   style = "font-size: 26px;
                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),

                    box(title = tags$p(tagList(icon("arrow-circle-right"),
                                               "Chromatogram Input"),
                                       style = "font-size: 24px;
                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(3,
                               sliderInput("ChromatogramBasePerRow",
                                           label =h4("Base Number Per Row"),
                                           min = 5,
                                           max = 200,
                                           value = ChromatogramParam[["baseNumPerRow"]]),
                               sliderInput("ChromatogramHeightPerRow",
                                           label = h4("Height Per Row"),
                                           min = 50,
                                           max = 600,
                                           value = ChromatogramParam[["heightPerRow"]]),
                        ),
                        column(3,
                               tags$hr(
                                   style =
                                       ("border-top: 4px hidden #A9A9A9;")),
                               numericInput(
                                   "ChromatogramSignalRatioCutoff",
                                   h3("Signal Ratio Cutoff"),
                                   value = ChromatogramParam[["signalRatioCutoff"]]),
                               checkboxInput(
                                   "ChromatogramCheckShowTrimmed",
                                   "Show trimmed region",
                                   value =
                                       ChromatogramParam[["showTrimmed"]])
                        ),
                        column(3,
                               tags$hr(
                                   style=("border-top: 4px hidden #A9A9A9;")),
                               uiOutput("ChromatogramtrimmedStartPos"),
                        ),
                        column(3,
                               tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                               uiOutput("ChromatogramtrimmedFinishPos"),
                        ),
                        actionBttn("saveChromatogramParam",
                                   "Apply Chromatogram Parameters",
                                   style = "simple", color = "danger",
                                   block = TRUE, size = "lg")
                    ),
                    box(title = tags$p(tagList(icon("arrow-circle-left"),
                                               "Chromatogram Output"),
                                       style = "font-size: 24px;
                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("chromatogramUIOutput"),
                        )
                    ),
                )
            )
        }
    })

    ############################################################################
    ### All other features (dynamic header / buttons)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic rightHeader text
    ### ------------------------------------------------------------------------
    #!!!!!!!! Fix
    # observeEventDynamicHeaderSC(input, output, session, trimmedRV)

    observeEvent(input$sidebar_menu, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (input$sidebar_menu == "Sanger Contig Overview") {
            html("rightHeader", "SangerContig Overview Page")
        } else if (!is.na(strtoi(readIndex)) &&
                   (directionParam == "Forward" ||
                    directionParam == "Reverse")) {
            html("rightHeader", paste(readIndex,
                                      directionParam, "SangerRead", "Page"))
        }
    })


















    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        log_info("@@@@@@@ 'close button' has been clicked")
        btn <- input$closeUI
        stopApp()
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        shinyjs::disable("closeUI")
        shinyjs::disable("saveS4")
        log_info("@@@@@@@ 'save button' has been clicked")
        newS4Object <- file.path(shinyDirectory, "SangerContig.Rda")
        showNotification(
            ui = column(12,
                       tags$p(tagList(icon("dot-circle"),
                                      "New S4 object is store as: "),
                              style = "font-size: 28px;
                                   font-weight: bold;"),
                       tags$p(paste0("'", newS4Object, "'"),
                              style = "font-size: 26px;
                                   font-style: italic"),
                       tags$br(),
                       tags$p(">> Run",
                              style = "font-size: 18px;
                                   font-style: italic"),
                       tags$p(paste0("readRDS(\"",newS4Object,
                                     "\")"),
                              style = "font-size: 18px;
                                      font-weight: bold;
                                      font-family: Courier, Monospace;"),
                       tags$p("   to load saved S4 object into R
                                     environment",
                              style = "font-size: 18px;
                                      font-style: italic")),
            closeButton = TRUE, id = "saveNotification",
            type = "message", duration = 10)
        ### --------------------------------------------------------------------
        ### Save SangerContig quality S4 object
        ### --------------------------------------------------------------------
        saveRDS(SangerContig, file=newS4Object)
        log_info("New S4 object is store as: ", newS4Object)
        NEW_SANGER_CONTIG <<- readRDS(file=newS4Object)
        shinyjs::enable("saveS4")
        shinyjs::enable("closeUI")
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Contig re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButton, {
        shinyjs::disable("closeUI")
        shinyjs::disable("recalculateButton")
        log_info("@@@@@@@ 'Reactive button' has been clicked")
        log_info("######## Start to recalculate contig")
        CSResult<-
            calculateContigSeq (SangerContig@inputSource,
                                SangerContig@forwardReadList,
                                SangerContig@reverseReadList,
                                SangerContig@refAminoAcidSeq,
                                SangerContig@minFractionCall,
                                SangerContig@maxFractionLost,
                                SangerContig@geneticCode,
                                SangerContig@acceptStopCodons,
                                SangerContig@readingFrame)
        SangerContig@contigSeq <<- CSResult$consensusGapfree
        SangerContig@differencesDF <<- CSResult$diffsDf
        SangerContig@alignment <<- CSResult$aln2
        SangerContig@distanceMatrix <<- CSResult$dist
        SangerContig@dendrogram <<- CSResult$dend
        SangerContig@indelsDF <<- CSResult$indels
        SangerContig@stopCodonsDF <<- CSResult$stopsDf
        SangerContig@secondaryPeakDF <<- CSResult$spDf
        contigParam[["contigSeq"]] <<- SangerContig@contigSeq
        contigParam[["differencesDF"]] <<- SangerContig@differencesDF
        contigParam[["alignment"]] <<- as.character(SangerContig@alignment)
        contigParam[["distanceMatrix"]] <<-SangerContig@distanceMatrix
        contigParam[["dendrogram"]] <<- SangerContig@dendrogram
        contigParam[["indelsDF"]] <<- SangerContig@indelsDF
        contigParam[["stopCodonsDF"]] <<- SangerContig@stopCodonsDF
        contigParam[["secondaryPeakDF"]] <<- SangerContig@secondaryPeakDF
        log_info("######## Finish recalculating contig")
        shinyjs::enable("recalculateButton")
        shinyjs::enable("closeUI")
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button chromatogram parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$saveChromatogramParam, {
        log_info("@@@@@@@ 'Reactive button' has been clicked")
        log_info("######## Start recalculating chromatogram")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        ### --------------------------------------------------------------------
        ### Update Chromatogram parameters
        ### --------------------------------------------------------------------
        if (!is.na(strtoi(readIndex))) {
            if (directionParam == "Forward") {
                SangerContig@forwardReadList[[readIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerContig@forwardReadList[[readIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerContig@forwardReadList[[readIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerContig@forwardReadList[[readIndex]]@ChromatogramParam@
                    showTrimmed <<- input$ChromatogramCheckShowTrimmed
            } else if (directionParam == "Reverse") {
                SangerContig@reverseReadList[[readIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerContig@reverseReadList[[readIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerContig@reverseReadList[[readIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerContig@reverseReadList[[readIndex]]@ChromatogramParam@
                    showTrimmed <<- input$ChromatogramCheckShowTrimmed
            }
            ### ------------------------------------------------------------
            ### Save 'ChromatogramParam' dynamic value
            ### ------------------------------------------------------------
            ChromatogramParam[["baseNumPerRow"]] <<- input$ChromatogramBasePerRow
            ChromatogramParam[["heightPerRow"]] <<- input$ChromatogramHeightPerRow
            ChromatogramParam[["signalRatioCutoff"]] <<-
                input$ChromatogramSignalRatioCutoff
            ChromatogramParam[["showTrimmed"]] <<-
                input$ChromatogramCheckShowTrimmed
        }
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button trimming parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$startTrimmingButton, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            if (SangerContig@trimmingMethodSC == "M1") {
                if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
                    as.numeric(input$M1TrimmingCutoffText) > 0 &&
                    as.numeric(input$M1TrimmingCutoffText) <= 1) {
                    inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
                } else {
                    inputM1TrimmingCutoffText <- 0.0001
                }
                if (directionParam == "Forward") {
                    ### --------------------------------------------------------
                    ### Start M1 trimming calculation
                    ### --------------------------------------------------------
                    trimmingPos <- M1inside_calculate_trimming(
                            SangerContig@forwardReadList[[readIndex]]@
                                QualityReport@qualityPhredScores,
                            SangerContig@forwardReadList[[readIndex]]@
                                QualityReport@qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerContig@forwardReadList[[readIndex]]@
                        QualityReport@M1TrimmingCutoff <<-
                        as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerContig@forwardReadList[[readIndex]]@
                        QualityReport@M1TrimmingCutoff
                } else if (directionParam == "Reverse") {
                    ### --------------------------------------------------------
                    ### Start M1 trimming calculation
                    ### --------------------------------------------------------
                    trimmingPos <- M1inside_calculate_trimming(
                            SangerContig@reverseReadList[[readIndex]]@
                                QualityReport@qualityPhredScores,
                            SangerContig@reverseReadList[[readIndex]]@
                                QualityReport@qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerContig@reverseReadList[[readIndex]]@
                        QualityReport@M1TrimmingCutoff <<-
                        as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerContig@reverseReadList[[readIndex]]@
                        QualityReport@M1TrimmingCutoff
                }
            } else if (SangerContig@trimmingMethodSC == "M2") {
                if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
                    strtoi(input$M2CutoffQualityScoreText) > 0 &&
                    strtoi(input$M2CutoffQualityScoreText) <= 60 &&
                    strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
                    inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
                } else {
                    inputM2CutoffQualityScoreText <- 20
                }
                if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
                    strtoi(input$M2SlidingWindowSizeText) > 0 &&
                    strtoi(input$M2SlidingWindowSizeText) <= 40 &&
                    strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
                    inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
                } else {
                    inputM2SlidingWindowSizeText <- 10
                }
                if (directionParam == "Forward") {
                    ### ------------------------------------------------------------
                    ### Start M2 trimming calculation
                    ### ------------------------------------------------------------
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerContig@forwardReadList[[readIndex]]@
                                QualityReport@qualityPhredScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText))
                    SangerContig@forwardReadList[[readIndex]]@
                        QualityReport@M2CutoffQualityScore <<-
                        strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerContig@forwardReadList[[readIndex]]@
                        QualityReport@M2CutoffQualityScore
                    SangerContig@forwardReadList[[readIndex]]@
                        QualityReport@M2SlidingWindowSize <<-
                        strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerContig@forwardReadList[[readIndex]]@
                        QualityReport@M2SlidingWindowSize
                } else if (directionParam == "Reverse") {
                    ### ------------------------------------------------------------
                    ### Start M2 trimming calculation
                    ### ------------------------------------------------------------
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerContig@reverseReadList[[readIndex]]@
                                QualityReport@qualityPhredScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText))
                    SangerContig@reverseReadList[[readIndex]]@
                        QualityReport@M2CutoffQualityScore <<-
                        strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerContig@reverseReadList[[readIndex]]@
                        QualityReport@M2CutoffQualityScore
                    SangerContig@reverseReadList[[readIndex]]@
                        QualityReport@M2SlidingWindowSize <<-
                        strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerContig@reverseReadList[[readIndex]]@
                        QualityReport@M2SlidingWindowSize
                }
            }
            if (directionParam == "Forward") {
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerContig@forwardReadList[[readIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]
                trimmedRV[["rawSeqLength"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerContig@forwardReadList[[readIndex]]@
                              QualityReport@remainingRatio * 100, 2)
            } else if (directionParam == "Reverse") {
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerContig@reverseReadList[[readIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]
                trimmedRV[["rawSeqLength"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerContig@reverseReadList[[readIndex]]@
                              QualityReport@remainingRatio * 100, 2)
            }
        }
    })

# <============================================================================>

    ############################################################################
    ### SangerContig (Function for Sanger Contig Overview)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Valuebox for basic information
    ### ------------------------------------------------------------------------
    valueBoxSCMinReadsNum(input, output,
                          SangerContig@minReadsNum, session)
    valueBoxSCMinReadLength(input, output,
                            SangerContig@minReadLength, session)
    valueBoxSCMinFractionCall(input, output,
                              SangerContig@minFractionCall, session)
    valueBoxSCMaxFractionLost(input, output,
                              SangerContig@maxFractionLost, session)
    valueBoxSCAcceptStopCodons(input, output,
                               SangerContig@acceptStopCodons, session)
    valueBoxSCReadingFrame(input, output,
                           SangerContig@readingFrame, session)
    ### ------------------------------------------------------------------------
    ### geneticCodeDF
    ### ------------------------------------------------------------------------
    output$geneticCodeDF <- renderExcel({
        suppressMessages(
            excelTable(data = t(data.frame(SangerContig@geneticCode)),
                       defaultColWidth = 50, editable = FALSE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })
    ### ------------------------------------------------------------------------
    ### refAminoAcidSeq
    ### ------------------------------------------------------------------------
    output$SCrefAminoAcidSeq <- renderUI({
        if (SangerContig@refAminoAcidSeq == "") {
            box(title = tags$p("Reference Amino Acid Sequence",
                               style = "font-size: 24px;
                                    font-weight: bold;"),
                collapsible = TRUE,
                status = "success", width = 12,
                column(width = 1),
                column(width = 11,
                       h4("Reference Amino Acid Sequence is not provided."))
            )
        } else {
            box(title = tags$p("Reference Amino Acid Sequence",
                               style = "font-size: 24px;
                                    font-weight: bold;"),
                collapsible = TRUE,
                status = "success", width = 12,
                column(width = 2,
                       tags$br(),
                       tags$p("AA Sequence:",
                              style = "font-size: 15px;
                                       font-weight: bold;"),
                ),
                column(width = 10,
                       excelOutput("SCrefAminoAcidSeqDF",
                                   width = "100%", height = "50"),
                       style = paste("height:100%; ",
                                     "overflow-y: hidden;",
                                     "overflow-x: scroll;")
                ),
            )
        }
    })
    output$SCrefAminoAcidSeqDF <- renderExcel({
        refAminoAcidSeqVec <- strsplit(SangerContig@refAminoAcidSeq, "")[[1]]
        names(refAminoAcidSeqVec) <- c(1:length(refAminoAcidSeqVec))
        suppressMessages(
            excelTable(data =
                           t(data.frame(refAminoAcidSeqVec)),
                       defaultColWidth = 50, editable = FALSE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })
    ### ------------------------------------------------------------------------
    ### Alignment
    ### ------------------------------------------------------------------------
    output$contigAlignmentHTML <- renderUI({
        browseSeqHTML <-
            file.path(shinyDirectory,
                      paste0(SangerContig@contigName,
                             "_Alignment_BrowseSeqs.html"))
        BrowseSeqs(DNAStringSet(contigParam[["alignment"]]),
                   openURL=FALSE, htmlFile=browseSeqHTML)
        column(width = 12,
               includeHTML(
                   file.path(shinyDirectory,
                             paste0(SangerContig@contigName,
                                    "_Alignment_BrowseSeqs.html"))),
               style = paste("height:100%; ",
                             "overflow-y: hidden;",
                             "overflow-x: scroll;")
        )
    })
    ### ------------------------------------------------------------------------
    ### Difference
    ### ------------------------------------------------------------------------
    output$SCDifferencesDFUI <- renderUI({
        if (all(dim(contigParam[["differencesDF"]]) == c(0,0))) {
            h4("*** 'Differences' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDifferencesDF")
        }
    })
    output$SCDifferencesDF = renderDataTable({
        contigParam[["differencesDF"]]
    })
    ### ------------------------------------------------------------------------
    ### dendrogram
    ### ------------------------------------------------------------------------
    output$dendrogramPlot <- renderPlot({
        # plot(contigParam[["dendrogram"]][[2]])
        ggdendrogram(contigParam[["dendrogram"]][[2]], rotate = TRUE)
    })
    output$dendrogramDF <- renderDataTable({
        contigParam[["dendrogram"]][[1]]
    })
    ### ------------------------------------------------------------------------
    ### Distance
    ### ------------------------------------------------------------------------
    output$SCDistanceMatrixPlotUI <- renderUI({
        if (all(dim(contigParam[["distanceMatrix"]]) == c(0,0))) {
            h4("*** 'Distance' dataframe is empty. (Cannot plot)***",
               style="font-weight: bold; font-style: italic;")
        } else {
            plotlyOutput("SCDistanceMatrixPlot")
        }
    })
    output$SCDistanceMatrixPlot <- renderPlotly({
        plot_ly(x = SangerReadBFN,
                y = SangerReadBFN,
                z = contigParam[["distanceMatrix"]],
                colors = colorRamp(c("white", "#32a852")),
                type = "heatmap")
    })
    output$SCDistanceMatrixUI <- renderUI({
        if (all(dim(contigParam[["distanceMatrix"]]) == c(0,0))) {
            h4("*** 'Distance' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDistanceMatrix")
        }
    })
    output$SCDistanceMatrix = renderDataTable({
        contigParam[["distanceMatrix"]]
    })
    ### ------------------------------------------------------------------------
    ### Indels
    ### ------------------------------------------------------------------------
    output$SCIndelsDFUI <- renderUI({
        if (all(dim(contigParam[["indelsDF"]]) == c(0,0))) {
            h4("*** 'Indels' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCIndelsDF")

        }
    })
    output$SCIndelsDF <- renderDataTable({
        contigParam[["indelsDF"]]
    })
    ### ------------------------------------------------------------------------
    ### StopCodons
    ### ------------------------------------------------------------------------
    output$SCStopCodonsDFUI <- renderUI({
        if (all(dim(contigParam[["stopCodonsDF"]]) == c(0,0))) {
            h4("*** 'Stop Codons' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCStopCodonsDF")
        }
    })
    output$SCStopCodonsDF <- renderDataTable({
        contigParam[["stopCodonsDF"]]
    })

# <============================================================================>

    ############################################################################
    ### SangerRead (Function for singel read)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Primary dataframe
    ### ------------------------------------------------------------------------
    output$primarySeqDF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(readIndex)) {
            primarySeqDisplay (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Secondary dataframe
    ### ------------------------------------------------------------------------
    output$secondSeqDF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(readIndex)) {
            secondarySeqDisplay (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Quality Score dataframe
    ### ------------------------------------------------------------------------
    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            if (directionParam == "Forward") {
                PhredScore <- SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                PhredScore <- SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityScoreDisplay (PhredScore)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acid dataframe (1)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS1DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            PrimAASeqS1Display (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acid dataframe (2)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS2DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            PrimAASeqS2Display (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acid dataframe (3)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS3DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            PrimAASeqS3Display (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Trimmed dataframe
    ### ------------------------------------------------------------------------
    output$primarySeqTrimmedDF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(readIndex)) {
            primarySeqTrimmedDisplay (input, output, session,
                                      sequenceParam, trimmedRV)
        }
    })
    ### ------------------------------------------------------------------------
    ### Secondary Trimmed dataframe
    ### ------------------------------------------------------------------------
    output$secondSeqTrimmedDF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(readIndex)) {
            secondSeqTrimmedDisplay (input, output, session,
                                     sequenceParam, trimmedRV)
        }
    })
    ### ------------------------------------------------------------------------
    ### Quality Score Trimmed dataframe
    ### ------------------------------------------------------------------------
    output$qualityScoreTrimmedDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(readIndex)) {
            if (directionParam == "Forward") {
                PhredScore <- SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@qualityPhredScores[
                        (trimmedRV[["trimmedStartPos"]]+1):
                            trimmedRV[["trimmedFinishPos"]]]
            } else if (directionParam == "Reverse") {
                PhredScore <- SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@qualityPhredScores[
                        (trimmedRV[["trimmedStartPos"]]+1):
                            trimmedRV[["trimmedFinishPos"]]]
            }
            qualityScoreDisplay (PhredScore)
        }
    })
    ############################################################################
    ### Trimming method selection functions
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Trimming Method Selection (Just check forward reads => not dynamic)
    ### ------------------------------------------------------------------------
    output$TrimmingMethodSelectionOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            if (SangerContig@forwardReadList[[1]]@
                QualityReport@TrimmingMethod == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Modified Mott Trimming'")
            } else if (SangerContig@forwardReadList[[1]]@
                       QualityReport@TrimmingMethod == "M2") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Trimmomatics Sliding Window Trimming'")
            }
        }
    })
    output$TrimmingMethodUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            ## For method, everyone is same, so just pick forward one.
            if (SangerContig@trimmingMethodSC == "M1") {
                if (directionParam == "Forward") {
                    if (is.null(SangerContig@forwardReadList[[readIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerContig@forwardReadList[[readIndex]]@
                            QualityReport@M1TrimmingCutoff <<-  0.0001
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerContig@reverseReadList[[readIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerContig@reverseReadList[[readIndex]]@
                            QualityReport@M1TrimmingCutoff <<-  0.0001
                    }
                }
                fluidRow(
                    column(6,
                           uiOutput("M1TrimmingCutoff") ,
                           tags$ul(
                               textInput("M1TrimmingCutoffText",
                                         label = p("Input Value"),
                                         value = toString(
                                             trimmedParam[["M1TrimmingCutoff"]]),
                                         width = '70%')
                           ),
                    ),
                )
            } else if (SangerContig@trimmingMethodSC == "M2") {
                if (directionParam == "Forward") {
                    if (is.null(SangerContig@forwardReadList[[readIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerContig@forwardReadList[[readIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerContig@forwardReadList[[readIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerContig@forwardReadList[[readIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  10
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerContig@reverseReadList[[readIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerContig@reverseReadList[[readIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerContig@reverseReadList[[readIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerContig@reverseReadList[[readIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  10
                    }
                }
                fluidRow(
                    column(6,
                           uiOutput("M2CutoffQualityScore") ,
                           tags$ul(
                               textInput("M2CutoffQualityScoreText",
                                         label = p("Input Value"),
                                         value = toString(
                                             trimmedParam[["M2CutoffQualityScore"]]),
                                         width = '70%')
                           ),
                    ),
                    column(6,
                           uiOutput("M2SlidingWindowSize") ,
                           tags$ul(
                               textInput("M2SlidingWindowSizeText",
                                         label = p("Input Value"),
                                         value = toString(
                                             trimmedParam[["M2SlidingWindowSize"]]),
                                         width = '70%')
                           ),
                    ),
                )
            }
        }
    })
    ### ------------------------------------------------------------------------
    ### Quality trimming related (value box)
    ### ------------------------------------------------------------------------
    valueBoxM1TrimmingCutoff (input, output, session)
    valueBoxM2CutoffQualityScore (input, output, session)
    valueBoxM2SlidingWindowSize (input, output, session)

    valueBoxRawSeqLength (input, output, session, trimmedRV)
    valueBoxRawMeanQualityScore (input, output, session, trimmedRV)
    valueBoxRawMinQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxTrimmedFinishPos (input, output, session, trimmedRV)
    valueBoxTrimmedSeqLength (input, output, session, trimmedRV)
    valueBoxTrimmedMeanQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedMinQualityScore (input, output, session, trimmedRV)
    valueBoxRemainingRatio (input, output, session, trimmedRV)

    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            if (directionParam == "Forward") {
                qualityPhredScores <-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                qualityPhredScores <-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityQualityBasePlotDisplay(input, output,session,
                                          trimmedRV, qualityPhredScores)
        }
    })

    valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)
    ### ------------------------------------------------------------------------
    ### chromatogram related feature
    ### ------------------------------------------------------------------------
    output$chromatogramUIOutput <- renderUI({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        trimmedRV[["trimmedSeqLength"]]
        if (!is.na(strtoi(readIndex))) {
            if (directionParam == "Forward") {
                chromatogramRowNumAns <-
                    chromatogramRowNum (
                        strtoi(ChromatogramParam[["baseNumPerRow"]]),
                        SangerContig@forwardReadList[[readIndex]]@
                            QualityReport@rawSeqLength,
                        SangerContig@forwardReadList[[readIndex]]@
                            QualityReport@trimmedSeqLength,
                        ChromatogramParam[["showTrimmed"]]) *
                    strtoi(ChromatogramParam[["heightPerRow"]])

            } else if (directionParam == "Reverse") {
                chromatogramRowNumAns <-
                    chromatogramRowNum (
                        strtoi(ChromatogramParam[["baseNumPerRow"]]),
                        SangerContig@reverseReadList[[readIndex]]@
                            QualityReport@rawSeqLength,
                        SangerContig@reverseReadList[[readIndex]]@
                            QualityReport@trimmedSeqLength,
                        ChromatogramParam[["showTrimmed"]]) *
                    strtoi(ChromatogramParam[["heightPerRow"]])
            }
            plotOutput("chromatogram", height = chromatogramRowNumAns)
        }
    })

    output$chromatogram <- renderPlot({
        ## !!!!! Update !!!!
        shinyjs::disable("startTrimmingButton")
        shinyjs::disable("saveChromatogramParam")
        shinyjs::disable("ChromatogramBasePerRow")
        shinyjs::disable("ChromatogramHeightPerRow")
        shinyjs::disable("ChromatogramSignalRatioCutoff")
        shinyjs::disable("ChromatogramCheckShowTrimmed")
        shinyjs::disable("M1TrimmingCutoffText")
        shinyjs::disable("M2CutoffQualityScoreText")
        shinyjs::disable("M2SlidingWindowSizeText")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(readIndex))) {
            removeNotification(id = "saveNotification")
            showNotification(
                ui = fluidRow(
                    column(1),
                    column(11,
                           tags$p(tagList(icon("hourglass-start"),
                                          "Plotting Chromatogram ...  "),
                                  style = "font-size: 28px;
                                   font-weight: bold;"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Base pairs per row : ",
                                          ChromatogramParam[["baseNumPerRow"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Height per row : ",
                                          ChromatogramParam[["heightPerRow"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Signal Ratio Cutoff : ",
                                          ChromatogramParam[["signalRatioCutoff"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Show trimmed : ",
                                          ChromatogramParam[["showTrimmed"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p("( If 'Signal Ratio Cutoff' is too small,
                                  it would need more time to replot
                                  the chromatogram)",
                                  style = "font-size: 16px;
                                   font-style: italic"),
                           )
                ), action = NULL, duration = NULL, closeButton = FALSE,
                id = "chromatogramNotification", type = "message")
            if (directionParam == "Forward") {
                rawSeqLength <-
                    SangerContig@forwardReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                log_info(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (forward)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerContig@forwardReadList[[readIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ----------------------------------------------------------------
                ### Update 'SangerContig'!
                ### ----------------------------------------------------------------
                SangerContig@forwardReadList[[readIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerContig@forwardReadList[[readIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerContig@forwardReadList[[readIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerContig@forwardReadList[[readIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ----------------------------------------------------------------
                ### Updating AASeqs
                ### ----------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerContig@forwardReadList[[readIndex]]@QualityReport@trimmedStartPos,
                                               SangerContig@forwardReadList[[readIndex]]@QualityReport@trimmedFinishPos,
                                               SangerContig@geneticCode)
                SangerContig@forwardReadList[[readIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerContig@forwardReadList[[readIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerContig@forwardReadList[[readIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ----------------------------------------------------------------
                ### Updating reactive values
                ### ----------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerContig@forwardReadList[[
                        readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerContig@forwardReadList[[
                        readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerContig@forwardReadList[[
                        readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerContig@forwardReadList[[
                        readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerContig@forwardReadList[[
                        readIndex]]@primaryAASeqS3)
            } else if (directionParam == "Reverse") {
                rawSeqLength <-
                    SangerContig@reverseReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                log_info(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (reverse)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerContig@reverseReadList[[readIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ----------------------------------------------------------------
                ### Update 'SangerContig'!
                ### ----------------------------------------------------------------
                SangerContig@reverseReadList[[readIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerContig@reverseReadList[[readIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerContig@reverseReadList[[readIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerContig@reverseReadList[[readIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ----------------------------------------------------------------
                ### Updating AASeqs
                ### ----------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerContig@reverseReadList[[readIndex]]@QualityReport@trimmedStartPos,
                                               SangerContig@reverseReadList[[readIndex]]@QualityReport@trimmedFinishPos,
                                               SangerContig@geneticCode)
                SangerContig@reverseReadList[[readIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerContig@reverseReadList[[readIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerContig@reverseReadList[[readIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ----------------------------------------------------------------
                ### Updating reactive values
                ### ----------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerContig@reverseReadList[[
                        readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerContig@reverseReadList[[
                        readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerContig@reverseReadList[[
                        readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerContig@reverseReadList[[
                        readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerContig@reverseReadList[[
                        readIndex]]@primaryAASeqS3)
            }

            # log_info(">>>>>>>>>>>> 'MakeBaseCalls' finished")
            chromatogram(hetcalls,
                         width = strtoi(
                             ChromatogramParam[["baseNumPerRow"]]),
                         height = 2,
                         trim5 = trimmedRV[["trimmedStartPos"]],
                         trim3 = rawSeqLength -
                             trimmedRV[["trimmedFinishPos"]],
                         showtrim = (ChromatogramParam[["showTrimmed"]]),
                         showcalls = "both")
            removeNotification(id = "chromatogramNotification")
            shinyjs::enable("ChromatogramBasePerRow")
            shinyjs::enable("ChromatogramHeightPerRow")
            shinyjs::enable("ChromatogramSignalRatioCutoff")
            shinyjs::enable("ChromatogramCheckShowTrimmed")
            shinyjs::enable("M1TrimmingCutoffText")
            shinyjs::enable("M2CutoffQualityScoreText")
            shinyjs::enable("M2SlidingWindowSizeText")
            shinyjs::enable("startTrimmingButton")
            shinyjs::enable("saveChromatogramParam")
        }
    })
}
