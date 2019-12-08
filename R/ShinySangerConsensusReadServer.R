### ============================================================================
### R shiny consensusRead server function
### ============================================================================
consensusReadServer <- function(input, output, session) {
    # Suppress Warning
    options(warn = -1)

    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusRead <- getShinyOption("SangerConsensusRead")
    shinyDirectory <- getShinyOption("shinyDirectory")
    SangerConsensus <- SangerConsensusRead[[1]]

    ### ------------------------------------------------------------------------
    ### ConsensusRead-related parameters initialization.
    ### ------------------------------------------------------------------------
    SCTrimmingMethod <- SangerConsensus@forwardReadsList[[1]]@
        QualityReport@TrimmingMethod
    if (SCTrimmingMethod == "M1") {
        SCTrimmingMethodName = "Method 1:
                                'Logarithmic Scale Trimming'"
    } else if (SCTrimmingMethod == "M2") {
        SCTrimmingMethodName = "Method 2:
                                'Logarithmic Scale Sliding Window Trimming'"
    }

    ### ------------------------------------------------------------------------
    ### Reads-related parameters initialization.
    ### ------------------------------------------------------------------------
    forwardReadNum <- length((SangerConsensus)@forwardReadsList)
    reverseReadNum <- length((SangerConsensus)@reverseReadsList)
    SangerSingleReadNum <- forwardReadNum + reverseReadNum
    # readFeature
    forwardReadFeature <- sapply(1:forwardReadNum, function(i)
        paste0(i, " ",
               SangerConsensus@forwardReadsList[[i]]@readFeature))
    reverseReadFeature <- sapply(1:reverseReadNum, function(i)
        paste0(i, " ",
               SangerConsensus@reverseReadsList[[i]]@readFeature))
    # readFileName (basename) (Fixed)
    forwardReadBFN <- sapply(1:forwardReadNum, function(i)
        basename(SangerConsensus@forwardReadsList[[i]]@readFileName))
    reverseReadBFN <- sapply(1:reverseReadNum, function(i)
        basename(SangerConsensus@reverseReadsList[[i]]@readFileName))
    SangerSingleReadBFN <- c(forwardReadBFN, reverseReadBFN)

    ### ------------------------------------------------------------------------
    ### ConsensusRead reactiveValue
    ### ------------------------------------------------------------------------
    consensusParam <-
        reactiveValues(
            consensusRead   = SangerConsensus@consensusRead,
            differencesDF   = SangerConsensus@differencesDF,
            alignment       = as.character(SangerConsensus@alignment),
            distanceMatrix  = SangerConsensus@distanceMatrix,
            dendrogram      = SangerConsensus@dendrogram,
            indelsDF        = SangerConsensus@indelsDF,
            stopCodonsDF    = SangerConsensus@stopCodonsDF,
            secondaryPeakDF = SangerConsensus@secondaryPeakDF)

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
    ### Functions for all UI page
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
    output$consensusRead_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (input$sidebar_menu == "Sanger Consensus Read Overview") {
            message(">>>>>>>> Inside '", input$sidebar_menu, "'")
            ### ------------------------------------------------------------
            ### First assign the ChromatogramParam parameter
            ### ------------------------------------------------------------
            ### ----------------------------------------------------------------
            ### Dynamic page navigation: consensusRead content overview
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
                                          "Re-calculate consensus read",
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
                                      h4(SangerConsensus@parentDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Consenesus Read Name: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerConsensus@consensusReadName),
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
                                      h4(SCTrimmingMethodName),
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
                                      h4(SangerConsensus@suffixForwardRegExp),
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
                                      h4(SangerConsensus@suffixReverseRegExp),
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
                    ################################################
                    #### Add this after having reference sample ####
                    ################################################
                    # If it is null
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Consensus Read Parameters",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(3,
                               uiOutput("SCMinReadsNum") ,
                        ),
                        column(3,
                               uiOutput("SCMinReadLength")  ,
                        ),
                        column(3,
                               uiOutput("SCMinFractionCall") ,
                        ),
                        column(3,
                               uiOutput("SCMaxFractionLost") ,
                        ),
                        column(3,
                               uiOutput("SCAcceptStopCodons") ,
                        ),
                        column(3,
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
                                           "Consensus Read Results: "),
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
                               htmlOutput("consensusAlignmentHTML"),
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
        } else if (!is.na(strtoi(singleReadIndex)) &&
                   (directionParam == "Forward" ||
                    directionParam == "Reverse")) {
            message(">>>>>>>> Inside '", input$sidebar_menu, "'")
            if (directionParam == "Forward") {
                sequenceParam[["primarySeq"]] <<-
                    as.character(
                        SangerConsensus@
                            forwardReadsList[[singleReadIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(
                        SangerConsensus@
                            forwardReadsList[[singleReadIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(
                        SangerConsensus@
                            forwardReadsList[[singleReadIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(
                        SangerConsensus@
                            forwardReadsList[[singleReadIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(
                        SangerConsensus@
                            forwardReadsList[[singleReadIndex]]@primaryAASeqS3)


                ChromatogramParam[["baseNumPerRow"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    ChromatogramParam@baseNumPerRow
                ChromatogramParam[["heightPerRow"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    ChromatogramParam@heightPerRow
                ChromatogramParam[["signalRatioCutoff"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    ChromatogramParam@signalRatioCutoff
                ChromatogramParam[["showTrimmed"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    ChromatogramParam@showTrimmed

                trimmedParam[["M1TrimmingCutoff"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@M1TrimmingCutoff
                trimmedParam[["M2CutoffQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@M2CutoffQualityScore
                trimmedParam[["M2SlidingWindowSize"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@M2SlidingWindowSize

                trimmedRV[["rawSeqLength"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerConsensus@forwardReadsList[[singleReadIndex]]@
                              QualityReport@remainingRatio * 100, 2)
                SSReadBFN <- basename(
                    SangerConsensus@
                        forwardReadsList[[singleReadIndex]]@readFileName)
                SSReadAFN <- SangerConsensus@
                    forwardReadsList[[singleReadIndex]]@readFileName
            } else if (directionParam == "Reverse") {
                sequenceParam[["primarySeq"]] <<-
                    as.character(
                        SangerConsensus@
                            reverseReadsList[[singleReadIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(
                        SangerConsensus@
                            reverseReadsList[[singleReadIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(
                        SangerConsensus@
                            reverseReadsList[[singleReadIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(
                        SangerConsensus@
                            reverseReadsList[[singleReadIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(
                        SangerConsensus@
                            reverseReadsList[[singleReadIndex]]@primaryAASeqS3)

                ChromatogramParam[["baseNumPerRow"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    ChromatogramParam@baseNumPerRow
                ChromatogramParam[["heightPerRow"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    ChromatogramParam@heightPerRow
                ChromatogramParam[["signalRatioCutoff"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    ChromatogramParam@signalRatioCutoff
                ChromatogramParam[["showTrimmed"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    ChromatogramParam@showTrimmed

                trimmedParam[["M1TrimmingCutoff"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    M1TrimmingCutoff
                trimmedParam[["M2CutoffQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    M2CutoffQualityScore
                trimmedParam[["M2SlidingWindowSize"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    M2SlidingWindowSize

                trimmedRV[["rawSeqLength"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerConsensus@reverseReadsList[[singleReadIndex]]@
                              QualityReport@remainingRatio * 100, 2)
                SSReadBFN <- basename(
                    SangerConsensus@
                        reverseReadsList[[singleReadIndex]]@readFileName)
                SSReadAFN <- SangerConsensus@
                    reverseReadsList[[singleReadIndex]]@readFileName
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
                            box(title = tags$p("Cumulative Ratio Plot",
                                               style = "font-size: 21px;
                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                plotlyOutput("qualityTrimmingRatioPlot") %>%
                                    withSpinner()),
                            box(title = tags$p("Quality BP Plot",
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
                                   "Whether show trimmed region",
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
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic rightHeader text
    ### ------------------------------------------------------------------------
    #!!!!!!!! Fix
    observeEventDynamicHeaderSC(input, output, session, trimmedRV)
    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        message("@@@@@@@ 'close button' has been clicked")
        btn <- input$closeUI
        stopApp()
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        shinyjs::disable("saveS4")
        message("@@@@@@@ 'save button' has been clicked")
        newS4Object <- file.path(shinyDirectory, "SangerConsensus.Rda")
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
        ### Save SangerConsensus quality S4 object
        ### --------------------------------------------------------------------
        saveRDS(SangerConsensus, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        NEW_SANGER_CONSENSUS_READ <<- readRDS(file=newS4Object)
        shinyjs::enable("saveS4")
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus read re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButton, {
        shinyjs::disable("recalculateButton")
        message("@@@@@@@ 'Reactive button' has been clicked")
        message("######## Start recalculating consensus read (SC")
        CSResult<-
            calculateConsensusRead (SangerConsensus@forwardReadsList,
                                    SangerConsensus@reverseReadsList,
                                    SangerConsensus@refAminoAcidSeq,
                                    SangerConsensus@minFractionCall,
                                    SangerConsensus@maxFractionLost,
                                    SangerConsensus@geneticCode,
                                    SangerConsensus@acceptStopCodons,
                                    SangerConsensus@readingFrame)

        SangerConsensus@consensusRead <<- CSResult$consensusGapfree
        SangerConsensus@differencesDF <<- CSResult$diffsDf
        SangerConsensus@alignment <<- CSResult$aln2
        SangerConsensus@distanceMatrix <<- CSResult$dist
        SangerConsensus@dendrogram <<- CSResult$dend
        SangerConsensus@indelsDF <<- CSResult$indels
        SangerConsensus@stopCodonsDF <<- CSResult$stopsDf
        SangerConsensus@secondaryPeakDF <<- CSResult$spDf

        consensusParam[["consensusRead"]] <<- SangerConsensus@consensusRead
        consensusParam[["differencesDF"]] <<- SangerConsensus@differencesDF
        consensusParam[["alignment"]] <<- as.character(SangerConsensus@alignment)
        consensusParam[["distanceMatrix"]] <<-SangerConsensus@distanceMatrix
        consensusParam[["dendrogram"]] <<- SangerConsensus@dendrogram
        consensusParam[["indelsDF"]] <<- SangerConsensus@indelsDF
        consensusParam[["stopCodonsDF"]] <<- SangerConsensus@stopCodonsDF
        consensusParam[["secondaryPeakDF"]] <<- SangerConsensus@secondaryPeakDF
        shinyjs::enable("recalculateButton")
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$saveChromatogramParam, {
        ## !!!!! Update !!!!
        message("@@@@@@@ 'Reactive button' has been clicked")
        message("######## Start recalculating chromatogram")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        ### --------------------------------------------------------------------
        ### Update ChromatogramBasePerRow
        ### --------------------------------------------------------------------
        if (!is.na(strtoi(singleReadIndex))) {
            if (directionParam == "Forward") {
                SangerConsensus@forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerConsensus@forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerConsensus@forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerConsensus@forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    showTrimmed <<- input$ChromatogramCheckShowTrimmed
            } else if (directionParam == "Reverse") {
                SangerConsensus@reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerConsensus@reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerConsensus@reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerConsensus@reverseReadsList[[singleReadIndex]]@ChromatogramParam@
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
    ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$startTrimmingButton, {
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            if (SangerConsensus@forwardReadsList[[1]]@
                QualityReport@TrimmingMethod == "M1") {
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
                            SangerConsensus@forwardReadsList[[singleReadIndex]]@
                                QualityReport@qualityPhredScores,
                            SangerConsensus@forwardReadsList[[singleReadIndex]]@
                                QualityReport@qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                        QualityReport@M1TrimmingCutoff <<-
                        as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                        QualityReport@M1TrimmingCutoff
                } else if (directionParam == "Reverse") {
                    ### --------------------------------------------------------
                    ### Start M1 trimming calculation
                    ### --------------------------------------------------------
                    trimmingPos <- M1inside_calculate_trimming(
                            SangerConsensus@reverseReadsList[[singleReadIndex]]@
                                QualityReport@qualityPhredScores,
                            SangerConsensus@reverseReadsList[[singleReadIndex]]@
                                QualityReport@qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                        QualityReport@M1TrimmingCutoff <<-
                        as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
                        QualityReport@M1TrimmingCutoff
                }
            } else if (SangerConsensus@forwardReadsList[[1]]@
                       QualityReport@TrimmingMethod == "M2") {
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
                    strtoi(input$M2SlidingWindowSizeText) <= 20 &&
                    strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
                    inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
                } else {
                    inputM2SlidingWindowSizeText <- 5
                }
                if (directionParam == "Forward") {
                    ### ------------------------------------------------------------
                    ### Start M2 trimming calculation
                    ### ------------------------------------------------------------
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerConsensus@forwardReadsList[[singleReadIndex]]@
                                QualityReport@qualityPhredScores,
                            SangerConsensus@forwardReadsList[[singleReadIndex]]@
                                QualityReport@qualityBaseScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText))
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                        QualityReport@M2CutoffQualityScore <<-
                        strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                        QualityReport@M2CutoffQualityScore
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                        QualityReport@M2SlidingWindowSize <<-
                        strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                        QualityReport@M2SlidingWindowSize
                } else if (directionParam == "Reverse") {
                    ### ------------------------------------------------------------
                    ### Start M2 trimming calculation
                    ### ------------------------------------------------------------
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerConsensus@reverseReadsList[[singleReadIndex]]@
                                QualityReport@qualityPhredScores,
                            SangerConsensus@reverseReadsList[[singleReadIndex]]@
                                QualityReport@qualityBaseScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText))
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                        QualityReport@M2CutoffQualityScore <<-
                        strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
                        QualityReport@M2CutoffQualityScore
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                        QualityReport@M2SlidingWindowSize <<-
                        strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
                        QualityReport@M2SlidingWindowSize
                }
            }
            if (directionParam == "Forward") {
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]
                trimmedRV[["rawSeqLength"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerConsensus@forwardReadsList[[singleReadIndex]]@
                              QualityReport@remainingRatio * 100, 2)
            } else if (directionParam == "Reverse") {
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]
                trimmedRV[["rawSeqLength"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerConsensus@reverseReadsList[[singleReadIndex]]@
                              QualityReport@remainingRatio * 100, 2)
            }
        }
    })

    ############################################################################
    ### ConsensusRead (Function for Sanger Consensus Read Overview)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Valuebox for basic information
    ### ------------------------------------------------------------------------
    valueBoxSCMinReadsNum(input, output,
                          SangerConsensus@minReadsNum, session)
    valueBoxSCMinReadLength(input, output,
                            SangerConsensus@minReadLength, session)
    valueBoxSCMinFractionCall(input, output,
                              SangerConsensus@minFractionCall, session)
    valueBoxSCMaxFractionLost(input, output,
                              SangerConsensus@maxFractionLost, session)
    valueBoxSCAcceptStopCodons(input, output,
                               SangerConsensus@acceptStopCodons, session)
    valueBoxSCReadingFrame(input, output,
                           SangerConsensus@readingFrame, session)
    ### ------------------------------------------------------------------------
    ### geneticCodeDF
    ### ------------------------------------------------------------------------
    output$geneticCodeDF <- renderExcel({
        suppressMessages(
            excelTable(data =  t(data.frame(SangerConsensus@geneticCode)),
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
        if (SangerConsensus@refAminoAcidSeq == "") {
            box(title = tags$p("Reference Amino Acids Sequence",
                               style = "font-size: 24px;
                                    font-weight: bold;"),
                collapsible = TRUE,
                status = "success", width = 12,
                column(width = 1),
                column(width = 11,
                       h4("Reference Amino Acid Sequence is not provided."))
            )
        } else {
            box(title = tags$p("Reference Amino Acids Sequence",
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
        refAminoAcidSeqVec <- strsplit(SangerConsensus@refAminoAcidSeq, "")[[1]]
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
    output$consensusAlignmentHTML <- renderUI({
        browseSeqHTML <-
            file.path(shinyDirectory,
                      paste0(SangerConsensus@consensusReadName,
                             "_Alignment_BrowseSeqs.html"))
        BrowseSeqs(DNAStringSet(consensusParam[["alignment"]]),
                   openURL=FALSE, htmlFile=browseSeqHTML)
        includeHTML(
            file.path(shinyDirectory,
                      paste0(SangerConsensus@consensusReadName,
                             "_Alignment_BrowseSeqs.html")))
    })
    ### ------------------------------------------------------------------------
    ### Difference
    ### ------------------------------------------------------------------------
    output$SCDifferencesDFUI <- renderUI({
        if (all(dim(consensusParam[["differencesDF"]]) == c(0,0))) {
            h4("*** 'Differences' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDifferencesDF")
        }
    })
    output$SCDifferencesDF = renderDataTable({
        consensusParam[["differencesDF"]]
    })
    ### ------------------------------------------------------------------------
    ### dendrogram
    ### ------------------------------------------------------------------------
    output$dendrogramPlot <- renderPlot({
        # plot(consensusParam[["dendrogram"]][[2]])
        ggdendrogram(consensusParam[["dendrogram"]][[2]], rotate = TRUE)
    })
    output$dendrogramDF <- renderDataTable({
        consensusParam[["dendrogram"]][[1]]
    })
    ### ------------------------------------------------------------------------
    ### Distance
    ### ------------------------------------------------------------------------
    output$SCDistanceMatrixPlotUI <- renderUI({
        if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
            h4("*** 'Distance' dataframe is empty. (Cannot plot)***",
               style="font-weight: bold; font-style: italic;")
        } else {
            plotlyOutput("SCDistanceMatrixPlot")
        }
    })
    output$SCDistanceMatrixPlot <- renderPlotly({
        plot_ly(x = SangerSingleReadBFN,
                y = SangerSingleReadBFN,
                z = consensusParam[["distanceMatrix"]],
                colors = colorRamp(c("white", "#32a852")),
                type = "heatmap")
    })
    output$SCDistanceMatrixUI <- renderUI({
        if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
            h4("*** 'Distance' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDistanceMatrix")
        }
    })
    output$SCDistanceMatrix = renderDataTable({
        consensusParam[["distanceMatrix"]]
    })
    ### ------------------------------------------------------------------------
    ### Indels
    ### ------------------------------------------------------------------------
    output$SCIndelsDFUI <- renderUI({
        if (all(dim(consensusParam[["indelsDF"]]) == c(0,0))) {
            h4("*** 'Indels' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCIndelsDF")

        }
    })
    output$SCIndelsDF <- renderDataTable({
        consensusParam[["indelsDF"]]
    })
    ### ------------------------------------------------------------------------
    ### StopCodons
    ### ------------------------------------------------------------------------
    output$SCStopCodonsDFUI <- renderUI({
        if (all(dim(consensusParam[["stopCodonsDF"]]) == c(0,0))) {
            h4("*** 'Stop Codons' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCStopCodonsDF")
        }
    })
    output$SCStopCodonsDF <- renderDataTable({
        consensusParam[["stopCodonsDF"]]
    })
    ############################################################################
    ### SangerSingleRead (Function for singel read in consensusRead)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Primary dataframe
    ### ------------------------------------------------------------------------
    output$primarySeqDF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(singleReadIndex)) {
            primarySeqDisplay (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Secondary dataframe
    ### ------------------------------------------------------------------------
    output$secondSeqDF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(singleReadIndex)) {
            secondarySeqDisplay (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Quality Score dataframe
    ### ------------------------------------------------------------------------
    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            if (directionParam == "Forward") {
                PhredScore <- SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                PhredScore <- SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityScoreDisplay (PhredScore)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acids dataframe (1)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS1DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            PrimAASeqS1Display (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acids dataframe (2)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS2DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            PrimAASeqS2Display (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acids dataframe (3)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS3DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            PrimAASeqS3Display (sequenceParam)
        }
    })
    ### ------------------------------------------------------------------------
    ### Primary Trimmed dataframe
    ### ------------------------------------------------------------------------
    output$primarySeqTrimmedDF <- renderExcel({
        ## !!!!! Update !!!!
        message(">>>>>>>>>>>> Update primarySeqTrimmedDF! ")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(singleReadIndex)) {
            primarySeqTrimmedDisplay (input, output, session,
                                      sequenceParam, trimmedRV)
        }
    })
    ### ------------------------------------------------------------------------
    ### Secondary Trimmed dataframe
    ### ------------------------------------------------------------------------
    output$secondSeqTrimmedDF <- renderExcel({
        ## !!!!! Update !!!!
        message(">>>>>>>>>>>> Update secondSeqTrimmedDF! ")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(singleReadIndex)) {
            secondSeqTrimmedDisplay (input, output, session,
                                     sequenceParam, trimmedRV)
        }
    })
    ### ------------------------------------------------------------------------
    ### Quality Score Trimmed dataframe
    ### ------------------------------------------------------------------------
    output$qualityScoreTrimmedDF <- renderExcel({
        message(">>>>>>>>>>>> Update qualityScoreTrimmedDF! ")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(singleReadIndex)) {
            if (directionParam == "Forward") {
                PhredScore <- SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores[
                        (trimmedRV[["trimmedStartPos"]]+1):
                            trimmedRV[["trimmedFinishPos"]]]
            } else if (directionParam == "Reverse") {
                PhredScore <- SangerConsensus@reverseReadsList[[singleReadIndex]]@
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
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            if (SangerConsensus@forwardReadsList[[1]]@
                QualityReport@TrimmingMethod == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Logarithmic Scale Trimming'")
            } else if (SangerConsensus@forwardReadsList[[1]]@
                       QualityReport@TrimmingMethod == "M2") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Logarithmic Scale Sliding Window Trimming'")
            }
        }
    })
    output$TrimmingMethodUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            ## For method, everyone is same, so just pick forward one.
            if (SangerConsensus@forwardReadsList[[1]]@
                QualityReport@TrimmingMethod == "M1") {
                if (directionParam == "Forward") {
                    if (is.null(SangerConsensus@forwardReadsList[[singleReadIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                            QualityReport@M1TrimmingCutoff <<-  0.0001
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerConsensus@reverseReadsList[[singleReadIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
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
            } else if (SangerConsensus@forwardReadsList[[1]]@
                       QualityReport@TrimmingMethod == "M2") {
                if (directionParam == "Forward") {
                    if (is.null(SangerConsensus@forwardReadsList[[singleReadIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerConsensus@forwardReadsList[[singleReadIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  5
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerConsensus@reverseReadsList[[singleReadIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerConsensus@reverseReadsList[[singleReadIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  5
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

    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]

        if (!is.na(strtoi(singleReadIndex))) {
            if (directionParam == "Forward") {
                qualityPhredScores <- SangerConsensus@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                qualityPhredScores <- SangerConsensus@
                    reverseReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityTrimmingRatioPlotDisplay(input, output,session,
                                            trimmedRV, qualityPhredScores)
        }
    })
    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            if (directionParam == "Forward") {
                qualityPhredScores <-
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                qualityPhredScores <-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
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
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        trimmedRV[["trimmedSeqLength"]]
        if (!is.na(strtoi(singleReadIndex))) {
            if (directionParam == "Forward") {
                chromatogramRowNumAns <-
                    chromatogramRowNum (
                        strtoi(ChromatogramParam[["baseNumPerRow"]]),
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                            QualityReport@rawSeqLength,
                        SangerConsensus@forwardReadsList[[singleReadIndex]]@
                            QualityReport@trimmedSeqLength,
                        ChromatogramParam[["showTrimmed"]]) *
                    strtoi(ChromatogramParam[["heightPerRow"]])

            } else if (directionParam == "Reverse") {
                chromatogramRowNumAns <-
                    chromatogramRowNum (
                        strtoi(ChromatogramParam[["baseNumPerRow"]]),
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
                            QualityReport@rawSeqLength,
                        SangerConsensus@reverseReadsList[[singleReadIndex]]@
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
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
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
                    SangerConsensus@forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                message(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (forward)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerConsensus@forwardReadsList[[singleReadIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ----------------------------------------------------------------
                ### Update 'SangerConsensus'!
                ### ----------------------------------------------------------------
                SangerConsensus@forwardReadsList[[singleReadIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerConsensus@forwardReadsList[[singleReadIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerConsensus@forwardReadsList[[singleReadIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerConsensus@forwardReadsList[[singleReadIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ----------------------------------------------------------------
                ### Updating AASeqs
                ### ----------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerConsensus@geneticCode)
                SangerConsensus@forwardReadsList[[singleReadIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerConsensus@forwardReadsList[[singleReadIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ----------------------------------------------------------------
                ### Updating reactive values
                ### ----------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerConsensus@forwardReadsList[[
                        singleReadIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerConsensus@forwardReadsList[[
                        singleReadIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerConsensus@forwardReadsList[[
                        singleReadIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerConsensus@forwardReadsList[[
                        singleReadIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerConsensus@forwardReadsList[[
                        singleReadIndex]]@primaryAASeqS3)
            } else if (directionParam == "Reverse") {
                rawSeqLength <-
                    SangerConsensus@reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                message(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (reverse)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerConsensus@reverseReadsList[[singleReadIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ----------------------------------------------------------------
                ### Update 'SangerConsensus'!
                ### ----------------------------------------------------------------
                SangerConsensus@reverseReadsList[[singleReadIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerConsensus@reverseReadsList[[singleReadIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerConsensus@reverseReadsList[[singleReadIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerConsensus@reverseReadsList[[singleReadIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ----------------------------------------------------------------
                ### Updating AASeqs
                ### ----------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerConsensus@geneticCode)
                SangerConsensus@reverseReadsList[[singleReadIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerConsensus@reverseReadsList[[singleReadIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ----------------------------------------------------------------
                ### Updating reactive values
                ### ----------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerConsensus@reverseReadsList[[
                        singleReadIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerConsensus@reverseReadsList[[
                        singleReadIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerConsensus@reverseReadsList[[
                        singleReadIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerConsensus@reverseReadsList[[
                        singleReadIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerConsensus@reverseReadsList[[
                        singleReadIndex]]@primaryAASeqS3)
            }

            # message(">>>>>>>>>>>> 'MakeBaseCalls' finished")
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
