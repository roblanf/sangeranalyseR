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

    # readFileName (basename) (Fixed)
    forwardReadBFN <- sapply(1:forwardReadNum, function(i)
        basename(SangerConsensus@forwardReadsList[[i]]@readFileName))
    reverseReadBFN <- sapply(1:reverseReadNum, function(i)
        basename(SangerConsensus@reverseReadsList[[i]]@readFileName))
    SangerSingleReadBFN <- c(forwardReadBFN, reverseReadBFN)




    # readFileName (absolute)
    forwardReadAFN <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@readFileName)
    reverseReadAFN <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@readFileName)
    SangerSingleReadAFN <- c(forwardReadAFN, reverseReadAFN)

    # readFeature
    forwardReadFeature <- sapply(1:forwardReadNum, function(i)
        paste0(i, " ",
               SangerConsensus@forwardReadsList[[i]]@readFeature))
    reverseReadFeature <- sapply(1:reverseReadNum, function(i)
        paste0(i, " ",
               SangerConsensus@reverseReadsList[[i]]@readFeature))
    SangerSingleReadFeature <- c(forwardReadFeature, reverseReadFeature)

    # QualityReport
    forwardReadQualReport <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@QualityReport)
    reverseReadQualReport <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@QualityReport)
    SangerSingleReadQualReport <- c(forwardReadQualReport,
                                    reverseReadQualReport)

    # ChromatogramParam
    forwardReadChromatogramParam <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@ChromatogramParam)
    reverseReadChromatogramParam <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@ChromatogramParam)
    SangerSingleReadChromatogramParam <- c(forwardReadChromatogramParam,
                                           reverseReadChromatogramParam)

    # primaryAASeqS1DF
    forwardReadPrimAASeqS1DF <- lapply(1:forwardReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@forwardReadsList[[i]]@primaryAASeqS1)
        AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    reverseReadPrimAASeqS1DF <- lapply(1:reverseReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@reverseReadsList[[i]]@primaryAASeqS1)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    SangerSingleReadPrimAASeqS1DF <- c(forwardReadPrimAASeqS1DF,
                                       reverseReadPrimAASeqS1DF)

    # primaryAASeqS2DF
    forwardReadPrimAASeqS2DF <- lapply(1:forwardReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@forwardReadsList[[i]]@primaryAASeqS2)
        AAString <- rbind(NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    reverseReadPrimAASeqS2DF <- lapply(1:reverseReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@reverseReadsList[[i]]@primaryAASeqS2)
        AAString <- rbind(NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    SangerSingleReadPrimAASeqS2DF <- c(forwardReadPrimAASeqS2DF,
                                       reverseReadPrimAASeqS2DF)

    # primaryAASeqS3DF
    forwardReadPrimAASeqS3DF <- lapply(1:forwardReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@forwardReadsList[[i]]@primaryAASeqS3)
        AAString <- rbind(NA, NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    reverseReadPrimAASeqS3DF <- lapply(1:reverseReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@reverseReadsList[[i]]@primaryAASeqS3)
        AAString <- rbind(NA, NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    SangerSingleReadPrimAASeqS3DF <- c(forwardReadPrimAASeqS3DF,
                                       reverseReadPrimAASeqS3DF)

    ### ------------------------------------------------------------------------
    ### ConsensusRead reactiveValue
    ### ------------------------------------------------------------------------
    consensusParam <-
        reactiveValues(
            consensusRead   = SangerConsensus@consensusRead,
            differencesDF   = SangerConsensus@differencesDF,
            alignment       = SangerConsensus@alignment,
            distanceMatrix  = SangerConsensus@distanceMatrix,
            dendrogram      = SangerConsensus@dendrogram,
            indelsDF        = SangerConsensus@indelsDF,
            stopCodonsDF    = SangerConsensus@stopCodonsDF,
            secondaryPeakDF = SangerConsensus@secondaryPeakDF)

    ### ------------------------------------------------------------------------
    ### SingleRead reactiveValue
    ### ------------------------------------------------------------------------
    sequenceParam <- reactiveValues(primarySeq = 0,
                                    secondarySeq = 0,
                                    primaryAASeqS1 = 0,
                                    primaryAASeqS2 = 0,
                                    primaryAASeqS3 = 0)

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
            h1(input$sidebar_menu)
            if (directionParam == "Forward") {
                SSReadBFN <- basename(
                    SangerConsensus@
                        forwardReadsList[[singleReadIndex]]@readFileName)
                SSReadAFN <- SangerConsensus@
                    forwardReadsList[[singleReadIndex]]@readFileName
            } else if (directionParam == "Reverse") {
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
                                       "Primary, Secondary DNA Sequences &
                                       Amino Acid Sequence (Before Trimming):"),
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
                                   style = "simple", color = "success",
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
                            box(title = tags$p("Cumulative Ratio Plot",
                                               style = "font-size: 21px;
                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                plotlyOutput("qualityQualityBasePlot") %>%
                                    withSpinner()),
                        ),
                    ),
                ),
                # box(title = tags$p(tagList(icon("dot-circle"),
                #                            "Chromatogram: "),
                #                    style = "font-size: 26px;
                #                            font-weight: bold;"),
                #     solidHeader = TRUE, collapsible = TRUE,
                #     status = "success", width = 12,
                #     tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                #
                #     box(title = tags$p(tagList(icon("arrow-circle-right"),
                #                                "Chromatogram Input"),
                #                        style = "font-size: 24px;
                #                            font-weight: bold;"),
                #         collapsible = TRUE,
                #         status = "success", width = 12,
                #         column(3,
                #                sliderInput("ChromatogramBasePerRow",
                #                            label =h4("Base Number Per Row"),
                #                            min = 5,
                #                            max = 200,
                #                            value = ChromatogramParam[["baseNumPerRow"]]),
                #                sliderInput("ChromatogramHeightPerRow",
                #                            label = h4("Height Per Row"),
                #                            min = 50,
                #                            max = 600,
                #                            value = ChromatogramParam[["heightPerRow"]]),
                #         ),
                #         column(3,
                #                tags$hr(
                #                    style =
                #                        ("border-top: 4px hidden #A9A9A9;")),
                #                numericInput(
                #                    "ChromatogramSignalRatioCutoff",
                #                    h3("Signal Ratio Cutoff"),
                #                    value = ChromatogramParam[["signalRatioCutoff"]]),
                #                checkboxInput(
                #                    "ChromatogramCheckShowTrimmed",
                #                    "Whether show trimmed region",
                #                    value =
                #                        ChromatogramParam[["showTrimmed"]])
                #         ),
                #         column(3,
                #                tags$hr(
                #                    style=("border-top: 4px hidden #A9A9A9;")),
                #                uiOutput("ChromatogramtrimmedStartPos"),
                #         ),
                #         column(3,
                #                tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                #                uiOutput("ChromatogramtrimmedFinishPos"),
                #         ),
                #         actionBttn("saveChromatogramParam",
                #                    "Apply Chromatogram Parameters",
                #                    style = "simple", color = "success",
                #                    block = TRUE, size = "lg")
                #     ),
                #     box(title = tags$p(tagList(icon("arrow-circle-left"),
                #                                "Chromatogram Output"),
                #                        style = "font-size: 24px;
                #                            font-weight: bold;"),
                #         collapsible = TRUE,
                #         status = "success", width = 12,
                #         column(width = 12,
                #                uiOutput("chromatogramUIOutput"),
                #         )
                #     ),
                # )
            )
        } else if (!is.na(strtoi(singleReadIndex)) &&
                   directionParam == "Reverse") {
            message(">>>>>>>> Inside '", input$sidebar_menu, "'")
            h1(input$sidebar_menu)
        }
    })


    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic rightHeader text
    ### ------------------------------------------------------------------------
    #!!!!!!!! Fix
    observeEventDynamicHeaderSC(input, output, session, trimmedRV,
                              SangerSingleReadQualReport)




    ############################################################################
    ### ConsensusRead (Function for Sanger Consensus Read Overview)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus read re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButton, {
        message("######## Reactive button clicked !!!")
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
        consensusParam[["alignment"]] <<- SangerConsensus@alignment
        consensusParam[["distanceMatrix"]] <<-SangerConsensus@distanceMatrix
        consensusParam[["dendrogram"]] <<- SangerConsensus@dendrogram
        consensusParam[["indelsDF"]] <<- SangerConsensus@indelsDF
        consensusParam[["stopCodonsDF"]] <<- SangerConsensus@stopCodonsDF
        consensusParam[["secondaryPeakDF"]] <<- SangerConsensus@secondaryPeakDF
        message("######## Finish recalculation")
    })
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
                       defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
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
                       defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
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
        BrowseSeqs(consensusParam[["alignment"]],
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
            if (directionParam == "Forward") {
                primarySeq <-unlist(strsplit(toString(
                    SangerConsensus@
                        forwardReadsList[[singleReadIndex]]@primarySeq), ""))

            } else if (directionParam == "Reverse") {
                primarySeq <- unlist(strsplit(toString(
                    SangerConsensus@
                        reverseReadsList[[singleReadIndex]]@primarySeq), ""))
            }
            primarySeqDF <- data.frame(
                t(data.frame(primarySeq)), stringsAsFactors = FALSE)
            colnames(primarySeqDF) <- substr(colnames(primarySeqDF), 2, 100)
            rownames(primarySeqDF) <- NULL
            AstyleList <- SetCharStyleList(primarySeqDF, "A", "#1eff00")
            TstyleList <- SetCharStyleList(primarySeqDF, "T", "#ff7a7a")
            CstyleList <- SetCharStyleList(primarySeqDF, "C", "#7ac3ff")
            GstyleList <- SetCharStyleList(primarySeqDF, "G", "#c9c9c9")
            styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
            suppressMessages(
                excelTable(data = primarySeqDF, defaultColWidth = 30,
                           editable = TRUE, rowResize = FALSE,
                           columnResize = FALSE, allowInsertRow = FALSE,
                           allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                           allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                           style = styleList, loadingSpin = TRUE)
            )
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
            if (directionParam == "Forward") {
                secondarySeq <-unlist(strsplit(toString(
                    SangerConsensus@
                        forwardReadsList[[singleReadIndex]]@secondarySeq), ""))
            } else if (directionParam == "Reverse") {
                secondarySeq <-unlist(strsplit(toString(
                    SangerConsensus@
                        reverseReadsList[[singleReadIndex]]@secondarySeq), ""))
            }
            secondarySeqDF <- data.frame(
                t(data.frame(secondarySeq)), stringsAsFactors = FALSE)
            colnames(secondarySeqDF) <- substr(colnames(secondarySeqDF), 2, 100)
            rownames(secondarySeqDF) <- NULL
            AstyleList <- SetCharStyleList(secondarySeqDF, "A", "#1eff00")
            TstyleList <- SetCharStyleList(secondarySeqDF, "T", "#ff7a7a")
            CstyleList <- SetCharStyleList(secondarySeqDF, "C", "#7ac3ff")
            GstyleList <- SetCharStyleList(secondarySeqDF, "G", "#c9c9c9")
            styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
            suppressMessages(
                excelTable(data = secondarySeqDF, defaultColWidth = 30,
                           editable = TRUE, rowResize = FALSE,
                           columnResize = FALSE, allowInsertRow = FALSE,
                           allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                           allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                           style = styleList, loadingSpin = TRUE)
            )
        }
    })
    ### ------------------------------------------------------------------------
    ### Quality Score dataframe
    ### ------------------------------------------------------------------------
    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex)) &&
            directionParam == "Forward") {
            PhredScore <- SangerConsensus@forwardReadsList[[singleReadIndex]]@
                QualityReport@qualityPhredScores
        } else if (!is.na(strtoi(singleReadIndex)) &&
                   directionParam == "Reverse") {
            PhredScore <- SangerConsensus@reverseReadsList[[singleReadIndex]]@
                QualityReport@qualityPhredScores
        }
        PhredScoreDF <- data.frame(
            t(data.frame(PhredScore)), stringsAsFactors = FALSE)
        colnames(PhredScoreDF) <- substr(colnames(PhredScoreDF), 2, 100)
        rownames(PhredScoreDF) <- NULL
        styleList <- SetAllStyleList(PhredScoreDF, "#ecffd9")
        suppressMessages(
            excelTable(data =
                           PhredScoreDF, defaultColWidth = 30, editable = TRUE,
                       rowResize = FALSE, columnResize = FALSE,
                       allowInsertRow = FALSE, allowInsertColumn = FALSE,
                       allowDeleteRow = FALSE, allowDeleteColumn = FALSE,
                       style = styleList, allowRenameColumn = FALSE,
                       loadingSpin = TRUE)
        )
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acids dataframe (1)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS1DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex)) &&
            directionParam == "Forward") {
            AAString <- data.frame(SangerConsensus@
                                       forwardReadsList[[singleReadIndex]]@
                                       primaryAASeqS1)
        } else if (!is.na(strtoi(singleReadIndex)) &&
                   directionParam == "Reverse") {
            AAString <- data.frame(SangerConsensus@
                                       reverseReadsList[[singleReadIndex]]@
                                       primaryAASeqS1)
        }
        AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        width <- rep(90, length(AAStringDF))
        styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
        styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
        styleList <- c(styleList1, styleList2)
        suppressMessages(
            excelTable(data = AAStringDF, columns = data.frame(width = width),
                       defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                       style = styleList, loadingSpin = TRUE)
        )
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acids dataframe (2)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS2DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex)) &&
            directionParam == "Forward") {
            AAString <- data.frame(SangerConsensus@
                                       forwardReadsList[[singleReadIndex]]@
                                       primaryAASeqS2)
        } else if (!is.na(strtoi(singleReadIndex)) &&
                   directionParam == "Reverse") {
            AAString <- data.frame(SangerConsensus@
                                       reverseReadsList[[singleReadIndex]]@
                                       primaryAASeqS2)
        }
        AAString <- rbind(NA, AAString)
        AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        width <- rep(90, length(AAStringDF) - 1)
        width <- c(30, width)
        styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
        styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
        styleList <- c(styleList1, styleList2)
        styleList[['A1']] <- 'background-color: black;'
        suppressMessages(
            excelTable(data = AAStringDF, columns = data.frame(width = width),
                       defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                       style = styleList, loadingSpin = TRUE)
        )
    })
    ### ------------------------------------------------------------------------
    ### Primary Amino Acids dataframe (3)
    ### ------------------------------------------------------------------------
    output$PrimAASeqS3DF <- renderExcel({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex)) &&
            directionParam == "Forward") {
            AAString <- data.frame(SangerConsensus@
                                       forwardReadsList[[singleReadIndex]]@
                                       primaryAASeqS3)
        } else if (!is.na(strtoi(singleReadIndex)) &&
                   directionParam == "Reverse") {
            AAString <- data.frame(SangerConsensus@
                                       reverseReadsList[[singleReadIndex]]@
                                       primaryAASeqS3)
        }
        AAString <- rbind(NA, NA, AAString)
        AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        width <- rep(90, length(AAStringDF) - 2)
        width <- c(30, 30, width)
        styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
        styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
        styleList <- c(styleList1, styleList2)
        styleList[['A1']] <- 'background-color: black;'
        styleList[['B1']] <- 'background-color: black;'
        suppressMessages(
            excelTable(data = AAStringDF, columns = data.frame(width = width),
                       defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                       style = styleList, loadingSpin = TRUE)
        )
    })
    ############################################################################
    ### Trimming method selection functions
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Trimming Method Selection
    ### ------------------------------------------------------------------------
    output$TrimmingMethodSelectionOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (!is.na(strtoi(singleReadIndex))) {
            if (SangerConsensus@forwardReadsList[[singleReadIndex]]@
                QualityReport@TrimmingMethod == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Logarithmic Scale Trimming'")
            } else if (SangerConsensus@forwardReadsList[[singleReadIndex]]@
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
            if (SangerConsensus@forwardReadsList[[singleReadIndex]]@
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
            } else if (SangerConsensus@forwardReadsList[[singleReadIndex]]@
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
}





## !!!!! Update !!!!
# sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
# singleReadIndex <- strtoi(sidebar_menu[[1]])
# directionParam <- sidebar_menu[[2]]
# if (!is.na(strtoi(singleReadIndex)) &&
#     directionParam == "Forward") {
#     SSReadBFN <- basename(SangerConsensus@forwardReadsList[[singleReadIndex]]@readFileName)
# } else if (!is.na(strtoi(singleReadIndex)) &&
#            directionParam == "Reverse") {
# }
