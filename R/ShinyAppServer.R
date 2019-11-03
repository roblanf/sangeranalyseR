### ============================================================================
### R shiny consensus read server function
### ============================================================================
consensusServer <- function(input, output, session) {
    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusRead <- getShinyOption("SangerConsensusReadSet")
    SangerConsensus <- SangerConsensusRead[[1]]


    ### ------------------------------------------------------------------------
    ### ConsensusRead-related parameters initialization.
    ### ------------------------------------------------------------------------
    SCMinReadsNum <- SangerConsensus@minReadsNum
    SCMinReadLength <- SangerConsensus@minReadLength


    SCRefAminoAcidSeq <- SangerConsensus@refAminoAcidSeq


    SCMinFractionCall <- SangerConsensus@minFractionCall
    SCMaxFractionLost <- SangerConsensus@maxFractionLost


    SCGeneticCode <- SangerConsensus@geneticCode


    SCAcceptStopCodons <- SangerConsensus@acceptStopCodons
    SCReadingFrame <- SangerConsensus@readingFrame

    as.character(SangerConsensus@consensusRead)

    # aln = AlignSeqs(c(DNAStringSet(SangerConsensus@consensusRead), DNAStringSet(SangerConsensus@consensusRead)),
    #                 processors = 1, verbose = FALSE)


    SCAlignment<- SangerConsensus@alignment
    SCDifferencesDF<- SangerConsensus@differencesDF
    SCDistanceMatrix <- SangerConsensus@distanceMatrix

    SCDendrogram <- SangerConsensus@dendrogram

    SCIndelsDF <- SangerConsensus@indelsDF
    SCStopCodonsDF <- SangerConsensus@stopCodonsDF

    SCSecondaryPeakDF <- SangerConsensus@secondaryPeakDF




    SangerConsensusForRegExp <- SangerConsensus@consenesusReadName
    SangerConsensusForRegExp <- SangerConsensus@suffixForwardRegExp
    SangerConsensusRevRegExp <- SangerConsensus@suffixReverseRegExp

    ### ------------------------------------------------------------------------
    ### Reads-related parameters initialization.
    ### ------------------------------------------------------------------------
    forwardReadNum <- length((SangerConsensus)@forwardReadsList)
    reverseReadNum <- length((SangerConsensus)@reverseReadsList)
    SangerSingleReadNum <- forwardReadNum + reverseReadNum

    # Forward & reverse reads list
    SangerConsensusFRReadsList <- c(SangerConsensus@forwardReadsList,
                                    SangerConsensus@reverseReadsList)

    # readFileName (basename)
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
        paste0(i, "_",
               SangerConsensus@forwardReadsList[[i]]@readFeature))
    reverseReadFeature <- sapply(1:reverseReadNum, function(i)
        paste0(i+forwardReadNum, "_",
               SangerConsensus@reverseReadsList[[i]]@readFeature))
    SangerSingleReadFeature <- c(forwardReadFeature, reverseReadFeature)

    # abifRawData
    forwardReadAbifRawData <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@abifRawData)
    reverseReadAbifRawData <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@abifRawData)
    SangerSingleReadAbifRawData <- c(forwardReadAbifRawData,
                                     reverseReadAbifRawData)

    # QualityReport
    forwardReadQualReport <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@QualityReport)
    reverseReadQualReport <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@QualityReport)
    SangerSingleReadQualReport <- c(forwardReadQualReport,
                                    reverseReadQualReport)

    # primarySeqID
    forwardReadPrimSeqID <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@primarySeqID)
    reverseReadPrimSeqID <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@primarySeqID)
    SangerSingleReadPrimSeqID <- c(forwardReadPrimSeqID, reverseReadPrimSeqID)
    # primarySeq
    forwardReadPrimSeq <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@primarySeq)
    reverseReadPrimSeq <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@primarySeq)
    SangerSingleReadPrimSeq <- c(forwardReadPrimSeq, reverseReadPrimSeq)

    # secondarySeqID
    forwardReadSecoSeqID <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@secondarySeqID)
    reverseReadSecoSeqID <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@secondarySeqID)
    SangerSingleReadSecoSeqID <- c(forwardReadSecoSeqID, reverseReadSecoSeqID)
    # secondarySeq
    forwardReadSecoSeq <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@secondarySeq)
    reverseReadSecoSeq <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@secondarySeq)
    SangerSingleReadSecoSeq <- c(forwardReadSecoSeq, reverseReadSecoSeq)
    # traceMatrix
    forwardReadTraceMat <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@traceMatrix)
    reverseReadTraceMat <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@traceMatrix)
    SangerSingleReadTraceMat <- c(forwardReadTraceMat, reverseReadTraceMat)
    # peakPosMatrix
    forwardReadReadPeakPosMat <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@peakPosMatrix)
    reverseReadReadPeakPosMat <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@peakPosMatrix)
    SangerSingleReadPeakPosMat <- c(forwardReadReadPeakPosMat,
                                    reverseReadReadPeakPosMat)
    # peakAmpMatrix
    # forwardReadPeakAmpMat <- sapply(1:forwardReadNum, function(i)
    #     SangerConsensus@forwardReadsList[[i]]@peakAmpMatrix)
    # reverseReadPeakAmpMat <- sapply(1:reverseReadNum, function(i)
    #     SangerConsensus@reverseReadsList[[i]]@peakAmpMatrix)
    # SangerSingleReadPeakAmpMat <- c(forwardReadPeakAmpMat, reverseReadPeakAmpMat)
    trimmedRV <- reactiveValues(trimmedStart = 0, trimmedEnd = 0,
                                remainingBP = 0, trimmedRatio = 0)

    ############################################################################
    ### output$ID
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar.
    ### ------------------------------------------------------------------------
    output$singleReadMenu <- renderMenu({
        menu_list <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(SangerSingleReadFeature[i],
                          tabName = SangerSingleReadFeature[i],
                          selected = TRUE, icon = icon("angle-right")))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Overview")})

    ### ------------------------------------------------------------------------
    ### Dynamic page navigation: consensusReadMenu_content
    ### ------------------------------------------------------------------------
    output$consensusReadMenu_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (input$sidebar_menu == "Overview") {
            fluidRow(
                useShinyjs(),
                box(title = tags$p("Parameters: ",
                                   style = "font-size: 26px;
                                       font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                    fluidRow(
                        column(3,
                               uiOutput("SCMinReadsNum"),
                               tags$ul(
                                   textInput("SCMinReadsNumText",
                                             label = p("Change Value"),
                                             value = toString(SCMinReadsNum),
                                             width = '90%')
                               ),
                        ),
                        column(3,
                               uiOutput("SCMinReadLength"),
                               tags$ul(
                                   textInput("SCMinReadLengthText",
                                             label = p("Change Value"),
                                             value = toString(SCMinReadLength),
                                             width = '90%')
                               ),
                        ),
                        column(3,
                               uiOutput("SCMinFractionCall"),
                               tags$ul(
                                   textInput("SCMinFractionCallText",
                                             label = p("Change Value"),
                                             value = toString(SCMinFractionCall),
                                             width = '90%')
                               ),
                        ),
                        column(3,
                               uiOutput("SCMaxFractionLost"),
                               tags$ul(
                                   textInput("SCMaxFractionLostText",
                                             label = p("Change Value"),
                                             value = toString(SCMaxFractionLost),
                                             width = '90%')
                               ),
                        ),
                        column(3,
                               uiOutput("SCAcceptStopCodons"),
                               tags$ul(
                                   textInput("SCAcceptStopCodonsText",
                                             label = p("Change Value"),
                                             value = toString(SCAcceptStopCodons),
                                             width = '90%')
                               ),
                        ),
                        column(3,
                               uiOutput("SCReadingFrame"),
                               tags$ul(
                                   textInput("SCReadingFrameText",
                                             label = p("Change Value"),
                                             value = toString(SCReadingFrame),
                                             width = '90%')
                               ),
                        ),
                    ),
                ),
                box(title = tags$p("Results: ",
                                   style = "font-size: 26px;
                                       font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                    box(title = tags$p("Differences Dataframe",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               dataTableOutput("differencesDF"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                    ),

                    box(title = tags$p("Secondary Peak Dataframe",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               dataTableOutput("secondaryPeakDF"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                    ),
                    # differencesDF             = "data.frame",
                    # distanceMatrix            = "matrix",
                    # dendrogram                = "list",
                    # indelsDF                  = "data.frame",
                    # stopCodonsDF              = "data.frame",
                    # secondaryPeakDF           = "data.frame"
                )
            )
        }
    })

    ### ------------------------------------------------------------------------
    ### Dynamic page navigation: singelReadMenu_content
    ### ------------------------------------------------------------------------
    output$singelReadMenu_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        data(mtcars)
        if (input$sidebar_menu != "Overview") {
            if (!is.na(as.numeric(sidebar_menu[[1]]))) {
                fluidRow(
                    useShinyjs(),
                    # box(sidebar_menu[[1]])
                    box(title = tags$p("Raw File: ",
                                style = "font-size: 26px;
                                         font-weight: bold;"),
                        solidHeader = TRUE,
                        status = "success", width = 12,
                        h1(paste0(
                            SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]])),
                        tags$h5(paste("( full path:",
                                      SangerSingleReadAFN[[
                                          strtoi(sidebar_menu[[1]])]],
                                      ")"), style = "font-style:italic")),
                    box(title = tags$p("Quality Report: ",
                                       style = "font-size: 26px;
                                       font-weight: bold;"),
                        solidHeader = TRUE, collapsible = TRUE,
                        status = "success", width = 12,
                        tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                        fluidRow(
                            column(3,
                                   uiOutput("cutoffQualityScore"),
                                   tags$ul(
                                       textInput("cutoffQualityScoreText",
                                                 label = p("Change Value"),
                                                 value = toString(
                                                     SangerSingleReadQualReport[[
                                                         strtoi(sidebar_menu[[1]])]]@
                                                         cutoffQualityScore),
                                                 width = '90%')
                                   ),
                            ),
                            column(3,
                                   uiOutput("slidingWindowSize"),
                                   tags$ul(
                                       textInput("slidingWindowSizeText",
                                                 label = p("Change Value"),
                                                 value =
                                                     toString(SangerSingleReadQualReport[[
                                                         strtoi(sidebar_menu[[1]])]]@
                                                             slidingWindowSize),
                                                 width = '90%')
                                   ),
                            ),
                            column(3,
                                   uiOutput("trimmingStartPos"),
                            ),
                            column(3,
                                   uiOutput("trimmingFinishPos"),
                            )
                        ),
                        tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                        fluidRow(
                            column(6,
                                   uiOutput("trimmedRatio")
                            ),
                            column(6,
                                   uiOutput("remainingBP")
                            )
                        ),
                        box(title = tags$p("Cumulative Ratio Plot",
                                           style = "font-size: 24px;
                                       font-weight: bold;"),
                            collapsible = TRUE,
                            status = "success", width = 6,
                            plotlyOutput("qualityTrimmingRatioPlot")),
                        box(title = tags$p("Cumulative Ratio Plot",
                                           style = "font-size: 24px;
                                       font-weight: bold;"),
                            collapsible = TRUE,
                            status = "success", width = 6,
                            plotlyOutput("qualityQualityBasePlot")),
                    )
                )
            }
        }
    })

    ############################################################################
    ### observeEvent(input$
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic rightHeader text
    ### ------------------------------------------------------------------------
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        # message("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
        if (!is.na(as.numeric(sidebar_menu[[1]]))) {
            trimmedRV[["trimmedStart"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmingStartPos
            trimmedRV[["trimmedEnd"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmingFinishPos
            qualityPhredScores = SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityPhredScores

            readLen = length(qualityPhredScores)
            trimmedRV[["remainingBP"]] <- trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1
            trimmedRV[["trimmedRatio"]] <- round(((trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1) / readLen) * 100, digits = 2)
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        btn <- input$saveS4
        id <- paste0('txt', btn)
        newS4Object <- file.path(tempdir(), "SangerConsensus.Rda")
        saveRDS(SangerConsensus, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        showNotification(paste("New S4 object is store as:", newS4Object),
                         type = "message", duration = 10)
        NEW_SANGER_CONSENSUS_READ <<- readRDS(file=newS4Object)
        # shinyOptions(NewSangerConsensusSet = newS4)
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        btn <- input$closeUI
        stopApp()
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Change cutoffQualityScoreText
    ### ------------------------------------------------------------------------
    observeEvent(input$cutoffQualityScoreText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (!is.na(strtoi(input$cutoffQualityScoreText)) &&
            strtoi(input$cutoffQualityScoreText) > 0 &&
            strtoi(input$cutoffQualityScoreText) <= 60 &&
            strtoi(input$cutoffQualityScoreText) %% 1 ==0) {
            inputCutoffQualityScoreText <- input$cutoffQualityScoreText
        } else {
            inputCutoffQualityScoreText <- 20
        }
        trimmingPos <-
            inside_calculate_trimming(
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    qualityBaseScore, strtoi(inputCutoffQualityScoreText),
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    slidingWindowSize)
        if (!is.null(trimmingPos[1]) && !is.null(trimmingPos[2])) {
            SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                cutoffQualityScore <<- strtoi(inputCutoffQualityScoreText)
            SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                trimmingStartPos <<- trimmingPos[1]
            SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                trimmingFinishPos <<- trimmingPos[2]
            trimmedRV[["trimmedStart"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmingStartPos
            trimmedRV[["trimmedEnd"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmingFinishPos
            qualityPhredScores = SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityPhredScores
            readLen = length(qualityPhredScores)
            trimmedRV[["remainingBP"]] <-
                trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1
            trimmedRV[["trimmedRatio"]] <-
                round(((trimmedRV[["trimmedEnd"]] -
                            trimmedRV[["trimmedStart"]] + 1) / readLen)*100,
                      digits = 2)
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Change slidingWindowSizeText
    ### ------------------------------------------------------------------------
    observeEvent(input$slidingWindowSizeText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (!is.na(strtoi(input$slidingWindowSizeText)) &&
            strtoi(input$slidingWindowSizeText) > 0 &&
            strtoi(input$slidingWindowSizeText) <= 20 &&
            strtoi(input$slidingWindowSizeText) %% 1 ==0) {
            inputSlidingWindowSizeText <- input$slidingWindowSizeText
        } else {
            inputSlidingWindowSizeText <- 5
        }
        trimmingPos <-
            inside_calculate_trimming(SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityBaseScore,
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@cutoffQualityScore,
                strtoi(inputSlidingWindowSizeText))
        if (!is.null(trimmingPos[1]) && !is.null(trimmingPos[2])) {
            SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                slidingWindowSize <<- strtoi(inputSlidingWindowSizeText)
            SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                trimmingStartPos <<- trimmingPos[1]
            SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                trimmingFinishPos <<- trimmingPos[2]
            trimmedRV[["trimmedStart"]] <-
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                trimmingStartPos
            trimmedRV[["trimmedEnd"]] <-
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                trimmingFinishPos
            qualityPhredScores = SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityPhredScores
            readLen = length(qualityPhredScores)
            trimmedRV[["remainingBP"]] <-
                trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1
            trimmedRV[["trimmedRatio"]] <-
                round(((trimmedRV[["trimmedEnd"]] -
                            trimmedRV[["trimmedStart"]] + 1) / readLen)*100,
                      digits = 2)
        }
    })

    ############################################################################
    ### ConsensusRead (Overview)
    ############################################################################
    # chromatogram
    # output$chromatogram <- renderUI({
    #     sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
    #     chromatogram(SangerConsensusFRReadsList[[strtoi(sidebar_menu[[1]])]])
    #
    #     chromatogramIn(SangerConsensusFRReadsList[[strtoi(sidebar_menu[[1]])]], trim5=0, trim3=0,
    #                     showcalls=c("primary", "secondary", "both", "none"),
    #                     width=100, height=2, cex.mtext=1, cex.base=1, ylim=3,
    #                     filename=NULL, showtrim=FALSE, showhets=TRUE)
    #     # valueBox(
    #     #     subtitle = tags$p("MinReadsNum",
    #     #                       style = "font-size: 15px;
    #     #                                     font-weight: bold;"),
    #     #     value = tags$p(strtoi(input$SCMinReadsNumText),
    #     #                    style = "font-size: 29px;"),
    #     #     icon = icon("cut", "fa-sm"),
    #     #     color = "olive",
    #     #     width = 12,
    #     # )
    # })

    ### ------------------------------------------------------------------------
    ### valueBox: Change cutoffQualityScore
    ### ------------------------------------------------------------------------
    output$SCMinReadsNum <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("MinReadsNum",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(input$SCMinReadsNumText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })

    output$SCMinReadLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("MinReadLength",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(input$SCMinReadLengthText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })

    output$SCMinFractionCall <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(input$SCMinFractionCallText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })

    output$SCMaxFractionLost <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(input$SCMaxFractionLostText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })

    output$SCAcceptStopCodons <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (input$SCAcceptStopCodonsText == TRUE) {
            inputSCAcceptStopCodonsText = "TRUE"
        } else if (input$SCAcceptStopCodonsText == FALSE) {
            inputSCAcceptStopCodonsText = "FALSE"
        } else {
            inputSCAcceptStopCodonsText = "NA"
        }
        valueBox(
            subtitle = tags$p("AcceptStopCodons",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(inputSCAcceptStopCodonsText,
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })

    output$SCReadingFrame <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("ReadingFrame",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(input$SCReadingFrameText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })

    ############################################################################
    ### Each Read
    ############################################################################
    ### ------------------------------------------------------------------------
    ### valueBox: Change cutoffQualityScore
    ### ------------------------------------------------------------------------
    output$cutoffQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (!is.na(strtoi(input$cutoffQualityScoreText)) &&
            strtoi(input$cutoffQualityScoreText) > 0 &&
            strtoi(input$cutoffQualityScoreText) <= 60 &&
            strtoi(input$cutoffQualityScoreText) %% 1 ==0) {
            inputCutoffQualityScoreText <- input$cutoffQualityScoreText
        } else {
            inputCutoffQualityScoreText <- 20
        }
        valueBox(
            subtitle = tags$p("Cut Off Quality Score",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(inputCutoffQualityScoreText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change slidingWindowSize
    ### ------------------------------------------------------------------------
    output$slidingWindowSize <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (!is.na(strtoi(input$slidingWindowSizeText)) &&
            strtoi(input$slidingWindowSizeText) > 0 &&
            strtoi(input$slidingWindowSizeText) <= 20 &&
            strtoi(input$slidingWindowSizeText) %% 1 ==0) {
            inputSlidingWindowSizeText <- input$slidingWindowSizeText
        } else {
            inputSlidingWindowSizeText <- 5
        }
        valueBox(
            # strtoi(input$cutoffQualityScoreText
            subtitle = tags$p("Sliding Window Size ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(inputSlidingWindowSizeText),
                           style = "font-size: 29px;"),
            icon = icon("expand", "fa-sm"),
            color = "olive", width = 12,
        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change trimmingStartPos
    ### ------------------------------------------------------------------------
    output$trimmingStartPos <- renderUI({
        input$cutoffQualityScoreText
        input$slidingWindowSizeText
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("Trimming Start Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedStart"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change trimmingFinishPos
    ### ------------------------------------------------------------------------
    output$trimmingFinishPos <- renderUI({
        input$cutoffQualityScoreText
        input$slidingWindowSizeText
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("Trimming End Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedEnd"]]),
                           style = "font-size: 29px;"),
            icon = icon("times-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change remainingBP
    ### ------------------------------------------------------------------------
    output$remainingBP <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(8,
               valueBox(
                   subtitle = tags$p("Remaining Read Length",
                                     style = "font-size: 15px;
                                            font-weight: bold;"),
                   value = tags$p(paste(strtoi(trimmedRV[["remainingBP"]]), "BPs"),
                                  style = "font-size: 32px;"),
                   icon = icon("dna", "fa-sm"),
                   color = "olive",
                   width = 12,
               ),
        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change trimmedRatio
    ### ------------------------------------------------------------------------
    output$trimmedRatio <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(8,
               valueBox(
                   subtitle = tags$p("Remaining Ratio",
                                     style = "font-size: 15px;
                                            font-weight: bold;"),
                   value = tags$p(paste(trimmedRV[["trimmedRatio"]], "%"),
                                  style = "font-size: 32px;"),
                   icon = icon("divide", "fa-sm"),
                   color = "olive",
                   width = 12,
               ),
        )
    })

    output$clientdataText <- renderText({
        SangerSingleReadNum
    })

    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        readFeature <- SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]]
        trimmingStartPos = trimmedRV[["trimmedStart"]]
        trimmingFinishPos = trimmedRV[["trimmedEnd"]]
        qualityPhredScores = SangerSingleReadQualReport[[
            strtoi(sidebar_menu[[1]])]]@qualityPhredScores
        readLen = length(qualityPhredScores)

        stepRatio = 1 / readLen
        trimmingStartPos / readLen
        trimmingFinishPos / readLen

        trimmedPer <- c()
        remainingPer <- c()

        for (i in 1:trimmingStartPos) {
            if (i != trimmingStartPos) {
                trimmedPer <- c(trimmedPer, stepRatio)
            }
        }

        for (i in trimmingStartPos:trimmingFinishPos) {
            trimmedPer <- c(trimmedPer, 0)
        }


        for (i in trimmingFinishPos:readLen) {
            if (i != trimmingFinishPos) {
                trimmedPer <- c(trimmedPer, stepRatio)
            }
        }

        trimmedPer <- cumsum(trimmedPer)
        remainingPer = 1 - trimmedPer

        PerData <- data.frame(1:length(trimmedPer),
                              trimmedPer, remainingPer)

        colnames(PerData) <- c("Base",
                               "Trimmed Ratio",
                               "Remaining Ratio")
        # Change font setting
        # f <- list(
        #     family = "Courier New, monospace",
        #     size = 18,
        #     color = "#7f7f7f"
        # )
        x <- list(
            title = "Base Pair Index"
            # titlefont = f
        )
        y <- list(
            title = "Read Ratio"
            # titlefont = f
        )
        PerDataPlot <- melt(PerData, id.vars = c("Base"))
        plot_ly(data=PerDataPlot,
                x=~Base,
                y=~value,
                mode="markers",
                color = ~variable,
                text = ~paste("BP Index : ",
                              Base, '<sup>th</sup><br>Read Ratio :',
                              round(value*100, digits = 2), '%')) %>%
            layout(xaxis = x, yaxis = y, legend = list(orientation = 'h',
                                                       xanchor = "center",  # use center of legend as anchor
                                                       x = 0.5, y = 1.1)) %>%
            add_annotations(
                text = "Trimmed Ratio (Each BP)",
                x = (trimmingStartPos + trimmingFinishPos) / 2,
                y = ((trimmedPer[1] + trimmedPer[length(trimmedPer)]) / 2)
                    + 0.06,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Remaining Ratio (Each BP)",
                x = (trimmingStartPos+trimmingFinishPos) / 2,
                y = ((remainingPer[1] + remainingPer[length(remainingPer)]) / 2)
                    - 0.06,
                showarrow=FALSE
            )
    })

    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        readFeature <- SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]]
        trimmingStartPos = trimmedRV[["trimmedStart"]]
        trimmingFinishPos = trimmedRV[["trimmedEnd"]]
        qualityPhredScores = SangerSingleReadQualReport[[
            strtoi(sidebar_menu[[1]])]]@qualityPhredScores
        readLen = length(qualityPhredScores)

        qualityPlotDf<- data.frame(1:length(qualityPhredScores),
                                   qualityPhredScores)
        colnames(qualityPlotDf) <- c("Index", "Score")
        x <- list(
            title = "Base Pair Index"
            # titlefont = f
        )
        y <- list(
            title = "Phred Quality Score"
            # titlefont = f
        )

        plot_ly(data=qualityPlotDf,
                x=~Index) %>%
            add_markers(y=~Score,
                        text = ~paste("BP Index : ",
                                      Index,
                                      '<sup>th</sup><br>Phred Quality Score :',
                                      Score),
                        name = 'Quality Each BP') %>%
            add_trace(x=seq(trimmingStartPos,
                            trimmingFinishPos,
                            len=trimmingFinishPos-trimmingStartPos+1),
                      y=rep(70, trimmingFinishPos-trimmingStartPos+1),
                      mode="lines", hoverinfo="text",
                      text=paste("Trimmed Reads BP length:",
                                 trimmingFinishPos-trimmingStartPos+1,
                                 "BPs <br>",
                                 "Trimmed Reads BP ratio:",
                                 round((trimmingFinishPos - trimmingStartPos+1)/
                                           readLen * 100,
                                       digits=2),
                                 "%"),
                      line = list(width = 12),
                      name = 'Trimmed Read') %>%
            add_trace(x=seq(0,readLen,len=readLen),
                      y=rep(80, readLen), mode="lines", hoverinfo="text",
                      text=paste("Whole Reads BP length:",
                                 readLen,
                                 "BPs <br>",
                                 "Trimmed Reads BP ratio: 100 %"),
                      line = list(width = 12),
                      name = 'Whole Read') %>%
            layout(xaxis = x, yaxis = y,
                   shapes = list(vline(trimmingStartPos),
                                 vline(trimmingFinishPos)),
                   legend = list(orientation = 'h',
                                 xanchor = "center",  # use center of legend as anchor
                                 x = 0.5, y = 1.1)) %>%
            # add_segments(x = trimmingStartPos, xend = trimmingFinishPos, y = 70, yend = 70, inherit = TRUE, width = 10, line = list(width = 8)) %>%
            # add_segments(x = 0, xend = readLen, y = 75, yend = 75, inherit = TRUE, width = 4, line = list(width = 8)) %>%
            add_annotations(
                text = "Trimming Strat <br> BP Index",
                x = trimmingStartPos + 40,
                y = 15,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Trimming End <br> BP Index",
                x = trimmingFinishPos - 40,
                y = 15,
                showarrow=FALSE
            )
            # add_markers(qualityPlotDf, x=~Index, y=~Score)
            # add_segments(x = trimmingStartPos, xend = trimmingFinishPos, y = 70, yend = 70, inherit = TRUE)
    })


    output$differencesDF = renderDataTable({

        SCDifferencesDF
        # SCDistanceMatrix <- SangerConsensus@distanceMatrix
        #
        # SCDendrogram <- SangerConsensus@dendrogram
        #
        # SCIndelsDF <- SangerConsensus@indelsDF
        # SCStopCodonsDF <- SangerConsensus@stopCodonsDF
        #
        # SCSecondaryPeakDF <- SangerConsensus@secondaryPeakDF
        # dataTableOutput("secondaryPeakDF")

        # differencesDF             = "data.frame",
        # distanceMatrix            = "matrix",
        # dendrogram                = "list",
        # indelsDF                  = "data.frame",
        # stopCodonsDF              = "data.frame",
        # secondaryPeakDF           = "data.frame"

    })


    output$secondaryPeakDF = renderDataTable({

        SCDistanceMatrix
    })

    output$info <- renderText({
        paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
    })
}

vline <- function(x = 0, color = "red") {
    list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = x,
        x1 = x,
        line = list(color = color)
    )
}




inside_calculate_trimming <- function(qualityBaseScore,
                                      cutoffQualityScore,
                                      slidingWindowSize) {
    readLen <- length(qualityBaseScore)
    qualityPbCutoff <- 10** (cutoffQualityScore / (-10.0))
    remainingIndex <- c()
    if (slidingWindowSize > 20 || slidingWindowSize < 0 ||
        slidingWindowSize%%1!=0 ||
        cutoffQualityScore > 60 || cutoffQualityScore < 0 ||
        cutoffQualityScore%%1!=0) {
        trimmingStartPos = NULL
        trimmingFinishPos = NULL
    } else {
        for (i in 1:(readLen-slidingWindowSize+1)) {
            meanSLidingWindow <-
                mean(qualityBaseScore[i:(i+slidingWindowSize-1)])
            if (meanSLidingWindow < qualityPbCutoff) {
                remainingIndex <- c(remainingIndex, i)
                # or ==> i + floor(slidingWindowSize/3)
            }
        }
        trimmingStartPos = remainingIndex[1]
        trimmingFinishPos = remainingIndex[length(remainingIndex)]
    }
    return(c(trimmingStartPos, trimmingFinishPos))
}

