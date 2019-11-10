### ============================================================================
### R shiny consensusRead server function
### ============================================================================
consensusReadServer <- function(input, output, session) {
    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusRead <- getShinyOption("SangerConsensusRead")
    SangerConsensus <- SangerConsensusRead[[1]]


    # ### ------------------------------------------------------------------------
    # ### ConsensusRead-related parameters initialization.
    # ### ------------------------------------------------------------------------
    SCMinReadsNum <- SangerConsensus@minReadsNum
    SCMinReadLength <- SangerConsensus@minReadLength


    SCRefAminoAcidSeq <- SangerConsensus@refAminoAcidSeq


    SCMinFractionCall <- SangerConsensus@minFractionCall
    SCMaxFractionLost <- SangerConsensus@maxFractionLost


    SCGeneticCode <- SangerConsensus@geneticCode


    SCAcceptStopCodons <- SangerConsensus@acceptStopCodons
    SCReadingFrame <- SangerConsensus@readingFrame

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
    peakAmpMatrix
    forwardReadPeakAmpMat <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@peakAmpMatrix)
    reverseReadPeakAmpMat <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@peakAmpMatrix)
    SangerSingleReadPeakAmpMat <- c(forwardReadPeakAmpMat, reverseReadPeakAmpMat)
    trimmedRV <- reactiveValues(trimmedStart = 0, trimmedEnd = 0,
                                remainingBP = 0, trimmedRatio = 0)

    ############################################################################
    ### output$ID
    ############################################################################
    dynamicMenuSideBarSC(input, output, session, SangerSingleReadNum, SangerSingleReadFeature)

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
    observeEventDynamicRightHeader(input, output, session, trimmedRV,
                                   SangerSingleReadQualReport)

    observeEventButtonSaveSC(input, output, session, SangerConsensus)
    observeEventButtonClose(input, output, session)


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

    valueBoxSCMinReadsNum(input, output, session)
    valueBoxSCMinReadLength(input, output, session)
    valueBoxSCMinFractionCall(input, output, session)
    valueBoxSCMaxFractionLost(input, output, session)
    valueBoxSCAcceptStopCodons(input, output, session)
    valueBoxSCReadingFrame(input, output, session)

    valueBoxCutoffQualityScore (input, output, session)
    valueBoxSlidingWindowSize (input, output, session)
    valueBoxTrimmingStartPos (input, output, session, trimmedRV)
    valueBoxTrimmingFinishPos (input, output, session, trimmedRV)
    valueBoxRemainingBP (input, output, session, trimmedRV)
    valueBoxTrimmedRatio (input, output, session, trimmedRV)
    clientdataText (input, output, session)

    qualityTrimmingRatioPlot (input, output, session, trimmedRV,
                              SangerSingleReadQualReport,
                              SangerSingleReadFeature)
    qualityQualityBasePlot (input, output, session, trimmedRV,
                            SangerSingleReadQualReport, SangerSingleReadFeature)
    overViewDifferencesDataFrame (input, output, session, SCDifferencesDF)
    overViewSecondaryPeakDataFrame (input, output, session, SCDistanceMatrix)

    output$info <- renderText({
        paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
    })
}
