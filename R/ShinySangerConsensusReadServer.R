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
        paste0(i, " ",
               SangerConsensus@forwardReadsList[[i]]@readFeature))
    reverseReadFeature <- sapply(1:reverseReadNum, function(i)
        paste0(i+forwardReadNum, " ",
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


    # basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    # basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    # aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    # basecalls1 <- basecalls1[1:length(aveposition)]
    # basecalls2 <- basecalls2[1:length(aveposition)]

    # primarySeq
    forwardReadPrimSeq <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@primarySeq)
    reverseReadPrimSeq <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@primarySeq)
    SangerSingleReadPrimSeq <- c(forwardReadPrimSeq, reverseReadPrimSeq)
    SangerSingleReadPrimSeqChara <- as.character(SangerSingleReadPrimSeq)











    # primarySeqDF
    forwardReadPrimSeqDF <- lapply(1:forwardReadNum, function(i) {
        basecalls1 <- unlist(strsplit(toString(SangerConsensus@forwardReadsList[[i]]@primarySeq), ""))
        aveposition <- rowMeans(SangerConsensus@forwardReadsList[[i]]@peakPosMatrix, na.rm=TRUE)
        basecalls1 <- basecalls1[1:length(aveposition)]
        basecalls1DF <- data.frame(t(data.frame(basecalls1)))
        return(basecalls1DF)
        }
    )
    reverseReadPrimSeqDF <- lapply(1:reverseReadNum, function(i) {
        basecalls1 <- unlist(strsplit(toString(SangerConsensus@reverseReadsList[[i]]@primarySeq), ""))
        aveposition <- rowMeans(SangerConsensus@reverseReadsList[[i]]@peakPosMatrix, na.rm=TRUE)
        basecalls1 <- basecalls1[1:length(aveposition)]
        basecalls1DF <- data.frame(t(data.frame(basecalls1)))
        return(basecalls1DF)
        }
    )
    SangerSingleReadPrimSeqDF <- c(forwardReadPrimSeqDF, reverseReadPrimSeqDF)


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
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (input$sidebar_menu == "Sanger Consensus Read Overview") {
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
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        data(mtcars)
        if (input$sidebar_menu != "Sanger Consensus Read Overview") {
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
                    ),
                    box(title = tags$p("Base Call: ",
                                       style = "font-size: 26px;
                                       font-weight: bold;"),
                        solidHeader = TRUE, collapsible = TRUE,
                        status = "success", width = 12,
                        tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                        box(title = tags$p("BP (ATCG)",
                                           style = "font-size: 24px;
                                       font-weight: bold;"),
                            collapsible = TRUE,
                            status = "success", width = 12,
                            column(width = 12,
                                   column(width = 12,
                                          rHandsontableOutput("primarySeqDF"),
                                          rHandsontableOutput("primarySeqDF2")
                                    ),
                            )
                                   # style = "height:100px; overflow-y: scroll;overflow-x: scroll;")
                        ),
                    )
                )
            }
        }
    })










    ############################################################################
    ### observeEvent(input$
    ############################################################################
    observeEventDynamicHeaderSC(input, output, session, trimmedRV,
                              SangerSingleReadQualReport)
    observeEventButtonSaveSC(input, output, session, SangerConsensus)
    observeEventButtonClose(input, output, session)


    observeEvent(input$cutoffQualityScoreText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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
    ### ConsensusRead (Sanger Consensus Read Overview)
    ############################################################################
    # chromatogram
    # output$chromatogram <- renderUI({
    #     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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









    # obj <- A_chloroticaSingleRead
    #
    # chromatogram(homosangerseq)
    #
    # # A_chloroticaSingleRead
    # chromatogram <- function(obj, trim5=0, trim3=0,
    #          showcalls="both",
    #          width=100, height=2, cex.mtext=1, cex.base=1, ylim=3,
    #          filename=NULL, showtrim=FALSE, showhets=TRUE) {
    #     originalpar <- par(no.readonly=TRUE)
    #     showcalls <- showcalls[1]
    #     traces <- obj@traceMatrix
    #     basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    #     basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    #     aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    #     basecalls1 <- basecalls1[1:length(aveposition)]
    #     basecalls2 <- basecalls2[1:length(aveposition)]
    #
    #
    #
    #
    #
    #     if(showtrim == FALSE) {
    #         if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
    #         else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
    #         if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
    #         else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
    #         aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)]
    #     }
    #     indexes <- 1:length(basecalls1)
    #     trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all
    #     #false if not trimmed
    #     if (!is.null(trim3)) {
    #         traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10,
    #                                 nrow(traces))), ]
    #     }
    #     if (!is.null(trim5)) {
    #         offset <- max(c(1, aveposition[1] - 10))
    #         traces <- traces[offset:nrow(traces),]
    #         aveposition <- aveposition - (offset-1)
    #     }
    #     maxsignal <- apply(traces, 1, max)
    #     ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))
    #     p <- c(0, aveposition, nrow(traces))
    #     midp <- diff(p)/2
    #     starts <- aveposition - midp[1:(length(midp)-1)]
    #     starthets <- starts
    #     starthets[basecalls1 == basecalls2] <- NA
    #     ends <- aveposition + midp[2:(length(midp))]
    #     endhets <- ends
    #     endhets[basecalls1 == basecalls2] <- NA
    #     starttrims <- starts
    #     starttrims[!trimmed] <- NA
    #     endtrims <- ends
    #     endtrims[!trimmed] <- NA
    #
    #     colortranslate <- c(A="green", C="blue", G="black", T="red")
    #     colorvector1 <- unname(colortranslate[basecalls1])
    #     colorvector1[is.na(colorvector1)] <- "purple"
    #     colorvector2 <- unname(colortranslate[basecalls2])
    #     colorvector2[is.na(colorvector2)] <- "purple"
    #
    #     valuesperbase <- nrow(traces)/length(basecalls1)
    #     tracewidth <- width*valuesperbase
    #     breaks <- seq(1,nrow(traces), by=tracewidth)
    #     numplots <- length(breaks)
    #     if(!is.null(filename)) pdf(filename, width=8.5, height=height*numplots)
    #     par(mar=c(2,2,2,1), mfrow=c(numplots, 1))
    #     basecallwarning1 = 0
    #     basecallwarning2 = 0
    #     j = 1
    #
    #     for(i in breaks) {
    #         range <- aveposition >= i & aveposition < (i+tracewidth)
    #         starthet <- starthets[range] - tracewidth*(j-1)
    #         starthet[starthet < 0] <- 0
    #         endhet <- endhets[range] - tracewidth*(j-1)
    #         endhet[endhet > tracewidth] <- tracewidth
    #         lab1 <- basecalls1[range]
    #         lab2 <- basecalls2[range]
    #         pos <- aveposition[range] - tracewidth*(j-1)
    #         colors1 <- colorvector1[range]
    #         colors2 <- colorvector2[range]
    #         starttrim <- starttrims[range] - tracewidth*(j-1)
    #         endtrim <- endtrims[range] - tracewidth*(j-1)
    #         plotrange <- i:min(i+tracewidth, nrow(traces))
    #         plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n",
    #              bty="n", xlab="", yaxt="n", , xlim=c(1,tracewidth))
    #         if (showhets==TRUE) {
    #             rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')
    #         }
    #         if (showtrim==TRUE) {
    #             rect(starttrim, 0, endtrim, ylims[2], col='red', border='transparent',
    #                  density=15)
    #         }
    #         lines(traces[plotrange,1], col="green")
    #         lines(traces[plotrange,2], col="blue")
    #         lines(traces[plotrange,3], col="black")
    #         lines(traces[plotrange,4], col="red")
    #         mtext(as.character(which(range)[1]), side=2, line=0, cex=cex.mtext)
    #
    #         for(k in 1:length(lab1)) {
    #             if (showcalls=="primary" | showcalls=="both") {
    #                 if (is.na(basecalls1[1]) & basecallwarning1==0) {
    #                     warning("Primary basecalls missing")
    #                     basecallwarning1 = 1
    #                 }
    #                 else if (length(lab1) > 0) {
    #                     axis(side=3, at=pos[k], labels=lab1[k], col.axis=colors1[k],
    #                          family="mono", cex=cex.base, line=ifelse(showcalls=="both", 0,
    #                                                                   -1), tick=FALSE)
    #                 }
    #             }
    #             if (showcalls=="secondary" | showcalls=="both") {
    #                 if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
    #                     warning("Secondary basecalls missing")
    #                     basecallwarning2 = 1
    #                 }
    #                 else if (length(lab2) > 0) {
    #                     axis(side=3, at=pos[k], labels=lab2[k], col.axis=colors2[k],
    #                          family="mono", cex=cex.base, line=-1, tick=FALSE)
    #                 }
    #             }
    #         }
    #         j = j + 1
    #     }
    #     if(!is.null(filename)) {
    #         dev.off()
    #         cat(paste("Chromatogram saved to", filename,
    #                   "in the current working directory"))
    #     }
    #     else par(originalpar)
    # }
















    output$primarySeqDF <- renderRHandsontable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        strtoi(sidebar_menu[[1]])
        rhandsontable(SangerSingleReadPrimSeqDF[[1]], rowHeaders = NULL, overflow = 'visible', mergeCells = TRUE)
    })


    output$primarySeqDF2 <- renderRHandsontable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        rhandsontable(SangerSingleReadPrimSeqDF[[strtoi(sidebar_menu[[1]])]], rowHeaders = NULL, overflow = 'visible', width = 100, height = 100)
    })























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
