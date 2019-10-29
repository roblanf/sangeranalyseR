### ============================================================================
### R shiny consensus read server function
### ============================================================================
consensusServer <- function(input, output, session) {
    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusRead <- getShinyOption("SangerConsensusReadSet")
    SangerConsensus <- SangerConsensusRead[[1]]
    SangerConsensusRegularExp <- SangerConsensus@readsRegularExp

    SangerSingleReadNum <- length((SangerConsensus)@SangerReadsList)
    SangerSingleReadBFN <- sapply(1:SangerSingleReadNum, function(i)
               basename(SangerConsensus@SangerReadsList[[i]]@readFileName))
    # readFeature
    SangerSingleReadFeature <- sapply(1:SangerSingleReadNum, function(i)
        paste0(i, "_",
               SangerConsensus@SangerReadsList[[i]]@readFeature))
    # readFileName
    SangerSingleReadAFN <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@readFileName)

    # abifRawData
    SangerSingleReadAbifRawData <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@abifRawData)

    # QualityReport
    SangerSingleReadQualReport <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@QualityReport)

    # primarySeqID
    SangerSingleReadPrimSeqID <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@primarySeqID)

    # primarySeq
    SangerSingleReadPrimSeq <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@primarySeq)

    # secondarySeqID
    SangerSingleReadSecoSeqID <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@secondarySeqID)

    # secondarySeq
    SangerSingleReadSecoSeq <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@secondarySeq)

    # traceMatrix
    SangerSingleReadTraceMat <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@traceMatrix)

    # peakPosMatrix
    SangerSingleReadPeakPosMat <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@peakPosMatrix)

    # peakAmpMatrix
    SangerSingleReadPeakAmpMat <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensus@SangerReadsList[[i]]@peakAmpMatrix)

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
                box("Overview")
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
                                      ")")), style = "font-style:italic"),
                    box(title = tags$p("Quality Report: ",
                                       style = "font-size: 26px;
                                       font-weight: bold;"),
                        solidHeader = TRUE, collapsible = TRUE,
                        status = "success", width = 12,
                        tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                        fluidRow(
                            uiOutput("cutoffQualityScore"),
                            uiOutput("slidingWindowSize"),
                            uiOutput("trimmingStartPos"),
                            uiOutput("trimmingFinishPos"),
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
                        column(6,
                               plotlyOutput("qualityTrimmingRatioPlot")),
                        column(6,
                               plotlyOutput("qualityQualityBasePlot")),
                    ),
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
        if (!is.na(strtoi(input$cutoffQualityScoreText))) {
            trimmingPos <-
                inside_calculate_trimming(
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        qualityBaseScore, strtoi(input$cutoffQualityScoreText),
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        slidingWindowSize)
            message("Outside!!!: ", trimmingPos[1], "  ", trimmingPos[1])
            if (!is.null(trimmingPos[1]) && !is.null(trimmingPos[2])) {
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    cutoffQualityScore <<- strtoi(input$cutoffQualityScoreText)
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
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Change slidingWindowSizeText
    ### ------------------------------------------------------------------------
    observeEvent(input$slidingWindowSizeText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (!is.na(strtoi(input$slidingWindowSizeText))) {
            trimmingPos <-
                inside_calculate_trimming(SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@qualityBaseScore,
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@cutoffQualityScore,
                    strtoi(input$slidingWindowSizeText))
            if (!is.null(trimmingPos[1]) && !is.null(trimmingPos[2])) {
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    slidingWindowSize <<- strtoi(input$slidingWindowSizeText)
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
        }
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change cutoffQualityScore
    ### ------------------------------------------------------------------------
    output$cutoffQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   subtitle = tags$p("Cut Off Quality Score",
                                     style = "font-size: 15px;
                                            font-weight: bold;"),
                   value = tags$p(strtoi(input$cutoffQualityScoreText),
                                  style = "font-size: 29px;"),
                   icon = icon("cut", "fa-sm"),
                   color = "olive",
                   width = 12,
               ),
               tags$ul(
                   textInput("cutoffQualityScoreText",
                             label = p("Change Value"),
                             value = toString(
                                 SangerSingleReadQualReport[[
                                     strtoi(sidebar_menu[[1]])]]@
                                     cutoffQualityScore),
                             width = '95%')
               ),
        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change slidingWindowSize
    ### ------------------------------------------------------------------------
    output$slidingWindowSize <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   # strtoi(input$cutoffQualityScoreText
                   subtitle = tags$p("Sliding Window Size ",
                                     style = "font-size: 15px;
                                            font-weight: bold;"),
                   value = tags$p(strtoi(input$slidingWindowSizeText),
                                  style = "font-size: 29px;"),
                   icon = icon("expand", "fa-sm"),
                   color = "olive", width = 12,
               ),
               tags$ul(
                   textInput("slidingWindowSizeText",
                             label = p("Change Value"),
                             value =
                                 toString(SangerSingleReadQualReport[[
                                     strtoi(sidebar_menu[[1]])]]@
                                         slidingWindowSize),
                             width = '95%')
               ),

        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change trimmingStartPos
    ### ------------------------------------------------------------------------
    output$trimmingStartPos <- renderUI({
        input$cutoffQualityScoreText
        input$slidingWindowSizeText
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   subtitle = tags$p("Trimming Start Pos ",
                                     style = "font-size: 15px;
                                            font-weight: bold;"),
                   value = tags$p(toString(trimmedRV[["trimmedStart"]]),
                                  style = "font-size: 29px;"),
                   icon = icon("check-circle", "fa-sm"),
                   color = "olive", width = 12,
               ),
        )
    })

    ### ------------------------------------------------------------------------
    ### valueBox: Change trimmingFinishPos
    ### ------------------------------------------------------------------------
    output$trimmingFinishPos <- renderUI({
        input$cutoffQualityScoreText
        input$slidingWindowSizeText
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   subtitle = tags$p("Trimming End Pos ",
                                     style = "font-size: 15px;
                                            font-weight: bold;"),
                   value = tags$p(toString(trimmedRV[["trimmedEnd"]]),
                                  style = "font-size: 40px;"),
                   icon = icon("times-circle", "fa-sm"),
                   color = "olive", width = 12,
               ),
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
                   subtitle = tags$p("Trimmed Ratio",
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
                              value, '%')) %>%
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
    if (slidingWindowSize > 20 || qualityBaseScore > 60) {
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

