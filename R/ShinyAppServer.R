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

    trimmedRV <- reactiveValues(trimmedStart=0, trimmedEnd=0)

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
                    box(title = "Raw File: ", solidHeader = TRUE,
                        status = "success", width = 12,
                        h1(paste0(
                            SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]])),
                        h5(
                            paste("( full path:",
                                  SangerSingleReadAFN[[
                                      strtoi(sidebar_menu[[1]])]],
                                  ")"))),
                    box(title = tags$p("Quality Visualization: ",
                                       style = "font-size: 26px;
                                       font-weight: bold;"),
                        status = "success", width = 12,
                        # data("mtcars"),
                        fluidRow(
                            uiOutput("cutoffQualityScore"),
                            uiOutput("slidingWindowSize"),
                            uiOutput("trimmingStartPos"),
                            uiOutput("trimmingFinishPos"),
                        ),
                        tags$hr(style = ("border-top: 5px double #A9A9A9;")),
                        plotlyOutput("qualityTrimmingRatioPlot")
                    ),








                    # readFeature <- SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]]
                    # trimmingStartPos = SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]
                    # trimmingFinishPos = object@trimmingFinishPos
                    # readLen = length(object@qualityPhredScores)
                    #
                    # stepRatio = 1 / readLen
                    # trimmingStartPos / readLen
                    # trimmingFinishPos / readLen
                    #
                    # trimmedPer <- c()
                    # remainingPer <- c()
                    #
                    # for (i in 1:trimmingStartPos) {
                    #     if (i != trimmingStartPos) {
                    #         trimmedPer <- c(trimmedPer, stepRatio)
                    #     }
                    # }
                    #
                    # for (i in trimmingStartPos:trimmingFinishPos) {
                    #     trimmedPer <- c(trimmedPer, 0)
                    # }
                    #
                    #
                    # for (i in trimmingFinishPos:readLen) {
                    #     if (i != trimmingFinishPos) {
                    #         trimmedPer <- c(trimmedPer, stepRatio)
                    #     }
                    # }
                    #
                    # trimmedPer <- cumsum(trimmedPer)
                    # remainingPer = 1 - trimmedPer
                    #
                    # PerData <- data.frame(1:length(trimmedPer),
                    #                       trimmedPer, remainingPer)
                    #
                    # colnames(PerData) <- c("Base",
                    #                        "Trimmed Percent",
                    #                        "Remaining Percent")
                    #
                    # PerDataPlot <- melt(PerData, id.vars = c("Base"))
                    # plot_ly(data=mtcars, x=~wt, y=~mpg, mode="markers")

                    # plot_ly(data=mtcars, x=~wt, y=~mpg, mode="markers"))
                    # h5(paste0("Absolute File Path: ",  SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]]))

                    # box(SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]]),
                    # box(SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]]),
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
            message("sidebar_menu[[1]]: ",
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@trimmingStartPos)
            message("sidebar_menu[[1]]: ",
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@trimmingFinishPos)
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
            message("Outside!!!: ", trimmingPos[1], "  ", trimmingPos[1])
            if (!is.null(trimmingPos[1]) && !is.null(trimmingPos[2])) {
                message("Iutside!!!")
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
                                  style = "font-size: 29px;"),
                   icon = icon("times-circle", "fa-sm"),
                   color = "olive", width = 12,
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
                               "Trimmed Percent",
                               "Remaining Percent")
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
            layout(xaxis = x, yaxis = y) %>%
            add_annotations(
                text = "Trimmed",
                x = (trimmingStartPos + trimmingFinishPos) / 2,
                y = ((trimmedPer[1] + trimmedPer[length(trimmedPer)]) / 2)
                    + 0.05,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Remaining",
                x = (trimmingStartPos+trimmingFinishPos) / 2,
                y = ((remainingPer[1] + remainingPer[length(remainingPer)]) / 2)
                    - 0.05,
                showarrow=FALSE
            )
    })

    output$info <- renderText({
        paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
    })
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

