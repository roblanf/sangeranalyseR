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

### ============================================================================
### Adding dynamic menu to sidebar.
### ============================================================================
dynamicMenuSideBarSC <- function(input, output, session,
                               SangerSingleReadNum, SangerSingleReadFeature) {
    output$singleReadMenu <- renderMenu({
        menu_list <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(text = SangerSingleReadFeature[i],
                          tabName = SangerSingleReadFeature[i],
                          selected = TRUE, icon = icon("minus")))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Overview")})
}

dynamicMenuSideBarSCSet <- function(input, output, session, SangerCSetParam) {
    output$singleReadMenu <- renderMenu({
        SangerCSNum <- length(SangerCSetParam)
        menu_list <- sapply(1:SangerCSNum, function(i) {
            SangerSingleReadNum <- SangerCSetParam[[i]]$SangerSingleReadNum
            SangerCSMenuSubItem <- sapply(1:SangerSingleReadNum, function(j) {
                list(menuSubItem(text = SangerCSetParam[[i]]$SangerSingleReadFeature[[j]],
                                 tabName = SangerCSetParam[[i]]$SangerSingleReadFeature[[j]]))
            })
            SangerCSetParam[[i]]$SCName
            list(menuItem(text = SangerCSetParam[[i]]$SCName,
                          tabName = SangerCSetParam[[i]]$SCName,
                          selected = TRUE, icon = icon("minus"),
                          SangerCSMenuSubItem))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Overview")})
}

### ============================================================================
### observeEvent: Adding dynamic rightHeader text
### ============================================================================
observeEventDynamicHeaderSC <- function(input, output, session, trimmedRV,
                                           SangerSingleReadQualReport) {
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        # message("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
        if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
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
}

observeEventDynamicHeaderSCSet <- function(input, output, session, trimmedRV,
                                           SangerCSetParam) {
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        # message("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
        if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
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
}

### ============================================================================
### observeEvent: Button Save S4 object
### ============================================================================
observeEventButtonSaveSC <- function(input, output, session, SangerConsensus) {
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
}

observeEventButtonSaveSCSet <- function(input, output, session, SangerCSetParam) {
    observeEvent(input$saveS4, {
        btn <- input$saveS4
        id <- paste0('txt', btn)
        newS4Object <- file.path(tempdir(), "SangerAlignedConsensusSet.Rda")
        saveRDS(SangerCSetParam, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        showNotification(paste("New S4 object is store as:", newS4Object),
                         type = "message", duration = 10)
        NEW_SANGER_CONSENSUS_READ <<- readRDS(file=newS4Object)
        # shinyOptions(NewSangerConsensusSet = newS4)
    })
}

### ============================================================================
### observeEvent: Button Close UI
### ============================================================================
observeEventButtonClose <- function(input, output, session) {
    observeEvent(input$closeUI, {
        btn <- input$closeUI
        stopApp()
    })
}

### ============================================================================
### valueBox: SCMinReadsNum
### ============================================================================
valueBoxSCMinReadsNum <- function(input, output, session) {
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
}

### ============================================================================
### valueBox: SCMinReadLength
### ============================================================================
valueBoxSCMinReadLength <- function(input, output, session) {
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
}

### ============================================================================
### valueBox: SCMinFractionCall
### ============================================================================
valueBoxSCMinFractionCall <- function(input, output, session) {
    output$SCMinFractionCall <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(suppressWarnings(as.numeric(input$SCMinFractionCallText)),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMaxFractionLost
### ============================================================================
valueBoxSCMaxFractionLost <- function(input, output, session) {
    output$SCMaxFractionLost <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(suppressWarnings(as.numeric(input$SCMaxFractionLostText)),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCAcceptStopCodons
### ============================================================================
valueBoxSCAcceptStopCodons <- function(input, output, session) {
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
}

### ============================================================================
### valueBox: SCReadingFrame
### ============================================================================
valueBoxSCReadingFrame <- function(input, output, session) {
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
}


################################################################################
### Each Read
################################################################################
### ============================================================================
### valueBox: Change cutoffQualityScore
### ============================================================================
valueBoxCutoffQualityScore <- function(input, output, session) {
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
}

### ============================================================================
### valueBox: Change slidingWindowSize
### ============================================================================
valueBoxSlidingWindowSize <- function(input, output, session) {
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
}

### ============================================================================
### valueBox: Change trimmingStartPos
### ============================================================================
valueBoxTrimmingStartPos <- function(input, output, session, trimmedRV) {
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
}

### ============================================================================
### valueBox: Change trimmingFinishPos
### ============================================================================
valueBoxTrimmingFinishPos <- function(input, output, session, trimmedRV) {
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
}

### ============================================================================
### valueBox: Change remainingBP
### ============================================================================
valueBoxRemainingBP <- function(input, output, session, trimmedRV) {
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
}

### ============================================================================
### valueBox: Change trimmedRatio
### ============================================================================
valueBoxTrimmedRatio <- function(input, output, session, trimmedRV) {
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
}

### ============================================================================
### clientdataText
### ============================================================================
clientdataText <- function(input, output, session) {
    output$clientdataText <- renderText({
        SangerSingleReadNum
    })
}


### ============================================================================
### qualityTrimmingRatioPlot
### ============================================================================
qualityTrimmingRatioPlot <- function(input, output, session, trimmedRV,
                                     SangerSingleReadQualReport,
                                     SangerSingleReadFeature) {
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
}

### ============================================================================
### qualityQualityBasePlot
### ============================================================================
qualityQualityBasePlot <- function(input, output, session, trimmedRV,
                                   SangerSingleReadQualReport,
                                   SangerSingleReadFeature) {
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
}


### ============================================================================
### differencesDataFrame
### ============================================================================
overViewDifferencesDataFrame <- function(input, output, session, SCDifferencesDF) {
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
}

### ============================================================================
### secondaryPeakDataFrame
### ============================================================================
overViewSecondaryPeakDataFrame <- function(input, output, session, SCDistanceMatrix) {
    output$secondaryPeakDF = renderDataTable({

        SCDistanceMatrix
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

