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
    ### Adding dynamic rightHeader text
    ### ------------------------------------------------------------------------
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
    })

    output$clientdataText <- renderText({
        SangerSingleReadNum
    })

    ### ------------------------------------------------------------------------
    ### Dynamic page navigation
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

    output$singelReadMenu_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
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
                    box(title = "Quality Visualization: ",
                        status = "success", width = 12,
                        # data("mtcars"),
                        fluidRow(
                            uiOutput("cutoffQualityScore"),
                            uiOutput("slidingWindowSize"),
                            uiOutput("trimmingStartPos"),
                            uiOutput("trimmingFinishPos"),
                        ),
                        tags$hr(),
                        fluidRow(
                            # Clicking this will increment the progress amount
                            box(width = 4, actionButton("count", "Increment progress"))
                        ),
                        tags$hr(),
                        h4(paste("Cut Off Quality Score\t: ", toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore))),
                        p("fdsffd"),
                        h4(paste("Sliding Window Size\t: ", toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@slidingWindowSize))),
                        h4(toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos)),
                        h4(toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingFinishPos)),
                        h4(toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@qualityPhredScores)),
                        h4(toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@qualityBaseScore)),
                    )
                        # plot_ly(data=mtcars, x=~wt, y=~mpg, mode="markers"))
                    # h5(paste0("Absolute File Path: ",  SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]]))

                    # box(SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]]),
                    # box(SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]]),
                )
            }
        }
    })

    ### ------------------------------------------------------------------------
    ### Button: Save S4 object
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
    ### Button: Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        btn <- input$closeUI
        stopApp()
    })


    observeEvent(input$cutoffQualityScoreText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        message("Before: SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos: ", SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos)
        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore <<- strtoi(input$cutoffQualityScoreText)
        trimmingPos <- inside_calculate_trimming(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@qualityBaseScore, SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore, SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@slidingWindowSize)
        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos <<- trimmingPos[1]
        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingFinishPos <<- trimmingPos[2]
        message("After: SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos: ", SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos)
        trimmedRV[["trimmedStart"]] <-  SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos
        trimmedRV[["trimmedEnd"]] <-  SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingFinishPos

    })

    observeEvent(input$slidingWindowSizeText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@slidingWindowSize <<- strtoi(input$slidingWindowSizeText)
        trimmingPos <- inside_calculate_trimming(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@qualityBaseScore, SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore, SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@slidingWindowSize)
        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos <<- trimmingPos[1]
        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingFinishPos <<- trimmingPos[2]
        trimmedRV[["trimmedStart"]] <-  SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos
        trimmedRV[["trimmedEnd"]] <-  SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingFinishPos
    })

    # eventReactive(input$trimmingStartPosText, {
    #     sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
    # })
    #
    # eventReactive(input$trimmingFinishPos, {
    #     sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
    # })

    output$cutoffQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   subtitle = tags$p("Cut Off Quality Score: ", style = "font-size: 15px; font-weight: bold;"), value = tags$p(strtoi(input$cutoffQualityScoreText), style = "font-size: 29px;"), icon = icon("list"),
                   color = "light-blue", width = 12,
               ),
               tags$ul(
                   textInput("cutoffQualityScoreText", label = p("Change Value"), value = toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore), width = '95%')
               ),
        )
    })


    output$slidingWindowSize <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   # strtoi(input$cutoffQualityScoreText
                   subtitle = tags$p("Sliding Window Size: ", style = "font-size: 15px; font-weight: bold;"), value = tags$p(strtoi(input$slidingWindowSizeText), style = "font-size: 29px;"), icon = icon("list"),
                   color = "light-blue", width = 12,
               ),
               tags$ul(
                   textInput("slidingWindowSizeText", label = p("Change Value"), value = toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@slidingWindowSize), width = '95%')
               ),

        )
    })

    output$trimmingStartPos <- renderUI({
        input$cutoffQualityScoreText
        input$slidingWindowSizeText
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   subtitle = tags$p("Cut Off Quality Score: ", style = "font-size: 15px; font-weight: bold;"), value = tags$p(toString(trimmedRV[["trimmedStart"]]), style = "font-size: 29px;"), icon = icon("list"),
                   color = "light-blue", width = 12,
               ),
               # tags$ul(
               #     # hidden(
               #      textInput("trimmingStartPosText", label = p("Change Value"), value = toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos), width = '95%')
               #     # )
               # ),
               # trimmedRV[["trimmedStart"]] <-  SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingStartPos
               # trimmedRV[["trimmedEnd"]] <-  SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingFinishPos
        )
    })

    output$trimmingFinishPos <- renderUI({
        input$cutoffQualityScoreText
        input$slidingWindowSizeText
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        column(3,
               valueBox(
                   subtitle = tags$p("Cut Off Quality Score: ", style = "font-size: 15px; font-weight: bold;"), value = tags$p(toString(trimmedRV[["trimmedEnd"]]), style = "font-size: 29px;"), icon = icon("list"),
                   color = "light-blue", width = 12,
               ),
               # tags$ul(
               #     # hidden(
               #      textInput("trimmingFinishPosText", label = p("Change Value"), value = toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@trimmingFinishPos), width = '95%')
               #     # )
               # ),
        )
    })
}





inside_calculate_trimming <- function(qualityBaseScore,
                                      cutoffQualityScore,
                                      slidingWindowSize) {
    readLen <- length(qualityBaseScore)
    qualityPbCutoff <- 10** (cutoffQualityScore / (-10.0))
    remainingIndex <- c()
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

    return(c(trimmingStartPos, trimmingFinishPos))
}

