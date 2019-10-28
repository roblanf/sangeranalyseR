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
                            # infoBoxOutput("cutoffQualityScore"),
                            uiOutput("cutoffQualityScore"),
                                # infoBox(
                                #     title = tags$p("Cut Off Quality Score: ", style = "font-size: 15px; font-weight: bold;"), value = tags$p(strtoi(input$cutoffQualityScoreText), style = "font-size: 29px;"), icon = icon("list"),
                                #     color = "purple", width = 12, fill = TRUE,
                                # ),
                            # uiOutput("cutoffQualityScore"),
                        ),
                        fluidRow(
                            # Clicking this will increment the progress amount
                            box(width = 4, actionButton("count", "Increment progress"))
                        ),
                        tags$hr(),
                        infoBox("Cut Off Quality Score", 5 * 2, icon = icon("credit-card"), fill = TRUE),
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
        message("Inside 'observeEvent(input$cutoffQualityScoreText': ", sidebar_menu[[1]])
        message("Original one ", sidebar_menu[[1]], ": ", SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore)
        message("Side bar: ", (sidebar_menu[[1]]), ". Update SangerSingleReadQualReport value !: ", input$cutoffQualityScoreText)

        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore <<- strtoi(input$cutoffQualityScoreText)
    })

    # observeEvent(input$do, {
    #     sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
    #     message("Inside 'observeEvent(input$cutoffQualityScoreText': ", sidebar_menu[[1]])
    #     message("Original one ", sidebar_menu[[1]], ": ", SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore)
    #     message("Side bar: ", (sidebar_menu[[1]]), ". Update SangerSingleReadQualReport value !: ", input$cutoffQualityScoreText)
    #     SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore <<- strtoi(input$cutoffQualityScoreText)
    #     message("Update one ", sidebar_menu[[1]], ": ", SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore)
    # })

    output$cutoffQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        message("Inside 'output$cutoffQualityScore': ", sidebar_menu[[1]])
        column(4,

               infoBox(
                   # strtoi(input$cutoffQualityScoreText
                   title = tags$p("Cut Off Quality Score: ", style = "font-size: 15px; font-weight: bold;"), value = tags$p(strtoi(input$cutoffQualityScoreText), style = "font-size: 29px;"), icon = icon("list"),
                   color = "purple", width = 12, fill = TRUE,
               ),
               tags$ul(
                   textInput("cutoffQualityScoreText", label = p("Change Value"), value = toString(SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@cutoffQualityScore), width = '95%')
               ),

        )
    })



    output$approvalBox2 <- renderInfoBox({
        infoBox(
            "Approval", "80%", icon = icon("thumbs-up", lib = "glyphicon"),
            color = "yellow", fill = TRUE
        )
    })
}

