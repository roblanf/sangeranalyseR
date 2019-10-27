### ============================================================================
### R shiny consensus read server function
### ============================================================================
consensusServer <- function(input, output, session) {
    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusRead <- getShinyOption("SangerConsensusReadSet")
    SangerConsensusRead <- SangerConsensusRead[[1]]

    SangerConsensusReadRegularExp <- SangerConsensusRead@readsRegularExp

    SangerSingleReadNum <- length((SangerConsensusRead)@SangerReadsList)
    SangerSingleReadBFN <- sapply(1:SangerSingleReadNum, function(i)
               basename(SangerConsensusRead@SangerReadsList[[i]]@readFileName))
    # readFeature
    SangerSingleReadFeature <- sapply(1:SangerSingleReadNum, function(i)
        paste0(i, "_",
               SangerConsensusRead@SangerReadsList[[i]]@readFeature))
    # readFileName
    SangerSingleReadAFN <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@readFileName)

    # abifRawData
    SangerSingleReadAbifRawData <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@abifRawData)

    # QualityReport
    SangerSingleReadQualReport <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@QualityReport)

    # primarySeqID
    SangerSingleReadPrimSeqID <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@primarySeqID)

    # primarySeq
    SangerSingleReadPrimSeq <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@primarySeq)

    # secondarySeqID
    SangerSingleReadSecoSeqID <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@secondarySeqID)

    # secondarySeq
    SangerSingleReadSecoSeq <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@secondarySeq)

    # traceMatrix
    SangerSingleReadTraceMat <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@traceMatrix)

    # peakPosMatrix
    SangerSingleReadPeakPosMat <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@peakPosMatrix)

    # peakAmpMatrix
    SangerSingleReadPeakAmpMat <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@peakAmpMatrix)

    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar.
    ### ------------------------------------------------------------------------
    output$singleReadMenu <- renderMenu({
        menu_list <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(SangerSingleReadFeature[i], tabName = SangerSingleReadFeature[i], selected = TRUE, icon = icon("angle-right")))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "consensusReadMenu")})

    ### ------------------------------------------------------------------------
    ### Other features to add.
    ### ------------------------------------------------------------------------
    output$res <- renderText({
        paste("You've selected:", input$sidebar_menu)
    })


    output$clientdataText <- renderText({
        SangerSingleReadNum
    })


    output$consensusReadMenu_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (input$sidebar_menu == "consensusReadMenu") box("consensusReadMenu")
    })


    output$singelReadMenu_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, "_")
        if (!is.na(as.numeric(sidebar_menu[[1]]))) {
            fluidRow(
                # box(sidebar_menu[[1]])
                box(title = "Row File: ", solidHeader = TRUE, status = "primary", width = 12, h1(paste0(SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]]))),
                # h5(paste0("Absolute File Path: ",  SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]]))

                box(SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]]),
                box(SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]]),
            )
        }
    })
}
