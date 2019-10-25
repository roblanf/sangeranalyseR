### ============================================================================
### R shiny consensus read server function
### ============================================================================
consensusServer <- function(input, output, session) {
    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusRead <- getShinyOption("SangerConsensusReadSet")
    SangerConsensusRead <- SangerConsensusRead[[1]]
    SangerSingleReadNum <- length((SangerConsensusRead)@SangerReadsList)
    SangerSingleReadBFN <- sapply(1:SangerSingleReadNum, function(i)
               basename(SangerConsensusRead@SangerReadsList[[i]]@readFileName))
    SangerSingleReadAFN <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@readFileName)

    SangerSingleReadAFN <- sapply(1:SangerSingleReadNum, function(i)
        SangerConsensusRead@SangerReadsList[[i]]@readFileName)

    SangerSingleReadFeature <- sapply(1:SangerSingleReadNum, function(i)
        paste0(i, "_",
                   SangerConsensusRead@SangerReadsList[[i]]@readFeature))

    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar.
    ### ------------------------------------------------------------------------
    output$singleReadMenu <- renderMenu({
        menu_list <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(SangerSingleReadFeature[i], tabName = SangerSingleReadFeature[i], selected = TRUE, icon = icon("angle-right")))
        })
        sidebarMenu(.list = menu_list)
    })

    ### ------------------------------------------------------------------------
    ### Other features to add.
    ### ------------------------------------------------------------------------
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
            # box(sidebar_menu[[1]])
            h1(paste0("File: ",  SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]]))
            # h5(paste0("Absolute File Path: ",  SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]]))

            # box(SangerSingleReadAFN[[strtoi(sidebar_menu[[1]])]])
            # box(SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]])
        }
    })
}
