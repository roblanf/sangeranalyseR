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
    SangerSingleReadFiles <- sapply(1:SangerSingleReadNum, function(i)
        paste0(i, "_",
               basename(
                   SangerConsensusRead@SangerReadsList[[i]]@readFileName)))

    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar.
    ### ------------------------------------------------------------------------
    output$menu <- renderMenu({
        menu_list2 <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(SangerSingleReadFiles[i], tabName = "main", selected = TRUE))
        })
        sidebarMenu(.list = menu_list2)
    })

    ### ------------------------------------------------------------------------
    ### Other features to add.
    ### ------------------------------------------------------------------------
    output$clientdataText <- renderText({
        SangerSingleReadNum
    })
}
