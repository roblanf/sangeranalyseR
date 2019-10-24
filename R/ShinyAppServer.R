consensusServer <- function(input, output, session) {
    SangerConsensusRead <- getShinyOption("SangerConsensusReadSet")
    SangerSingleReadNum <- length((SangerConsensusRead[[1]])@SangerReadsList)

    output$clientdataText <- renderText({
        SangerSingleReadNum
    })

    # output$menu <- renderMenu({
    #     for (iterater in c(1:SangerSingleReadNum)) {
    #         sidebarMenu(
    #             menuItem("Menu item", icon = icon("calendar"))
    #         )
    #     }
    # })
    output$menuitem <- renderMenu({
        menuItem("Menu item", icon = icon("calendar"))
    })

    # # A histogram
    # output$myplot <- renderPlot({
    #     hist(rnorm(input$obs), main = "Generated in renderPlot()")
    # })

    set.seed(122)
    histdata <- rnorm(500)

    # output$selected_var <- renderText(cdata)

    output$plot1 <- renderPlot({
        data <- histdata[seq_len(input$slider)]
        hist(data)
    })
}
