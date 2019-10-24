consensusServer <- function(input, output, session) {
    SangerConsensusRead <- getShinyOption("SangerConsensusReadSet")
    SangerConsensusRead <- SangerConsensusRead[[1]]
    SangerSingleReadNum <- length((SangerConsensusRead)@SangerReadsList)
    SangerSingleReadFiles <- sapply(1:SangerSingleReadNum, function(i)
        paste0(i, "_",
               basename(
                   SangerConsensusRead@SangerReadsList[[i]]@readFileName)))

    output$menu <- renderMenu({
        menu_list2 <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(SangerSingleReadFiles[i], tabName = "main", selected = TRUE))
        })
        sidebarMenu(.list = menu_list2)
    })
    output$clientdataText <- renderText({
        SangerSingleReadNum
    })
}




    #     output$res <- renderText({
    #         paste("You've selected:", input$tabs)
    #     })
    #     output$menu <- renderMenu({
    #         sidebarMenu(
    #             # Setting id makes input$tabs give the tabName of currently-selected tab
    #             id = "tabs",
    #
    #             menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    #             menuItem("Widgets", icon = icon("th"), tabName = "widgets"),
    #             menuItem("Charts", icon = icon("bar-chart-o"),
    #                      menuSubItem("Sub-item 1", tabName = "subitem1"),
    #                      menuSubItem("Sub-item 2", tabName = "subitem2")
    #             )
    #         )
    #     })
    # }







#     SangerConsensusRead <- getShinyOption("SangerConsensusReadSet")
#     SangerSingleReadNum <- length((SangerConsensusRead[[1]])@SangerReadsList)
#
#     output$clientdataText <- renderText({
#         SangerSingleReadNum
#     })
#
#
#     inserted <- c()
#
#     # for (iterater in c(1:SangerSingleReadNum)) {
#     #     insertUI(
#     #         selector = '#placeholder',
#     #         ## wrap element in a div with id for ease of removal
#     #         ui = menuItem("Menu item", icon = icon("calendar"))
#     #     )
#     #     inserted <<- c(iterater, inserted)
#     # }
#
#
#
#
#     # output$menuitem <- renderMenu({
#     #     # Code to generate each of the messageItems here, in a list. This assumes
#     #     # that messageData is a data frame with two columns, 'from' and 'message'.
#     #     msgs <-
#     #         c(1:SangerSingleReadNum)
#     #
#     #     # This is equivalent to calling:
#     #     #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
#     #     dropdownMenu(type = "messages", .list = msgs)
#     # })
#
# #
# #     output$messageMenu <- renderMenu({
# #         # Code to generate each of the messageItems here, in a list. This assumes
# #         # that messageData is a data frame with two columns, 'from' and 'message'.
# #         msgs <- apply(messageData, 1, function(row) {
# #             messageItem(from = row[["from"]], message = row[["message"]])
# #         })
# #
# #         # This is equivalent to calling:
# #         #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
# #         dropdownMenu(type = "messages", .list = msgs)
# #     })
#
#
#
#
#
#
#     # # A histogram
#     # output$myplot <- renderPlot({
#     #     hist(rnorm(input$obs), main = "Generated in renderPlot()")
#     # })
#
#     set.seed(122)
#     histdata <- rnorm(500)
#
#     # output$selected_var <- renderText(cdata)
#
#     output$plot1 <- renderPlot({
#         data <- histdata[seq_len(input$slider)]
#         hist(data)
#     })
# }

#
#
# m1 <- matrix(C<-(1:10),nrow=5, ncol=6)
# m1
# a_m1 <- apply(m1, 2, sum)
# a_m1
