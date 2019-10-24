# Define UI for app that draws a histogram ----
consensusUI <- dashboardPage(
    dashboardHeader(title = "sangeranalyseR"),
    dashboardSidebar(
        sidebarMenu(
            menuItemOutput("menuitem")
        )
    ),
    dashboardBody(
        # Boxes need to be put in a row (or column)
        textOutput("selected_var"),
        h3("clientData values"),
        verbatimTextOutput("clientdataText"),
        plotOutput("myplot"),
        fluidRow(
            box(plotOutput("plot1", height = 250)),

            box(
                title = "Controls",
                sliderInput("slider", "Number of observations:", 1, 100, 50)
            )
        )
    )
)
