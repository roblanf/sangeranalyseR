# Define UI for app that draws a histogram ----
ui <- fluidPage(

    # App title ----
    titlePanel("Welcome to sangeranalyseR !"),

    ### This is just example code for now !!
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        # Sidebar panel for inputs ----
        sidebarPanel(
            # Input: Slider for the number of bins ----
            sliderInput(inputId = "bins",
                        label = "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Main panel for displaying outputs ----
        mainPanel(
            # Output: Histogram ----
            plotOutput(outputId = "distPlot")
        )
    )
)
