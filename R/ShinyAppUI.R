# Define UI for app that draws a histogram ----
consensusUI <- dashboardPage(
    dashboardHeader(title = "Dynamic sidebar"),
    dashboardSidebar(
        sidebarMenuOutput("menu")
    ),
    dashboardBody(
        textOutput("selected_var"),
        h3("clientData values"),
        verbatimTextOutput("clientdataText")
    )
)

#
#         dashboardHeader(),
#         dashboardSidebar(
#             sidebarMenuOutput("menu"),
#             textOutput("res")
#         ),
#         dashboardBody(
#             tabItems(
#                 tabItem("dashboard", "Dashboard tab content"),
#                 tabItem("widgets", "Widgets tab content"),
#                 tabItem("subitem1", "Sub-item 1 tab content"),
#                 tabItem("subitem2", "Sub-item 2 tab content")
#             )
#         )
#     )





#     dashboardHeader(title = "sangeranalyseR"),
#     dashboardSidebar(
#         # tags$div(id = 'placeholder'),
#         sidebarMenu(
#             id = 'placeholder'),
#         sidebarMenu(
#             menuItemOutput("menuitem")
#         )
#     ),
#     dashboardBody(
#         # Boxes need to be put in a row (or column)
#         dashboardHeader(dropdownMenuOutput("messageMenu")),
#
#
#
#
#
#         textOutput("selected_var"),
#         h3("clientData values"),
#         verbatimTextOutput("clientdataText"),
#         plotOutput("myplot"),
#         fluidRow(
#             box(plotOutput("plot1", height = 250)),
#
#             box(
#                 title = "Controls",
#                 sliderInput("slider", "Number of observations:", 1, 100, 50)
#             )
#         )
#     )
# )
