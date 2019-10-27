### ============================================================================
### R shiny consensus read UI function
### ============================================================================
consensusUI <- dashboardPage(
    skin = "green",
    dashboardHeader(title = "Dynamic sidebar"),
    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar. (Two layers)
    ### ------------------------------------------------------------------------
    dashboardSidebar(
        sidebarMenu(
            id = "sidebar_menu",
            menuItem("Consensus Read", tabName = "consensusReadMenu", icon=icon("dashboard")),
            sidebarMenuOutput("singleReadMenu")
        )
    ),
    ### ------------------------------------------------------------------------
    ### Others
    ### ------------------------------------------------------------------------
    dashboardBody(
        tags$style(
            HTML(".shiny-notification {
             position:fixed;
             top: calc(50% - 150px);
             left: calc(50% - 150px);
             width: 500px;
             font-size: 24px;
             font-weight: bold;
             overflow-wrap: break-word;
             }
             "
            )
        ),
        textOutput("res"),
        actionButton('saveS4', 'Save S4 object'),
        tags$div(id = 'placeholder'),
        uiOutput("consensusReadMenu_content"),
        uiOutput("singelReadMenu_content"),
        tags$head(tags$style(".sidebar-menu li { margin-bottom: 10px; }")),
        textOutput("selected_var"),
        h3("clientData values"),
        verbatimTextOutput("clientdataText")
    )
)
