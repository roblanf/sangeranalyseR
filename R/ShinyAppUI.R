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
        textOutput("res"),
        uiOutput("consensusReadMenu_content"),
        uiOutput("singelReadMenu_content"),
        tags$head(tags$style(".sidebar-menu li { margin-bottom: 10px; }")),
        textOutput("selected_var"),
        h3("clientData values"),
        verbatimTextOutput("clientdataText")
    )
)
