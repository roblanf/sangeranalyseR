### ============================================================================
### R shiny consensus read UI function
### ============================================================================
consensusUI <- dashboardPage(
    dashboardHeader(title = "Dynamic sidebar"),
    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar. (Two layers)
    ### ------------------------------------------------------------------------
    dashboardSidebar(
        sidebarMenuOutput("menu")
    ),

    ### ------------------------------------------------------------------------
    ### Others
    ### ------------------------------------------------------------------------
    dashboardBody(
        textOutput("selected_var"),
        h3("clientData values"),
        verbatimTextOutput("clientdataText")
    )
)
