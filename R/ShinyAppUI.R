### ============================================================================
### R shiny consensus read UI function
### ============================================================================
consensusUI <- dashboardPage(
    skin = "green",
    dashboardHeader(
        title = "sangeranalyseR"
    ),
    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar. (Two layers)
    ### ------------------------------------------------------------------------
    dashboardSidebar(
        useShinyjs(debug = TRUE),
        sidebarMenu(
            id = "sidebar_menu",
            menuItem("Consensus Read", tabName = "Overview", icon=icon("dashboard")),
            sidebarMenuOutput("singleReadMenu")
        )
    ),
    ### ------------------------------------------------------------------------
    ### Others
    ### ------------------------------------------------------------------------
    dashboardBody(
        useShinyjs(debug = TRUE),

        # tabItems(
        #     id = "tab_items",
        #     # menuItem("Consensus Read", tabName = "Overview", icon=icon("dashboard")),
        #     sidebarMenuOutput("singleReadMenu")
        # )

        tags$script(HTML('
            $(document).ready(function() {
            $("header").find("nav").append(\'<span id="rightHeader" class="myClass"> Overview </span>\');
            })
        ')),
        tags$head(tags$style(HTML(
            '.myClass {
            font-size: 25px;
            line-height: 50px;
            text-align: left;
            font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
            padding: 0 15px;
            overflow: hidden;
            color: white;
            font-weight: bold;
        }'))),
        tags$head(
            tags$style(HTML('#closeUI{background-color:red; color:white; font-size: 23; font-weight: bold;}')),
            tags$style(HTML('#saveS4{background-color: #b3d9ff; font-size: 23; font-weight: bold;}'))
        ),
        tags$style(
            HTML(".shiny-notification {
             color: white;
             background-color:#0E8C3A;
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
        # textOutput("res"),

        box(status = "success", width = 12, htmlOutput("res"), actionButton('saveS4', 'Save S4 object'), actionButton('closeUI', 'Close UI')),
        uiOutput("consensusReadMenu_content"),
        uiOutput("singelReadMenu_content"),
        tags$head(tags$style(".sidebar-menu li { margin-bottom: 10px; }")),
        textOutput("selected_var"),
        h3("clientData values"),
        verbatimTextOutput("clientdataText")
    )
)
