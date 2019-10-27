### ============================================================================
### R shiny consensus read UI function
### ============================================================================
consensusUI <- dashboardPage(
    skin = "green",
    dashboardHeader(
        title = "Dynamic sidebar"
    ),
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
        tags$script(HTML('
            $(document).ready(function() {
            $("header").find("nav").append(\'<span id="pageHeader" class="myClass"> Dynamic Text </span>\');
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
        verbatimTextOutput("clientdataText"),
        useShinyjs()
    )
)
