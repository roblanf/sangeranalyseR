### ============================================================================
### R shiny alignedConsensusSet UI function
### ============================================================================
alignedConsensusSetUI <- dashboardPage(
    skin = "green",
    dashboardHeader(
        title = "sangeranalyseR",
        tags$li(class = "dropdown",
                actionButton('saveS4', 'Save S4 object'),
                actionButton('closeUI', 'Close UI'))
    ),
    ### ------------------------------------------------------------------------
    ### Adding dynamic menu to sidebar. (Two layers)
    ### ------------------------------------------------------------------------
    dashboardSidebar(
        useShinyjs(debug = TRUE),
        sidebarMenu(
            id = "sidebar_menu",
            menuItem("Aligned Consensus Set", tabName = "Sanger Aligned Consensus Set Overview", icon=icon("dashboard")),
            sidebarMenuOutput("singleReadMenu")
        )
    ),
    ### ------------------------------------------------------------------------
    ### Others
    ### ------------------------------------------------------------------------
    dashboardBody(
        useShinyjs(debug = TRUE),
        ### --------------------------------------------------------------------
        ### Box style changing
        ### --------------------------------------------------------------------
        tags$style(HTML("
                    .box.box-solid.box-success>.box-header {
                    }
                    .box.box-solid.box-success{
                    border-bottom-color:#f3ffe6;
                    border-left-color:#f3ffe6;
                    border-right-color:#f3ffe6;
                    border-top-color:#f3ffe6;
                    background:#f3ffe6
                    }
                    ")),
        tags$script(HTML('
            $(document).ready(function() {
            $("header").find("nav").append(\'<span id="rightHeader" class="myClass"> Sanger Aligned Consensus Set Overview </span>\');
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
            tags$style(HTML('#closeUI{background-color:red; color:white; padding:15px; font-size: 35; font-weight: bold;}')),
            tags$style(HTML('#saveS4{background-color:#0083FF; color:white; padding:15px; font-size: 35; font-weight: bold;}')),
            # tags$style(HTML(".fa { font-size: 30px; }"))
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
        tags$style(HTML(".sidebar-menu li a[data-value='Sanger Aligned Consensus Set Overview']
                        { font-size: 18px; font-weight: bold }")),
        tags$style(HTML(".sidebar-menu ul li a
                        { font-size: 15px}")),
        textOutput("res")
    )
)
