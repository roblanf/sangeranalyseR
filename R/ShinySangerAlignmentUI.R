### ============================================================================
### R shiny SangerAlignment UI function
### ============================================================================
SangerAlignmentUI <- dashboardPage(
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
            menuItem("Contigs Alignment",
                     tabName = "Contigs Alignment Overview Page _",
                     icon=icon("dashboard")),
            sidebarMenuOutput("singleReadMenu")
        )
    ),
    ### ------------------------------------------------------------------------
    ### Others
    ### ------------------------------------------------------------------------
    dashboardBody(
        useShinyjs(debug = TRUE),
        ### --------------------------------------------------------------------
        ### Main Box style changing
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
        ### --------------------------------------------------------------------
        ### Second layer Box style changing
        ### --------------------------------------------------------------------
        tags$style(HTML("
                    .box.box-success{
                    border-bottom-color:#f4faf0;
                    border-left-color:#f4faf0;
                    border-right-color:#f4faf0;
                    background:#f4faf0
                    }
                    ")),

        ### --------------------------------------------------------------------
        ### Right navigation bar style changing
        ### --------------------------------------------------------------------
        tags$script(HTML('
            $(document).ready(function() {
            $("header").find("nav").append(\'<span id="rightHeader" class="myClass"> Contigs Alignment Overview </span>\');
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
            ### ----------------------------------------------------------------
            ### Close button style changing
            ### ----------------------------------------------------------------
            tags$style(HTML('#closeUI{background-color:red; color:white;
                            padding:15px; font-size: 35; font-weight: bold;}')),
            ### ----------------------------------------------------------------
            ### Save button style changing
            ### ----------------------------------------------------------------
            tags$style(HTML('#saveS4{background-color:#0083FF; color:white;
                            padding:15px; font-size: 35; font-weight: bold;}')),
            ### ----------------------------------------------------------------
            ### Overwrite .jexcel_content Height!!
            ### ----------------------------------------------------------------
            tags$style(HTML(".jexcel_content { overflow-y: auto;
                            height: 100% !important;}")),
            ### ----------------------------------------------------------------
            ### Trimming method selection text style changing
            ### ----------------------------------------------------------------
            tags$style(HTML("#TrimmingMethodSelectionOutput{font-size: 18px;
                                 margin-bottom: 30px;
                                 }")),
            ### ----------------------------------------------------------------
            ### Suppress error message in Shiny app
            ### ----------------------------------------------------------------
            tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
            )
        ),
        ### --------------------------------------------------------------------
        ### Pop-up notification style changing
        ### --------------------------------------------------------------------
        tags$style(
            HTML(".shiny-notification {
             color: white;
             background-color:#1b6940;
             opacity:1.0;
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

        ### --------------------------------------------------------------------
        ### Pop-up notification style changing
        ### --------------------------------------------------------------------
        tags$style(HTML(".sidebar-menu li a
                        [data-value='Sanger Consensus Read Overview']
                        { font-size: 18px; font-weight: bold }")),
        tags$style(HTML(".sidebar-menu ul li a
                        { font-size: 15px}")),

        tags$style(HTML(".nav-tabs-custom .nav-tabs li.active:hover a,
                        .nav-tabs-custom .nav-tabs li.active a
                        {background-color: transparent;
                        border-color: transparent;}")),
        tags$style(HTML(".nav-tabs-custom .nav-tabs li.active
                        {border-top-color: #5cb85c;}")),


        tags$style(HTML(".sidebar-menu li a[data-value='Contigs Alignment Overview']
                        { font-size: 18px; font-weight: bold }")),
        tags$style(HTML(".sidebar-menu ul li a
                        { font-size: 15px}")),
        uiOutput("aligned_contigSeq_content"),
    )
)
