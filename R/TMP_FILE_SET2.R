# consensusReadIndex <- strtoi(sidebar_menu[[1]])
# forwardReadNum <-
#     SangerCSetParam[[consensusReadIndex]]$forwardReadNum
# reverseReadNum <-
#     SangerCSetParam[[consensusReadIndex]]$reverseReadNum
# SCTrimmingMethod <-
#     SangerCSetParam[[consensusReadIndex]]$
#     SangerSingleReadFReadsList[[1]]@
#     QualityReport@TrimmingMethod
#
# if (SCTrimmingMethod == "M1") {
#     SCTrimmingMethodName =
#         "Method 1: 'Logarithmic Scale Trimming'"
# } else if (SCTrimmingMethod == "M2") {
#     SCTrimmingMethodName =
#         "Method 2: 'Logarithmic Scale Sliding Window Trimming'"
# }
# consensusParam[["consensusReadName"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@consensusReadName
# consensusParam[["consensusRead"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@consensusRead
# consensusParam[["differencesDF"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@differencesDF
# consensusParam[["alignment"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@alignment
# consensusParam[["distanceMatrix"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@distanceMatrix
# consensusParam[["dendrogram"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@dendrogram
# consensusParam[["indelsDF"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@indelsDF
# consensusParam[["stopCodonsDF"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@stopCodonsDF
# consensusParam[["secondaryPeakDF"]] <<-
#     SangerConsensusSet@
#     consensusReadsList[[consensusReadIndex]]@secondaryPeakDF
#
# fluidRow(
#     useShinyjs(),
#     box(title = tags$p(tagList(icon("dot-circle"),
#                                "Basic Information: "),
#                        style = "font-size: 26px;
#                                        font-weight: bold;"),
#         solidHeader = TRUE, collapsible = TRUE,
#         status = "success", width = 12,
#         tags$hr(
#             style = ("border-top: 0.2px hidden #A9A9A9;")),
#         fluidRow(
#             column(width = 12,
#                    actionBttn("recalculateButton",
#                               "Re-calculate consensus read",
#                               icon = icon("calculator"),
#                               style = "simple",
#                               color = "danger",
#                               block = TRUE, size = "lg")
#             ),
#             column(12,
#                    tags$hr(
#                        style = ("border-top: 2px
#                                                     hidden #A9A9A9;")),
#             ),
#             column(12,
#                    column(3,
#                           tags$p(tagList(icon("caret-right"),
#                                          "Output Directory:"),
#                                  style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(shinyDirectory),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(
#                               tagList(icon("caret-right"),
#                                       "Raw ABI Parent Directory: "),
#                               style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SangerConsensusSet@
#                                  parentDirectory),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(
#                               tagList(icon("caret-right"),
#                                       "Consenesus Read Name:"),
#                               style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(consensusParam[["consensusReadName"]]),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(
#                               tagList(icon("caret-right"),
#                                       "Trimming Method: "),
#                               style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SCTrimmingMethodName),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(
#                               tagList(icon("caret-right"),
#                                       "Forward Suffix RegExp:"),
#                               style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SangerConsensusSet@
#                                  suffixForwardRegExp),
#                    )
#
#             ),
#             column(12,
#                    column(3,
#                           tags$p(
#                               tagList(icon("caret-right"),
#                                       "Forward Read Number:"),
#                               style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(forwardReadNum),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(
#                               tagList(icon("caret-right"),
#                                       "Reverse Suffix RegExp:"),
#                               style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SangerConsensusSet@
#                                  suffixReverseRegExp),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(
#                               tagList(icon("caret-right"),
#                                       "Reverse Read Number:"),
#                               style = "font-size: 20px;
#                                               font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(reverseReadNum),
#                    )
#             ),
#         ),
#         ################################################
#         #### Add this after having reference sample ####
#         ################################################
#         tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
#         box(title = tags$p("Consensus Read Parameters",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(3,
#                    uiOutput("SCMinReadsNum") ,
#             ),
#             column(3,
#                    uiOutput("SCMinReadLength")  ,
#             ),
#             column(3,
#                    uiOutput("SCMinFractionCall") ,
#             ),
#             column(3,
#                    uiOutput("SCMaxFractionLost") ,
#             ),
#             column(3,
#                    uiOutput("SCAcceptStopCodons") ,
#             ),
#             column(3,
#                    uiOutput("SCReadingFrame") ,
#             ),
#         ),
#         box(title = tags$p("Genetic Code Data Frame",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 2,
#                    tags$p("Tri-nucleotide:",
#                           style = "font-size: 15px;
#                                        font-weight: bold;"),
#                    tags$p("Amino Acid : ",
#                           style = "font-size: 15px;
#                                        font-weight: bold;"),
#                    tags$p("('*' : stop codon) ",
#                           style = "font-size: 12px;
#                                        font-weight: italic;"),
#             ),
#             column(width = 10,
#                    excelOutput("geneticCodeDF",
#                                width = "100%",
#                                height = "50"),
#                    style = paste("height:100%; ",
#                                  "overflow-y: hidden;",
#                                  "overflow-x: scroll;")
#             ),
#         ),
#         uiOutput("SCrefAminoAcidSeq") ,
#     ),
#
#     box(title = tags$p(tagList(icon("dot-circle"),
#                                "Consensus Read Results: "),
#                        style = "font-size: 26px;
#                                        font-weight: bold;"),
#         solidHeader = TRUE, collapsible = TRUE,
#         status = "success", width = 12,
#         tags$hr(style =("border-top: 4px hidden #A9A9A9;")),
#         box(title = tags$p("Alignment",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    htmlOutput("consensusAlignmentHTML"),
#             ),
#         ),
#         box(title = tags$p("Differences Data frame",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    uiOutput("SCDifferencesDFUI"),
#                    style =
#                        paste("height:100%; overflow-y:",
#                              "scroll;overflow-x: scroll;")
#             )
#         ),
#         box(title = tags$p("Dendrogram",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    plotOutput("dendrogramPlot"),
#                    style =
#                        paste("height:100%; overflow-y:",
#                              "scroll;overflow-x: scroll;")
#             ),
#             column(width = 12,
#                    tags$hr(
#                        style =("border-top: 4px hidden #A9A9A9;")),
#             ),
#             column(width = 12,
#                    dataTableOutput("dendrogramDF"),
#                    style =
#                        paste("height:100%; overflow-y:",
#                              "scroll;overflow-x: scroll;")
#             )
#         ),
#         box(title = tags$p("Samples Distance",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    # plot()
#                    uiOutput("SCDistanceMatrixPlotUI"),
#                    style =
#                        paste("height:100%; overflow-y:",
#                              "scroll;overflow-x: scroll;")
#             ),
#             column(width = 12,
#                    tags$hr(
#                        style =("border-top: 4px hidden #A9A9A9;")),
#             ),
#             column(width = 12,
#                    uiOutput("SCDistanceMatrixUI"),
#                    style =
#                        paste("height:100%; overflow-y:",
#                              "scroll;overflow-x: scroll;")
#             )
#         ),
#         box(title = tags$p("Indels Data frame",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    uiOutput("SCIndelsDFUI"),
#                    style =
#                        paste("height:100%; overflow-y:",
#                              "scroll;overflow-x: scroll;")
#             )
#         ),
#         box(title = tags$p("Stop Codons Data frame",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    uiOutput("SCStopCodonsDFUI"),
#                    style =
#                        paste("height:100%; overflow-y:",
#                              "scroll;overflow-x: scroll;")
#             )
#         )
#     ),
# )
