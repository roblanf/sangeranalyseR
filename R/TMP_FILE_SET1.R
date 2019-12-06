#
# consensusParamSet[["consensusReadSCSet"]] <<-
#     SangerConsensusSet@consensusReadSCSet
# consensusParamSet[["alignmentSCSet"]] <<-
#     SangerConsensusSet@alignmentSCSet
# consensusParamSet[["alignmentTreeSCSet"]] <<-
#     SangerConsensusSet@alignmentTreeSCSet
#
# message(">>>>>>>> Inside 'Sanger Aligned Consensus Set Overview'")
# SCTrimmingMethod <-
#     SangerCSetParam[[1]]$SangerSingleReadFReadsList[[1]]@
#     QualityReport@TrimmingMethod
# if (SCTrimmingMethod == "M1") {
#     SCTrimmingMethodName = "Method 1: 'Logarithmic Scale Trimming'"
# } else if (SCTrimmingMethod == "M2") {
#     SCTrimmingMethodName =
#         "Method 2: 'Logarithmic Scale Sliding Window Trimming'"
# }
# fluidRow(
#     useShinyjs(),
#     box(title = tags$p("Input Parameters: ",
#                        style = "font-size: 26px;
#                                        font-weight: bold;"),
#         solidHeader = TRUE, collapsible = TRUE,
#         status = "success", width = 12,
#         tags$hr(style = ("border-top: 0.2px hidden #A9A9A9;")),
#         fluidRow(
#             column(width = 12,
#                    actionBttn("recalculateButtonSCSet",
#                               "Re-calculate
#                                           consensusread (read set)",
#                               icon = icon("calculator"),
#                               style = "simple", color = "danger",
#                               block = TRUE, size = "lg")
#             ),
#             column(12,
#                    tags$hr(
#                        style = ("border-top: 2px hidden #A9A9A9;"))
#             ),
#             column(12,
#                    column(3,
#                           tags$p(tagList(icon("caret-right"),
#                                          "Output Directory: "),
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
#                                       "Raw ABI Parent Directory:"),
#                               style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SangerConsensusSet@parentDirectory),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(tagList(icon("caret-right"),
#                                          "Trimming Method: "),
#                                  style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SCTrimmingMethodName),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(tagList(icon("caret-right"),
#                                          "Forward Suffix RegExp: "),
#                                  style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SangerConsensusSet@suffixForwardRegExp),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(tagList(icon("caret-right"),
#                                          "Reverse Suffix RegExp: "),
#                                  style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SangerConsensusSet@suffixReverseRegExp),
#                    )
#             ),
#             column(12,
#                    column(3,
#                           tags$p(tagList(icon("caret-right"),
#                                          "Consensus Read Number: "),
#                                  style = "font-size: 20px;
#                                        font-weight: bold;"),
#                    ),
#                    column(9,
#                           h4(SangerConsensusSetNum),
#                    )
#             ),
#             box(title = tags$p("Genetic Code Data Frame",
#                                style = "font-size: 24px;
#                                        font-weight: bold;"),
#                 collapsible = TRUE,
#                 status = "success", width = 12,
#                 column(width = 2,
#                        tags$p("Tri-nucleotide:",
#                               style = "font-size: 15px;
#                                        font-weight: bold;"),
#                        tags$p("Amino Acid : ",
#                               style = "font-size: 15px;
#                                        font-weight: bold;"),
#                        tags$p("('*' : stop codon) ",
#                               style = "font-size: 12px;
#                                        font-weight: italic;"),
#                 ),
#                 column(width = 10,
#                        excelOutput("geneticCodeDF",
#                                    width = "100%", height = "50"),
#                        style = paste("height:100%; ",
#                                      "overflow-y: hidden;",
#                                      "overflow-x: scroll;")
#                 ),
#             ),
#             uiOutput("SCrefAminoAcidSeq") ,
#         ),
#     ),
#
#
#     box(title = tags$p(tagList(icon("dot-circle"),
#                                "Consensus Readset Results: "),
#                        style = "font-size: 26px;
#                                        font-weight: bold;"),
#         solidHeader = TRUE, collapsible = TRUE,
#         status = "success", width = 12,
#         tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
#         box(title = tags$p("Consensus Reads Alignment",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    htmlOutput("consensusSetAlignmentHTML"),
#             ),
#         ),
#         box(title = tags$p("Consensus Reads Tree",
#                            style = "font-size: 24px;
#                                        font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 12,
#                    plotOutput("SCSetConsensusReadTreePlot"),
#                    style = paste("height:100%; overflow-y:",
#                                  "scroll;overflow-x: scroll;")
#             )
#         ),
#     ),
# )
