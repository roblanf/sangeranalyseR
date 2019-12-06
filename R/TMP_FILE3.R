#         singleReadIndex <- strtoi(sidebar_menu[[1]])
#
#         sequenceParam[["primarySeq"]] <<-
#             SangerSingleReadPrimSeqDF[[singleReadIndex]]
#         sequenceParam[["secondarySeq"]] <<-
#             SangerSingleReadSecoSeqDF[[singleReadIndex]]
#         sequenceParam[["primaryAASeqS1"]] <<-
#             SangerSingleReadPrimAASeqS1DF[[singleReadIndex]]
#         sequenceParam[["primaryAASeqS2"]] <<-
#             SangerSingleReadPrimAASeqS2DF[[singleReadIndex]]
#         sequenceParam[["primaryAASeqS3"]] <<-
#             SangerSingleReadPrimAASeqS3DF[[singleReadIndex]]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#         ChromatogramParam[["baseNumPerRow"]] <<-
#             SangerSingleReadChromatogramParam[[singleReadIndex]]@
#             baseNumPerRow
#         ChromatogramParam[["heightPerRow"]] <<-
#             SangerSingleReadChromatogramParam[[singleReadIndex]]@
#             heightPerRow
#         ChromatogramParam[["signalRatioCutoff"]] <<-
#             SangerSingleReadChromatogramParam[[singleReadIndex]]@
#             signalRatioCutoff
#         ChromatogramParam[["showTrimmed"]] <<-
#             SangerSingleReadChromatogramParam[[singleReadIndex]]@
#             showTrimmed
#
#         trimmedParam[["M1TrimmingCutoff"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             M1TrimmingCutoff
#         trimmedParam[["M2CutoffQualityScore"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             M2CutoffQualityScore
#         trimmedParam[["M2SlidingWindowSize"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             M2SlidingWindowSize
#
#
#         trimmedRV[["rawSeqLength"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             rawSeqLength
#         trimmedRV[["rawMeanQualityScore"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             rawMeanQualityScore
#         trimmedRV[["rawMinQualityScore"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             rawMinQualityScore
#         trimmedRV[["trimmedStartPos"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedStartPos
#         trimmedRV[["trimmedFinishPos"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedFinishPos
#         trimmedRV[["trimmedSeqLength"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedSeqLength
#         trimmedRV[["trimmedMeanQualityScore"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMeanQualityScore
#         trimmedRV[["trimmedMinQualityScore"]] <<-
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMinQualityScore
#         trimmedRV[["remainingRatio"]] <<-
#             round(SangerSingleReadQualReport[[singleReadIndex]]@
#                       remainingRatio * 100, 2)
#         fluidRow(
#             useShinyjs(),
#             box(title = tags$p(tagList(icon("dot-circle"),
#                                        "Raw File: "),
#                                style = "font-size: 26px;
#                                  font-weight: bold;"),
#                 solidHeader = TRUE,
#                 status = "success", width = 12,
#                 h1(paste0(
#                     SangerSingleReadBFN[[strtoi(singleReadIndex)]])),
#                 tags$h5(paste("( full path:",
#                               SangerSingleReadAFN[[
#                                   strtoi(singleReadIndex)]],
#                               ")"), style = "font-style:italic")),
#             box(title = tags$p(tagList(icon("dot-circle"),
#                                        "Primary, Secondary DNA Sequences & Amino Acid Sequence (Before Trimming):"),
#                                style = "font-size: 26px;
#                                font-weight: bold;"),
#                 solidHeader = TRUE, collapsible = TRUE,
#                 status = "success", width = 12,
#                 tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
#                 column(width = 12,
#                        tags$p(tagList(icon("bars"),
#                                       "Primary Sequence"),
#                               style = "font-size: 22px;
#                                font-weight: bold;"),
#                        excelOutput("primarySeqDF",
#                                    width = "100%", height = "50"),
#                        tags$br(),
#                        tags$br(),
#                        tags$p(tagList(icon("bars"),
#                                       "Secondary Sequence"),
#                               style = "font-size: 22px;
#                                font-weight: bold;"),
#                        excelOutput("secondSeqDF",
#                                    width = "100%", height = "50"),
#                        tags$br(),
#                        tags$br(),
#                        tags$p(tagList(icon("bars"),
#                                       "Quality Phred Score"),
#                               style = "font-size: 22px;
#                                font-weight: bold;"),
#                        excelOutput("qualityScoreDF",
#                                    width = "100%", height = "50"),
#                        tags$br(),
#                        tags$br(),
#                        tags$p(tagList(icon("bars"),
#                                       "AA Sequence 1"),
#                               style = "font-size: 22px;
#                                font-weight: bold;"),
#                        excelOutput("PrimAASeqS1DF",
#                                    width = "100%", height = "50"),
#                        tags$br(),
#                        tags$br(),
#                        tags$p(tagList(icon("bars"),
#                                       "AA Sequence 2"),
#                               style = "font-size: 22px;
#                                font-weight: bold;"),
#                        excelOutput("PrimAASeqS2DF",
#                                    width = "100%", height = "50"),
#                        tags$br(),
#                        tags$br(),
#                        tags$p(tagList(icon("bars"),
#                                       "AA Sequence 3"),
#                               style = "font-size: 22px;
#                                font-weight: bold;"),
#                        excelOutput("PrimAASeqS3DF",
#                                    width = "100%", height = "50"),
#                        style = paste("overflow-y: hidden;",
#                                      "overflow-x: scroll;")
#                 ),
#             ),
#             box(title = tags$p(tagList(icon("dot-circle"),
#                                        "Quality Report: "),
#                                style = "font-size: 26px;
#                                font-weight: bold;"),
#                 solidHeader = TRUE, collapsible = TRUE,
#                 status = "success", width = 12,
#                 tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
#                 box(title = tags$p(tagList(icon("arrow-circle-right"),
#                                            "Trimming Parameters Input"),
#                                    style = "font-size: 24px;
#                                font-weight: bold;"),
#                     collapsible = TRUE,
#                     status = "success", width = 12,
#                     fluidRow(
#                         column(width = 12,
#                               uiOutput("TrimmingMethodSelectionOutput"),
#                             ),
#                         ),
#                         column(width = 12,
#                                uiOutput("TrimmingMethodUI"),
#                         ),
#                     actionBttn("startTrimmingButton",
#                                "Apply Trimming Parameters",
#                                style = "simple", color = "success",
#                                block = TRUE, size = "lg")
#                     ),
#                 box(title = tags$p(tagList(icon("arrow-circle-left"),
#                                            "Trimmed Result Output"),
#                                    style = "font-size: 24px;
#                                font-weight: bold;"),
#                     collapsible = TRUE,
#                     status = "success", width = 12,
#                     fluidRow(
#                         box(title = tags$p("Before Trimming",
#                                            style = "font-size: 21px;
#                                font-weight: bold;"),
#                             collapsible = TRUE,
#                             status = "success", width = 12,
#                             column(width = 12,
#                                    column(4,
#                                           uiOutput("rawSeqLength") ,
#                                    ),
#                                    column(4,
#                                         uiOutput("rawMeanQualityScore"),
#                                    ),
#                                    column(4,
#                                         uiOutput("rawMinQualityScore"),
#                                    ),
#                             ),
#                         ),
#                     ),
#                     fluidRow(
#                         box(title = tags$p("After Trimming",
#                                            style = "font-size: 21px;
#                                font-weight: bold;"),
#                             collapsible = TRUE,
#                             status = "success", width = 12,
#                             column(width = 12,
#                                    column(4,
#                                           uiOutput("trimmedSeqLength"),
#                                    ),
#                                    column(4,
#                                           uiOutput("trimmedMeanQualityScore"),
#                                    ),
#                                    column(4,
#                                           uiOutput("trimmedMinQualityScore"),
#                                    ),
#                             ),
#
#                             column(width = 12,
#                                    column(4,
#                                           uiOutput("trimmedStartPos") ,
#                                    ),
#                                    column(4,
#                                           uiOutput("trimmedFinishPos") ,
#                                    ),
#                                    column(4,
#                                           uiOutput("remainingRatio") ,
#                                    )
#                             ),
#                         ),
#                     ),
#                     tags$hr(
#                         style = ("border-top: 6px double #A9A9A9;")),
#                     fluidRow(
#                         box(title = tags$p("Cumulative Ratio Plot",
#                                            style = "font-size: 21px;
#                                font-weight: bold;"),
#                             collapsible = TRUE,
#                             status = "success", width = 12,
#                             plotlyOutput("qualityTrimmingRatioPlot") %>%
#                                 withSpinner()),
#                         box(title = tags$p("Cumulative Ratio Plot",
#                                            style = "font-size: 21px;
#                                font-weight: bold;"),
#                             collapsible = TRUE,
#                             status = "success", width = 12,
#                             plotlyOutput("qualityQualityBasePlot") %>%
#                                 withSpinner()),
#                     ),
#                 ),
#             ),
#             box(title = tags$p(tagList(icon("dot-circle"),
#                                        "Chromatogram: "),
#                                style = "font-size: 26px;
#                                font-weight: bold;"),
#                 solidHeader = TRUE, collapsible = TRUE,
#                 status = "success", width = 12,
#                 tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
#
#                 box(title = tags$p(tagList(icon("arrow-circle-right"),
#                                            "Chromatogram Input"),
#                                    style = "font-size: 24px;
#                                font-weight: bold;"),
#                     collapsible = TRUE,
#                     status = "success", width = 12,
#                     column(3,
#                            sliderInput("ChromatogramBasePerRow",
#                                        label =h4("Base Number Per Row"),
#                                        min = 5,
#                                        max = 200,
#                                        value = ChromatogramParam[["baseNumPerRow"]]),
#                            sliderInput("ChromatogramHeightPerRow",
#                                        label = h4("Height Per Row"),
#                                        min = 50,
#                                        max = 600,
#                                        value = ChromatogramParam[["heightPerRow"]]),
#                     ),
#                     column(3,
#                            tags$hr(
#                                style =
#                                    ("border-top: 4px hidden #A9A9A9;")),
#                            numericInput(
#                                "ChromatogramSignalRatioCutoff",
#                                h3("Signal Ratio Cutoff"),
#                                value = ChromatogramParam[["signalRatioCutoff"]]),
#                            checkboxInput(
#                                "ChromatogramCheckShowTrimmed",
#                                "Whether show trimmed region",
#                                value =
#                                    ChromatogramParam[["showTrimmed"]])
#                     ),
#                     column(3,
#                            tags$hr(
#                               style=("border-top: 4px hidden #A9A9A9;")),
#                            uiOutput("ChromatogramtrimmedStartPos"),
#                     ),
#                     column(3,
#                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
#                            uiOutput("ChromatogramtrimmedFinishPos"),
#                     ),
#                     actionBttn("saveChromatogramParam",
#                                "Apply Chromatogram Parameters",
#                                style = "simple", color = "success",
#                                block = TRUE, size = "lg")
#                 ),
#                 box(title = tags$p(tagList(icon("arrow-circle-left"),
#                                            "Chromatogram Output"),
#                                    style = "font-size: 24px;
#                                font-weight: bold;"),
#                     collapsible = TRUE,
#                     status = "success", width = 12,
#                     column(width = 12,
#                            uiOutput("chromatogramUIOutput"),
#                     )
#                 ),
#             )
#         )
