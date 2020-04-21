### ============================================================================
### Adding dynamic menu to sidebar.
### ============================================================================
dynamicMenuSideBarSC <- function(input, output, session,
                                  forwardReadNum, reverseReadNum,
                                  forwardReadFeature, reverseReadFeature) {
    output$singleReadMenu <- renderMenu({
        fmenuSub_list <- lapply(seq_len(forwardReadNum), function(i) {
            list(menuSubItem(text = paste(strsplit(forwardReadFeature[i], " ")[[1]][1],
                                          "Forward SangerRead"),
                          tabName = forwardReadFeature[i], icon = icon("minus")))
        })
        rmenuSub_list <- lapply(seq_len(reverseReadNum), function(i) {
            list(menuSubItem(text = paste(strsplit(forwardReadFeature[i], " ")[[1]][1],
                                          "Reverse SangerRead"),
                          tabName = reverseReadFeature[i], icon = icon("minus")))
        })
        fmenu_list <- menuItem(text = tags$p(tagList(icon("circle"),
                                                     HTML('&nbsp;'),
                                                     "Forward SangerReads"),
                                             style = "font-size: 15px;
                                           font-weight: bold;"),
                               tabName = "forwardReads", fmenuSub_list)

        rmenu_list <- menuItem(text = tags$p(tagList(icon("circle"),
                                                     HTML('&nbsp;'),
                                                     "Reverse SangerReads"),
                                             style = "font-size: 15px;
                                           font-weight: bold;"),
                               tabName = "reverseReads", rmenuSub_list)
        sidebarMenu(.list = list(fmenu_list, rmenu_list))
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu",
                            "Sanger Contig Overview")})
}

dynamicMenuSideBarSA <- function(input, output, session, SangerAlignmentParam) {
    output$singleReadMenu <- renderMenu({
        SangerCSNum <- length(SangerAlignmentParam)
        menu_list <- lapply(seq_len(SangerCSNum), function(i) {
            forwardReadNum <- SangerAlignmentParam[[i]]$forwardReadNum
            reverseReadNum <- SangerAlignmentParam[[i]]$reverseReadNum

            forwardReadFeature <- SangerAlignmentParam[[i]]$forwardReadFeature
            reverseReadFeature <- SangerAlignmentParam[[i]]$reverseReadFeature
            fmenuSub_list <- lapply(seq_len(forwardReadNum), function(j) {
                list(menuSubItem(text = paste(i, "-",
                                              strsplit(forwardReadFeature[j], " ")[[1]][1],
                                              "Forward SangerRead"),
                                 tabName = paste(i, "Contig -", forwardReadFeature[j]),
                                 icon = icon("minus")))
            })
            rmenuSub_list <- lapply(seq_len(reverseReadNum), function(j) {
                list(menuSubItem(text = paste(i, "-",
                                              strsplit(reverseReadFeature[j], " ")[[1]][1],
                                              "Reverse SangerRead"),
                                 tabName = paste(i, "Contig -", reverseReadFeature[j]),
                                 icon = icon("minus")))
            })
            fmenu_list <- menuItem(text = tags$p(tagList(HTML('&nbsp;'),
                                                         HTML('&nbsp;'),
                                                         icon("circle"),
                                                         HTML('&nbsp;'),
                                                         "Forward SangerReads"),
                                                         style = "font-size: 14px;
                                                 font-weight: bold;"),
                                   tabName = "forwardReads", fmenuSub_list)

            rmenu_list <- menuItem(text = tags$p(tagList(HTML('&nbsp;'),
                                                         HTML('&nbsp;'),
                                                         icon("circle"),
                                                         HTML('&nbsp;'),
                                                         "Reverse SangerReads"),
                                                         style = "font-size: 14px;
                                                 font-weight: bold;"),
                                   tabName = "reverseReads",
                                   rmenuSub_list)
            SangerCSMenuSubItem <- list(fmenu_list, rmenu_list)

            SangerCSMenuSubItem <- c(list(menuSubItem(text = paste(i, "SangerContig Overview"),
                                                      tabName = paste(i, "Sanger Contig Overview Page"),
                                                      icon = icon("align-left"))),
                                     SangerCSMenuSubItem)
            SangerAlignmentParam[[i]]$SCName
            list(menuItem(text = tags$p(tagList(icon("align-left"), i, "SangerContig"),
                                        style = "font-size: 17px;
                                           font-weight: bold;"),
                          tabName = SangerAlignmentParam[[i]]$SCName,
                          SangerCSMenuSubItem))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu",
                            "Contigs Alignment Overview Page .")})
}

### ============================================================================
### observeEvent: Adding dynamic rightHeader text
### ============================================================================
observeEventDynamicHeaderSC <- function(input, output, session, trimmedRV) {
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        log_info("menuItem: ", menuItem)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
    })
}

# observeEventDynamicHeaderSA <- function(input, output, session, trimmedRV,
#                                            SangerAlignmentParam) {
#
#     output$res <- renderText({
#         paste("You've selected:", input$sidebar_menu)
#     })
#
#     observeEvent(input$sidebar_menu, {
#         menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
#         html("rightHeader", menuItem)
#         sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#         # log_info("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
#         if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
#         #     trimmedRV[["trimmedStartPos"]] <-
#         #         SangerReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedStartPos
#         #     trimmedRV[["trimmedFinishPos"]] <-
#         #         SangerReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
#         #     qualityPhredScores = SangerReadQualReport[[
#         #         strtoi(sidebar_menu[[1]])]]@qualityPhredScores
#         #
#         #     readLen = length(qualityPhredScores)
#         #     trimmedRV[["trimmedSeqLength"]] <- trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1
#         #     trimmedRV[["remainingRatio"]] <- round(((trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1) / readLen) * 100, digits = 2)
#         }
#     })
# }


### ============================================================================
### valueBox for SangerAlignment
### ============================================================================
valueBoxSAMinFractionCallSA <- function(input, output,
                                      SAMinFractionCall, session) {
    output$SAMinFractionCallSA <- renderUI({
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SAMinFractionCall),
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}
valueBoxSAMaxFractionLostSA <- function(input, output,
                                      SAMaxFractionLost, session) {
    output$SAMaxFractionLostSA <- renderUI({
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                      font-weight: bold;"),
            value = tags$p(as.numeric(SAMaxFractionLost),
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

valueBoxSAMinReadsNum <- function(input, output, SangerConsensusSet, session) {
    output$SAMinReadsNum <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMinReadsNum <- SangerConsensusSet@contigList[[contigIndex]]@minReadsNum
            valueBox(
                subtitle = tags$p("MinReadsNum",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(SAMinReadsNum),
                               style = "font-size: 29px;"),
                icon = icon("sliders-h", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}
valueBoxSAMinReadLength <- function(input, output, SangerConsensusSet, session) {
    output$SAMinReadLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMinReadLength <- SangerConsensusSet@contigList[[contigIndex]]@minReadLength
            valueBox(
                subtitle = tags$p("MinReadLength",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(SAMinReadLength),
                               style = "font-size: 29px;"),
                icon = icon("sliders-h", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}
valueBoxSAMinFractionCall <- function(input, output, SangerConsensusSet, session) {
    output$SAMinFractionCall <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMinFractionCall <- SangerConsensusSet@contigList[[contigIndex]]@minFractionCall
            valueBox(
                subtitle = tags$p("MinFractionCall",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(as.numeric(SAMinFractionCall),
                               style = "font-size: 29px;"),
                icon = icon("sliders-h", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}
valueBoxSAMaxFractionLost <- function(input, output, SangerConsensusSet, session) {
    output$SAMaxFractionLost <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMaxFractionLost <- SangerConsensusSet@contigList[[contigIndex]]@maxFractionLost
            valueBox(
                subtitle = tags$p("MaxFractionLost",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(as.numeric(SAMaxFractionLost),
                               style = "font-size: 29px;"),
                icon = icon("sliders-h", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}
valueBoxSAAcceptStopCodons <- function(input, output, SangerConsensusSet, session) {
    output$SAAcceptStopCodons <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAAcceptStopCodons <- SangerConsensusSet@contigList[[contigIndex]]@acceptStopCodons
            valueBox(
                subtitle = tags$p("AcceptStopCodons",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(SAAcceptStopCodons,
                               style = "font-size: 29px;"),
                icon = icon("sliders-h", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}
valueBoxSAReadingFrame <- function(input, output, SangerConsensusSet, session) {
    output$SAReadingFrame <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAReadingFrame <-
                SangerConsensusSet@
                contigList[[contigIndex]]@readingFrame
            valueBox(
                subtitle = tags$p("ReadingFrame",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(SAReadingFrame),
                               style = "font-size: 29px;"),
                icon = icon("sliders-h", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}

### ============================================================================
### valueBox for SangerContig
### ============================================================================
valueBoxSCMinReadsNum <- function(input, output, SCMinReadsNum, session) {
    output$SCMinReadsNum <- renderUI({
        valueBox(
            subtitle = tags$p("MinReadsNum",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCMinReadsNum),
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}
valueBoxSCMinReadLength <- function(input, output, SCMinReadLength, session) {
    output$SCMinReadLength <- renderUI({
        valueBox(
            subtitle = tags$p("MinReadLength",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCMinReadLength),
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}
valueBoxSCMinFractionCall <- function(input, output,
                                      SCMinFractionCall, session) {
    output$SCMinFractionCall <- renderUI({
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SCMinFractionCall),
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}
valueBoxSCMaxFractionLost <- function(input, output,
                                      SCMaxFractionLost, session) {
    output$SCMaxFractionLost <- renderUI({
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                      font-weight: bold;"),
            value = tags$p(as.numeric(SCMaxFractionLost),
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}
valueBoxSCAcceptStopCodons <- function(input, output,
                                       SCAcceptStopCodons, session) {
    output$SCAcceptStopCodons <- renderUI({
        valueBox(
            subtitle = tags$p("AcceptStopCodons",
                              style = "font-size: 15px;
                                       font-weight: bold;"),
            value = tags$p(SCAcceptStopCodons,
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}
valueBoxSCReadingFrame <- function(input, output, SCReadingFrame, session) {
    output$SCReadingFrame <- renderUI({
        valueBox(
            subtitle = tags$p("ReadingFrame",
                              style = "font-size: 15px;
                                      font-weight: bold;"),
            value = tags$p(strtoi(SCReadingFrame),
                           style = "font-size: 29px;"),
            icon = icon("sliders-h", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox for SangerRead
### ============================================================================
valueBoxM1TrimmingCutoff <- function(input, output, session) {
    output$M1TrimmingCutoff <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(readIndex))) {
            if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
                as.numeric(input$M1TrimmingCutoffText) > 0 &&
                as.numeric(input$M1TrimmingCutoffText) <= 1) {
                inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
            } else {
                inputM1TrimmingCutoffText <- 0.0001
            }
            valueBox(
                subtitle = tags$p("Cutoff Score",
                                  style = "font-size: 15px;
                                       font-weight: bold;"),
                value = tags$p(as.numeric(inputM1TrimmingCutoffText),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 10,
            )
        }
    })
}
valueBoxM2CutoffQualityScore <- function(input, output, session) {
    output$M2CutoffQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(readIndex))) {
            if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
                strtoi(input$M2CutoffQualityScoreText) > 0 &&
                strtoi(input$M2CutoffQualityScoreText) <= 60 &&
                strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
                inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
            } else {
                inputM2CutoffQualityScoreText <- 20
            }
            valueBox(
                subtitle = tags$p("Cutoff Quality Score",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(inputM2CutoffQualityScoreText),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 10,
            )
        }
    })
}
valueBoxM2SlidingWindowSize <- function(input, output, session) {
    output$M2SlidingWindowSize <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(readIndex))) {
            if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
                strtoi(input$M2SlidingWindowSizeText) > 0 &&
                strtoi(input$M2SlidingWindowSizeText) <= 40 &&
                strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
                inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
            } else {
                inputM2SlidingWindowSizeText <- 10
            }
            valueBox(
                # strtoi(input$M2CutoffQualityScoreText
                subtitle = tags$p("Sliding Window Size ",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(inputM2SlidingWindowSizeText),
                               style = "font-size: 29px;"),
                icon = icon("expand", "fa-sm"),
                color = "olive", width = 10,
            )
        }
    })
}
valueBoxRawSeqLength <- function(input, output, session, trimmedRV) {
    output$rawSeqLength <- renderUI({
        valueBox(
            subtitle = tags$p("Raw Seqence Len ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["rawSeqLength"]]),
                           style = "font-size: 29px;"),
            icon = icon("ruler", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxRawMeanQualityScore <- function(input, output, session, trimmedRV) {
    output$rawMeanQualityScore <- renderUI({
        valueBox(
            subtitle = tags$p("Raw Mean Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(round(trimmedRV[["rawMeanQualityScore"]], 2)),
                           style = "font-size: 29px;"),
            icon = icon("cogs", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxRawMinQualityScore <- function(input, output, session, trimmedRV) {
    output$rawMinQualityScore <- renderUI({
        valueBox(
            subtitle = tags$p("Raw Min Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["rawMinQualityScore"]]),
                           style = "font-size: 29px;"),
            icon = icon("cogs", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxTrimmedStartPos <- function(input, output, session, trimmedRV) {
    output$trimmedStartPos <- renderUI({
        valueBox(
            subtitle = tags$p("Trimming Start Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedStartPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxTrimmedFinishPos <- function(input, output, session, trimmedRV) {
    output$trimmedFinishPos <- renderUI({
        valueBox(
            subtitle = tags$p("Trimming End Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedFinishPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("times-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxTrimmedSeqLength <- function(input, output, session, trimmedRV) {
    output$trimmedSeqLength <- renderUI({
        valueBox(
            subtitle = tags$p("Trimmed Seqence Length ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedSeqLength"]]),
                           style = "font-size: 29px;"),
            icon = icon("ruler", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxTrimmedMeanQualityScore <- function(input, output, session, trimmedRV) {
    output$trimmedMeanQualityScore <- renderUI({
        valueBox(
            subtitle = tags$p("Trimmed Mean Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(round(trimmedRV[["trimmedMeanQualityScore"]], 2)),
                           style = "font-size: 29px;"),
            icon = icon("cogs", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxTrimmedMinQualityScore <- function(input, output, session, trimmedRV) {
    output$trimmedMinQualityScore <- renderUI({
        valueBox(
            subtitle = tags$p("Trimmed Min Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedMinQualityScore"]]),
                           style = "font-size: 29px;"),
            icon = icon("cogs", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxRemainingRatio <- function(input, output, session, trimmedRV) {
    output$remainingRatio <- renderUI({
        valueBox(
            subtitle = tags$p("Remaining Ratio",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(paste(trimmedRV[["remainingRatio"]], "%"),
                           style = "font-size: 32px;"),
            icon = icon("divide", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}
valueBoxChromTrimmedStartPos <- function(input, output, session, trimmedRV) {
    output$ChromatogramtrimmedStartPos <- renderUI({
        valueBox(
            subtitle = tags$p("Trimming Start Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedStartPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}
valueBoxChromTrimmedFinishPos <- function(input, output, session, trimmedRV) {
    output$ChromatogramtrimmedFinishPos <- renderUI({
        valueBox(
            subtitle = tags$p("Trimming End Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedFinishPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("times-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### Plotting : qualityQualityBasePlot
### ============================================================================
qualityQualityBasePlotDisplay <- function(input, output, session,
                                          trimmedRV, qualityPhredScores) {
    trimmedStartPos <- trimmedRV[["trimmedStartPos"]]
    trimmedFinishPos <- trimmedRV[["trimmedFinishPos"]]
    readLen = length(qualityPhredScores)
    qualityPlotDf<- data.frame(seq_len(length(qualityPhredScores)),
                               qualityPhredScores)
    colnames(qualityPlotDf) <- c("Index", "Score")
    x <- list(
        title = "Base Pair Index"
        # titlefont = f
    )
    y <- list(
        title = "Phred Quality Score"
        # titlefont = f
    )
    p <- QualityBasePlotly(trimmedStartPos, trimmedFinishPos,
                           readLen, qualityPlotDf, x,  y)
    p
}

### ============================================================================
### SangerRead Sequence Render Function (DNA / AA) dynamic
### ============================================================================
primarySeqDisplay <- function(sequenceParam) {
    primarySeq <- unlist(strsplit(
        sequenceParam[["primarySeq"]], ""))
    primarySeqDF <- data.frame(
        t(data.frame(primarySeq)), stringsAsFactors = FALSE)
    colnames(primarySeqDF) <- substr(colnames(primarySeqDF), 2, 100)
    rownames(primarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(primarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(primarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(primarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(primarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = primarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
secondarySeqDisplay <- function(sequenceParam) {
    secondarySeq <-unlist(strsplit(
        sequenceParam[["secondarySeq"]], ""))
    secondarySeqDF <- data.frame(
        t(data.frame(secondarySeq)), stringsAsFactors = FALSE)
    colnames(secondarySeqDF) <- substr(colnames(secondarySeqDF), 2, 100)
    rownames(secondarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(secondarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(secondarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(secondarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(secondarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = secondarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
qualityScoreDisplay <- function(PhredScore) {
    PhredScoreDF <- data.frame(
        t(data.frame(PhredScore)), stringsAsFactors = FALSE)
    colnames(PhredScoreDF) <- substr(colnames(PhredScoreDF), 2, 100)
    rownames(PhredScoreDF) <- NULL
    styleList <- SetAllStyleList(PhredScoreDF, "#ecffd9")
    suppressMessages(
        excelTable(data =
                       PhredScoreDF, defaultColWidth = 30, editable = FALSE,
                   rowResize = FALSE, columnResize = FALSE,
                   allowInsertRow = FALSE, allowInsertColumn = FALSE,
                   allowDeleteRow = FALSE, allowDeleteColumn = FALSE,
                   style = styleList, allowRenameColumn = FALSE,
                   loadingSpin = TRUE)
    )
}
PrimAASeqS1Display <- function(sequenceParam) {
    AAString <- data.frame(AAString(sequenceParam[["primaryAASeqS1"]]))
    if (nchar(sequenceParam[["primaryAASeqS1"]]) == 0) {
        AAString <- rbind(NA, AAString)
        AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        width <- c(30)
        styleList <- list()
        styleList[['A1']] <- 'background-color: black;'
    } else {
        AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        width <- rep(90, length(AAStringDF))
        styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
        styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
        styleList <- c(styleList1, styleList2)
    }
    suppressMessages(
        excelTable(data = AAStringDF, columns = data.frame(width = width),
                   defaultColWidth = 90, editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
PrimAASeqS2Display <- function(sequenceParam) {
    AAString <- data.frame(AAString(sequenceParam[["primaryAASeqS2"]]))
    AAString <- rbind(NA, AAString)
    AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
    colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
    rownames(AAStringDF) <- NULL
    width <- rep(90, length(AAStringDF) - 1)
    width <- c(30, width)
    styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
    styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
    styleList <- c(styleList1, styleList2)
    styleList[['A1']] <- 'background-color: black;'
    suppressMessages(
        excelTable(data = AAStringDF, columns = data.frame(width = width),
                   defaultColWidth = 90, editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
PrimAASeqS3Display <- function(sequenceParam) {
    AAString <- data.frame(AAString(sequenceParam[["primaryAASeqS3"]]))
    AAString <- rbind(NA, NA, AAString)
    AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
    colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
    rownames(AAStringDF) <- NULL
    width <- rep(90, length(AAStringDF) - 2)
    width <- c(30, 30, width)
    styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
    styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
    styleList <- c(styleList1, styleList2)
    styleList[['A1']] <- 'background-color: black;'
    styleList[['B1']] <- 'background-color: black;'
    suppressMessages(
        excelTable(data = AAStringDF, columns = data.frame(width = width),
                   defaultColWidth = 90, editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
primarySeqTrimmedDisplay <- function(input, output, session,
                                     sequenceParam, trimmedRV) {
    primarySeq <- unlist(strsplit(
        substr(sequenceParam[["primarySeq"]],
               trimmedRV[["trimmedStartPos"]] + 1,
               trimmedRV[["trimmedFinishPos"]])
        , ""))
    primarySeqDF <- data.frame(
        t(data.frame(primarySeq)), stringsAsFactors = FALSE)
    if ((trimmedRV[["trimmedFinishPos"]]-trimmedRV[["trimmedStartPos"]]) == 1) {
        colnames(primarySeqDF) <- "1"
    } else {
        colnames(primarySeqDF) <- substr(colnames(primarySeqDF), 2, 100)
    }
    rownames(primarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(primarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(primarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(primarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(primarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = primarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}

secondSeqTrimmedDisplay <- function(input, output, session,
                                    sequenceParam, trimmedRV) {
    secondarySeq <- unlist(strsplit(
        substr(sequenceParam[["secondarySeq"]],
               trimmedRV[["trimmedStartPos"]] + 1,
               trimmedRV[["trimmedFinishPos"]])
        , ""))
    secondarySeqDF <- data.frame(
        t(data.frame(secondarySeq)), stringsAsFactors = FALSE)
    if ((trimmedRV[["trimmedFinishPos"]]-trimmedRV[["trimmedStartPos"]]) == 1) {
        colnames(secondarySeqDF) <- "1"
    } else {
        colnames(secondarySeqDF) <- substr(colnames(secondarySeqDF), 2, 100)
    }
    rownames(secondarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(secondarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(secondarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(secondarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(secondarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = secondarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}

