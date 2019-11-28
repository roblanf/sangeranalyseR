M1inside_calculate_trimming <- function(qualityPhredScores,
                                        qualityBaseScore,
                                        M1TrimmingCutoff) {
    rawSeqLength <- length(qualityBaseScore)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    start = FALSE
    trimmedStartPos = 0
    qualityBaseScoreCutOff = M1TrimmingCutoff - qualityBaseScore
    ### ------------------------------------------------------------------------
    ### calculate cummulative score
    ### if cumulative value < 0, set it to 0
    ### the BioPython implementation always trims the first base,
    ### this implementation does not.
    ### ------------------------------------------------------------------------
    score = qualityBaseScoreCutOff[1]
    if(score < 0){
        score = 0
    }else{
        trimmedStartPos = 1
        start = TRUE
    }
    cummul_score = c(score)
    ### ------------------------------------------------------------------------
    ### trimmedStartPos = value when cummulative score is first > 0
    ### ------------------------------------------------------------------------
    ### ------------------------------------------------------------------------
    ### trimmedFinishPos = index of highest cummulative score,
    ### marking the end of sequence segment with highest cummulative score
    ### ------------------------------------------------------------------------
    for(i in 2:length(qualityBaseScoreCutOff)){
        score = cummul_score[length(cummul_score)] + qualityBaseScoreCutOff[i]
        if (score <= 0) {
            cummul_score = c(cummul_score, 0)
        }else{
            cummul_score = c(cummul_score, score)
            if(start == FALSE){
                trimmedStartPos = i
                start = TRUE
            }
        }
        trimmedFinishPos = which.max(cummul_score)
    }
    ### ------------------------------------------------------------------------
    ### fix an edge case, where all scores are worse than the cutoff
    ### in this case you wouldn't want to keep any bases at all
    ### ------------------------------------------------------------------------
    if(sum(cummul_score)==0){trimmedFinishPos = 0}
    trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
    trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
    trimmedMinQualityScore <- min(trimmedQualityPhredScore)
    trimmedSeqLength = trimmedFinishPos - trimmedStartPos + 1
    remainingRatio = trimmedSeqLength / rawSeqLength

    return(c(rawSeqLength, rawMeanQualityScore, rawMinQualityScore,
             trimmedStartPos, trimmedFinishPos, trimmedSeqLength,
             trimmedMeanQualityScore, trimmedMinQualityScore, remainingRatio))
}

M2inside_calculate_trimming <- function(qualityPhredScores,
                                      qualityBaseScore,
                                      M2CutoffQualityScore,
                                      M2SlidingWindowSize) {
    rawSeqLength <- length(qualityBaseScore)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    qualityPbCutoff <- 10** (M2CutoffQualityScore / (-10.0))
    remainingIndex <- c()
    if (M2SlidingWindowSize > 20 || M2SlidingWindowSize < 0 ||
        M2SlidingWindowSize%%1!=0 ||
        M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
        M2CutoffQualityScore%%1!=0) {
        trimmedStartPos = NULL
        trimmedFinishPos = NULL
    } else {
        for (i in 1:(rawSeqLength-M2SlidingWindowSize+1)) {
            meanSLidingWindow <-
                mean(qualityBaseScore[i:(i+M2SlidingWindowSize-1)])
            if (meanSLidingWindow < qualityPbCutoff) {
                remainingIndex <- c(remainingIndex, i)
                # or ==> i + floor(M2SlidingWindowSize/3)
            }
        }
        trimmedStartPos = remainingIndex[1]
        trimmedFinishPos = remainingIndex[length(remainingIndex)]
        trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
        trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
        trimmedMinQualityScore <- min(trimmedQualityPhredScore)
        trimmedSeqLength = trimmedFinishPos - trimmedStartPos + 1
        remainingRatio = trimmedSeqLength / rawSeqLength
    }
    return(c(rawSeqLength, rawMeanQualityScore, rawMinQualityScore,
             trimmedStartPos, trimmedFinishPos, trimmedSeqLength,
             trimmedMeanQualityScore, trimmedMinQualityScore, remainingRatio))
}

### ============================================================================
### Adding dynamic menu to sidebar.
### ============================================================================
dynamicMenuSideBarSC <- function(input, output, session,
                               SangerSingleReadNum, SangerSingleReadFeature) {
    output$singleReadMenu <- renderMenu({
        menu_list <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(text = SangerSingleReadFeature[i],
                          tabName = SangerSingleReadFeature[i],
                          selected = TRUE, icon = icon("minus")))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Sanger Consensus Read Overview")})
}

dynamicMenuSideBarSCSet <- function(input, output, session, SangerCSetParam) {
    output$singleReadMenu <- renderMenu({
        SangerCSNum <- length(SangerCSetParam)
        menu_list <- sapply(1:SangerCSNum, function(i) {
            SangerSingleReadNum <- SangerCSetParam[[i]]$SangerSingleReadNum
            SangerCSMenuSubItem <- sapply(1:SangerSingleReadNum, function(j) {
                list(menuSubItem(text = SangerCSetParam[[i]]$SangerSingleReadFeature[[j]],
                                 tabName = paste0(i, " Consensus Read - ", SangerCSetParam[[i]]$SangerSingleReadFeature[[j]])))
            })
            SangerCSMenuSubItem <- c(list(menuSubItem(text = paste(SangerCSetParam[[i]]$SCName, "Overview"),
                                                      tabName = paste0(i, " Sanger Consensus Read Overview"), icon = icon("align-left"))),
                                     SangerCSMenuSubItem)
            SangerCSetParam[[i]]$SCName
            list(menuItem(text = SangerCSetParam[[i]]$SCName,
                          tabName = SangerCSetParam[[i]]$SCName,
                          icon = icon("minus"), SangerCSMenuSubItem))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Sanger Aligned Consensus Set Overview")})
}

### ============================================================================
### observeEvent: Adding dynamic rightHeader text
### ============================================================================
observeEventDynamicHeaderSC <- function(input, output, session, trimmedRV,
                                           SangerSingleReadQualReport) {
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        # message("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
        if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
            trimmedRV[["trimmedStartPos"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmedStartPos
            trimmedRV[["trimmedFinishPos"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
            qualityPhredScores = SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityPhredScores

            readLen = length(qualityPhredScores)
            trimmedRV[["trimmedSeqLength"]] <- trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1
            trimmedRV[["remainingRatio"]] <- round(((trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1) / readLen) * 100, digits = 2)
        }
    })
}






















# observeEventDynamicHeaderSCSet <- function(input, output, session, trimmedRV,
#                                            SangerCSetParam) {
#
#     output$res <- renderText({
#         paste("You've selected:", input$sidebar_menu)
#     })
#
#     observeEvent(input$sidebar_menu, {
#         menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
#         html("rightHeader", menuItem)
#         sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#         # message("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
#         if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
#         #     trimmedRV[["trimmedStartPos"]] <-
#         #         SangerSingleReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedStartPos
#         #     trimmedRV[["trimmedFinishPos"]] <-
#         #         SangerSingleReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
#         #     qualityPhredScores = SangerSingleReadQualReport[[
#         #         strtoi(sidebar_menu[[1]])]]@qualityPhredScores
#         #
#         #     readLen = length(qualityPhredScores)
#         #     trimmedRV[["trimmedSeqLength"]] <- trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1
#         #     trimmedRV[["remainingRatio"]] <- round(((trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1) / readLen) * 100, digits = 2)
#         }
#     })
# }

### ============================================================================
### valueBox: SCMinReadsNum
### ============================================================================


# SCRefAminoAcidSeq <- SangerConsensus@refAminoAcidSeq
# SCGeneticCode <- SangerConsensus@geneticCode
# SCAlignment<- SangerConsensus@alignment
# SCDifferencesDF<- SangerConsensus@differencesDF
# SCDistanceMatrix <- SangerConsensus@distanceMatrix
# SCDendrogram <- SangerConsensus@dendrogram
# SCIndelsDF <- SangerConsensus@indelsDF
# SCStopCodonsDF <- SangerConsensus@stopCodonsDF
# SCSecondaryPeakDF <- SangerConsensus@secondaryPeakDF
# SangerConsensusForRegExp <- SangerConsensus@consenesusReadName
# SangerConsensusForRegExp <- SangerConsensus@suffixForwardRegExp
# SangerConsensusRevRegExp <- SangerConsensus@suffixReverseRegExp








valueBoxSCMinReadsNum <- function(input, output, SCMinReadsNum, session) {
    output$SCMinReadsNum <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("MinReadsNum",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCMinReadsNum),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMinReadLength
### ============================================================================
valueBoxSCMinReadLength <- function(input, output, SCMinReadLength, session) {
    output$SCMinReadLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("MinReadLength",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCMinReadLength),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMinFractionCall
### ============================================================================
valueBoxSCMinFractionCall <- function(input, output, SCMinFractionCall, session) {
    output$SCMinFractionCall <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SCMinFractionCall),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMaxFractionLost
### ============================================================================
valueBoxSCMaxFractionLost <- function(input, output, SCMaxFractionLost, session) {
    output$SCMaxFractionLost <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SCMaxFractionLost),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCAcceptStopCodons
### ============================================================================
valueBoxSCAcceptStopCodons <- function(input, output, SCAcceptStopCodons, session) {
    output$SCAcceptStopCodons <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("AcceptStopCodons",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(SCAcceptStopCodons,
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCReadingFrame
### ============================================================================
valueBoxSCReadingFrame <- function(input, output, SCReadingFrame, session) {
    output$SCReadingFrame <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("ReadingFrame",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCReadingFrame),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

# trimmedQS <- reactiveValues(cuffOffQuality = 0, M2SlidingWindowSize = 0)


################################################################################
### Each Read
################################################################################
### ============================================================================
### valueBox: Change M2CutoffQualityScore
### ============================================================================
valueBoxM1TrimmingCutoff <- function(input, output, session, SangerSingleReadQualReport) {
    output$M1TrimmingCutoff <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
            as.numeric(input$M1TrimmingCutoffText) > 0 &&
            as.numeric(input$M1TrimmingCutoffText) <= 1) {
            inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
        } else {
            inputM1TrimmingCutoffText <- 0.0001
        }


        if (input$TrimmingMethodSelection == "M1") {
            # message("&&&& Dynamic M1")
            inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
            trimmingPos <-
                M1inside_calculate_trimming(
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        qualityPhredScores,
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@qualityBaseScore,
                    inputM1TrimmingCutoffText)
            rawSeqLength <- trimmingPos[1]
            rawMeanQualityScore <- trimmingPos[2]
            rawMinQualityScore <- trimmingPos[3]
            trimmedStartPos <- trimmingPos[4]
            trimmedFinishPos <- trimmingPos[5]
            trimmedSeqLength <- trimmingPos[6]
            trimmedMeanQualityScore <- trimmingPos[7]
            trimmedMinQualityScore <- trimmingPos[8]
            remainingRatio <- trimmingPos[9]

            if (!is.null(rawSeqLength) && !is.null(rawMeanQualityScore) &&
                !is.null(rawMinQualityScore ) && !is.null(trimmedStartPos) &&
                !is.null(trimmedFinishPos) && !is.null(trimmedSeqLength) &&
                !is.null(trimmedMeanQualityScore) &&
                !is.null(trimmedMinQualityScore)) {

                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)

                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    rawSeqLength <<- rawSeqLength
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    rawMeanQualityScore <<- rawMeanQualityScore
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    rawMinQualityScore <<- rawMinQualityScore
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    trimmedStartPos <<- trimmedStartPos
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    trimmedFinishPos <<- trimmedFinishPos
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    trimmedSeqLength <<- trimmedSeqLength
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    trimmedMeanQualityScore <<- trimmedMeanQualityScore
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    trimmedMinQualityScore <<- trimmedMinQualityScore
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    remainingRatio <<- remainingRatio

                trimmedRV[["rawSeqLength"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@remainingRatio * 100, 2)
            }
        }














        valueBox(
            subtitle = tags$p("Cut Off Log Score",
                              style = "font-size: 15px;
                                       font-weight: bold;"),
            value = tags$p(as.numeric(inputM1TrimmingCutoffText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 10,
        )
    })
}

### ============================================================================
### valueBox: Change M2CutoffQualityScore
### ============================================================================
valueBoxM2CutoffQualityScore <- function(input, output, session, SangerSingleReadQualReport) {
    output$M2CutoffQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        # trimmedQS[["cuffOffQuality"]],
        # trimmedQS[["M2SlidingWindowSize"]])
        if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
            strtoi(input$M2CutoffQualityScoreText) > 0 &&
            strtoi(input$M2CutoffQualityScoreText) <= 60 &&
            strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
            inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
        } else {
            inputM2CutoffQualityScoreText <- 20
        }



        if (input$TrimmingMethodSelection == "M2") {
            # message("&&&& Dynamic M2")
            if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
                strtoi(input$M2CutoffQualityScoreText) > 0 &&
                strtoi(input$M2CutoffQualityScoreText) <= 60 &&
                strtoi(input$M2CutoffQualityScoreText) %% 1 ==0 &&
                !is.na(strtoi(input$M2SlidingWindowSizeText)) &&
                strtoi(input$M2SlidingWindowSizeText) > 0 &&
                strtoi(input$M2SlidingWindowSizeText) <= 20 &&
                strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
                inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
                inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
                trimmingPos <-
                    M2inside_calculate_trimming(
                        SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                            qualityPhredScores,
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@qualityBaseScore,
                        inputM2CutoffQualityScoreText,
                        inputM2SlidingWindowSizeText)
                rawSeqLength <- trimmingPos[1]
                rawMeanQualityScore <- trimmingPos[2]
                rawMinQualityScore <- trimmingPos[3]
                trimmedStartPos <- trimmingPos[4]
                trimmedFinishPos <- trimmingPos[5]
                trimmedSeqLength <- trimmingPos[6]
                trimmedMeanQualityScore <- trimmingPos[7]
                trimmedMinQualityScore <- trimmingPos[8]
                remainingRatio <- trimmingPos[9]

                if (!is.null(rawSeqLength) && !is.null(rawMeanQualityScore) &&
                    !is.null(rawMinQualityScore ) && !is.null(trimmedStartPos) &&
                    !is.null(trimmedFinishPos) && !is.null(trimmedSeqLength) &&
                    !is.null(trimmedMeanQualityScore) &&
                    !is.null(trimmedMinQualityScore)) {

                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        M1TrimmingCutoff <<- as.numeric(inputM2CutoffQualityScoreText)
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)

                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        rawSeqLength <<- rawSeqLength
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        rawMeanQualityScore <<- rawMeanQualityScore
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        rawMinQualityScore <<- rawMinQualityScore
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        trimmedStartPos <<- trimmedStartPos
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        trimmedFinishPos <<- trimmedFinishPos
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        trimmedSeqLength <<- trimmedSeqLength
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        trimmedMeanQualityScore <<- trimmedMeanQualityScore
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        trimmedMinQualityScore <<- trimmedMinQualityScore
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        remainingRatio <<- remainingRatio

                    trimmedRV[["rawSeqLength"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@rawSeqLength
                    trimmedRV[["rawMeanQualityScore"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@rawMeanQualityScore
                    trimmedRV[["rawMinQualityScore"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@rawMinQualityScore
                    trimmedRV[["trimmedStartPos"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@trimmedStartPos
                    trimmedRV[["trimmedFinishPos"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
                    trimmedRV[["trimmedSeqLength"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@trimmedSeqLength
                    trimmedRV[["trimmedMeanQualityScore"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@trimmedMeanQualityScore
                    trimmedRV[["trimmedMinQualityScore"]] <<-
                        SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@trimmedMinQualityScore
                    trimmedRV[["remainingRatio"]] <<-
                        round(SangerSingleReadQualReport[[
                            strtoi(sidebar_menu[[1]])]]@remainingRatio * 100, 2)
                }
            }
        }





        valueBox(
            subtitle = tags$p("Cut Off Quality Score",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(inputM2CutoffQualityScoreText),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 10,
        )
    })
}

### ============================================================================
### valueBox: Change M2SlidingWindowSize
### ============================================================================
valueBoxM2SlidingWindowSize <- function(input, output, session) {
    output$M2SlidingWindowSize <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
            strtoi(input$M2SlidingWindowSizeText) > 0 &&
            strtoi(input$M2SlidingWindowSizeText) <= 20 &&
            strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
            inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
        } else {
            inputM2SlidingWindowSizeText <- 5
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
    })
}









































### ============================================================================
### valueBox: Change rawSeqLength
### ============================================================================
valueBoxRawSeqLength <- function(input, output, session, trimmedRV) {
    output$rawSeqLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Raw Seqence Len ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["rawSeqLength"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}


### ============================================================================
### valueBox: Change rawMeanQualityScore
### ============================================================================
valueBoxRawMeanQualityScore <- function(input, output, session, trimmedRV) {
    output$rawMeanQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Raw Mean Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(round(trimmedRV[["rawMeanQualityScore"]], 2)),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}


### ============================================================================
### valueBox: Change rawMinQualityScore
### ============================================================================
valueBoxRawMinQualityScore <- function(input, output, session, trimmedRV) {
    output$rawMinQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Raw Min Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["rawMinQualityScore"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}


### ============================================================================
### valueBox: Change trimmedStartPos
### ============================================================================
valueBoxTrimmedStartPos <- function(input, output, session, trimmedRV) {
    output$trimmedStartPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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

### ============================================================================
### valueBox: Change trimmedFinishPos
### ============================================================================
valueBoxTrimmedFinishPos <- function(input, output, session, trimmedRV) {
    output$trimmedFinishPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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
### valueBox: Change trimmedSeqLength
### ============================================================================
valueBoxTrimmedSeqLength <- function(input, output, session, trimmedRV) {
    output$trimmedSeqLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimmed Seqence Length ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedSeqLength"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedMeanQualityScore
### ============================================================================
valueBoxTrimmedMeanQualityScore <- function(input, output, session, trimmedRV) {
    output$trimmedMeanQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimmed Mean Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(round(trimmedRV[["trimmedMeanQualityScore"]], 2)),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedMinQualityScore
### ============================================================================
valueBoxTrimmedMinQualityScore <- function(input, output, session, trimmedRV) {
    output$trimmedMinQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimmed Min Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedMinQualityScore"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change remainingRatio
### ============================================================================
valueBoxRemainingRatio <- function(input, output, session, trimmedRV) {
    output$remainingRatio <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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



### ============================================================================
### valueBox: Change trimmedStartPos
### ============================================================================
valueBoxChromTrimmedStartPos <- function(input, output, session, trimmedRV) {
    output$ChromatogramtrimmedStartPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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

### ============================================================================
### valueBox: Change trimmedFinishPos
### ============================================================================
valueBoxChromTrimmedFinishPos <- function(input, output, session, trimmedRV) {
    output$ChromatogramtrimmedFinishPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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
### qualityTrimmingRatioPlot
### ============================================================================
qualityTrimmingRatioPlot <- function(input, output, session, trimmedRV,
                                     SangerSingleReadQualReport,
                                     SangerSingleReadFeature) {
    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readFeature <- SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]]
        trimmedStartPos = trimmedRV[["trimmedStartPos"]]
        trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
        qualityPhredScores = SangerSingleReadQualReport[[
            strtoi(sidebar_menu[[1]])]]@qualityPhredScores
        readLen = length(qualityPhredScores)

        stepRatio = 1 / readLen
        trimmedStartPos / readLen
        trimmedFinishPos / readLen

        trimmedPer <- c()
        remainingPer <- c()

        for (i in 1:trimmedStartPos) {
            if (i != trimmedStartPos) {
                trimmedPer <- c(trimmedPer, stepRatio)
            }
        }

        for (i in trimmedStartPos:trimmedFinishPos) {
            trimmedPer <- c(trimmedPer, 0)
        }


        for (i in trimmedFinishPos:readLen) {
            if (i != trimmedFinishPos) {
                trimmedPer <- c(trimmedPer, stepRatio)
            }
        }

        trimmedPer <- cumsum(trimmedPer)
        remainingPer = 1 - trimmedPer

        PerData <- data.frame(1:length(trimmedPer),
                              trimmedPer, remainingPer)

        colnames(PerData) <- c("Base",
                               "Trimmed Ratio",
                               "Remaining Ratio")
        # Change font setting
        # f <- list(
        #     family = "Courier New, monospace",
        #     size = 18,
        #     color = "#7f7f7f"
        # )
        x <- list(
            title = "Base Pair Index"
            # titlefont = f
        )
        y <- list(
            title = "Read Ratio"
            # titlefont = f
        )
        PerDataPlot <- melt(PerData, id.vars = c("Base"))
        plot_ly(data=PerDataPlot,
                x=~Base,
                y=~value,
                mode="markers",
                color = ~variable,
                text = ~paste("BP Index : ",
                              Base, '<sup>th</sup><br>Read Ratio :',
                              round(value*100, digits = 2), '%')) %>%
            layout(xaxis = x, yaxis = y, legend = list(orientation = 'h',
                                                       xanchor = "center",  # use center of legend as anchor
                                                       x = 0.5, y = 1.1)) %>%
            add_annotations(
                text = "Trimmed Ratio (Each BP)",
                x = (trimmedStartPos + trimmedFinishPos) / 2,
                y = ((trimmedPer[1] + trimmedPer[length(trimmedPer)]) / 2)
                + 0.06,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Remaining Ratio (Each BP)",
                x = (trimmedStartPos+trimmedFinishPos) / 2,
                y = ((remainingPer[1] + remainingPer[length(remainingPer)]) / 2)
                - 0.06,
                showarrow=FALSE
            )
    })
}

### ============================================================================
### qualityQualityBasePlot
### ============================================================================
qualityQualityBasePlot <- function(input, output, session, trimmedRV,
                                   SangerSingleReadQualReport,
                                   SangerSingleReadFeature) {
    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readFeature <- SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]]
        trimmedStartPos = trimmedRV[["trimmedStartPos"]]
        trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
        qualityPhredScores = SangerSingleReadQualReport[[
            strtoi(sidebar_menu[[1]])]]@qualityPhredScores
        readLen = length(qualityPhredScores)

        qualityPlotDf<- data.frame(1:length(qualityPhredScores),
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

        plot_ly(data=qualityPlotDf,
                x=~Index) %>%
            add_markers(y=~Score,
                        text = ~paste("BP Index : ",
                                      Index,
                                      '<sup>th</sup><br>Phred Quality Score :',
                                      Score),
                        name = 'Quality Each BP') %>%
            add_trace(x=seq(trimmedStartPos,
                            trimmedFinishPos,
                            len=trimmedFinishPos-trimmedStartPos+1),
                      y=rep(70, trimmedFinishPos-trimmedStartPos+1),
                      mode="lines", hoverinfo="text",
                      text=paste("Trimmed Reads BP length:",
                                 trimmedFinishPos-trimmedStartPos+1,
                                 "BPs <br>",
                                 "Trimmed Reads BP ratio:",
                                 round((trimmedFinishPos - trimmedStartPos+1)/
                                           readLen * 100,
                                       digits=2),
                                 "%"),
                      line = list(width = 12),
                      name = 'Trimmed Read') %>%
            add_trace(x=seq(0,readLen,len=readLen),
                      y=rep(80, readLen), mode="lines", hoverinfo="text",
                      text=paste("Whole Reads BP length:",
                                 readLen,
                                 "BPs <br>",
                                 "Trimmed Reads BP ratio: 100 %"),
                      line = list(width = 12),
                      name = 'Whole Read') %>%
            layout(xaxis = x, yaxis = y,
                   shapes = list(vline(trimmedStartPos),
                                 vline(trimmedFinishPos)),
                   legend = list(orientation = 'h',
                                 xanchor = "center",  # use center of legend as anchor
                                 x = 0.5, y = 1.1)) %>%
            # add_segments(x = trimmedStartPos, xend = trimmedFinishPos, y = 70, yend = 70, inherit = TRUE, width = 10, line = list(width = 8)) %>%
            # add_segments(x = 0, xend = readLen, y = 75, yend = 75, inherit = TRUE, width = 4, line = list(width = 8)) %>%
            add_annotations(
                text = "Trimming Strat <br> BP Index",
                x = trimmedStartPos + 40,
                y = 15,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Trimming End <br> BP Index",
                x = trimmedFinishPos - 40,
                y = 15,
                showarrow=FALSE
            )
        # add_markers(qualityPlotDf, x=~Index, y=~Score)
        # add_segments(x = trimmedStartPos, xend = trimmedFinishPos, y = 70, yend = 70, inherit = TRUE)
    })
}





### ============================================================================
### chromatogram row number counting
### ============================================================================
chromatogramRowNum <- function(obj, width) {
    traces <- obj@traceMatrix
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    basecalls1 <- basecalls1[1:length(aveposition)]
    valuesperbase <- nrow(traces)/length(basecalls1)
    tracewidth <- width*valuesperbase
    breaks <- seq(1,nrow(traces), by=tracewidth)
    numplots <- length(breaks)
    return(numplots)
}


vline <- function(x = 0, color = "red") {
    list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = x,
        x1 = x,
        line = list(color = color)
    )
}

