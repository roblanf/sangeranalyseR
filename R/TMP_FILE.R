# ### ------------------------------------------------------------------------
# ### observeEvent: Button Save S4 object
# ### ------------------------------------------------------------------------
# observeEvent(input$saveS4, {
#     newS4Object <- file.path(shinyDirectory, "SangerConsensus.Rda")
#     showNotification(paste("New S4 object is store as:", newS4Object),
#                      type = "message", duration = 10)
#     ### --------------------------------------------------------------------
#     ### Save SangerConsensus quality S4 object
#     ### --------------------------------------------------------------------
#     forwardReadNum <- length((SangerConsensus)@forwardReadsList)
#     reverseReadNum <- length((SangerConsensus)@reverseReadsList)
#     SangerSingleReadNum <- forwardReadNum + reverseReadNum
#     sapply(1:forwardReadNum, function(i) {
#         SangerConsensus@forwardReadsList[[i]]@QualityReport <<-
#             SangerSingleReadQualReport[[i]]
#         SangerConsensus@forwardReadsList[[i]]@ChromatogramParam <<-
#             SangerSingleReadChromatogramParam[[i]]
#         message("save SangerConsensus quality S4 object Forward")
#     })
#     sapply(1:reverseReadNum, function(i) {
#         SangerConsensus@reverseReadsList[[i]]@QualityReport <<-
#             SangerSingleReadQualReport[[forwardReadNum + i]]
#         SangerConsensus@reverseReadsList[[i]]@ChromatogramParam <<-
#             SangerSingleReadChromatogramParam[[forwardReadNum + i]]
#         message("save SangerConsensus quality S4 object Reverse")
#     })
#     saveRDS(SangerConsensus, file=newS4Object)
#     message("New S4 object is store as: ", newS4Object)
#     NEW_SANGER_CONSENSUS_READ <<- readRDS(file=newS4Object)
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Close UI
# ### ------------------------------------------------------------------------
# observeEvent(input$closeUI, {
#     btn <- input$closeUI
#     stopApp()
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Consensus read re-calculating UI
# ### ------------------------------------------------------------------------
# observeEvent(input$recalculateButton, {
#     message("######## Reactive button clicked !!!")
#     message("######## Start recalculating consensus read (SC")
#     CSResult<-
#         calculateConsensusRead (SangerConsensus@forwardReadsList,
#                                 SangerConsensus@reverseReadsList,
#                                 SangerConsensus@refAminoAcidSeq,
#                                 SangerConsensus@minFractionCall,
#                                 SangerConsensus@maxFractionLost,
#                                 SangerConsensus@geneticCode,
#                                 SangerConsensus@acceptStopCodons,
#                                 SangerConsensus@readingFrame)
#
#     SangerConsensus@consensusRead <<- CSResult$consensusGapfree
#     SangerConsensus@differencesDF <<- CSResult$diffsDf
#     SangerConsensus@alignment <<- CSResult$aln2
#     SangerConsensus@distanceMatrix <<- CSResult$dist
#     SangerConsensus@dendrogram <<- CSResult$dend
#     SangerConsensus@indelsDF <<- CSResult$indels
#     SangerConsensus@stopCodonsDF <<- CSResult$stopsDf
#     SangerConsensus@secondaryPeakDF <<- CSResult$spDf
#
#     consensusParam[["consensusRead"]] <<- SangerConsensus@consensusRead
#     consensusParam[["differencesDF"]] <<- SangerConsensus@differencesDF
#     consensusParam[["alignment"]] <<- SangerConsensus@alignment
#     consensusParam[["distanceMatrix"]] <<-SangerConsensus@distanceMatrix
#     consensusParam[["dendrogram"]] <<- SangerConsensus@dendrogram
#     consensusParam[["indelsDF"]] <<- SangerConsensus@indelsDF
#     consensusParam[["stopCodonsDF"]] <<- SangerConsensus@stopCodonsDF
#     consensusParam[["secondaryPeakDF"]] <<-
#         SangerConsensus@secondaryPeakDF
#     message("######## Finish recalculation")
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
# ### ------------------------------------------------------------------------
# observeEvent(input$startTrimmingButton, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.null(SangerSingleReadQualReport[[singleReadIndex]])) {
#         if (SangerSingleReadQualReport[[singleReadIndex]]@
#             TrimmingMethod == "M1") {
#             if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
#                 as.numeric(input$M1TrimmingCutoffText) > 0 &&
#                 as.numeric(input$M1TrimmingCutoffText) <= 1) {
#                 inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
#             } else {
#                 inputM1TrimmingCutoffText <- 0.0001
#             }
#             ### ------------------------------------------------------------
#             ### Start M1 trimming calculation
#             ### ------------------------------------------------------------
#             trimmingPos <-
#                 M1inside_calculate_trimming(
#                     SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityPhredScores,
#                     SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityBaseScores,
#                     as.numeric(inputM1TrimmingCutoffText))
#
#             SangerSingleReadQualReport[[singleReadIndex]]@
#                 M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)
#             trimmedParam[["M1TrimmingCutoff"]] <<-
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M1TrimmingCutoff
#
#         } else if (SangerSingleReadQualReport[[singleReadIndex]]@
#                    TrimmingMethod == "M2") {
#             if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
#                 strtoi(input$M2CutoffQualityScoreText) > 0 &&
#                 strtoi(input$M2CutoffQualityScoreText) <= 60 &&
#                 strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
#                 inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
#             } else {
#                 inputM2CutoffQualityScoreText <- 20
#             }
#             if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
#                 strtoi(input$M2SlidingWindowSizeText) > 0 &&
#                 strtoi(input$M2SlidingWindowSizeText) <= 20 &&
#                 strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
#                 inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
#             } else {
#                 inputM2SlidingWindowSizeText <- 5
#             }
#             ### ------------------------------------------------------------
#             ### Start M2 trimming calculation
#             ### ------------------------------------------------------------
#             trimmingPos <-
#                 M2inside_calculate_trimming(
#                     SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityPhredScores,
#                     SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityBaseScores,
#                     strtoi(inputM2CutoffQualityScoreText),
#                     SangerSingleReadQualReport[[singleReadIndex]]@
#                         M2SlidingWindowSize)
#
#             SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2CutoffQualityScore <<-
#                 strtoi(inputM2CutoffQualityScoreText)
#             trimmedParam[["M2CutoffQualityScore"]] <<-
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2CutoffQualityScore
#             SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)
#             trimmedParam[["M2SlidingWindowSize"]] <<-
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2SlidingWindowSize
#         }
#
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             rawSeqLength <<- trimmingPos[["rawSeqLength"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
#         SangerSingleReadQualReport[[singleReadIndex]]@
#             remainingRatio <<- trimmingPos[["remainingRatio"]]
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
#         ### ------------------------------------------------------------
#         ### Save SangerConsensus quality S4 object
#         ### ------------------------------------------------------------
#         forwardReadNum <- length((SangerConsensus)@forwardReadsList)
#         reverseReadNum <- length((SangerConsensus)@reverseReadsList)
#         SangerSingleReadNum <- forwardReadNum + reverseReadNum
#         if (singleReadIndex <= forwardReadNum) {
#             # This is forward list
#             SangerConsensus@
#                 forwardReadsList[[singleReadIndex]]@
#                 QualityReport <<-
#                 SangerSingleReadQualReport[[singleReadIndex]]
#         } else {
#             # This is reverse list
#             SangerConsensus@
#                 reverseReadsList[[singleReadIndex-forwardReadNum]]@
#                 QualityReport <<-
#                 SangerSingleReadQualReport[[singleReadIndex]]
#         }
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
# ### ------------------------------------------------------------------------
# observeEvent(input$saveChromatogramParam, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     ### ------------------------------------------------------------
#     ### Update ChromatogramBasePerRow
#     ### ------------------------------------------------------------
#     SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         baseNumPerRow <<- input$ChromatogramBasePerRow
#     SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         heightPerRow <<- input$ChromatogramHeightPerRow
#     SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
#     SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         showTrimmed <<- input$ChromatogramCheckShowTrimmed
#
#     ### ------------------------------------------------------------
#     ### Save SangerConsensus quality S4 object
#     ### ------------------------------------------------------------
#     forwardReadNum <- length((SangerConsensus)@forwardReadsList)
#     reverseReadNum <- length((SangerConsensus)@reverseReadsList)
#     SangerSingleReadNum <- forwardReadNum + reverseReadNum
#
#     if (singleReadIndex <= forwardReadNum) {
#         # This is forward list
#         SangerConsensus@
#             forwardReadsList[[singleReadIndex]]@ChromatogramParam <<-
#             SangerSingleReadChromatogramParam[[singleReadIndex]]
#     } else {
#         # This is reverse list
#         SangerConsensus@reverseReadsList[[singleReadIndex-forwardReadNum]]@
#             ChromatogramParam <<-
#             SangerSingleReadChromatogramParam[[singleReadIndex]]
#     }
#     ChromatogramParam[["baseNumPerRow"]] <<-
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         baseNumPerRow
#     ChromatogramParam[["heightPerRow"]] <<-
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         heightPerRow
#     ChromatogramParam[["signalRatioCutoff"]] <<-
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         signalRatioCutoff
#     ChromatogramParam[["showTrimmed"]] <<-
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         showTrimmed
# })
#
#
# ############################################################################
# ### ConsensusRead (Function for Sanger Consensus Read Overview)
# ############################################################################
# output$geneticCodeDF <- renderExcel({
#     suppressMessages(
#         excelTable(data =  t(data.frame(SangerConsensus@geneticCode)),
#                    defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
#     )
# })
#
# output$SCrefAminoAcidSeqDF <- renderExcel({
#     refAminoAcidSeqVec <- strsplit(SangerConsensus@refAminoAcidSeq, "")[[1]]
#     names(refAminoAcidSeqVec) <- c(1:length(refAminoAcidSeqVec))
#     suppressMessages(
#         excelTable(data =
#                        t(data.frame(refAminoAcidSeqVec)),
#                    defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
#     )
# })
#
# ### ------------------------------------------------------------------------
# ### refAminoAcidSeq
# ### ------------------------------------------------------------------------
# output$SCrefAminoAcidSeq <- renderUI({
#     if (SangerConsensus@refAminoAcidSeq == "") {
#         box(title = tags$p("Reference Amino Acids Sequence",
#                            style = "font-size: 24px;
#                                 font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 1),
#             column(width = 11,
#                    h4("Reference Amino Acid Sequence is not provided."))
#         )
#     } else {
#         box(title = tags$p("Reference Amino Acids Sequence",
#                            style = "font-size: 24px;
#                                 font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 2,
#                    tags$br(),
#                    tags$p("AA Sequence:",
#                           style = "font-size: 15px;
#                                    font-weight: bold;"),
#             ),
#             column(width = 10,
#                    excelOutput("SCrefAminoAcidSeqDF",
#                                width = "100%", height = "50"),
#                    style = paste("height:100%; ",
#                                  "overflow-y: hidden;",
#                                  "overflow-x: scroll;")
#             ),
#         )
#     }
# })
# ### ------------------------------------------------------------------------
# ### Alignment
# ### ------------------------------------------------------------------------
# output$consensusAlignmentHTML <- renderUI({
#     browseSeqHTML <-
#         file.path(shinyDirectory,
#                   paste0(SangerConsensus@consensusReadName,
#                          "_Alignment_BrowseSeqs.html"))
#     BrowseSeqs(consensusParam[["alignment"]],
#                openURL=FALSE, htmlFile=browseSeqHTML)
#     includeHTML(
#         file.path(shinyDirectory,
#                   paste0(SangerConsensus@consensusReadName,
#                          "_Alignment_BrowseSeqs.html")))
# })
#
# ### ------------------------------------------------------------------------
# ### Difference
# ### ------------------------------------------------------------------------
# output$SCDifferencesDFUI <- renderUI({
#     if (all(dim(consensusParam[["differencesDF"]]) == c(0,0))) {
#         h4("*** 'Differences' dataframe is empty. ***",
#            style="font-weight: bold; font-style: italic;")
#     } else {
#         dataTableOutput("SCDifferencesDF")
#
#     }
# })
# output$SCDifferencesDF = renderDataTable({
#     consensusParam[["differencesDF"]]
# })
#
# ### ------------------------------------------------------------------------
# ### dendrogram
# ### ------------------------------------------------------------------------
# output$dendrogramPlot <- renderPlot({
#     # plot(consensusParam[["dendrogram"]][[2]])
#     ggdendrogram(consensusParam[["dendrogram"]][[2]], rotate = TRUE)
# })
# output$dendrogramDF <- renderDataTable({
#     consensusParam[["dendrogram"]][[1]]
# })
#
# ### ------------------------------------------------------------------------
# ### Distance
# ### ------------------------------------------------------------------------
# output$SCDistanceMatrixPlotUI <- renderUI({
#     if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
#         h4("*** 'Distance' dataframe is empty. (Cannot plot)***",
#            style="font-weight: bold; font-style: italic;")
#     } else {
#         plotlyOutput("SCDistanceMatrixPlot")
#     }
# })
# output$SCDistanceMatrixPlot <- renderPlotly({
#     plot_ly(x = SangerSingleReadBFN,
#             y = SangerSingleReadBFN,
#             z = consensusParam[["distanceMatrix"]],
#             colors = colorRamp(c("white", "#32a852")),
#             type = "heatmap")
# })
# output$SCDistanceMatrixUI <- renderUI({
#     if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
#         h4("*** 'Distance' dataframe is empty. ***",
#            style="font-weight: bold; font-style: italic;")
#     } else {
#         dataTableOutput("SCDistanceMatrix")
#     }
# })
# output$SCDistanceMatrix = renderDataTable({
#     consensusParam[["distanceMatrix"]]
# })
#
#
# ### ------------------------------------------------------------------------
# ### Indels
# ### ------------------------------------------------------------------------
# output$SCIndelsDFUI <- renderUI({
#     if (all(dim(consensusParam[["indelsDF"]]) == c(0,0))) {
#         h4("*** 'Indels' dataframe is empty. ***",
#            style="font-weight: bold; font-style: italic;")
#     } else {
#         dataTableOutput("SCIndelsDF")
#
#     }
# })
# output$SCIndelsDF <- renderDataTable({
#     consensusParam[["indelsDF"]]
# })
#
# ### ------------------------------------------------------------------------
# ### StopCodons
# ### ------------------------------------------------------------------------
# output$SCStopCodonsDFUI <- renderUI({
#     if (all(dim(consensusParam[["stopCodonsDF"]]) == c(0,0))) {
#         h4("*** 'Stop Codons' dataframe is empty. ***",
#            style="font-weight: bold; font-style: italic;")
#     } else {
#         dataTableOutput("SCStopCodonsDF")
#     }
# })
# output$SCStopCodonsDF <- renderDataTable({
#     consensusParam[["stopCodonsDF"]]
# })
#
# ### ------------------------------------------------------------------------
# ### Valuebox for basic information
# ### ------------------------------------------------------------------------
# valueBoxSCMinReadsNum(input, output,
#                       SangerConsensus@minReadsNum, session)
# valueBoxSCMinReadLength(input, output,
#                         SangerConsensus@minReadLength, session)
# valueBoxSCMinFractionCall(input, output,
#                           SangerConsensus@minFractionCall, session)
# valueBoxSCMaxFractionLost(input, output,
#                           SangerConsensus@maxFractionLost, session)
# valueBoxSCAcceptStopCodons(input, output,
#                            SangerConsensus@acceptStopCodons, session)
# valueBoxSCReadingFrame(input, output,
#                        SangerConsensus@readingFrame, session)
#
# ############################################################################
# ### SangerSingleRead (Function for singel read in consensusRead)
# ############################################################################
# ### ------------------------------------------------------------------------
# ### primarySeq & secondSeq related
# ### ------------------------------------------------------------------------
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
# # sequenceParam[["primarySeq"]] <<-
# #     SangerSingleReadPrimSeqDF[[singleReadIndex]]
# # sequenceParam[["secondarySeq"]] <<-
# #     SangerSingleReadSecoSeqDF[[singleReadIndex]]
# # sequenceParam[["primaryAASeqS1"]] <<-
# #     SangerSingleReadPrimAASeqS1DF[[singleReadIndex]]
# # sequenceParam[["primaryAASeqS2"]] <<-
# #     SangerSingleReadPrimAASeqS2DF[[singleReadIndex]]
# # sequenceParam[["primaryAASeqS3"]] <<-
# #     SangerSingleReadPrimAASeqS3DF[[singleReadIndex]]
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
# output$primarySeqDF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     AstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "A", "#1eff00")
#     TstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "T", "#ff7a7a")
#     CstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "C", "#7ac3ff")
#     GstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "G", "#c9c9c9")
#     styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
#     suppressMessages(
#         excelTable(data =
#                        SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                    defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                    style = styleList, loadingSpin = TRUE)
#     )
# })
#
# output$secondSeqDF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     AstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "A", "#1eff00")
#     TstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "T", "#ff7a7a")
#     CstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "C", "#7ac3ff")
#     GstyleList <-
#         getStopList(SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                     "G", "#c9c9c9")
#     styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
#     suppressMessages(
#         excelTable(data =
#                        SangerSingleReadSecoSeqDF[[singleReadIndex]],
#                    defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                    style = styleList, loadingSpin = TRUE)
#     )
# })
#
# output$qualityScoreDF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     suppressMessages(
#         excelTable(data =
#                        SangerSingleReadQSDF[[singleReadIndex]],
#                    defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                    loadingSpin = TRUE)
#     )
# })
#
# output$PrimAASeqS1DF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     width <- rep(90, length(SangerSingleReadPrimAASeqS1DF[[singleReadIndex]]))
#     styleList <-
#         getStopList (SangerSingleReadPrimAASeqS1DF[[singleReadIndex]],
#                      "*", "#cf0000")
#     suppressMessages(
#         excelTable(data =
#                        SangerSingleReadPrimAASeqS1DF[[singleReadIndex]],
#                    columns = data.frame(width = width),
#                    defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                    style = styleList, loadingSpin = TRUE)
#     )
# })
#
# output$PrimAASeqS2DF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     width <- rep(90, length(SangerSingleReadPrimAASeqS2DF[[singleReadIndex]])-1)
#     widthFinal <- c(30, width)
#     styleList <-
#         getStopList (SangerSingleReadPrimAASeqS2DF[[singleReadIndex]],
#                      "*", "#cf0000")
#     styleList[['A1']] <- 'background-color: black;'
#     suppressMessages(
#         excelTable(data =
#                        SangerSingleReadPrimAASeqS2DF[[singleReadIndex]],
#                    columns = data.frame(width = widthFinal),
#                    defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                    style = styleList, loadingSpin = TRUE)
#     )
# })
#
# output$PrimAASeqS3DF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     width <- rep(90, length(SangerSingleReadPrimAASeqS3DF[[singleReadIndex]])-2)
#     widthFinal <- c(30, 30, width)
#     styleList <-
#         getStopList (SangerSingleReadPrimAASeqS3DF[[singleReadIndex]],
#                      "*", "#cf0000")
#     styleList[['A1']] <- 'background-color: black;'
#     styleList[['B1']] <- 'background-color: black;'
#     suppressMessages(
#         excelTable(data =
#                        SangerSingleReadPrimAASeqS3DF[[singleReadIndex]],
#                    columns = data.frame(width = widthFinal),
#                    defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
#                    columnResize = FALSE, allowInsertRow = FALSE,
#                    allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                    allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                    style = styleList, loadingSpin = TRUE)
#     )
# })
#
# ### ------------------------------------------------------------------------
# ### Quality trimming related (value box)
# ### ------------------------------------------------------------------------
# valueBoxM1TrimmingCutoff (input, output, session)
# valueBoxM2CutoffQualityScore (input, output, session)
# valueBoxM2SlidingWindowSize (input, output, session)
#
# valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
# valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)
#
# valueBoxRawSeqLength (input, output, session, trimmedRV)
# valueBoxRawMeanQualityScore (input, output, session, trimmedRV)
# valueBoxRawMinQualityScore (input, output, session, trimmedRV)
# valueBoxTrimmedStartPos (input, output, session, trimmedRV)
# valueBoxTrimmedFinishPos (input, output, session, trimmedRV)
# valueBoxTrimmedSeqLength (input, output, session, trimmedRV)
# valueBoxTrimmedMeanQualityScore (input, output, session, trimmedRV)
# valueBoxTrimmedMinQualityScore (input, output, session, trimmedRV)
# valueBoxRemainingRatio (input, output, session, trimmedRV)
#
#
# qualityTrimmingRatioPlot (input, output, session, trimmedRV,
#                           SangerSingleReadQualReport,
#                           SangerSingleReadFeature)
# qualityQualityBasePlot (input, output, session, trimmedRV,
#                         SangerSingleReadQualReport, SangerSingleReadFeature)
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Method 1 trimming parameter (M1TrimmingCutoff)
# ### ------------------------------------------------------------------------
# observeEvent(input$M1TrimmingCutoffText, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(singleReadIndex)) {
#         message("************ You have input ", input$M1TrimmingCutoffText,
#                 " in the 'M1 Trimming Cutoff' input box")
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Method 2 trimming parameter (M2CutoffQualityScore)
# ### ------------------------------------------------------------------------
# observeEvent(input$M2CutoffQualityScoreText, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(singleReadIndex)) {
#         message("************ You have input ",
#                 input$M2CutoffQualityScoreText,
#                 " in the 'M2 Cutoff Quality Score' input box")
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Method 2 trimming parameter (M2SlidingWindowSize)
# ### ------------------------------------------------------------------------
# observeEvent(input$M2SlidingWindowSizeText, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(singleReadIndex)) {
#         message("************ You have input ",
#                 input$M2SlidingWindowSizeText,
#                 " in the 'M2 Sliding Window Size' input box")
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### chromatogram related feature
# ### ------------------------------------------------------------------------
# output$chromatogramUIOutput <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (input$sidebar_menu != "Sanger Consensus Read Overview") {
#         if (!is.na(as.numeric(sidebar_menu[[1]]))) {
#             trimmedRV[["trimmedSeqLength"]]
#             chromatogramRowNumAns <-
#                 chromatogramRowNum (
#                     strtoi(ChromatogramParam[["baseNumPerRow"]]),
#                     SangerSingleReadQualReport[[singleReadIndex]]@
#                         rawSeqLength,
#                     SangerSingleReadQualReport[[singleReadIndex]]@
#                         trimmedSeqLength,
#                     ChromatogramParam[["showTrimmed"]]) *
#                 strtoi(ChromatogramParam[["heightPerRow"]])
#                 plotOutput("chromatogram", height = chromatogramRowNumAns)
#         }
#     }
# })
#
# output$chromatogram <- renderPlot({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (input$sidebar_menu != "Sanger Consensus Read Overview") {
#         if (!is.na(as.numeric(sidebar_menu[[1]]))) {
#             rawSeqLength <-
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 rawSeqLength
#             message(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function")
#             ### ------------------------------------------------------------
#             ### Re-run 'MakeBaseCall' function
#             ### ------------------------------------------------------------
#             forwardReadNum <- length((SangerConsensus)@forwardReadsList)
#             reverseReadNum <- length((SangerConsensus)@reverseReadsList)
#             SangerSingleReadNum <- forwardReadNum + reverseReadNum
#             if (singleReadIndex <= forwardReadNum) {
#                 # This is forward list
#                 index <- singleReadIndex
#                 hetcalls <-
#                     MakeBaseCalls(SangerConsensus@forwardReadsList[[
#                         index]], signalRatioCutoff = as.numeric(
#                             ChromatogramParam[["signalRatioCutoff"]]))
#
#                 SangerConsensus@forwardReadsList[[index]]@peakPosMatrix <-
#                     hetcalls@peakPosMatrix
#                 SangerConsensus@forwardReadsList[[index]]@peakAmpMatrix <-
#                     hetcalls@peakAmpMatrix
#                 SangerConsensus@forwardReadsList[[index]]@primarySeq <-
#                     hetcalls@primarySeq
#                 SangerConsensus@forwardReadsList[[index]]@secondarySeq <-
#                     hetcalls@secondarySeq
#
#                 ### --------------------------------------------------------
#                 ### Updating AASeqs
#                 ### --------------------------------------------------------
#                 AASeqResult <- calculateAASeq (hetcalls@primarySeq,
#                                                hetcalls@secondarySeq)
#                 SangerConsensus@forwardReadsList[[index]]@primaryAASeqS1 <-
#                     AASeqResult[["primaryAASeqS1"]]
#                 SangerConsensus@forwardReadsList[[index]]@primaryAASeqS2 <-
#                     AASeqResult[["primaryAASeqS2"]]
#                 SangerConsensus@forwardReadsList[[index]]@primaryAASeqS3 <-
#                     AASeqResult[["primaryAASeqS3"]]
#             } else {
#                 # This is reverse list
#                 index <- singleReadIndex-forwardReadNum
#                 hetcalls <-
#                     MakeBaseCalls(SangerConsensus@reverseReadsList[[
#                         index]], signalRatioCutoff = as.numeric(
#                             ChromatogramParam[["signalRatioCutoff"]]))
#
#                 SangerConsensus@reverseReadsList[[index]]@peakPosMatrix <-
#                     hetcalls@peakPosMatrix
#                 SangerConsensus@reverseReadsList[[index]]@peakAmpMatrix <-
#                     hetcalls@peakAmpMatrix
#                 SangerConsensus@reverseReadsList[[index]]@primarySeq <-
#                     hetcalls@primarySeq
#                 SangerConsensus@reverseReadsList[[index]]@secondarySeq <-
#                     hetcalls@secondarySeq
#
#                 ### --------------------------------------------------------
#                 ### Updating AASeqs
#                 ### --------------------------------------------------------
#                 AASeqResult <- calculateAASeq (hetcalls@primarySeq,
#                                                hetcalls@secondarySeq)
#                 SangerConsensus@reverseReadsList[[index]]@primaryAASeqS1 <-
#                     AASeqResult[["primaryAASeqS1"]]
#                 SangerConsensus@reverseReadsList[[index]]@primaryAASeqS2 <-
#                     AASeqResult[["primaryAASeqS2"]]
#                 SangerConsensus@reverseReadsList[[index]]@primaryAASeqS3 <-
#                     AASeqResult[["primaryAASeqS3"]]
#             }
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
#             message(">>>>>>>>>>>> 'MakeBaseCalls' finished")
#             chromatogram(hetcalls,
#                          width = strtoi(
#                              ChromatogramParam[["baseNumPerRow"]]),
#                          height = 2,
#                          trim5 = trimmedRV[["trimmedStartPos"]],
#                          trim3 = rawSeqLength -
#                              trimmedRV[["trimmedFinishPos"]],
#                          showtrim = (ChromatogramParam[["showTrimmed"]]),
#                          showcalls = "both")
#             }
#
#         }
# })
#
#
# ############################################################################
# ### Switch trimming method related function
# ############################################################################
# output$TrimmingMethodUI <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.null(SangerSingleReadQualReport[[singleReadIndex]])) {
#         if (SangerSingleReadQualReport[[singleReadIndex]]@
#             TrimmingMethod == "M1") {
#             if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@
#                         M1TrimmingCutoff)) {
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                     M1TrimmingCutoff <<-  0.0001
#             }
#             fluidRow(
#                 column(6,
#                        uiOutput("M1TrimmingCutoff") ,
#                        tags$ul(
#                            textInput("M1TrimmingCutoffText",
#                                      label = p("Input Value"),
#                                      value = toString(
#                                          trimmedParam[["M1TrimmingCutoff"]]),
#                                      width = '70%')
#                        ),
#                 ),
#             )
#         } else if (SangerSingleReadQualReport[[singleReadIndex]]@
#                    TrimmingMethod == "M2") {
#             message("Inside Method 2!!")
#             if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@
#                         M2CutoffQualityScore)) {
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                     M2CutoffQualityScore <<-  20
#             }
#             if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@
#                         M2SlidingWindowSize )) {
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                     M2SlidingWindowSize <<-  5
#             }
#
#             fluidRow(
#                 column(6,
#                        uiOutput("M2CutoffQualityScore") ,
#                        tags$ul(
#                            textInput("M2CutoffQualityScoreText",
#                                      label = p("Input Value"),
#                                      value = toString(
#                                          trimmedParam[["M2CutoffQualityScore"]]),
#                                      width = '70%')
#                        ),
#                 ),
#                 column(6,
#                        uiOutput("M2SlidingWindowSize") ,
#                        tags$ul(
#                            textInput("M2SlidingWindowSizeText",
#                                      label = p("Input Value"),
#                                      value = toString(
#                                          trimmedParam[["M2SlidingWindowSize"]]),
#                                      width = '70%')
#                        ),
#                 ),
#             )
#         }
#     }
# })
#
# output$TrimmingMethodSelectionOutput <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     singleReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.null(SangerSingleReadQualReport[[singleReadIndex]])) {
#         if (SangerSingleReadQualReport[[singleReadIndex]]@
#             TrimmingMethod == "M1") {
#             tagList(icon("check-circle"),
#                     "Your trimming method selection :
#                     'Logarithmic Scale Trimming'")
#         } else if (SangerSingleReadQualReport[[singleReadIndex]]@
#                    TrimmingMethod == "M2") {
#             tagList(icon("check-circle"),
#                     "Your trimming method selection :
#                     'Logarithmic Scale Sliding Window Trimming'")
#         }
#     }
# })
