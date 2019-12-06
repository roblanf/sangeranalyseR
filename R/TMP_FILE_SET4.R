#
# ############################################################################
# ### All other features (dynamic header / button save / button close)
# ############################################################################
# observeEvent(input$sidebar_menu, {
#     menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
#     html("rightHeader", menuItem)
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Save S4 object
# ### ------------------------------------------------------------------------
# observeEvent(input$saveS4, {
#     newS4Object <- file.path(shinyDirectory,
#                              "SangerAlignedConsensusSet.Rda")
#     showNotification(paste("New S4 object is store as:", newS4Object),
#                      type = "message", duration = 10)
#     sapply(1:SangerConsensusSetNum, function(i) {
#         forwardReadNum <-
#             length(SangerConsensusSet@
#                        consensusReadsList[[i]]@forwardReadsList)
#         reverseReadNum <-
#             length(SangerConsensusSet@
#                        consensusReadsList[[i]]@reverseReadsList)
#
#         sapply(1:forwardReadNum, function(j) {
#             SangerConsensusSet@consensusReadsList[[i]]@
#                 forwardReadsList[[j]]@QualityReport <<-
#                 SangerCSetParam[[i]]$SangerSingleReadQualReport[[j]]
#             SangerConsensusSet@consensusReadsList[[i]]@
#                 forwardReadsList[[j]]@ChromatogramParam <<-
#                 SangerCSetParam[[i]]$
#                 SangerSingleReadChromatogramParam[[j]]
#             message("save SangerConsensus quality S4 object Forward")
#         })
#         sapply(1:reverseReadNum, function(j) {
#             SangerConsensusSet@consensusReadsList[[i]]@
#                 reverseReadsList[[j]]@QualityReport <<-
#                 SangerCSetParam[[i]]$
#                 SangerSingleReadQualReport[[forwardReadNum + j]]
#             SangerConsensusSet@consensusReadsList[[i]]@
#                 reverseReadsList[[j]]@ChromatogramParam <<-
#                 SangerCSetParam[[i]]$
#                 SangerSingleReadChromatogramParam[[forwardReadNum+j]]
#             message("save SangerConsensus quality S4 object Reverse")
#         })
#     })
#     saveRDS(SangerConsensusSet, file=newS4Object)
#     message("New S4 object is store as: ", newS4Object)
#     NEW_SANGER_ALIGNED_CONSENSUS_READ_SET <<- readRDS(file=newS4Object)
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
#     message("######## Start recalculating consensus read (SCSet)")
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     CSResult<-
#         calculateConsensusRead (
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@forwardReadsList,
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@reverseReadsList,
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@refAminoAcidSeq,
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@minFractionCall,
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@maxFractionLost,
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@geneticCode,
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@acceptStopCodons,
#             SangerConsensusSet@
#                 consensusReadsList[[consensusReadIndex]]@readingFrame)
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@consensusRead <<-
#         CSResult$consensusGapfree
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@differencesDF <<-
#         CSResult$diffsDf
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@alignment <<-
#         CSResult$aln2
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@distanceMatrix <<-
#         CSResult$dist
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@dendrogram <<-
#         CSResult$dend
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@indelsDF <<-
#         CSResult$indels
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@stopCodonsDF <<-
#         CSResult$stopsDf
#     SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@secondaryPeakDF <<-
#         CSResult$spDf
#
#     consensusParam[["consensusReadName"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@consensusReadName
#     consensusParam[["consensusRead"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@consensusRead
#     consensusParam[["differencesDF"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@differencesDF
#     consensusParam[["alignment"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@alignment
#     consensusParam[["distanceMatrix"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@distanceMatrix
#     consensusParam[["dendrogram"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@dendrogram
#     consensusParam[["indelsDF"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@indelsDF
#     consensusParam[["stopCodonsDF"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@stopCodonsDF
#     consensusParam[["secondaryPeakDF"]] <<-
#         SangerConsensusSet@
#         consensusReadsList[[consensusReadIndex]]@secondaryPeakDF
#     message("######## Finish recalculation")
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Consensus read re-calculating (SCSet) UI
# ### ------------------------------------------------------------------------
# observeEvent(input$recalculateButtonSCSet, {
#     message("######## Reactive button clicked !!!")
#     message("######## Start recalculating consensus read (SC)")
#     if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
#         CSSetResult <-
#             alignConsensusReads (SangerConsensusSet@consensusReadsList,
#                                  SangerConsensusSet@geneticCode,
#                                  SangerConsensusSet@refAminoAcidSeq,
#                                  SangerConsensusSet@minFractionCallSCSet,
#                                  SangerConsensusSet@maxFractionLostSCSet,
#                                  1)
#
#         SangerConsensusSet@consensusReadSCSet <<- CSSetResult$consensus
#         SangerConsensusSet@alignmentSCSet <<- CSSetResult$aln
#         SangerConsensusSet@alignmentTreeSCSet <<- CSSetResult$aln.tree
#
#         consensusParamSet[["consensusReadSCSet"]] <<- SangerConsensusSet@consensusReadSCSet
#         consensusParamSet[["alignmentSCSet"]] <<- SangerConsensusSet@alignmentSCSet
#         consensusParamSet[["alignmentTreeSCSet"]] <<- SangerConsensusSet@alignmentTreeSCSet
#         message("######## Finish recalculation consensus read")
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
# ### ------------------------------------------------------------------------
# observeEvent(input$startTrimmingButton, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         if (SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             TrimmingMethod == "M1") {
#
#             if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
#                 as.numeric(input$M1TrimmingCutoffText) > 0 &&
#                 as.numeric(input$M1TrimmingCutoffText) <= 1) {
#                 inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
#             } else {
#                 inputM1TrimmingCutoffText <- 0.0001
#             }
#             trimmingPos <-
#                 M1inside_calculate_trimming(
#                     SangerCSetParam[[consensusReadIndex]]$
#                         SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityPhredScores,
#                     SangerCSetParam[[consensusReadIndex]]$
#                         SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityBaseScores,
#                     as.numeric(inputM1TrimmingCutoffText))
#             SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)
#             trimmedParam[["M1TrimmingCutoff"]] <<-
#                 SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M1TrimmingCutoff
#
#         } else if (SangerCSetParam[[consensusReadIndex]]$
#                    SangerSingleReadQualReport[[singleReadIndex]]@
#                    TrimmingMethod == "M2") {
#
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
#
#             trimmingPos <-
#                 M2inside_calculate_trimming(
#                     SangerCSetParam[[consensusReadIndex]]$
#                         SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityPhredScores,
#                     SangerCSetParam[[consensusReadIndex]]$
#                         SangerSingleReadQualReport[[singleReadIndex]]@
#                         qualityBaseScores,
#                     strtoi(inputM2CutoffQualityScoreText),
#                     strtoi(inputM2SlidingWindowSizeText)
#                 )
#             SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2CutoffQualityScore <<- strtoi(inputM2CutoffQualityScoreText)
#             trimmedParam[["M2CutoffQualityScore"]] <<-
#                 SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2CutoffQualityScore
#
#             SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)
#             trimmedParam[["M2SlidingWindowSize"]] <<-
#                 SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]@
#                 M2SlidingWindowSize
#         }
#
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             rawSeqLength <<- trimmingPos[["rawSeqLength"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
#         SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             remainingRatio <<- trimmingPos[["remainingRatio"]]
#
#         trimmedRV[["rawSeqLength"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
#         trimmedRV[["rawMeanQualityScore"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             rawMeanQualityScore
#         trimmedRV[["rawMinQualityScore"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@rawMinQualityScore
#         trimmedRV[["trimmedStartPos"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@trimmedStartPos
#         trimmedRV[["trimmedFinishPos"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@trimmedFinishPos
#         trimmedRV[["trimmedSeqLength"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
#         trimmedRV[["trimmedMeanQualityScore"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMeanQualityScore
#         trimmedRV[["trimmedMinQualityScore"]] <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@
#             trimmedMinQualityScore
#         trimmedRV[["remainingRatio"]] <<-
#             round(SangerCSetParam[[consensusReadIndex]]$
#                       SangerSingleReadQualReport[[singleReadIndex]]@
#                       remainingRatio * 100, 2)
#         ### ------------------------------------------------------------
#         ### Save SangerConsensus quality S4 object
#         ### ------------------------------------------------------------
#         forwardReadNum <-
#             length(SangerConsensusSet@
#                        consensusReadsList[[consensusReadIndex]]@forwardReadsList)
#         reverseReadNum <-
#             length(SangerConsensusSet@
#                        consensusReadsList[[consensusReadIndex]]@reverseReadsList)
#         SangerSingleReadNum <- forwardReadNum + reverseReadNum
#         if (singleReadIndex <= forwardReadNum) {
#             # This is forward list
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[singleReadIndex]]@QualityReport <<-
#                 SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]
#         } else {
#             # This is reverse list
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[singleReadIndex-forwardReadNum]]@QualityReport <<-
#                 SangerCSetParam[[consensusReadIndex]]$
#                 SangerSingleReadQualReport[[singleReadIndex]]
#         }
#     }
#
# })
#
# ### ------------------------------------------------------------------------
# ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
# ### ------------------------------------------------------------------------
# observeEvent(input$saveChromatogramParam, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     ### ------------------------------------------------------------
#     ### Update ChromatogramBasePerRow
#     ### ------------------------------------------------------------
#     SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         baseNumPerRow <<- input$ChromatogramBasePerRow
#     SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         heightPerRow <<- input$ChromatogramHeightPerRow
#     SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
#     SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         showTrimmed <<- input$ChromatogramCheckShowTrimmed
#
#     ### ------------------------------------------------------------
#     ### Save SangerConsensus quality S4 object
#     ### ------------------------------------------------------------
#     forwardReadNum <-
#         length(SangerConsensusSet@
#                    consensusReadsList[[consensusReadIndex]]@
#                    forwardReadsList)
#     reverseReadNum <-
#         length(SangerConsensusSet@
#                    consensusReadsList[[consensusReadIndex]]@
#                    reverseReadsList)
#     SangerSingleReadNum <- forwardReadNum + reverseReadNum
#
#     if (singleReadIndex <= forwardReadNum) {
#         # This is forward list
#         SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#             forwardReadsList[[singleReadIndex]]@ChromatogramParam <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadChromatogramParam[[singleReadIndex]]
#     } else {
#         # This is reverse list
#         SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#             reverseReadsList[[singleReadIndex - forwardReadNum]]@ChromatogramParam <<-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadChromatogramParam[[singleReadIndex]]
#     }
#
#     ChromatogramParam[["baseNumPerRow"]] <<-
#         SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         baseNumPerRow
#     ChromatogramParam[["heightPerRow"]] <<-
#         SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         heightPerRow
#     ChromatogramParam[["signalRatioCutoff"]] <<-
#         SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         signalRatioCutoff
#     ChromatogramParam[["showTrimmed"]] <<-
#         SangerCSetParam[[consensusReadIndex]]$
#         SangerSingleReadChromatogramParam[[singleReadIndex]]@
#         showTrimmed
# })
#
# ############################################################################
# ### ConsensusReadSet (Function for Sanger Consensus Read Set Overview)
# ############################################################################
# ### ------------------------------------------------------------------------
# ### Alignment
# ### ------------------------------------------------------------------------
# output$consensusSetAlignmentHTML<-renderUI({
#     if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
#         # consensusParamSet[["alignmentSCSet"]] <-
#         #     SangerConsensusSet@alignmentSCSet
#         browseSeqHTML <-
#             file.path(shinyDirectory,
#                       "Consensus_Readset_Alignment_BrowseSeqs.html")
#         BrowseSeqs(consensusParamSet[["alignmentSCSet"]] ,
#                    openURL=FALSE, htmlFile=browseSeqHTML)
#         includeHTML(browseSeqHTML)
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### Consensus Reads Tree
# ### ------------------------------------------------------------------------
# output$SCSetConsensusReadTreePlot <- renderPlot({
#     if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
#         # consensusParamSet[["alignmentTreeSCSet"]] <<- SangerConsensusSet@alignmentTreeSCSet
#         plot(consensusParamSet[["alignmentTreeSCSet"]])
#     }
# })
#
# ############################################################################
# ### ConsensusRead (Function for Sanger Consensus Read Overview)
# ############################################################################
# ### ------------------------------------------------------------------------
# ### genetic code
# ### ------------------------------------------------------------------------
# output$geneticCodeDF <- renderExcel({
#     SCGeneticCode <- SangerConsensusSet@geneticCode
#     suppressMessages(
#         excelTable(data = t(data.frame(SCGeneticCode)),
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
#     if (SangerConsensusSet@refAminoAcidSeq == "") {
#         box(title = tags$p("Reference Amino Acids Sequence",
#                            style = "font-size: 24px;
#                                     font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 1),
#             column(width = 11,
#                    h4("Reference Amino Acid Sequence is not provided."))
#         )
#     } else {
#         box(title = tags$p("Reference Amino Acids Sequence",
#                            style = "font-size: 24px;
#                                     font-weight: bold;"),
#             collapsible = TRUE,
#             status = "success", width = 12,
#             column(width = 2,
#                    tags$br(),
#                    tags$p("AA Sequence:",
#                           style = "font-size: 15px;
#                                        font-weight: bold;"),
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
#
# output$SCrefAminoAcidSeqDF <- renderExcel({
#     refAminoAcidSeqVec <-
#         strsplit(SangerConsensusSet@refAminoAcidSeq, "")[[1]]
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
# ### Alignment
# ### ------------------------------------------------------------------------
# output$consensusAlignmentHTML<-renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["consensusReadName"]] <<-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@consensusReadName
#         consensusParam[["alignment"]] <-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@alignment
#
#         browseSeqHTML <-
#             file.path(shinyDirectory, "BrowseSeqs_html",
#                       paste0(sidebar_menu[[1]], "_",
#                              consensusParam[["consensusReadName"]],
#                              "_Alignment_BrowseSeqs.html"))
#         if (!dir.exists(file.path(shinyDirectory, "BrowseSeqs_html"))) {
#             dir.create(file.path(shinyDirectory, "BrowseSeqs_html"))
#         }
#         BrowseSeqs(consensusParam[["alignment"]],
#                    openURL=FALSE, htmlFile=browseSeqHTML)
#         includeHTML(browseSeqHTML)
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### difference
# ### ------------------------------------------------------------------------
# output$SCDifferencesDFUI <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["differencesDF"]] <<-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@differencesDF
#         if (all(dim(consensusParam[["differencesDF"]]) == c(0,0))) {
#             h4("*** 'Differences' dataframe is empty. ***",
#                style="font-weight: bold; font-style: italic;")
#         } else {
#             dataTableOutput("SCDifferencesDF")
#         }
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
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["dendrogram"]] <<-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@dendrogram
#         plot(consensusParam[["dendrogram"]][[2]])
#         ggdendrogram(consensusParam[["dendrogram"]][[2]], rotate = TRUE)
#     }
# })
# output$dendrogramDF <- renderDataTable({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["dendrogram"]] <<-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@dendrogram
#         consensusParam[["dendrogram"]][[1]]
#     }
# })
#
# ### ------------------------------------------------------------------------
# ### distance
# ### ------------------------------------------------------------------------
# output$SCDistanceMatrixPlotUI <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["distanceMatrix"]] <<-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@distanceMatrix
#         if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
#             h4("*** 'Distance' dataframe is empty. (Cannot plot)***",
#                style="font-weight: bold; font-style: italic;")
#         } else {
#             plotlyOutput("SCDistanceMatrixPlot")
#         }
#     }
# })
#
# output$SCDistanceMatrixPlot <- renderPlotly({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         SangerSingleReadBFN <-
#             SangerCSetParam[[consensusReadIndex]]$SangerSingleReadBFN
#         suppressPlotlyMessage(
#             plot_ly(x = SangerSingleReadBFN,
#                     y = SangerSingleReadBFN,
#                     z = consensusParam[["distanceMatrix"]],
#                     colors = colorRamp(c("white", "#32a852")),
#                     type = "heatmap")
#         )
#     }
# })
# output$SCDistanceMatrixUI <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["distanceMatrix"]] <<-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@distanceMatrix
#         if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
#             h4("*** 'Distance' dataframe is empty. ***",
#                style="font-weight: bold; font-style: italic;")
#         } else {
#             dataTableOutput("SCDistanceMatrix")
#         }
#     }
# })
# output$SCDistanceMatrix = renderDataTable({
#     consensusParam[["distanceMatrix"]]
# })
#
# ### ------------------------------------------------------------------------
# ### SCIndelsDF
# ### ------------------------------------------------------------------------
# output$SCIndelsDFUI <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["indelsDF"]] <<-
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@indelsDF
#         if (all(dim(consensusParam[["indelsDF"]] ) == c(0,0))) {
#             h4("*** 'Indels' data frame is empty. ***",
#                style="font-weight: bold; font-style: italic;")
#         } else {
#             dataTableOutput("SCIndelsDF")
#         }
#     }
# })
# output$SCIndelsDF <- renderDataTable({
#     consensusParam[["indelsDF"]]
# })
#
# ### ------------------------------------------------------------------------
# ### SCStopCodons
# ### ------------------------------------------------------------------------
# output$SCStopCodonsDFUI <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     if (!is.na(consensusReadIndex)) {
#         consensusParam[["stopCodonsDF"]] <<-
#             SangerConsensusSet@
#             consensusReadsList[[consensusReadIndex]]@stopCodonsDF
#         if (all(dim(consensusParam[["stopCodonsDF"]]) == c(0,0))) {
#             h4("*** 'Stop Codons' dataframe is empty. ***",
#                style="font-weight: bold; font-style: italic;")
#         } else {
#             dataTableOutput("SCStopCodonsDF")
#         }
#     }
# })
# output$SCStopCodonsDF <- renderDataTable({
#     consensusParam[["stopCodonsDF"]]
# })
#
# ############################################################################
# ### SangerSingleRead (Function for singel read in consensusRead)
# ############################################################################
# output$primarySeqDF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         AstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "A", "#1eff00")
#         TstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "T", "#ff7a7a")
#         CstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "C", "#7ac3ff")
#         GstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "G", "#c9c9c9")
#         styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
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
#         suppressMessages(
#             excelTable(data =SangerCSetParam[[consensusReadIndex]]$
#                            SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                        defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
#                        columnResize = FALSE, allowInsertRow = FALSE,
#                        allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                        allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                        style = styleList, loadingSpin = TRUE)
#         )
#     }
# })
#
# output$secondSeqDF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         AstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "A", "#1eff00")
#         TstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "T", "#ff7a7a")
#         CstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "C", "#7ac3ff")
#         GstyleList <-
#             getStopList(SangerCSetParam[[consensusReadIndex]]$
#                             SangerSingleReadPrimSeqDF[[singleReadIndex]],
#                         "G", "#c9c9c9")
#         styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
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
#         suppressMessages(
#             excelTable(data =
#                            SangerCSetParam[[consensusReadIndex]]$
#                            SangerSingleReadSecoSeqDF[[singleReadIndex]],
#                        defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
#                        columnResize = FALSE, allowInsertRow = FALSE,
#                        allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                        allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                        style = styleList, loadingSpin = TRUE)
#         )
#     }
# })
#
# output$qualityScoreDF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         suppressMessages(
#             excelTable(data =
#                            SangerCSetParam[[consensusReadIndex]]$
#                            SangerSingleReadQSDF[[singleReadIndex]],
#                        defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
#                        columnResize = FALSE, allowInsertRow = FALSE,
#                        allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                        allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
#         )
#     }
# })
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
# output$PrimAASeqS1DF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         width <- rep(90, length(SangerCSetParam[[consensusReadIndex]]$
#                                     SangerSingleReadPrimAASeqS1DF[[singleReadIndex]]))
#         styleList <-
#             getStopList (SangerCSetParam[[consensusReadIndex]]$
#                              SangerSingleReadPrimAASeqS1DF[[singleReadIndex]],
#                          "*", "#cf0000")
#         suppressMessages(
#             excelTable(data =
#                            SangerCSetParam[[consensusReadIndex]]$
#                            SangerSingleReadPrimAASeqS1DF[[singleReadIndex]],
#                        columns = data.frame(width = width),
#                        defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
#                        columnResize = FALSE, allowInsertRow = FALSE,
#                        allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                        allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                        style = styleList, loadingSpin = TRUE)
#         )
#     }
# })
# output$PrimAASeqS2DF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         width <- rep(90, length(SangerCSetParam[[consensusReadIndex]]$
#                                     SangerSingleReadPrimAASeqS2DF[[singleReadIndex]])-1)
#         widthFinal <- c(30, width)
#         styleList <-
#             getStopList (SangerCSetParam[[consensusReadIndex]]$
#                              SangerSingleReadPrimAASeqS2DF[[singleReadIndex]],
#                          "*", "#cf0000")
#         styleList[['A1']] <- 'background-color: black;'
#         suppressMessages(
#             excelTable(data =
#                            SangerCSetParam[[consensusReadIndex]]$
#                            SangerSingleReadPrimAASeqS2DF[[singleReadIndex]],
#                        columns = data.frame(width = widthFinal),
#                        defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
#                        columnResize = FALSE, allowInsertRow = FALSE,
#                        allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                        allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                        style = styleList, loadingSpin = TRUE)
#         )
#     }
# })
# output$PrimAASeqS3DF <- renderExcel({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         width <- rep(90, length(SangerCSetParam[[consensusReadIndex]]$
#                                     SangerSingleReadPrimAASeqS3DF[[singleReadIndex]])-2)
#         widthFinal <- c(30, 30, width)
#         styleList <-
#             getStopList (SangerCSetParam[[consensusReadIndex]]$
#                              SangerSingleReadPrimAASeqS3DF[[singleReadIndex]],
#                          "*", "#cf0000")
#         styleList[['A1']] <- 'background-color: black;'
#         styleList[['B1']] <- 'background-color: black;'
#         suppressMessages(
#             excelTable(data =
#                            SangerCSetParam[[consensusReadIndex]]$
#                            SangerSingleReadPrimAASeqS3DF[[singleReadIndex]],
#                        columns = data.frame(width = widthFinal),
#                        defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
#                        columnResize = FALSE, allowInsertRow = FALSE,
#                        allowInsertColumn = FALSE, allowDeleteRow = FALSE,
#                        allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
#                        style = styleList, loadingSpin = TRUE)
#         )
#     }
# })
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
# valueBoxSCMinReadsNumCSSet (input, output, SangerConsensusSet, session)
# valueBoxSCMinReadLengthCSSet (input, output, SangerConsensusSet, session)
# valueBoxSCMinFractionCallCSSet (input, output, SangerConsensusSet, session)
# valueBoxSCMaxFractionLostCSSet (input, output, SangerConsensusSet, session)
# valueBoxSCAcceptStopCodonsCSSet (input, output, SangerConsensusSet, session)
# valueBoxSCReadingFrameCSSet (input, output, SangerConsensusSet, session)
#
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
# observeEvent(input$M1TrimmingCutoffText, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         message("************ You have input ",
#                 input$M1TrimmingCutoffText,
#                 " in the 'M1 Trimming Cutoff' input box")
#     }
# })
#
# observeEvent(input$M2CutoffQualityScoreText, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         message("************ You have input ",
#                 input$M2CutoffQualityScoreText,
#                 " in the 'M2 Cutoff Quality Score' input box")
#     }
# })
#
# observeEvent(input$M2SlidingWindowSizeText, {
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         message("************ You have input ",
#                 input$M2SlidingWindowSizeText,
#                 " in the 'M2 Sliding Window Size' input box")
#     }
# })
#
# output$qualityTrimmingRatioPlot <- renderPlotly({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         readFeature <- SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadFeature[[singleReadIndex]]
#         trimmedStartPos = trimmedRV[["trimmedStartPos"]]
#         trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
#         qualityPhredScores <-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
#         readLen = length(qualityPhredScores)
#
#         stepRatio = 1 / readLen
#         trimmedStartPos / readLen
#         trimmedFinishPos / readLen
#
#         trimmedPer <- c()
#         remainingPer <- c()
#
#         for (i in 1:trimmedStartPos) {
#             if (i != trimmedStartPos) {
#                 trimmedPer <- c(trimmedPer, stepRatio)
#             }
#         }
#
#         for (i in trimmedStartPos:trimmedFinishPos) {
#             trimmedPer <- c(trimmedPer, 0)
#         }
#
#
#         for (i in trimmedFinishPos:readLen) {
#             if (i != trimmedFinishPos) {
#                 trimmedPer <- c(trimmedPer, stepRatio)
#             }
#         }
#
#         trimmedPer <- cumsum(trimmedPer)
#         remainingPer = 1 - trimmedPer
#
#         PerData <- data.frame(1:length(trimmedPer),
#                               trimmedPer, remainingPer)
#
#         colnames(PerData) <- c("Base",
#                                "Trimmed Ratio",
#                                "Remaining Ratio")
#         x <- list(
#             title = "Base Pair Index"
#             # titlefont = f
#         )
#         y <- list(
#             title = "Read Ratio"
#             # titlefont = f
#         )
#         PerDataPlot <- melt(PerData, id.vars = c("Base"))
#         suppressPlotlyMessage(
#             plot_ly(data=PerDataPlot,
#                     x=~Base,
#                     y=~value,
#                     mode="markers",
#                     color = ~variable,
#                     text = ~paste("BP Index : ",
#                                   Base, '<sup>th</sup><br>Read Ratio :',
#                                   round(value*100, digits = 2), '%')) %>%
#                 layout(xaxis = x,
#                        yaxis = y,
#                        legend = list(orientation = 'h',
#                                      xanchor = "center",
#                                      x = 0.5, y = 1.1)) %>%
#                 add_annotations(
#                     text = "Trimmed Ratio (Each BP)",
#                     x = (trimmedStartPos + trimmedFinishPos) / 2,
#                     y = ((trimmedPer[1] + trimmedPer[length(trimmedPer)]) / 2)
#                     + 0.06,
#                     showarrow=FALSE
#                 ) %>%
#                 add_annotations(
#                     text = "Remaining Ratio (Each BP)",
#                     x = (trimmedStartPos+trimmedFinishPos) / 2,
#                     y = ((remainingPer[1]+remainingPer[length(remainingPer)])/2)
#                     - 0.06,
#                     showarrow=FALSE
#                 ))
#     }
# })
#
# output$qualityQualityBasePlot <- renderPlotly({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         readFeature <-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadFeature[[singleReadIndex]]
#         trimmedStartPos = trimmedRV[["trimmedStartPos"]]
#         trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
#         qualityPhredScores <-
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
#         readLen = length(qualityPhredScores)
#
#         qualityPlotDf<- data.frame(1:length(qualityPhredScores),
#                                    qualityPhredScores)
#         colnames(qualityPlotDf) <- c("Index", "Score")
#         x <- list(
#             title = "Base Pair Index"
#             # titlefont = f
#         )
#         y <- list(
#             title = "Phred Quality Score"
#             # titlefont = f
#         )
#         suppressPlotlyMessage(
#             plot_ly(data=qualityPlotDf,
#                     x=~Index) %>%
#                 add_markers(y=~Score,
#                             text =
#                                 ~paste("BP Index : ",
#                                        Index,
#                                        '<sup>th</sup><br>Phred Quality Score :',
#                                        Score),
#                             name = 'Quality Each BP') %>%
#                 add_trace(x=seq(trimmedStartPos,
#                                 trimmedFinishPos,
#                                 len=trimmedFinishPos-trimmedStartPos+1),
#                           y=rep(70, trimmedFinishPos-trimmedStartPos+1),
#                           mode="lines", hoverinfo="text",
#                           text=paste("Trimmed Reads BP length:",
#                                      trimmedFinishPos-trimmedStartPos+1,
#                                      "BPs <br>",
#                                      "Trimmed Reads BP ratio:",
#                                      round((trimmedFinishPos-trimmedStartPos+1)/
#                                                readLen * 100,
#                                            digits=2),
#                                      "%"),
#                           line = list(width = 12),
#                           name = 'Trimmed Read') %>%
#                 add_trace(x=seq(0,readLen,len=readLen),
#                           y=rep(80, readLen), mode="lines", hoverinfo="text",
#                           text=paste("Whole Reads BP length:",
#                                      readLen,
#                                      "BPs <br>",
#                                      "Trimmed Reads BP ratio: 100 %"),
#                           line = list(width = 12),
#                           name = 'Whole Read') %>%
#                 layout(xaxis = x, yaxis = y,
#                        shapes = list(vline(trimmedStartPos),
#                                      vline(trimmedFinishPos)),
#                        legend = list(orientation = 'h',
#                                      xanchor = "center",
#                                      x = 0.5, y = 1.1)) %>%
#                 add_annotations(
#                     text = "Trimming Strat <br> BP Index",
#                     x = trimmedStartPos + 40,
#                     y = 15,
#                     showarrow=FALSE
#                 ) %>%
#                 add_annotations(
#                     text = "Trimming End <br> BP Index",
#                     x = trimmedFinishPos - 40,
#                     y = 15,
#                     showarrow=FALSE
#                 ))
#     }
# })
#
# valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
# valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)
#
# # chromatogram
# output$chromatogramUIOutput <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         sidebar_menu[[2]] == "Consensus" &&
#         sidebar_menu[[3]] == "Read" &&
#         sidebar_menu[[4]] == "-" &&
#         !is.na(singleReadIndex) &&
#         (sidebar_menu[[6]] == "Forward" ||
#          sidebar_menu[[6]] == "Reverse") &&
#         sidebar_menu[[7]] == "Read") {
#         trimmedRV[["trimmedSeqLength"]]
#         rawSeqLength =
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
#         trimmedSeqLength =
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
#
#         chromatogramRowNumAns <- chromatogramRowNum (
#             strtoi(ChromatogramParam[["baseNumPerRow"]]),
#             rawSeqLength, trimmedSeqLength,
#             ChromatogramParam[["showTrimmed"]]) *
#             strtoi(ChromatogramParam[["heightPerRow"]])
#         plotOutput("chromatogram", height = chromatogramRowNumAns)
#     }
# })
#
# output$chromatogram <- renderPlot({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         sidebar_menu[[2]] == "Consensus" &&
#         sidebar_menu[[3]] == "Read" &&
#         sidebar_menu[[4]] == "-" &&
#         !is.na(singleReadIndex) &&
#         (sidebar_menu[[6]] == "Forward" || sidebar_menu[[6]] == "Reverse")&&
#         sidebar_menu[[7]] == "Read") {
#         rawSeqLength =
#             SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
#
#         #### Should be test when adding dynamic secondary peak
#         ### ------------------------------------------------------------
#         ### Save SangerConsensus quality S4 object
#         ### ------------------------------------------------------------
#         forwardReadNum <-
#             length(SangerConsensusSet@
#                        consensusReadsList[[consensusReadIndex]]@forwardReadsList)
#         reverseReadNum <-
#             length(SangerConsensusSet@
#                        consensusReadsList[[consensusReadIndex]]@reverseReadsList)
#         SangerSingleReadNum <- forwardReadNum + reverseReadNum
#         if (singleReadIndex <= forwardReadNum) {
#             # This is forward list
#             index <- singleReadIndex
#             hetcalls <-
#                 MakeBaseCalls(
#                     SangerConsensusSet@
#                         consensusReadsList[[consensusReadIndex]]@
#                         forwardReadsList[[index]],
#                     signalRatioCutoff = as.numeric(
#                         ChromatogramParam[["signalRatioCutoff"]]))
#
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@QualityReport@qualityScoresID <-
#                 hetcalls@QualityReport@qualityScoresID
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@
#                 QualityReport@qualityPhredScores <-
#                 hetcalls@QualityReport@qualityPhredScores
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@
#                 QualityReport@qualityBaseScoresRaw <-
#                 hetcalls@QualityReport@qualityBaseScoresRaw
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@peakPosMatrix <-
#                 hetcalls@peakPosMatrix
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@peakAmpMatrix <-
#                 hetcalls@peakAmpMatrix
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@primarySeqID <-
#                 hetcalls@primarySeqID
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@primarySeq <-
#                 hetcalls@primarySeq
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@secondarySeqID <-
#                 hetcalls@secondarySeqID
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 forwardReadsList[[index]]@secondarySeq <-
#                 hetcalls@secondarySeq
#         } else {
#             # This is reverse list
#             index <- singleReadIndex-forwardReadNum
#             hetcalls <-
#                 MakeBaseCalls(
#                     SangerConsensusSet@
#                         consensusReadsList[[consensusReadIndex]]@
#                         reverseReadsList[[index]],
#                     signalRatioCutoff = as.numeric(
#                         ChromatogramParam[["signalRatioCutoff"]]))
#
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@
#                 QualityReport@qualityScoresID <-
#                 hetcalls@QualityReport@qualityScoresID
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@
#                 QualityReport@qualityPhredScores <-
#                 hetcalls@QualityReport@qualityPhredScores
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@
#                 QualityReport@qualityBaseScoresRaw <-
#                 hetcalls@QualityReport@qualityBaseScoresRaw
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@peakPosMatrix <-
#                 hetcalls@peakPosMatrix
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@peakAmpMatrix <-
#                 hetcalls@peakAmpMatrix
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@primarySeqID <-
#                 hetcalls@primarySeqID
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@primarySeq <-
#                 hetcalls@primarySeq
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@secondarySeqID <-
#                 hetcalls@secondarySeqID
#             SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
#                 reverseReadsList[[index]]@secondarySeq <-
#                 hetcalls@secondarySeq
#         }
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
#         chromatogram(hetcalls,
#                      width = strtoi(ChromatogramParam[["baseNumPerRow"]]),
#                      height = 2, trim5 = trimmedRV[["trimmedStartPos"]],
#                      trim3 = rawSeqLength - trimmedRV[["trimmedFinishPos"]],
#                      showtrim = (ChromatogramParam[["showTrimmed"]]),
#                      showcalls = "both")
#     }
# })
#
# output$TrimmingMethodUI <- renderUI({
#     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         if ( SangerCSetParam[[consensusReadIndex]]$
#              SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod== "M1") {
#             if (is.null(SangerCSetParam[[consensusReadIndex]]$
#                         SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff)) {
#                 SangerCSetParam[[consensusReadIndex]]$
#                     SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff <<-  0.0001
#             }
#             fluidRow(
#                 column(6,
#                        uiOutput("M1TrimmingCutoff") ,
#                        tags$ul(
#                            textInput("M1TrimmingCutoffText",
#                                      label = p("Change Value"),
#                                      value = toString(
#                                          trimmedParam[["M1TrimmingCutoff"]]),
#                                      width = '70%')
#                        ),
#                 ),
#             )
#         } else if (SangerCSetParam[[consensusReadIndex]]$
#                    SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M2") {
#             if (is.null(SangerCSetParam[[consensusReadIndex]]$
#                         SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore)) {
#                 SangerCSetParam[[consensusReadIndex]]$
#                     SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore <<-  20
#             }
#             if (is.null(SangerCSetParam[[consensusReadIndex]]$
#                         SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize )) {
#                 SangerCSetParam[[consensusReadIndex]]$
#                     SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize <<-  5
#             }
#             fluidRow(
#                 column(6,
#                        uiOutput("M2CutoffQualityScore") ,
#                        tags$ul(
#                            textInput("M2CutoffQualityScoreText",
#                                      label = p("Change Value"),
#                                      value = toString(
#                                          trimmedParam[["M2CutoffQualityScore"]]),
#                                      width = '70%')
#                        ),
#                 ),
#                 column(6,
#                        uiOutput("M2SlidingWindowSize") ,
#                        tags$ul(
#                            textInput("M2SlidingWindowSizeText",
#                                      label = p("Change Value"),
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
#     consensusReadIndex <- strtoi(sidebar_menu[[1]])
#     singleReadIndex <- strtoi(sidebar_menu[[5]])
#     if (!is.na(consensusReadIndex) &&
#         !is.na(singleReadIndex)) {
#         if (SangerCSetParam[[consensusReadIndex]]$
#             SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M1") {
#             tagList(icon("check-circle"),
#                     "Your trimming method selection :
#                         'Logarithmic Scale Trimming'")
#         } else if (SangerCSetParam[[consensusReadIndex]]$
#                    SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M2") {
#             tagList(icon("check-circle"),
#                     "Your trimming method selection :
#                         'Logarithmic Scale Sliding Window Trimming'")
#         }
#     }
# })
