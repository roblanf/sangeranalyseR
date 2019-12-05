### ============================================================================
### R shiny alignedConsensusSet server function
### ============================================================================
alignedConsensusSetServer <- function(input, output, session) {
    # Suppress Warning
    options(warn = -1)
    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusSet <- getShinyOption("SangerAlignedConsensusSet")
    shinyDirectory <- getShinyOption("shinyDirectory")
    SangerConsensusSet <- SangerConsensusSet[[1]]

    SangerConsensusSetNum <- length(SangerConsensusSet@consensusReadsList)


    SangerCSetParam <- lapply(1:SangerConsensusSetNum, function(i) {
        ### --------------------------------------------------------------------
        ### ConsensusRead-related parameters initialization.
        ### --------------------------------------------------------------------
        # readFeature
        SCName <- paste0(i, " Consensus Read")

        # Forward & reverse reads list
        SangerSingleReadFReadsList <-
            SangerConsensusSet@consensusReadsList[[i]]@forwardReadsList
        SangerSingleReadRReadsList <-
            SangerConsensusSet@consensusReadsList[[i]]@reverseReadsList
        forwardReadNum <- length(SangerSingleReadFReadsList)
        reverseReadNum <- length(SangerSingleReadRReadsList)
        SangerSingleReadNum <- forwardReadNum + reverseReadNum

        # Forward + reverse reads list
        SangerConsensusFRReadsList <- c(SangerSingleReadFReadsList,
                                        SangerSingleReadRReadsList)

        # readFileName (basename)
        forwardReadBFN <- sapply(1:forwardReadNum, function(j)
            basename(SangerSingleReadFReadsList[[j]]@readFileName))
        reverseReadBFN <- sapply(1:reverseReadNum, function(j)
            basename(SangerSingleReadRReadsList[[j]]@readFileName))
        SangerSingleReadBFN <- c(forwardReadBFN, reverseReadBFN)

        # readFileName (absolute)
        forwardReadAFN <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@readFileName)
        reverseReadAFN <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@readFileName)
        SangerSingleReadAFN <- c(forwardReadAFN, reverseReadAFN)

        # readFeature
        forwardReadFeature <- sapply(1:forwardReadNum, function(j)
            paste0(j, " ",
                   SangerSingleReadFReadsList[[j]]@readFeature))
        reverseReadFeature <- sapply(1:reverseReadNum, function(j)
            paste0(j+forwardReadNum, " ",
                   SangerSingleReadRReadsList[[j]]@readFeature))
        SangerSingleReadFeature <- c(forwardReadFeature, reverseReadFeature)

        # QualityReport
        forwardReadQualReport <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@QualityReport)
        reverseReadQualReport <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@QualityReport)
        SangerSingleReadQualReport <- c(forwardReadQualReport,
                                        reverseReadQualReport)
        # Quality Score Dataframe
        forwardQualityScoreDF <- lapply(1:forwardReadNum, function(j) {
            PhredScore <- forwardReadQualReport[[j]]@qualityPhredScores
            PhredScoreDF <- data.frame(
                t(data.frame(PhredScore)), stringsAsFactors = FALSE)
            colnames(PhredScoreDF) <- substr(colnames(PhredScoreDF), 2, 100)
            rownames(PhredScoreDF) <- NULL
            return(PhredScoreDF)
            }
        )
        reverseQualityScoreDF <- lapply(1:reverseReadNum, function(j) {
            PhredScore <- reverseReadQualReport[[j]]@qualityPhredScores
            PhredScoreDF <- data.frame(
                t(data.frame(PhredScore)), stringsAsFactors = FALSE)
            colnames(PhredScoreDF) <- substr(colnames(PhredScoreDF), 2, 100)
            rownames(PhredScoreDF) <- NULL
            return(PhredScoreDF)
            }
        )
        SangerSingleReadQSDF <- c(forwardQualityScoreDF, reverseQualityScoreDF)

        # ChromatogramParam
        forwardReadChromatogramParam <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@ChromatogramParam)
        reverseReadChromatogramParam <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@ChromatogramParam)
        SangerSingleReadChromatogramParam <- c(forwardReadChromatogramParam,
                                               reverseReadChromatogramParam)

        # primarySeqDF
        forwardReadPrimSeqDF <- lapply(1:forwardReadNum, function(j) {
            basecalls1 <- unlist(strsplit(
                toString(SangerSingleReadFReadsList[[j]]@primarySeq), ""))
            aveposition <- rowMeans(
                SangerSingleReadFReadsList[[j]]@peakPosMatrix, na.rm=TRUE)
            basecalls1 <- basecalls1[1:length(aveposition)]
            basecalls1DF <- data.frame(
                t(data.frame(basecalls1)), stringsAsFactors = FALSE)
            colnames(basecalls1DF) <- substr(colnames(basecalls1DF), 2, 100)
            rownames(basecalls1DF) <- NULL
            return(basecalls1DF)
            }
        )
        reverseReadPrimSeqDF <- lapply(1:reverseReadNum, function(j) {
            basecalls1 <- unlist(strsplit(
                toString(SangerSingleReadRReadsList[[j]]@primarySeq), ""))
            aveposition <- rowMeans(
                SangerSingleReadRReadsList[[j]]@peakPosMatrix, na.rm=TRUE)
            basecalls1 <- basecalls1[1:length(aveposition)]
            basecalls1DF <- data.frame(
                t(data.frame(basecalls1)), stringsAsFactors = FALSE)
            colnames(basecalls1DF) <- substr(colnames(basecalls1DF), 2, 100)
            rownames(basecalls1DF) <- NULL
            return(basecalls1DF)
            }
        )
        SangerSingleReadPrimSeqDF <- c(forwardReadPrimSeqDF,
                                       reverseReadPrimSeqDF)

        # secondarySeqDF
        forwardReadSecoSeqDF <- lapply(1:forwardReadNum, function(j) {
            basecalls2 <- unlist(strsplit(
                toString(SangerSingleReadFReadsList[[j]]@secondarySeq), ""))
            aveposition <- rowMeans(
                SangerSingleReadFReadsList[[j]]@peakPosMatrix, na.rm=TRUE)
            basecalls2 <- basecalls2[1:length(aveposition)]
            basecalls2DF <- data.frame(
                t(data.frame(basecalls2)), stringsAsFactors = FALSE)
            colnames(basecalls2DF) <- substr(colnames(basecalls2DF), 2, 100)
            rownames(basecalls2DF) <- NULL
            return(basecalls2DF)
            }
        )
        reverseReadSecoSeqDF <- lapply(1:reverseReadNum, function(j) {
            basecalls2 <- unlist(strsplit(toString(
                SangerSingleReadRReadsList[[j]]@secondarySeq), ""))
            aveposition <- rowMeans(
                SangerSingleReadRReadsList[[j]]@peakPosMatrix, na.rm=TRUE)
            basecalls2 <- basecalls2[1:length(aveposition)]
            basecalls2DF <- data.frame(
                t(data.frame(basecalls2)), stringsAsFactors = FALSE)
            colnames(basecalls2DF) <- substr(colnames(basecalls2DF), 2, 100)
            rownames(basecalls2DF) <- NULL
            return(basecalls2DF)
            }
        )
        SangerSingleReadSecoSeqDF <- c(forwardReadSecoSeqDF,
                                       reverseReadSecoSeqDF)

        # primaryAASeqDF
        forwardReadPrimAASeqDF <- lapply(1:forwardReadNum, function(j) {
            AAString <- SangerSingleReadFReadsList[[j]]@primaryAASeq
            AAStringDF <- data.frame(
                t(data.frame(AAString)), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
            }
        )
        reverseReadPrimAASeqDF <- lapply(1:reverseReadNum, function(j) {
            AAString <- SangerSingleReadRReadsList[[j]]@primaryAASeq
            AAStringDF <- data.frame(
                t(data.frame(AAString)), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
            }
        )
        SangerSingleReadPrimAASeqDF <- c(forwardReadPrimAASeqDF,
                                         reverseReadPrimAASeqDF)
        return(list(SCName = SCName,
                    SangerSingleReadFReadsList = SangerSingleReadFReadsList,
                    SangerSingleReadRReadsList = SangerSingleReadRReadsList,
                    forwardReadNum = forwardReadNum,
                    reverseReadNum = reverseReadNum,
                    SangerSingleReadNum = SangerSingleReadNum,
                    SangerConsensusFRReadsList = SangerConsensusFRReadsList,
                    SangerSingleReadBFN = SangerSingleReadBFN,
                    SangerSingleReadAFN = SangerSingleReadAFN,
                    SangerSingleReadFeature = SangerSingleReadFeature,
                    SangerSingleReadQualReport = SangerSingleReadQualReport,
                    SangerSingleReadChromatogramParam = SangerSingleReadChromatogramParam,
                    forwardReadPrimSeqDF = forwardReadPrimSeqDF,
                    reverseReadPrimSeqDF = reverseReadPrimSeqDF,
                    SangerSingleReadPrimSeqDF = SangerSingleReadPrimSeqDF,
                    forwardReadSecoSeqDF = forwardReadSecoSeqDF,
                    reverseReadSecoSeqDF = reverseReadSecoSeqDF,
                    SangerSingleReadPrimAASeqDF = SangerSingleReadPrimAASeqDF,
                    SangerSingleReadQSDF = SangerSingleReadQSDF,
                    SangerSingleReadSecoSeqDF = SangerSingleReadSecoSeqDF))
    })


    consensusParamSet <- reactiveValues(consensusReadSCSet  = SangerConsensusSet@consensusReadSCSet,
                                        alignmentSCSet      = SangerConsensusSet@alignmentSCSet,
                                        alignmentTreeSCSet  = SangerConsensusSet@alignmentTreeSCSet)

    consensusParam <- reactiveValues(consensusRead      = NULL,
                                     consensusReadName = NULL,
                                     differencesDF      = NULL,
                                     alignment          = NULL,
                                     distanceMatrix     = NULL,
                                     dendrogram         = NULL,
                                     indelsDF           = NULL,
                                     stopCodonsDF       = NULL,
                                     secondaryPeakDF    = NULL)

    trimmedRV <- reactiveValues(rawSeqLength = 0,
                                rawMeanQualityScore = 0,
                                rawMinQualityScore = 0,
                                trimmedStartPos = 0,
                                trimmedFinishPos = 0,
                                trimmedSeqLength = 0,
                                trimmedMeanQualityScore = 0,
                                trimmedMinQualityScore = 0,
                                remainingRatio = 0)

    trimmedParam <- reactiveValues(M1TrimmingCutoff     = 0,
                                   M2CutoffQualityScore = 0,
                                   M2SlidingWindowSize  = 0)

    ChromatogramParam <- reactiveValues(baseNumPerRow     = 0,
                                        heightPerRow      = 0,
                                        signalRatioCutoff = 0,
                                        showTrimmed       = TRUE)

    ############################################################################
    ### output$ID
    ############################################################################
    dynamicMenuSideBarSCSet(input, output, session, SangerCSetParam)

    output$aligned_consensusRead_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
            message(">>>>>>>> Inside 'Sanger Aligned Consensus Set Overview'")
            SCTrimmingMethod <-
                SangerCSetParam[[1]]$SangerSingleReadFReadsList[[1]]@
                QualityReport@TrimmingMethod
            if (SCTrimmingMethod == "M1") {
                SCTrimmingMethodName = "Method 1: 'Logarithmic Scale Trimming'"
            } else if (SCTrimmingMethod == "M2") {
                SCTrimmingMethodName =
                    "Method 2: 'Logarithmic Scale Sliding Window Trimming'"
            }
            # consensusParamSet[["consensusReadSCSet"]] <<-
            #     SangerConsensusSet@consensusReadSCSet
            # consensusParamSet[["alignmentSCSet"]] <<-
            #     SangerConsensusSet@alignmentSCSet
            # consensusParamSet[["alignmentTreeSCSet"]] <<-
            #     SangerConsensusSet@alignmentTreeSCSet
            fluidRow(
                useShinyjs(),
                box(title = tags$p("Input Parameters: ",
                                   style = "font-size: 26px;
                                       font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 0.2px hidden #A9A9A9;")),
                    fluidRow(
                        column(width = 12,
                               actionBttn("recalculateButtonSCSet",
                                          "Re-calculate
                                          consensusread (read set)",
                                          icon = icon("calculator"),
                                          style = "simple", color = "danger",
                                          block = TRUE, size = "lg")
                        ),
                        column(12,
                               tags$hr(
                                   style = ("border-top: 2px hidden #A9A9A9;"))
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Output Directory: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(shinyDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Raw ABI Parent Directory:"),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerConsensusSet@parentDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Trimming Method: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SCTrimmingMethodName),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Forward Suffix RegExp: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerConsensusSet@suffixForwardRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Reverse Suffix RegExp: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerConsensusSet@suffixReverseRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Consensus Read Number: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerConsensusSetNum),
                               )
                        ),
                        box(title = tags$p("Genetic Code Data Frame",
                                           style = "font-size: 24px;
                                       font-weight: bold;"),
                            collapsible = TRUE,
                            status = "success", width = 12,
                            column(width = 2,
                                   tags$p("Tri-nucleotide:",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                                   tags$p("Amino Acid : ",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                                   tags$p("('*' : stop codon) ",
                                          style = "font-size: 12px;
                                       font-weight: italic;"),
                            ),
                            column(width = 10,
                                   excelOutput("geneticCodeDF",
                                               width = "100%", height = "50"),
                                   style = paste("height:100%; ",
                                                 "overflow-y: hidden;",
                                                 "overflow-x: scroll;")
                            ),
                        ),
                        uiOutput("SCrefAminoAcidSeq") ,
                    ),
                ),


                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Consensus Readset Results: "),
                                   style = "font-size: 26px;
                                       font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Consensus Reads Alignment",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               htmlOutput("consensusSetAlignmentHTML"),
                        ),
                    ),
                    box(title = tags$p("Consensus Reads Tree",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               plotOutput("SCSetConsensusReadTreePlot"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                ),
            )
        } else {
            if (!is.na(as.numeric(sidebar_menu[[1]]))) {
                if (sidebar_menu[[2]] == "Sanger" &&
                    sidebar_menu[[3]] == "Consensus" &&
                    sidebar_menu[[4]] == "Read" &&
                    sidebar_menu[[5]] == "Overview") {
                    message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                    consensusReadIndex <- strtoi(sidebar_menu[[1]])
                    forwardReadNum <-
                        SangerCSetParam[[consensusReadIndex]]$forwardReadNum
                    reverseReadNum <-
                        SangerCSetParam[[consensusReadIndex]]$reverseReadNum
                    SCTrimmingMethod <-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadFReadsList[[1]]@
                        QualityReport@TrimmingMethod

                    if (SCTrimmingMethod == "M1") {
                        SCTrimmingMethodName =
                            "Method 1: 'Logarithmic Scale Trimming'"
                    } else if (SCTrimmingMethod == "M2") {
                        SCTrimmingMethodName =
                            "Method 2: 'Logarithmic Scale Sliding Window Trimming'"
                    }
                    consensusParam[["consensusReadName"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@consensusReadName
                    consensusParam[["consensusRead"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@consensusRead
                    consensusParam[["differencesDF"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@differencesDF
                    consensusParam[["alignment"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@alignment
                    consensusParam[["distanceMatrix"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@distanceMatrix
                    consensusParam[["dendrogram"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@dendrogram
                    consensusParam[["indelsDF"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@indelsDF
                    consensusParam[["stopCodonsDF"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@stopCodonsDF
                    consensusParam[["secondaryPeakDF"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@secondaryPeakDF

                    fluidRow(
                        useShinyjs(),
                        box(title = tags$p(tagList(icon("dot-circle"),
                                                   "Basic Information: "),
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(
                                style = ("border-top: 0.2px hidden #A9A9A9;")),
                            fluidRow(
                                column(width = 12,
                                       actionBttn("recalculateButton",
                                                  "Re-calculate consensus read",
                                                  icon = icon("calculator"),
                                                  style = "simple",
                                                  color = "danger",
                                                  block = TRUE, size = "lg")
                                ),
                                column(12,
                                       tags$hr(
                                           style = ("border-top: 2px
                                                    hidden #A9A9A9;")),
                                ),
                                column(12,
                                       column(3,
                                              tags$p(tagList(icon("caret-right"),
                                                            "Output Directory:"),
                                                     style = "font-size: 20px;
                                       font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(shinyDirectory),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              tags$p(
                                                  tagList(icon("caret-right"),
                                                          "Raw ABI Parent Directory: "),
                                                     style = "font-size: 20px;
                                       font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SangerConsensusSet@
                                                     parentDirectory),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              tags$p(
                                                  tagList(icon("caret-right"),
                                                          "Consenesus Read Name:"),
                                                     style = "font-size: 20px;
                                       font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(consensusParam[["consensusReadName"]]),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              tags$p(
                                                  tagList(icon("caret-right"),
                                                          "Trimming Method: "),
                                                     style = "font-size: 20px;
                                       font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SCTrimmingMethodName),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              tags$p(
                                                  tagList(icon("caret-right"),
                                                          "Forward Suffix RegExp:"),
                                                     style = "font-size: 20px;
                                       font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SangerConsensusSet@
                                                     suffixForwardRegExp),
                                       )

                                ),
                                column(12,
                                       column(3,
                                              tags$p(
                                                  tagList(icon("caret-right"),
                                                          "Forward Read Number:"),
                                                     style = "font-size: 20px;
                                       font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(forwardReadNum),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              tags$p(
                                                  tagList(icon("caret-right"),
                                                          "Reverse Suffix RegExp:"),
                                                     style = "font-size: 20px;
                                       font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SangerConsensusSet@
                                                     suffixReverseRegExp),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              tags$p(
                                                  tagList(icon("caret-right"),
                                                          "Reverse Read Number:"),
                                                     style = "font-size: 20px;
                                              font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(reverseReadNum),
                                       )
                                ),
                            ),
                            ################################################
                            #### Add this after having reference sample ####
                            ################################################
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                            box(title = tags$p("Consensus Read Parameters",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(3,
                                       uiOutput("SCMinReadsNum") ,
                                ),
                                column(3,
                                       uiOutput("SCMinReadLength")  ,
                                ),
                                column(3,
                                       uiOutput("SCMinFractionCall") ,
                                ),
                                column(3,
                                       uiOutput("SCMaxFractionLost") ,
                                ),
                                column(3,
                                       uiOutput("SCAcceptStopCodons") ,
                                ),
                                column(3,
                                       uiOutput("SCReadingFrame") ,
                                ),
                            ),
                            box(title = tags$p("Genetic Code Data Frame",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 2,
                                       tags$p("Tri-nucleotide:",
                                              style = "font-size: 15px;
                                       font-weight: bold;"),
                                       tags$p("Amino Acid : ",
                                              style = "font-size: 15px;
                                       font-weight: bold;"),
                                       tags$p("('*' : stop codon) ",
                                              style = "font-size: 12px;
                                       font-weight: italic;"),
                                ),
                                column(width = 10,
                                       excelOutput("geneticCodeDF",
                                                   width = "100%",
                                                   height = "50"),
                                       style = paste("height:100%; ",
                                                     "overflow-y: hidden;",
                                                     "overflow-x: scroll;")
                                ),
                            ),
                            uiOutput("SCrefAminoAcidSeq") ,
                        ),

                        box(title = tags$p(tagList(icon("dot-circle"),
                                                   "Consensus Read Results: "),
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style =("border-top: 4px hidden #A9A9A9;")),
                            box(title = tags$p("Alignment",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       htmlOutput("consensusAlignmentHTML"),
                                ),
                            ),
                            box(title = tags$p("Differences Data frame",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       uiOutput("SCDifferencesDFUI"),
                                       style =
                                           paste("height:100%; overflow-y:",
                                                 "scroll;overflow-x: scroll;")
                                )
                            ),
                            box(title = tags$p("Dendrogram",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       plotOutput("dendrogramPlot"),
                                       style =
                                           paste("height:100%; overflow-y:",
                                                 "scroll;overflow-x: scroll;")
                                ),
                                column(width = 12,
                                       tags$hr(
                                           style =("border-top: 4px hidden #A9A9A9;")),
                                ),
                                column(width = 12,
                                       dataTableOutput("dendrogramDF"),
                                       style =
                                           paste("height:100%; overflow-y:",
                                                 "scroll;overflow-x: scroll;")
                                )
                            ),
                            box(title = tags$p("Samples Distance",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       # plot()
                                       uiOutput("SCDistanceMatrixPlotUI"),
                                       style =
                                           paste("height:100%; overflow-y:",
                                                 "scroll;overflow-x: scroll;")
                                ),
                                column(width = 12,
                                       tags$hr(
                                           style =("border-top: 4px hidden #A9A9A9;")),
                                ),
                                column(width = 12,
                                       uiOutput("SCDistanceMatrixUI"),
                                       style =
                                           paste("height:100%; overflow-y:",
                                                 "scroll;overflow-x: scroll;")
                                )
                            ),
                            box(title = tags$p("Indels Data frame",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       uiOutput("SCIndelsDFUI"),
                                       style =
                                           paste("height:100%; overflow-y:",
                                                 "scroll;overflow-x: scroll;")
                                )
                            ),
                            box(title = tags$p("Stop Codons Data frame",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       uiOutput("SCStopCodonsDFUI"),
                                       style =
                                           paste("height:100%; overflow-y:",
                                                 "scroll;overflow-x: scroll;")
                                )
                            )
                        ),
                    )
                } else if (sidebar_menu[[2]] == "Consensus" &&
                           sidebar_menu[[3]] == "Read" &&
                           sidebar_menu[[4]] == "-" &&
                           !is.na(as.numeric(sidebar_menu[[5]])) &&
                           (sidebar_menu[[6]] == "Forward" ||
                            sidebar_menu[[6]] == "Reverse") &&
                           sidebar_menu[[7]] == "Read") {
                    message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                    consensusReadIndex <- strtoi(sidebar_menu[[1]])
                    singleReadIndex <- strtoi(sidebar_menu[[5]])
                    SangerSingleReadBFN <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadBFN
                    SangerSingleReadAFN <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadAFN

                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        M1TrimmingCutoff
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        M2CutoffQualityScore
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        M2SlidingWindowSize

                    ChromatogramParam[["baseNumPerRow"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadChromatogramParam[[singleReadIndex]]@
                        baseNumPerRow
                    ChromatogramParam[["heightPerRow"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadChromatogramParam[[singleReadIndex]]@
                        heightPerRow
                    ChromatogramParam[["signalRatioCutoff"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadChromatogramParam[[singleReadIndex]]@
                        signalRatioCutoff
                    ChromatogramParam[["showTrimmed"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadChromatogramParam[[singleReadIndex]]@
                        showTrimmed

                    trimmedRV[["rawSeqLength"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
                    trimmedRV[["rawMeanQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        rawMeanQualityScore
                    trimmedRV[["rawMinQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@rawMinQualityScore
                    trimmedRV[["trimmedStartPos"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@trimmedStartPos
                    trimmedRV[["trimmedFinishPos"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@trimmedFinishPos
                    trimmedRV[["trimmedSeqLength"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
                    trimmedRV[["trimmedMeanQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        trimmedMeanQualityScore
                    trimmedRV[["trimmedMinQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        trimmedMinQualityScore
                    trimmedRV[["remainingRatio"]] <<-
                        round(SangerCSetParam[[consensusReadIndex]]$
                                  SangerSingleReadQualReport[[singleReadIndex]]@
                                  remainingRatio * 100, 2)
                    fluidRow(
                        useShinyjs(),
                        box(title = tags$p(tagList(icon("dot-circle"),
                                                   "Raw File: "),
                                           style = "font-size: 26px;
                                         font-weight: bold;"),
                            solidHeader = TRUE,
                            status = "success", width = 12,
                            h1(paste0(
                                SangerSingleReadBFN[[strtoi(singleReadIndex)]])),
                            tags$h5(paste("( full path:",
                                          SangerSingleReadAFN[[
                                              strtoi(singleReadIndex)]],
                                          ")"), style = "font-style:italic")),
                        box(title = tags$p(
                            tagList(icon("dot-circle"),
                                    "Primary & Secondary Sequences:"),
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                            column(width = 2,
                                   tags$p("Primary Sequence",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p("Primary AA Sequence",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p("Secondary Sequence",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p("Quality Phred Score",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                            ),
                            column(width = 10,
                                   excelOutput("primarySeqDF",
                                               width = "100%", height = "50"),
                                   excelOutput("PrimAASeqDF",
                                               width = "100%", height = "50"),
                                   excelOutput("secondSeqDF",
                                               width = "100%", height = "50"),
                                   excelOutput("qualityScoreDF",
                                               width = "100%", height = "50"),

                                   style = paste("overflow-y: hidden;",
                                                 "overflow-x: scroll;")
                            )
                        ),
                        box(title = tags$p(tagList(icon("dot-circle"),
                                                   "Quality Report: "),
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(
                                style = ("border-top: 4px hidden #A9A9A9;")),
                            box(title =
                                    tags$p(tagList(icon("arrow-circle-right"),
                                                   "Trimming Parameters Input"),
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,

                                fluidRow(
                                    column(width = 12,
                                           uiOutput("TrimmingMethodSelectionOutput"),
                                    ),
                                ),
                                column(width = 12,
                                       uiOutput("TrimmingMethodUI"),
                                ),
                            ),
                            box(title = tags$p(tagList(icon("arrow-circle-left"),
                                                       "Trimmed Result Output"),
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                fluidRow(
                                    box(title = tags$p("Before Trimming",
                                                       style = "font-size: 21px;
                                       font-weight: bold;"),
                                        collapsible = TRUE,
                                        status = "success", width = 12,
                                        column(width = 12,
                                               column(4,
                                                      uiOutput("rawSeqLength") ,
                                               ),
                                               column(4,
                                                      uiOutput("rawMeanQualityScore"),
                                               ),
                                               column(4,
                                                      uiOutput("rawMinQualityScore"),
                                               ),
                                        ),
                                    ),
                                ),
                                fluidRow(
                                    box(title = tags$p("After Trimming",
                                                       style = "font-size: 21px;
                                       font-weight: bold;"),
                                        collapsible = TRUE,
                                        status = "success", width = 12,
                                        column(width = 12,
                                               column(4,
                                                      uiOutput("trimmedSeqLength"),
                                               ),
                                               column(4,
                                                      uiOutput("trimmedMeanQualityScore"),
                                               ),
                                               column(4,
                                                      uiOutput("trimmedMinQualityScore"),
                                               ),
                                        ),

                                        column(width = 12,
                                               column(4,
                                                      uiOutput("trimmedStartPos") ,
                                               ),
                                               column(4,
                                                      uiOutput("trimmedFinishPos") ,
                                               ),
                                               column(4,
                                                      uiOutput("remainingRatio") ,
                                               )
                                        ),
                                    ),
                                ),
                                tags$hr(
                                    style = ("border-top: 6px double #A9A9A9;")),
                                fluidRow(
                                    box(title = tags$p("Cumulative Ratio Plot",
                                                       style = "font-size: 21px;
                                       font-weight: bold;"),
                                        collapsible = TRUE,
                                        status = "success", width = 12,
                                        plotlyOutput("qualityTrimmingRatioPlot") %>%
                                            withSpinner()),
                                    box(title = tags$p("Cumulative Ratio Plot",
                                                       style = "font-size: 21px;
                                       font-weight: bold;"),
                                        collapsible = TRUE,
                                        status = "success", width = 12,
                                        plotlyOutput("qualityQualityBasePlot") %>%
                                            withSpinner()),
                                ),
                            ),
                        ),
                        box(title = tags$p(tagList(icon("dot-circle"),
                                                   "Chromatogram: "),
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),

                            box(title = tags$p(tagList(icon("arrow-circle-right"),
                                                       "Chromatogram Input"),
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(3,
                                       sliderInput("ChromatogramBasePerRow",
                                                   label = h4("Base Number Per Row"),
                                                   min = 5,
                                                   max = 200,
                                                   value = ChromatogramParam[["baseNumPerRow"]]),
                                       sliderInput("ChromatogramHeightPerRow",
                                                   label = h4("Height Per Row"),
                                                   min = 50,
                                                   max = 600,
                                                   value =ChromatogramParam[["heightPerRow"]]),
                                ),
                                column(3,
                                       numericInput(
                                           "ChromatogramSignalRatioCutoff",
                                           h3("Signal Ratio Cutoff"),
                                           value = ChromatogramParam[["signalRatioCutoff"]]),
                                       checkboxInput(
                                           "ChromatogramCheckShowTrimmed",
                                           "Whether show trimmed region",
                                           value = ChromatogramParam[["showTrimmed"]])
                                ),
                                column(3,
                                       uiOutput("ChromatogramtrimmedStartPos"),
                                ),
                                column(3,
                                       uiOutput("ChromatogramtrimmedFinishPos"),
                                ),
                                actionBttn("saveChromatogramParam",
                                           "Apply Chromatogram Parameters",
                                           style = "simple", color = "success",
                                           block = TRUE, size = "lg")
                            ),
                            box(title = tags$p(tagList(icon("arrow-circle-left"),
                                                       "Chromatogram Output"),
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       uiOutput("chromatogramUIOutput"),
                                )
                            ),
                        )
                    )
                }
            }
        }
    })





    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        newS4Object <- file.path(shinyDirectory,
                                 "SangerAlignedConsensusSet.Rda")
        showNotification(paste("New S4 object is store as:", newS4Object),
                         type = "message", duration = 10)
        sapply(1:SangerConsensusSetNum, function(i) {
            forwardReadNum <-
                length(SangerConsensusSet@
                           consensusReadsList[[i]]@forwardReadsList)
            reverseReadNum <-
                length(SangerConsensusSet@
                           consensusReadsList[[i]]@reverseReadsList)

            sapply(1:forwardReadNum, function(j) {
                SangerConsensusSet@consensusReadsList[[i]]@
                    forwardReadsList[[j]]@QualityReport <<-
                    SangerCSetParam[[i]]$SangerSingleReadQualReport[[j]]
                SangerConsensusSet@consensusReadsList[[i]]@
                    forwardReadsList[[j]]@ChromatogramParam <<-
                    SangerCSetParam[[i]]$
                    SangerSingleReadChromatogramParam[[j]]
                message("save SangerConsensus quality S4 object Forward")
                })
            sapply(1:reverseReadNum, function(j) {
                SangerConsensusSet@consensusReadsList[[i]]@
                    reverseReadsList[[j]]@QualityReport <<-
                    SangerCSetParam[[i]]$
                    SangerSingleReadQualReport[[forwardReadNum + j]]
                SangerConsensusSet@consensusReadsList[[i]]@
                    reverseReadsList[[j]]@ChromatogramParam <<-
                    SangerCSetParam[[i]]$
                    SangerSingleReadChromatogramParam[[forwardReadNum+j]]
                message("save SangerConsensus quality S4 object Reverse")
            })
        })
        saveRDS(SangerConsensusSet, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        NEW_SANGER_ALIGNED_CONSENSUS_READ_SET <<- readRDS(file=newS4Object)
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        btn <- input$closeUI
        stopApp()
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus read re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButton, {
        message("######## Reactive button clicked !!!")
        message("######## Start recalculating consensus read (SCSet")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        CSResult<-
            calculateConsensusRead (
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@forwardReadsList,
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@reverseReadsList,
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@refAminoAcidSeq,
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@minFractionCall,
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@maxFractionLost,
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@geneticCode,
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@acceptStopCodons,
                SangerConsensusSet@
                    consensusReadsList[[consensusReadIndex]]@readingFrame)
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@consensusRead <<-
            CSResult$consensusGapfree
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@differencesDF <<-
            CSResult$diffsDf
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@alignment <<-
            CSResult$aln2
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@distanceMatrix <<-
            CSResult$dist
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@dendrogram <<-
            CSResult$dend
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@indelsDF <<-
            CSResult$indels
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@stopCodonsDF <<-
            CSResult$stopsDf
        SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@secondaryPeakDF <<-
            CSResult$spDf


        consensusParam[["consensusReadName"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@consensusReadName
        consensusParam[["consensusRead"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@consensusRead
        consensusParam[["differencesDF"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@differencesDF
        consensusParam[["alignment"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@alignment
        consensusParam[["distanceMatrix"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@distanceMatrix
        consensusParam[["dendrogram"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@dendrogram
        consensusParam[["indelsDF"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@indelsDF
        consensusParam[["stopCodonsDF"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@stopCodonsDF
        consensusParam[["secondaryPeakDF"]] <<-
            SangerConsensusSet@
            consensusReadsList[[consensusReadIndex]]@secondaryPeakDF
        message("######## Finish recalculation")
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus read re-calculating (SCSet) UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButtonSCSet, {
        message("######## Reactive button clicked !!!")
        message("######## Start recalculating consensus read (SC)")
        CSSetResult <-
            alignConsensusReads (SangerConsensusSet@consensusReadsList,
                                 SangerConsensusSet@geneticCode,
                                 SangerConsensusSet@refAminoAcidSeq,
                                 SangerConsensusSet@minFractionCallSCSet,
                                 SangerConsensusSet@maxFractionLostSCSet,
                                 1)

        SangerConsensusSet@consensusReadSCSet <<- CSSetResult$consensus
        SangerConsensusSet@alignmentSCSet <<- CSSetResult$aln
        SangerConsensusSet@alignmentTreeSCSet <<- CSSetResult$aln.tree

        consensusParamSet[["consensusReadSCSet"]] <<- SangerConsensusSet@consensusReadSCSet
        consensusParamSet[["alignmentSCSet"]] <<- SangerConsensusSet@alignmentSCSet
        consensusParamSet[["alignmentTreeSCSet"]] <<- SangerConsensusSet@alignmentTreeSCSet
        message("######## Finish recalculation consensus read")
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$saveChromatogramParam, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        ### ------------------------------------------------------------
        ### Update ChromatogramBasePerRow
        ### ------------------------------------------------------------
        SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            baseNumPerRow <<- input$ChromatogramBasePerRow
        SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            heightPerRow <<- input$ChromatogramHeightPerRow
        SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
        SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            showTrimmed <<- input$ChromatogramCheckShowTrimmed

        ### ------------------------------------------------------------
        ### Save SangerConsensus quality S4 object
        ### ------------------------------------------------------------
        forwardReadNum <-
            length(SangerConsensusSet@
                       consensusReadsList[[consensusReadIndex]]@
                       forwardReadsList)
        reverseReadNum <-
            length(SangerConsensusSet@
                       consensusReadsList[[consensusReadIndex]]@
                       reverseReadsList)
        SangerSingleReadNum <- forwardReadNum + reverseReadNum

        if (singleReadIndex <= forwardReadNum) {
            # This is forward list
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                forwardReadsList[[singleReadIndex]]@ChromatogramParam <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadChromatogramParam[[singleReadIndex]]
        } else {
            # This is reverse list
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                reverseReadsList[[singleReadIndex - forwardReadNum]]@ChromatogramParam <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadChromatogramParam[[singleReadIndex]]
        }

        ChromatogramParam[["baseNumPerRow"]] <<-
            SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            baseNumPerRow
        ChromatogramParam[["heightPerRow"]] <<-
            SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            heightPerRow
        ChromatogramParam[["signalRatioCutoff"]] <<-
            SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            signalRatioCutoff
        ChromatogramParam[["showTrimmed"]] <<-
            SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            showTrimmed
    })

    ############################################################################
    ### ConsensusReadSet (Function for Sanger Consensus Read Set Overview)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Alignment
    ### ------------------------------------------------------------------------
    output$consensusSetAlignmentHTML<-renderUI({
        consensusParamSet[["alignmentSCSet"]] <-
            SangerConsensusSet@alignmentSCSet
        browseSeqHTML <-
            file.path(shinyDirectory,
                      "Consensus_Readset_Alignment_BrowseSeqs.html")
        BrowseSeqs(consensusParamSet[["alignmentSCSet"]] ,
                   openURL=FALSE, htmlFile=browseSeqHTML)
        includeHTML(browseSeqHTML)
    })

    ### ------------------------------------------------------------------------
    ### Consensus Reads Tree
    ### ------------------------------------------------------------------------
    output$SCSetConsensusReadTreePlot <- renderPlot({
        consensusParamSet[["alignmentTreeSCSet"]] <<- SangerConsensusSet@alignmentTreeSCSet
        plot(consensusParamSet[["alignmentTreeSCSet"]])
    })

    ############################################################################
    ### ConsensusRead (Function for Sanger Consensus Read Overview)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### genetic code
    ### ------------------------------------------------------------------------
    output$geneticCodeDF <- renderExcel({
        SCGeneticCode <- SangerConsensusSet@geneticCode
        suppressMessages(
            excelTable(data = t(data.frame(SCGeneticCode)),
                       defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })

    ### ------------------------------------------------------------------------
    ### refAminoAcidSeq
    ### ------------------------------------------------------------------------
    output$SCrefAminoAcidSeq <- renderUI({
        if (SangerConsensusSet@refAminoAcidSeq == "") {
            box(title = tags$p("Reference Amino Acids Sequence",
                               style = "font-size: 24px;
                                    font-weight: bold;"),
                collapsible = TRUE,
                status = "success", width = 12,
                column(width = 1),
                column(width = 11,
                       h4("Reference Amino Acid Sequence is not provided."))
            )
        } else {
            box(title = tags$p("Reference Amino Acids Sequence",
                               style = "font-size: 24px;
                                    font-weight: bold;"),
                collapsible = TRUE,
                status = "success", width = 12,
                column(width = 2,
                       tags$br(),
                       tags$p("AA Sequence:",
                              style = "font-size: 15px;
                                       font-weight: bold;"),
                ),
                column(width = 10,
                       excelOutput("SCrefAminoAcidSeqDF",
                                   width = "100%", height = "50"),
                       style = paste("height:100%; ",
                                     "overflow-y: hidden;",
                                     "overflow-x: scroll;")
                ),
            )
        }
    })

    output$SCrefAminoAcidSeqDF <- renderExcel({
        refAminoAcidSeqVec <-
            strsplit(SangerConsensusSet@refAminoAcidSeq, "")[[1]]
        names(refAminoAcidSeqVec) <- c(1:length(refAminoAcidSeqVec))
        suppressMessages(
            excelTable(data =
                           t(data.frame(refAminoAcidSeqVec)),
                       defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })

    ### ------------------------------------------------------------------------
    ### Alignment
    ### ------------------------------------------------------------------------
    output$consensusAlignmentHTML<-renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["consensusReadName"]] <<-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@consensusReadName
            consensusParam[["alignment"]] <-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@alignment

            browseSeqHTML <-
                file.path(shinyDirectory, "BrowseSeqs_html",
                          paste0(sidebar_menu[[1]], "_",
                                 consensusParam[["consensusReadName"]],
                                 "_Alignment_BrowseSeqs.html"))
            if (!dir.exists(file.path(shinyDirectory, "BrowseSeqs_html"))) {
                dir.create(file.path(shinyDirectory, "BrowseSeqs_html"))
            }
            BrowseSeqs(consensusParam[["alignment"]],
                       openURL=FALSE, htmlFile=browseSeqHTML)
            includeHTML(browseSeqHTML)
        }
    })

    ### ------------------------------------------------------------------------
    ### difference
    ### ------------------------------------------------------------------------
    output$SCDifferencesDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["differencesDF"]] <<-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@differencesDF
            if (all(dim(consensusParam[["differencesDF"]]) == c(0,0))) {
                h4("*** 'Differences' dataframe is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SCDifferencesDF")
            }
        }
    })
    output$SCDifferencesDF = renderDataTable({
        consensusParam[["differencesDF"]]
    })

    ### ------------------------------------------------------------------------
    ### dendrogram
    ### ------------------------------------------------------------------------
    output$dendrogramPlot <- renderPlot({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["dendrogram"]] <<-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@dendrogram
            plot(consensusParam[["dendrogram"]][[2]])
            ggdendrogram(consensusParam[["dendrogram"]][[2]], rotate = TRUE)
        }
    })
    output$dendrogramDF <- renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["dendrogram"]] <<-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@dendrogram
            consensusParam[["dendrogram"]][[1]]
        }
    })

    ### ------------------------------------------------------------------------
    ### distance
    ### ------------------------------------------------------------------------
    output$SCDistanceMatrixPlotUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["distanceMatrix"]] <<-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@distanceMatrix
            if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
                h4("*** 'Distance' dataframe is empty. (Cannot plot)***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                plotlyOutput("SCDistanceMatrixPlot")
            }
        }
    })

    output$SCDistanceMatrixPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            SangerSingleReadBFN <-
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadBFN
            suppressPlotlyMessage(
                plot_ly(x = SangerSingleReadBFN,
                        y = SangerSingleReadBFN,
                        z = consensusParam[["distanceMatrix"]],
                        colors = colorRamp(c("white", "#32a852")),
                        type = "heatmap")
            )
        }
    })
    output$SCDistanceMatrixUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["distanceMatrix"]] <<-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@distanceMatrix
            if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
                h4("*** 'Distance' dataframe is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SCDistanceMatrix")
            }
        }
    })
    output$SCDistanceMatrix = renderDataTable({
        consensusParam[["distanceMatrix"]]
    })

    ### ------------------------------------------------------------------------
    ### SCIndelsDF
    ### ------------------------------------------------------------------------
    output$SCIndelsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["indelsDF"]] <<-
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@indelsDF
            if (all(dim(consensusParam[["indelsDF"]] ) == c(0,0))) {
                h4("*** 'Indels' data frame is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SCIndelsDF")
            }
        }
    })
    output$SCIndelsDF <- renderDataTable({
        consensusParam[["indelsDF"]]
    })

    ### ------------------------------------------------------------------------
    ### SCStopCodons
    ### ------------------------------------------------------------------------
    output$SCStopCodonsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            consensusParam[["stopCodonsDF"]] <<-
                SangerConsensusSet@
                consensusReadsList[[consensusReadIndex]]@stopCodonsDF
            if (all(dim(consensusParam[["stopCodonsDF"]]) == c(0,0))) {
                h4("*** 'Stop Codons' dataframe is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SCStopCodonsDF")
            }
        }
    })
    output$SCStopCodonsDF <- renderDataTable({
        consensusParam[["stopCodonsDF"]]
    })

    ############################################################################
    ### SangerSingleRead (Function for singel read in consensusRead)
    ############################################################################
    output$primarySeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            suppressMessages(
                excelTable(data =SangerCSetParam[[consensusReadIndex]]$
                               SangerSingleReadPrimSeqDF[[singleReadIndex]],
                           defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                           columnResize = FALSE, allowInsertRow = FALSE,
                           allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                           allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
            )
        }
    })

    output$secondSeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            suppressMessages(
                excelTable(data =
                               SangerCSetParam[[consensusReadIndex]]$
                               SangerSingleReadSecoSeqDF[[singleReadIndex]],
                           defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                           columnResize = FALSE, allowInsertRow = FALSE,
                           allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                           allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
            )
        }
    })

    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            suppressMessages(
                excelTable(data =
                               SangerCSetParam[[consensusReadIndex]]$
                               SangerSingleReadQSDF[[singleReadIndex]],
                           defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                           columnResize = FALSE, allowInsertRow = FALSE,
                           allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                           allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
            )
        }
    })

    output$PrimAASeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            suppressMessages(
                excelTable(data =
                               SangerCSetParam[[consensusReadIndex]]$
                               SangerSingleReadPrimAASeqDF[[singleReadIndex]],
                           defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
                           columnResize = FALSE, allowInsertRow = FALSE,
                           allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                           allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
            )
        }
    })

    valueBoxSCMinReadsNumCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCMinReadLengthCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCMinFractionCallCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCMaxFractionLostCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCAcceptStopCodonsCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCReadingFrameCSSet (input, output, SangerConsensusSet, session)

    valueBoxM1TrimmingCutoff (input, output, session)
    valueBoxM2CutoffQualityScore (input, output, session)
    valueBoxM2SlidingWindowSize (input, output, session)

    valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)

    valueBoxRawSeqLength (input, output, session, trimmedRV)
    valueBoxRawMeanQualityScore (input, output, session, trimmedRV)
    valueBoxRawMinQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxTrimmedFinishPos (input, output, session, trimmedRV)
    valueBoxTrimmedSeqLength (input, output, session, trimmedRV)
    valueBoxTrimmedMeanQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedMinQualityScore (input, output, session, trimmedRV)
    valueBoxRemainingRatio (input, output, session, trimmedRV)


    observeEvent(input$M1TrimmingCutoffText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
            as.numeric(input$M1TrimmingCutoffText) > 0 &&
            as.numeric(input$M1TrimmingCutoffText) <= 1) {
            inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
        } else {
            inputM1TrimmingCutoffText <- 0.0001
        }
        if (SangerCSetParam[[consensusReadIndex]]$
            SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M1") {
            trimmingPos <-
                M1inside_calculate_trimming(
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityPhredScores,
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityBaseScore,
                    as.numeric(inputM1TrimmingCutoffText))
            rawSeqLength <- trimmingPos[["rawSeqLength"]]
            rawMeanQualityScore <- trimmingPos[["rawMeanQualityScore"]]
            rawMinQualityScore <- trimmingPos[["rawMinQualityScore"]]
            trimmedStartPos <- trimmingPos[["trimmedStartPos"]]
            trimmedFinishPos <- trimmingPos[["trimmedFinishPos"]]
            trimmedSeqLength <- trimmingPos[["trimmedSeqLength"]]
            trimmedMeanQualityScore <- trimmingPos[["trimmedMeanQualityScore"]]
            trimmedMinQualityScore <- trimmingPos[["trimmedMinQualityScore"]]
            remainingRatio <- trimmingPos[["remainingRatio"]]

            if (!is.null(rawSeqLength) && !is.null(rawMeanQualityScore) &&
                !is.null(rawMinQualityScore ) && !is.null(trimmedStartPos) &&
                !is.null(trimmedFinishPos) && !is.null(trimmedSeqLength) &&
                !is.null(trimmedMeanQualityScore) &&
                !is.null(trimmedMinQualityScore)) {

                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)

                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawSeqLength <<- rawSeqLength
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawMeanQualityScore <<- rawMeanQualityScore
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawMinQualityScore <<- rawMinQualityScore
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedStartPos <<- trimmedStartPos
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedFinishPos <<- trimmedFinishPos
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedSeqLength <<- trimmedSeqLength
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMeanQualityScore <<- trimmedMeanQualityScore
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMinQualityScore <<- trimmedMinQualityScore
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    remainingRatio <<- remainingRatio

                trimmedRV[["rawSeqLength"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerCSetParam[[consensusReadIndex]]$
                              SangerSingleReadQualReport[[singleReadIndex]]@
                              remainingRatio * 100, 2)
                ### ------------------------------------------------------------
                ### Save SangerConsensus quality S4 object
                ### ------------------------------------------------------------
                forwardReadNum <-
                    length(SangerConsensusSet@
                               consensusReadsList[[consensusReadIndex]]@forwardReadsList)
                reverseReadNum <-
                    length(SangerConsensusSet@
                               consensusReadsList[[consensusReadIndex]]@reverseReadsList)
                SangerSingleReadNum <- forwardReadNum + reverseReadNum
                if (singleReadIndex <= forwardReadNum) {
                    # This is forward list
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        forwardReadsList[[singleReadIndex]]@QualityReport <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]
                } else {
                    # This is reverse list
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        reverseReadsList[[singleReadIndex-forwardReadNum]]@QualityReport <<-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]
                }
            }
        }
    })

    observeEvent(input$M2CutoffQualityScoreText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
            strtoi(input$M2CutoffQualityScoreText) > 0 &&
            strtoi(input$M2CutoffQualityScoreText) <= 60 &&
            strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
            inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
        } else {
            inputM2CutoffQualityScoreText <- 20
        }

        trimmingPos <-
            M2inside_calculate_trimming(
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    qualityPhredScores,
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    qualityBaseScore,
                strtoi(inputM2CutoffQualityScoreText),
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    M2SlidingWindowSize)

        rawSeqLength <- trimmingPos[["rawSeqLength"]]
        rawMeanQualityScore <- trimmingPos[["rawMeanQualityScore"]]
        rawMinQualityScore <- trimmingPos[["rawMinQualityScore"]]
        trimmedStartPos <- trimmingPos[["trimmedStartPos"]]
        trimmedFinishPos <- trimmingPos[["trimmedFinishPos"]]
        trimmedSeqLength <- trimmingPos[["trimmedSeqLength"]]
        trimmedMeanQualityScore <- trimmingPos[["trimmedMeanQualityScore"]]
        trimmedMinQualityScore <- trimmingPos[["trimmedMinQualityScore"]]
        remainingRatio <- trimmingPos[["remainingRatio"]]

        if (!is.null(rawSeqLength) && !is.null(rawMeanQualityScore) &&
            !is.null(rawMinQualityScore ) && !is.null(trimmedStartPos) &&
            !is.null(trimmedFinishPos) && !is.null(trimmedSeqLength) &&
            !is.null(trimmedMeanQualityScore) &&
            !is.null(trimmedMinQualityScore)) {
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                M2CutoffQualityScore <<- strtoi(inputM2CutoffQualityScoreText)

            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength <<- rawSeqLength
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore <<- rawMeanQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore <<- rawMinQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos <<- trimmedStartPos
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos <<- trimmedFinishPos
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength <<- trimmedSeqLength
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore <<- trimmedMeanQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore <<- trimmedMinQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                remainingRatio <<- remainingRatio

            trimmedRV[["rawSeqLength"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
            trimmedRV[["rawMeanQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore
            trimmedRV[["rawMinQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@rawMinQualityScore
            trimmedRV[["trimmedStartPos"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@trimmedStartPos
            trimmedRV[["trimmedFinishPos"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@trimmedFinishPos
            trimmedRV[["trimmedSeqLength"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
            trimmedRV[["trimmedMeanQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore
            trimmedRV[["trimmedMinQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore
            trimmedRV[["remainingRatio"]] <<-
                round(SangerCSetParam[[consensusReadIndex]]$
                          SangerSingleReadQualReport[[singleReadIndex]]@
                          remainingRatio * 100, 2)
            ### ------------------------------------------------------------
            ### Save SangerConsensus quality S4 object
            ### ------------------------------------------------------------
            forwardReadNum <-
                length(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@forwardReadsList)
            reverseReadNum <-
                length(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@reverseReadsList)
            SangerSingleReadNum <- forwardReadNum + reverseReadNum
            if (singleReadIndex <= forwardReadNum) {
                # This is forward list
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]
            } else {
                # This is reverse list
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex - forwardReadNum]]@QualityReport <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]
            }
        }
    })

    observeEvent(input$M2SlidingWindowSizeText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
            strtoi(input$M2SlidingWindowSizeText) > 0 &&
            strtoi(input$M2SlidingWindowSizeText) <= 20 &&
            strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
            inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
        } else {
            inputM2SlidingWindowSizeText <- 5
        }
        trimmingPos <-
            M2inside_calculate_trimming(
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    qualityPhredScores,
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    qualityBaseScore,
                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    M2CutoffQualityScore,
                strtoi(inputM2SlidingWindowSizeText))
        rawSeqLength <- trimmingPos[["rawSeqLength"]]
        rawMeanQualityScore <- trimmingPos[["rawMeanQualityScore"]]
        rawMinQualityScore <- trimmingPos[["rawMinQualityScore"]]
        trimmedStartPos <- trimmingPos[["trimmedStartPos"]]
        trimmedFinishPos <- trimmingPos[["trimmedFinishPos"]]
        trimmedSeqLength <- trimmingPos[["trimmedSeqLength"]]
        trimmedMeanQualityScore <- trimmingPos[["trimmedMeanQualityScore"]]
        trimmedMinQualityScore <- trimmingPos[["trimmedMinQualityScore"]]
        remainingRatio <- trimmingPos[["remainingRatio"]]

        if (!is.null(rawSeqLength) && !is.null(rawMeanQualityScore) &&
            !is.null(rawMinQualityScore ) && !is.null(trimmedStartPos) &&
            !is.null(trimmedFinishPos) && !is.null(trimmedSeqLength) &&
            !is.null(trimmedMeanQualityScore) &&
            !is.null(trimmedMinQualityScore)) {

            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)

            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength <<- rawSeqLength
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore <<- rawMeanQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore <<- rawMinQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos <<- trimmedStartPos
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos <<- trimmedFinishPos
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength <<- trimmedSeqLength
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore <<- trimmedMeanQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore <<- trimmedMinQualityScore
            SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                remainingRatio <<- remainingRatio

            trimmedRV[["rawSeqLength"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
            trimmedRV[["rawMeanQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore
            trimmedRV[["rawMinQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@rawMinQualityScore
            trimmedRV[["trimmedStartPos"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@trimmedStartPos
            trimmedRV[["trimmedFinishPos"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@trimmedFinishPos
            trimmedRV[["trimmedSeqLength"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
            trimmedRV[["trimmedMeanQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore
            trimmedRV[["trimmedMinQualityScore"]] <<-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore
            trimmedRV[["remainingRatio"]] <<-
                round(SangerCSetParam[[consensusReadIndex]]$
                          SangerSingleReadQualReport[[singleReadIndex]]@
                          remainingRatio * 100, 2)
            ### ------------------------------------------------------------
            ### Save SangerConsensus quality S4 object
            ### ------------------------------------------------------------
            forwardReadNum <-
                length(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@forwardReadsList)
            reverseReadNum <-
                length(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@reverseReadsList)
            SangerSingleReadNum <- forwardReadNum + reverseReadNum
            if (singleReadIndex <= forwardReadNum) {
                # This is forward list
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]
            } else {
                # This is reverse list
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex-forwardReadNum]]@QualityReport <<-
                    SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]
            }
        }
    })

    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            readFeature <- SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadFeature[[singleReadIndex]]
            trimmedStartPos = trimmedRV[["trimmedStartPos"]]
            trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
            qualityPhredScores <-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
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
            x <- list(
                title = "Base Pair Index"
                # titlefont = f
            )
            y <- list(
                title = "Read Ratio"
                # titlefont = f
            )
            PerDataPlot <- melt(PerData, id.vars = c("Base"))
            suppressPlotlyMessage(
                plot_ly(data=PerDataPlot,
                        x=~Base,
                        y=~value,
                        mode="markers",
                        color = ~variable,
                        text = ~paste("BP Index : ",
                                      Base, '<sup>th</sup><br>Read Ratio :',
                                      round(value*100, digits = 2), '%')) %>%
                    layout(xaxis = x,
                           yaxis = y,
                           legend = list(orientation = 'h',
                                         xanchor = "center",
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
                        y = ((remainingPer[1]+remainingPer[length(remainingPer)])/2)
                        - 0.06,
                        showarrow=FALSE
                    ))
        }
    })

    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            readFeature <-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadFeature[[singleReadIndex]]
            trimmedStartPos = trimmedRV[["trimmedStartPos"]]
            trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
            qualityPhredScores <-
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
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
            suppressPlotlyMessage(
                plot_ly(data=qualityPlotDf,
                        x=~Index) %>%
                    add_markers(y=~Score,
                                text =
                                    ~paste("BP Index : ",
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
                                         round((trimmedFinishPos-trimmedStartPos+1)/
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
                                         xanchor = "center",
                                         x = 0.5, y = 1.1)) %>%
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
                    ))
        }
    })

    valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)

    # chromatogram
    output$chromatogramUIOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            sidebar_menu[[2]] == "Consensus" &&
            sidebar_menu[[3]] == "Read" &&
            sidebar_menu[[4]] == "-" &&
            !is.na(singleReadIndex) &&
            (sidebar_menu[[6]] == "Forward" ||
             sidebar_menu[[6]] == "Reverse") &&
            sidebar_menu[[7]] == "Read") {
            trimmedRV[["trimmedSeqLength"]]
            rawSeqLength =
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
            trimmedSeqLength =
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength

            chromatogramRowNumAns <- chromatogramRowNum (
                strtoi(ChromatogramParam[["baseNumPerRow"]]),
                rawSeqLength, trimmedSeqLength,
                ChromatogramParam[["showTrimmed"]]) *
                strtoi(ChromatogramParam[["heightPerRow"]])
            plotOutput("chromatogram", height = chromatogramRowNumAns)
        }
    })

    output$chromatogram <- renderPlot({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            sidebar_menu[[2]] == "Consensus" &&
            sidebar_menu[[3]] == "Read" &&
            sidebar_menu[[4]] == "-" &&
            !is.na(singleReadIndex) &&
            (sidebar_menu[[6]] == "Forward" || sidebar_menu[[6]] == "Reverse")&&
            sidebar_menu[[7]] == "Read") {
            rawSeqLength =
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
            hetcalls <-
                MakeBaseCalls(
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerConsensusFRReadsList[[singleReadIndex]],
                    signalRatioCutoff = as.numeric(
                        ChromatogramParam[["signalRatioCutoff"]]))
            chromatogram(hetcalls,
                         width = strtoi(ChromatogramParam[["baseNumPerRow"]]),
                         height = 2, trim5 = trimmedRV[["trimmedStartPos"]],
                         trim3 = rawSeqLength - trimmedRV[["trimmedFinishPos"]],
                         showtrim = (ChromatogramParam[["showTrimmed"]]),
                         showcalls = "both")
        }
    })

    output$TrimmingMethodUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if ( SangerCSetParam[[consensusReadIndex]]$
                 SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod== "M1") {
                if (is.null(SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff)) {
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff <<-  0.0001
                }
                fluidRow(
                    column(6,
                           uiOutput("M1TrimmingCutoff") ,
                           tags$ul(
                               textInput("M1TrimmingCutoffText",
                                         label = p("Change Value"),
                                         value = toString(
                                             trimmedParam[["M1TrimmingCutoff"]]),
                                         width = '70%')
                           ),
                    ),
                )
            } else if (SangerCSetParam[[consensusReadIndex]]$
                       SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M2") {
                if (is.null(SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore)) {
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore <<-  20
                }
                if (is.null(SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize )) {
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize <<-  5
                }
                fluidRow(
                    column(6,
                           uiOutput("M2CutoffQualityScore") ,
                           tags$ul(
                               textInput("M2CutoffQualityScoreText",
                                         label = p("Change Value"),
                                         value = toString(
                                             trimmedParam[["M2CutoffQualityScore"]]),
                                         width = '70%')
                           ),
                    ),
                    column(6,
                           uiOutput("M2SlidingWindowSize") ,
                           tags$ul(
                               textInput("M2SlidingWindowSizeText",
                                         label = p("Change Value"),
                                         value = toString(
                                             trimmedParam[["M2SlidingWindowSize"]]),
                                         width = '70%')
                           ),
                    ),
                )
            }
        }
    })

    output$TrimmingMethodSelectionOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if (SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Logarithmic Scale Trimming'")
            } else if (SangerCSetParam[[consensusReadIndex]]$
                       SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M2") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Logarithmic Scale Sliding Window Trimming'")
            }
        }
    })
}

