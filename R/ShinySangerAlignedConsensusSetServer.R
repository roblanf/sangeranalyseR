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
            paste0(j, " ",
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

        # primaryAASeqS1DF
        forwardReadPrimAASeqS1DF <- lapply(1:forwardReadNum, function(j) {
            AAString <- data.frame(SangerSingleReadFReadsList[[j]]@primaryAASeqS1)
            AAStringDF <- data.frame(
                t(AAString), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
            }
        )
        reverseReadPrimAASeqS1DF <- lapply(1:reverseReadNum, function(j) {
            AAString <- data.frame(SangerSingleReadRReadsList[[j]]@primaryAASeqS1)
            AAStringDF <- data.frame(
                t(AAString), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
            }
        )
        SangerSingleReadPrimAASeqS1DF <- c(forwardReadPrimAASeqS1DF,
                                           reverseReadPrimAASeqS1DF)

        # primaryAASeqS2DF
        forwardReadPrimAASeqS2DF <- lapply(1:forwardReadNum, function(j) {
            AAString <- data.frame(SangerSingleReadFReadsList[[j]]@primaryAASeqS2)
            AAString <- rbind(NA, AAString)
            AAStringDF <- data.frame(
                t(AAString), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
        }
        )
        reverseReadPrimAASeqS2DF <- lapply(1:reverseReadNum, function(j) {
            AAString <- data.frame(SangerSingleReadRReadsList[[j]]@primaryAASeqS2)
            AAString <- rbind(NA, AAString)
            AAStringDF <- data.frame(
                t(AAString), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
        }
        )
        SangerSingleReadPrimAASeqS2DF <- c(forwardReadPrimAASeqS2DF,
                                           reverseReadPrimAASeqS2DF)

        # primaryAASeqS3DF
        forwardReadPrimAASeqS3DF <- lapply(1:forwardReadNum, function(j) {
            AAString <- data.frame(SangerSingleReadFReadsList[[j]]@primaryAASeqS3)
            AAString <- rbind(NA, AAString)
            AAString <- rbind(NA, AAString)
            AAStringDF <- data.frame(
                t(AAString), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
        }
        )
        reverseReadPrimAASeqS3DF <- lapply(1:reverseReadNum, function(j) {
            AAString <- data.frame(SangerSingleReadRReadsList[[j]]@primaryAASeqS3)
            AAString <- rbind(NA, AAString)
            AAString <- rbind(NA, AAString)
            AAStringDF <- data.frame(
                t(AAString), stringsAsFactors = FALSE)
            colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
            rownames(AAStringDF) <- NULL
            return(AAStringDF)
        }
        )
        SangerSingleReadPrimAASeqS3DF <- c(forwardReadPrimAASeqS3DF,
                                           reverseReadPrimAASeqS3DF)

        return(list(SCName = SCName,
                    SangerSingleReadFReadsList = SangerSingleReadFReadsList,
                    SangerSingleReadRReadsList = SangerSingleReadRReadsList,
                    forwardReadNum = forwardReadNum,
                    reverseReadNum = reverseReadNum,
                    SangerSingleReadNum = SangerSingleReadNum,
                    SangerConsensusFRReadsList = SangerConsensusFRReadsList,
                    SangerSingleReadBFN = SangerSingleReadBFN,
                    SangerSingleReadAFN = SangerSingleReadAFN,
                    forwardReadFeature = forwardReadFeature,
                    reverseReadFeature = reverseReadFeature,
                    SangerSingleReadFeature = SangerSingleReadFeature,
                    SangerSingleReadQualReport = SangerSingleReadQualReport,
                    SangerSingleReadChromatogramParam = SangerSingleReadChromatogramParam,
                    forwardReadPrimSeqDF = forwardReadPrimSeqDF,
                    reverseReadPrimSeqDF = reverseReadPrimSeqDF,
                    SangerSingleReadPrimSeqDF = SangerSingleReadPrimSeqDF,
                    forwardReadSecoSeqDF = forwardReadSecoSeqDF,
                    reverseReadSecoSeqDF = reverseReadSecoSeqDF,
                    SangerSingleReadPrimAASeqS1DF = SangerSingleReadPrimAASeqS1DF,
                    SangerSingleReadPrimAASeqS2DF = SangerSingleReadPrimAASeqS2DF,
                    SangerSingleReadPrimAASeqS3DF = SangerSingleReadPrimAASeqS3DF,


                    SangerSingleReadQSDF = SangerSingleReadQSDF,
                    SangerSingleReadSecoSeqDF = SangerSingleReadSecoSeqDF))
    })

    SCTrimmingMethod <-
        SangerConsensusSet@consensusReadsList[[1]]@
        forwardReadsList[[1]]@QualityReport@TrimmingMethod
    if (SCTrimmingMethod == "M1") {
        SCTrimmingMethodName = "Method 1: 'Logarithmic Scale Trimming'"
    } else if (SCTrimmingMethod == "M2") {
        SCTrimmingMethodName =
            "Method 2: 'Logarithmic Scale Sliding Window Trimming'"
    }

    ### ------------------------------------------------------------------------
    ### ConsensusRead Set reactiveValue
    ### ------------------------------------------------------------------------
    consensusParamSet <-
        reactiveValues(
            consensusReadSCSet = SangerConsensusSet@consensusReadSCSet,
            alignmentSCSet     = SangerConsensusSet@alignmentSCSet,
            alignmentTreeSCSet = SangerConsensusSet@alignmentTreeSCSet)

    ### ------------------------------------------------------------------------
    ### ConsensusRead reactiveValue
    ### ------------------------------------------------------------------------
    consensusParam <-
        reactiveValues(
            consensusRead   = NULL,
            differencesDF   = NULL,
            alignment       = NULL,
            distanceMatrix  = NULL,
            dendrogram      = NULL,
            indelsDF        = NULL,
            stopCodonsDF    = NULL,
            secondaryPeakDF = NULL)

    ### ------------------------------------------------------------------------
    ### SingleRead reactiveValue
    ### ------------------------------------------------------------------------
    sequenceParam <- reactiveValues(primarySeq = "",
                                    secondarySeq = "",
                                    primaryAASeqS1 = "",
                                    primaryAASeqS2 = "",
                                    primaryAASeqS3 = "")

    trimmedRV <- reactiveValues(rawSeqLength            = 0,
                                rawMeanQualityScore     = 0,
                                rawMinQualityScore      = 0,
                                trimmedStartPos         = 0,
                                trimmedFinishPos        = 0,
                                trimmedSeqLength        = 0,
                                trimmedMeanQualityScore = 0,
                                trimmedMinQualityScore  = 0,
                                remainingRatio          = 0)

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
    output$aligned_consensusRead_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        message(input$sidebar_menu)
        if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
            h1(input$sidebar_menu)
            message(">>>>>>>> Inside 'Sanger Aligned Consensus Set Overview'")
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
                                      h4(length(SangerConsensusSet@
                                                    consensusReadsList)),
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
            if (!is.na(consensusReadIndex)) {
                if (sidebar_menu[[2]] == "Sanger" &&
                    sidebar_menu[[3]] == "Consensus" &&
                    sidebar_menu[[4]] == "Read" &&
                    sidebar_menu[[5]] == "Overview") {
                    message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                    consensusParam[["consensusReadName"]] <<-
                        SangerConsensusSet@
                        consensusReadsList[[consensusReadIndex]]@consensusReadName
                    consensusParam[["consensusRead"]] <<-
                        as.character(SangerConsensusSet@
                                         consensusReadsList[[consensusReadIndex]]@consensusRead)
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
                                              h4(length(SangerConsensusSet@
                                                     consensusReadsList[[consensusReadIndex]]@
                                                     forwardReadsList)),
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
                                              h4(length(SangerConsensusSet@
                                                     consensusReadsList[[consensusReadIndex]]@
                                                     reverseReadsList)),
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
                    )
                )

                } else if (sidebar_menu[[2]] == "CR" &&
                           sidebar_menu[[3]] == "-" &&
                           !is.na(singleReadIndex) &&
                           (directionParam == "Forward" ||
                            directionParam == "Reverse") &&
                           sidebar_menu[[6]] == "Read") {
                    if (directionParam == "Forward") {
                        sequenceParam[["primarySeq"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             forwardReadsList[[singleReadIndex]]@primarySeq)
                        sequenceParam[["secondarySeq"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             forwardReadsList[[singleReadIndex]]@secondarySeq)
                        sequenceParam[["primaryAASeqS1"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             forwardReadsList[[singleReadIndex]]@primaryAASeqS1)
                        sequenceParam[["primaryAASeqS2"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             forwardReadsList[[singleReadIndex]]@primaryAASeqS2)
                        sequenceParam[["primaryAASeqS3"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             forwardReadsList[[singleReadIndex]]@primaryAASeqS3)
                        trimmedParam[["M1TrimmingCutoff"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@
                            M1TrimmingCutoff
                        trimmedParam[["M2CutoffQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@
                            M2CutoffQualityScore
                        trimmedParam[["M2SlidingWindowSize"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@
                            M2SlidingWindowSize

                        trimmedRV[["rawSeqLength"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@rawSeqLength
                        trimmedRV[["rawMeanQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@
                            rawMeanQualityScore
                        trimmedRV[["rawMinQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@rawMinQualityScore
                        trimmedRV[["trimmedStartPos"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@trimmedStartPos
                        trimmedRV[["trimmedFinishPos"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@trimmedFinishPos
                        trimmedRV[["trimmedSeqLength"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@trimmedSeqLength
                        trimmedRV[["trimmedMeanQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@
                            trimmedMeanQualityScore
                        trimmedRV[["trimmedMinQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@QualityReport@
                            trimmedMinQualityScore
                        trimmedRV[["remainingRatio"]] <<-
                            round(SangerConsensusSet@
                                      consensusReadsList[[consensusReadIndex]]@
                                      forwardReadsList[[singleReadIndex]]@QualityReport@
                                      remainingRatio * 100, 2)

                        ChromatogramParam[["baseNumPerRow"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                            baseNumPerRow
                        ChromatogramParam[["heightPerRow"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                            heightPerRow
                        ChromatogramParam[["signalRatioCutoff"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                            signalRatioCutoff
                        ChromatogramParam[["showTrimmed"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                            showTrimmed

                        SangerSingleReadBFN <-
                            basename(SangerConsensusSet@
                                         consensusReadsList[[consensusReadIndex]]@
                                         forwardReadsList[[singleReadIndex]]@readFileName)
                        SangerSingleReadAFN <-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@readFileName
                    } else if (directionParam == "Reverse") {
                        sequenceParam[["primarySeq"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             reverseReadsList[[singleReadIndex]]@primarySeq)
                        sequenceParam[["secondarySeq"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             reverseReadsList[[singleReadIndex]]@secondarySeq)
                        sequenceParam[["primaryAASeqS1"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             reverseReadsList[[singleReadIndex]]@primaryAASeqS1)
                        sequenceParam[["primaryAASeqS2"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             reverseReadsList[[singleReadIndex]]@primaryAASeqS2)
                        sequenceParam[["primaryAASeqS3"]] <<-
                            as.character(SangerConsensusSet@
                                             consensusReadsList[[consensusReadIndex]]@
                                             reverseReadsList[[singleReadIndex]]@primaryAASeqS3)

                        trimmedParam[["M1TrimmingCutoff"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@
                            M1TrimmingCutoff
                        trimmedParam[["M2CutoffQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@
                            M2CutoffQualityScore
                        trimmedParam[["M2SlidingWindowSize"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@
                            M2SlidingWindowSize

                        trimmedRV[["rawSeqLength"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@rawSeqLength
                        trimmedRV[["rawMeanQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@
                            rawMeanQualityScore
                        trimmedRV[["rawMinQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@rawMinQualityScore
                        trimmedRV[["trimmedStartPos"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@trimmedStartPos
                        trimmedRV[["trimmedFinishPos"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@trimmedFinishPos
                        trimmedRV[["trimmedSeqLength"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@trimmedSeqLength
                        trimmedRV[["trimmedMeanQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@
                            trimmedMeanQualityScore
                        trimmedRV[["trimmedMinQualityScore"]] <<-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@QualityReport@
                            trimmedMinQualityScore
                        trimmedRV[["remainingRatio"]] <<-
                            round(SangerConsensusSet@
                                      consensusReadsList[[consensusReadIndex]]@
                                      reverseReadsList[[singleReadIndex]]@QualityReport@
                                      remainingRatio * 100, 2)

                        ChromatogramParam[["baseNumPerRow"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                            baseNumPerRow
                        ChromatogramParam[["heightPerRow"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                            heightPerRow
                        ChromatogramParam[["signalRatioCutoff"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                            signalRatioCutoff
                        ChromatogramParam[["showTrimmed"]] <<-
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                            showTrimmed

                        SangerSingleReadBFN <-
                            basename(SangerConsensusSet@
                                         consensusReadsList[[consensusReadIndex]]@
                                         reverseReadsList[[singleReadIndex]]@readFileName)
                        SangerSingleReadAFN <-
                            SangerConsensusSet@
                            consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@readFileName
                    }

                    message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                    h1(input$sidebar_menu)
                    fluidRow(
                        useShinyjs(),
                        box(title = tags$p(tagList(icon("dot-circle"),
                                                   "Raw File: "),
                                           style = "font-size: 26px;
                                                             font-weight: bold;"),
                            solidHeader = TRUE,
                            status = "success", width = 12,
                            h1(SangerSingleReadBFN),
                            tags$h5(paste("( full path:",
                                          SangerSingleReadAFN,
                                          ")"), style = "font-style:italic")),
                        box(title = tags$p(
                            tagList(icon("dot-circle"),
                                    "Primary, Secondary DNA Sequences & Amino Acid Sequence (Before Trimming):"),
                            style = "font-size: 26px;
                                                           font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                            column(width = 12,
                                   tags$p(tagList(icon("bars"),
                                                  "Primary Sequence"),
                                          style = "font-size: 22px;
                                                           font-weight: bold;"),
                                   excelOutput("primarySeqDF",
                                               width = "100%", height = "50"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p(tagList(icon("bars"),
                                                  "Secondary Sequence"),
                                          style = "font-size: 22px;
                                                           font-weight: bold;"),
                                   excelOutput("secondSeqDF",
                                               width = "100%", height = "50"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p(tagList(icon("bars"),
                                                  "Quality Phred Score"),
                                          style = "font-size: 22px;
                                                           font-weight: bold;"),
                                   excelOutput("qualityScoreDF",
                                               width = "100%", height = "50"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p(tagList(icon("bars"),
                                                  "AA Sequence 1"),
                                          style = "font-size: 22px;
                                                           font-weight: bold;"),
                                   excelOutput("PrimAASeqS1DF",
                                               width = "100%", height = "50"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p(tagList(icon("bars"),
                                                  "AA Sequence 2"),
                                          style = "font-size: 22px;
                                                           font-weight: bold;"),
                                   excelOutput("PrimAASeqS2DF",
                                               width = "100%", height = "50"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p(tagList(icon("bars"),
                                                  "AA Sequence 3"),
                                          style = "font-size: 22px;
                                                           font-weight: bold;"),
                                   excelOutput("PrimAASeqS3DF",
                                               width = "100%", height = "50"),
                                   style = paste("overflow-y: hidden;",
                                                 "overflow-x: scroll;")
                            ),
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
                                actionBttn("startTrimmingButton",
                                           "Apply Trimming Parameters",
                                           style = "simple", color = "success",
                                           block = TRUE, size = "lg")
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
    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic sidebar
    ### ------------------------------------------------------------------------
    # !!! Fix !!!
    dynamicMenuSideBarSCSet(input, output, session, SangerCSetParam)
    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        message("@@@@@@@ 'close button' has been clicked")
        btn <- input$closeUI
        stopApp()
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        message("@@@@@@@ 'save button' has been clicked")
        newS4Object <- file.path(shinyDirectory,
                                 "SangerAlignedConsensusSet.Rda")
        showNotification(paste("New S4 object is store as:", newS4Object),
                         closeButton = TRUE, type = "message", duration = 10)
        ### --------------------------------------------------------------------
        ### Save SangerConsensus quality S4 object
        ### --------------------------------------------------------------------
        saveRDS(SangerConsensusSet, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        NEW_SANGER_ALIGNED_CONSENSUS_READ_SET <<- readRDS(file=newS4Object)
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus reads re-calculating (SCSet) UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButtonSCSet, {
        message("######## Reactive button clicked !!!")
        message("######## Start recalculating consensus read (SC)")
        if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
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
        }
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus read re-calculating (SC) UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButton, {
        message("@@@@@@@ 'Reactive button' has been clicked")
        message("######## Start recalculating consensus read (SC")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(consensusReadIndex)) {
            CSResult<-
                calculateConsensusRead (
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@forwardReadsList,
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@reverseReadsList,
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@refAminoAcidSeq,
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@minFractionCall,
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@maxFractionLost,
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@geneticCode,
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@acceptStopCodons,
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@readingFrame)
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@consensusRead <<- CSResult$consensusGapfree
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@differencesDF <<- CSResult$diffsDf
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@alignment <<- CSResult$aln2
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@distanceMatrix <<- CSResult$dist
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@dendrogram <<- CSResult$dend
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@indelsDF <<- CSResult$indels
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@stopCodonsDF <<- CSResult$stopsDf
            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@secondaryPeakDF <<- CSResult$spDf
            consensusParam[["consensusRead"]] <<- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@consensusRead
            consensusParam[["differencesDF"]] <<- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@differencesDF
            consensusParam[["alignment"]] <<- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@alignment
            consensusParam[["distanceMatrix"]] <<-SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@distanceMatrix
            consensusParam[["dendrogram"]] <<- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@dendrogram
            consensusParam[["indelsDF"]] <<- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@indelsDF
            consensusParam[["stopCodonsDF"]] <<- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@stopCodonsDF
            consensusParam[["secondaryPeakDF"]] <<- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@secondaryPeakDF
        }
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button apply trimming parameters
    ### ------------------------------------------------------------------------
    observeEvent(input$startTrimmingButton, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if (SangerConsensusSet@consensusReadsList[[1]]@
                forwardReadsList[[1]]@QualityReport@
                TrimmingMethod == "M1") {
                if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
                    as.numeric(input$M1TrimmingCutoffText) > 0 &&
                    as.numeric(input$M1TrimmingCutoffText) <= 1) {
                    inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
                } else {
                    inputM1TrimmingCutoffText <- 0.0001
                }
                if (directionParam == "Forward") {
                    trimmingPos <-
                        M1inside_calculate_trimming(
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                forwardReadsList[[singleReadIndex]]@QualityReport@
                                qualityPhredScores,
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                forwardReadsList[[singleReadIndex]]@QualityReport@
                                qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        forwardReadsList[[singleReadIndex]]@QualityReport@
                        M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        forwardReadsList[[singleReadIndex]]@QualityReport@M1TrimmingCutoff
                } else if (directionParam == "Reverse") {
                    trimmingPos <-
                        M1inside_calculate_trimming(
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                reverseReadsList[[singleReadIndex]]@QualityReport@
                                qualityPhredScores,
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                reverseReadsList[[singleReadIndex]]@QualityReport@
                                qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        reverseReadsList[[singleReadIndex]]@QualityReport@
                        M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        reverseReadsList[[singleReadIndex]]@QualityReport@M1TrimmingCutoff
                }
            } else if (SangerCSetParam[[1]]$
                       SangerSingleReadQualReport[[1]]@
                       TrimmingMethod == "M2") {
                if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
                    strtoi(input$M2CutoffQualityScoreText) > 0 &&
                    strtoi(input$M2CutoffQualityScoreText) <= 60 &&
                    strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
                    inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
                } else {
                    inputM2CutoffQualityScoreText <- 20
                }
                if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
                    strtoi(input$M2SlidingWindowSizeText) > 0 &&
                    strtoi(input$M2SlidingWindowSizeText) <= 20 &&
                    strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
                    inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
                } else {
                    inputM2SlidingWindowSizeText <- 5
                }
                if (directionParam == "Forward") {
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                forwardReadsList[[singleReadIndex]]@QualityReport@
                                qualityPhredScores,
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                forwardReadsList[[singleReadIndex]]@QualityReport@
                                qualityBaseScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText)
                        )
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        forwardReadsList[[singleReadIndex]]@QualityReport@
                        M2CutoffQualityScore <<- strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        forwardReadsList[[singleReadIndex]]@QualityReport@
                        M2CutoffQualityScore
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        forwardReadsList[[singleReadIndex]]@QualityReport@
                        M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        forwardReadsList[[singleReadIndex]]@QualityReport@
                        M2SlidingWindowSize
                } else if ("Reverse") {
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                reverseReadsList[[singleReadIndex]]@QualityReport@
                                qualityPhredScores,
                            SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                reverseReadsList[[singleReadIndex]]@QualityReport@
                                qualityBaseScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText)
                        )
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        reverseReadsList[[singleReadIndex]]@QualityReport@
                        M2CutoffQualityScore <<- strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        reverseReadsList[[singleReadIndex]]@QualityReport@
                        M2CutoffQualityScore
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        reverseReadsList[[singleReadIndex]]@QualityReport@
                        M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                        reverseReadsList[[singleReadIndex]]@QualityReport@
                        M2SlidingWindowSize
                }
            }
            if (directionParam == "Forward") {
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]

                trimmedRV[["rawSeqLength"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                              forwardReadsList[[singleReadIndex]]@QualityReport@
                              remainingRatio * 100, 2)
            } else if (directionParam == "Reverse") {
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]

                trimmedRV[["rawSeqLength"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                              reverseReadsList[[singleReadIndex]]@QualityReport@
                              remainingRatio * 100, 2)
            }
        }
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$saveChromatogramParam, {
        message("@@@@@@@ 'Reactive button' has been clicked")
        message("######## Start recalculating chromatogram")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        ### ------------------------------------------------------------
        ### Update ChromatogramBasePerRow
        ### ------------------------------------------------------------
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if (directionParam == "Forward") {
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@ChromatogramParam@
                    showTrimmed <<- input$ChromatogramCheckShowTrimmed
            } else if (directionParam == "Reverse") {
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@ChromatogramParam@
                    showTrimmed <<- input$ChromatogramCheckShowTrimmed
            }
            ### --------------------------------------------------------------------
            ### Save 'ChromatogramParam' dynamic value
            ### --------------------------------------------------------------------
            ChromatogramParam[["baseNumPerRow"]] <<-
                input$ChromatogramBasePerRow
            ChromatogramParam[["heightPerRow"]] <<-
                input$ChromatogramHeightPerRow
            ChromatogramParam[["signalRatioCutoff"]] <<-
                input$ChromatogramSignalRatioCutoff
            ChromatogramParam[["showTrimmed"]] <<-
                input$ChromatogramCheckShowTrimmed
        }
    })



    ############################################################################
    ### ConsensusReadSet (Function for Sanger Consensus Read Set Overview)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Alignment
    ### ------------------------------------------------------------------------
    output$consensusSetAlignmentHTML<-renderUI({
        if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
            # consensusParamSet[["alignmentSCSet"]] <-
            #     SangerConsensusSet@alignmentSCSet
            browseSeqHTML <-
                file.path(shinyDirectory,
                          "Consensus_Readset_Alignment_BrowseSeqs.html")
            BrowseSeqs(consensusParamSet[["alignmentSCSet"]] ,
                       openURL=FALSE, htmlFile=browseSeqHTML)
            includeHTML(browseSeqHTML)
        }
    })
    ### ------------------------------------------------------------------------
    ### Consensus Reads Tree
    ### ------------------------------------------------------------------------
    output$SCSetConsensusReadTreePlot <- renderPlot({
        if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
            # consensusParamSet[["alignmentTreeSCSet"]] <<- SangerConsensusSet@alignmentTreeSCSet
            plot(consensusParamSet[["alignmentTreeSCSet"]])
        }
    })



    ############################################################################
    ### ConsensusRead / ConsensusReadSet Share
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



    ############################################################################
    ### ConsensusRead (Function for Sanger Consensus Read Overview)
    ############################################################################
    valueBoxSCMinReadsNumCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCMinReadLengthCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCMinFractionCallCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCMaxFractionLostCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCAcceptStopCodonsCSSet (input, output, SangerConsensusSet, session)
    valueBoxSCReadingFrameCSSet (input, output, SangerConsensusSet, session)
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
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            primarySeqDisplay (sequenceParam)
        }
    })
    output$secondSeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            secondarySeqDisplay (sequenceParam)
        }
    })
    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if (directionParam == "Forward") {
                PhredScore <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                PhredScore <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityScoreDisplay (PhredScore)
        }
    })
    output$PrimAASeqS1DF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            PrimAASeqS1Display (sequenceParam)
        }
    })
    output$PrimAASeqS2DF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            PrimAASeqS2Display (sequenceParam)
        }
    })
    output$PrimAASeqS3DF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            PrimAASeqS3Display (sequenceParam)
        }
    })
    ############################################################################
    ### Trimming method selection functions
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Trimming Method Selection (Just check forward reads => not dynamic)
    ### ------------------------------------------------------------------------
    output$TrimmingMethodSelectionOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if (SangerConsensusSet@consensusReadsList[[1]]@
                forwardReadsList[[1]]@QualityReport@TrimmingMethod == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                            'Logarithmic Scale Trimming'")
            } else if (SangerConsensusSet@consensusReadsList[[1]]@
                       forwardReadsList[[1]]@QualityReport@TrimmingMethod == "M2") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                            'Logarithmic Scale Sliding Window Trimming'")
            }
        }
    })
    output$TrimmingMethodUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            ## For method, everyone is same, so just pick forward one.
            if (SangerConsensusSet@consensusReadsList[[1]]@
                forwardReadsList[[1]]@QualityReport@TrimmingMethod== "M1") {
                if (directionParam == "Forward") {
                    if (is.null(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                forwardReadsList[[singleReadIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@
                            QualityReport@M1TrimmingCutoff <<-  0.0001
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                reverseReadsList[[singleReadIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@
                            QualityReport@M1TrimmingCutoff <<-  0.0001
                    }
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
            } else if (SangerConsensusSet@consensusReadsList[[1]]@
                       forwardReadsList[[1]]@QualityReport@TrimmingMethod == "M2") {


                if (directionParam == "Forward") {
                    if (is.null(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                forwardReadsList[[singleReadIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                forwardReadsList[[singleReadIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            forwardReadsList[[singleReadIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  5
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                reverseReadsList[[singleReadIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                reverseReadsList[[singleReadIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                            reverseReadsList[[singleReadIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  5
                    }
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

    valueBoxM1TrimmingCutoff (input, output, session)
    valueBoxM2CutoffQualityScore (input, output, session)
    valueBoxM2SlidingWindowSize (input, output, session)

    valueBoxRawSeqLength (input, output, session, trimmedRV)
    valueBoxRawMeanQualityScore (input, output, session, trimmedRV)
    valueBoxRawMinQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxTrimmedFinishPos (input, output, session, trimmedRV)
    valueBoxTrimmedSeqLength (input, output, session, trimmedRV)
    valueBoxTrimmedMeanQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedMinQualityScore (input, output, session, trimmedRV)
    valueBoxRemainingRatio (input, output, session, trimmedRV)

    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if (directionParam == "Forward") {
                qualityPhredScores <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                qualityPhredScores <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityTrimmingRatioPlotDisplay(input, output,session,
                                            trimmedRV, qualityPhredScores)
        }
    })
    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(consensusReadIndex) &&
            !is.na(singleReadIndex)) {
            if (directionParam == "Forward") {
                qualityPhredScores <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                qualityPhredScores <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityQualityBasePlotDisplay(input, output,session,
                                          trimmedRV, qualityPhredScores)
        }
    })

    valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)

    # chromatogram
    output$chromatogramUIOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (sidebar_menu[[2]] == "CR" &&
            sidebar_menu[[3]] == "-" &&
            !is.na(singleReadIndex) &&
            (directionParam == "Forward" ||
             directionParam == "Reverse") &&
            sidebar_menu[[6]] == "Read") {
            trimmedRV[["trimmedSeqLength"]]

            if (directionParam == "Forward") {
                rawSeqLength <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                trimmedSeqLength <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport@trimmedSeqLength
            } else if (directionParam == "Reverse") {
                rawSeqLength <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@
                    QualityReport@rawSeqLength
                trimmedSeqLength <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@
                    QualityReport@trimmedSeqLength
            }
            chromatogramRowNumAns <- chromatogramRowNum (
                strtoi(ChromatogramParam[["baseNumPerRow"]]),
                rawSeqLength, trimmedSeqLength,
                ChromatogramParam[["showTrimmed"]]) *
                strtoi(ChromatogramParam[["heightPerRow"]])
            plotOutput("chromatogram", height = chromatogramRowNumAns)
        }
    })
    output$chromatogram <- renderPlot({
        ## !!!!! Update !!!!
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (sidebar_menu[[2]] == "CR" &&
            sidebar_menu[[3]] == "-" &&
            !is.na(singleReadIndex) &&
            (directionParam == "Forward" ||
             directionParam == "Reverse") &&
            sidebar_menu[[6]] == "Read") {
            showNotification(
                ui = fluidRow(
                    column(1),
                    column(11,
                           tags$p(tagList(icon("hourglass-start"),
                                          "Plotting Chromatogram ...  "),
                                  style = "font-size: 28px;
                                   font-weight: bold;"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Base pairs per row : ",
                                          ChromatogramParam[["baseNumPerRow"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Height per row : ",
                                          ChromatogramParam[["heightPerRow"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Signal Ratio Cutoff : ",
                                          ChromatogramParam[["signalRatioCutoff"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Show trimmed : ",
                                          ChromatogramParam[["showTrimmed"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),)
                ),
                action = NULL, duration = NULL, closeButton = FALSE,
                id = "chromatogramNotification",
                type = "message")
            if (directionParam == "Forward") {
                rawSeqLength <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@QualityReport@rawSeqLength
                message(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (forward)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                      forwardReadsList[[singleReadIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ----------------------------------------------------------------
                ### Update 'SangerConsensus'!
                ### ----------------------------------------------------------------
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ----------------------------------------------------------------
                ### Updating AASeqs
                ### ----------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerConsensusSet@geneticCode)
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    forwardReadsList[[singleReadIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ----------------------------------------------------------------
                ### Updating reactive values
                ### ----------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     forwardReadsList[[singleReadIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     forwardReadsList[[singleReadIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     forwardReadsList[[singleReadIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     forwardReadsList[[singleReadIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     forwardReadsList[[singleReadIndex]]@primaryAASeqS3)
            } else if (directionParam == "Reverse") {
                rawSeqLength <-
                    SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@QualityReport@rawSeqLength
                message(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (reverse)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                      reverseReadsList[[singleReadIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ----------------------------------------------------------------
                ### Update 'SangerConsensus'!
                ### ----------------------------------------------------------------
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ----------------------------------------------------------------
                ### Updating AASeqs
                ### ----------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerConsensusSet@geneticCode)
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                    reverseReadsList[[singleReadIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ----------------------------------------------------------------
                ### Updating reactive values
                ### ----------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     reverseReadsList[[singleReadIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     reverseReadsList[[singleReadIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     reverseReadsList[[singleReadIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     reverseReadsList[[singleReadIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@
                                     reverseReadsList[[singleReadIndex]]@primaryAASeqS3)
            }
            # message(">>>>>>>>>>>> 'MakeBaseCalls' finished")
            chromatogram(hetcalls,
                         width = strtoi(
                             ChromatogramParam[["baseNumPerRow"]]),
                         height = 2,
                         trim5 = trimmedRV[["trimmedStartPos"]],
                         trim3 = rawSeqLength -
                             trimmedRV[["trimmedFinishPos"]],
                         showtrim = (ChromatogramParam[["showTrimmed"]]),
                         showcalls = "both")
            removeNotification(id = "chromatogramNotification")
        }
    })
}

