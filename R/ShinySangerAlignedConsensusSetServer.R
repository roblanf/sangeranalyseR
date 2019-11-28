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
    SangerCSetParentDir <- SangerConsensusSet@parentDirectory

    SangerCSetSuffixForwardRegExp <- SangerConsensusSet@suffixForwardRegExp
    SangerCSetSuffixReverseRegExp <- SangerConsensusSet@suffixReverseRegExp
    SangerCSetList <- SangerConsensusSet@consensusReadsList
    SangerConsensusSetNum <- length(SangerCSetList)

    SangerCSetParam <- lapply(1:SangerConsensusSetNum, function(i) {
        ### --------------------------------------------------------------------
        ### ConsensusRead-related parameters initialization.
        ### --------------------------------------------------------------------
        # readFeature
        SCName <- paste0(i, " Consensus Read")
        SCMinReadsNum <- SangerCSetList[[i]]@minReadsNum
        SCMinReadLength <- SangerCSetList[[i]]@minReadLength
        SCRefAminoAcidSeq <- SangerCSetList[[i]]@refAminoAcidSeq
        SCMinFractionCall <- SangerCSetList[[i]]@minFractionCall
        SCMaxFractionLost <- SangerCSetList[[i]]@maxFractionLost
        SCGeneticCode <- SangerCSetList[[i]]@geneticCode
        SCAcceptStopCodons <- SangerCSetList[[i]]@acceptStopCodons
        SCReadingFrame <- SangerCSetList[[i]]@readingFrame
        SCConsensusRead <- as.character(SangerCSetList[[i]]@consensusRead)
        SCAlignment<- SangerCSetList[[i]]@alignment
        SCDifferencesDF<- SangerCSetList[[i]]@differencesDF
        SCDistanceMatrix <- SangerCSetList[[i]]@distanceMatrix
        SCDendrogram <- SangerCSetList[[i]]@dendrogram
        SCIndelsDF <- SangerCSetList[[i]]@indelsDF
        SCStopCodonsDF <- SangerCSetList[[i]]@stopCodonsDF
        SCSecondaryPeakDF <- SangerCSetList[[i]]@secondaryPeakDF
        SCconsenesusReadName <- SangerCSetList[[i]]@consenesusReadName
        SangerConsensusForRegExp <- SangerCSetList[[i]]@suffixForwardRegExp
        SangerConsensusRevRegExp <- SangerCSetList[[i]]@suffixReverseRegExp

        # Forward & reverse reads list
        SangerSingleReadFReadsList <- SangerCSetList[[i]]@forwardReadsList
        SangerSingleReadRReadsList <- SangerCSetList[[i]]@reverseReadsList
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

        # abifRawData
        forwardReadAbifRawData <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@abifRawData)
        reverseReadAbifRawData <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@abifRawData)
        SangerSingleReadAbifRawData <- c(forwardReadAbifRawData,
                                         reverseReadAbifRawData)

        # QualityReport
        forwardReadQualReport <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@QualityReport)
        reverseReadQualReport <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@QualityReport)
        SangerSingleReadQualReport <- c(forwardReadQualReport,
                                        reverseReadQualReport)

        # primarySeqID
        forwardReadPrimSeqID <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@primarySeqID)
        reverseReadPrimSeqID <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@primarySeqID)
        SangerSingleReadPrimSeqID <- c(forwardReadPrimSeqID,
                                       reverseReadPrimSeqID)

        # primarySeq
        forwardReadPrimSeq <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@primarySeq)
        reverseReadPrimSeq <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@primarySeq)
        SangerSingleReadPrimSeq <- c(forwardReadPrimSeq, reverseReadPrimSeq)


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
            rownames(basecalls1DF) <- "Primary"
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
            rownames(basecalls1DF) <- "Primary "
            return(basecalls1DF)
            }
        )
        SangerSingleReadPrimSeqDF <- c(forwardReadPrimSeqDF,
                                       reverseReadPrimSeqDF)


        # secondarySeqID
        forwardReadSecoSeqID <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@secondarySeqID)
        reverseReadSecoSeqID <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@secondarySeqID)
        SangerSingleReadSecoSeqID <- c(forwardReadSecoSeqID,
                                       reverseReadSecoSeqID)

        # secondarySeq
        forwardReadSecoSeq <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@secondarySeq)
        reverseReadSecoSeq <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@secondarySeq)
        SangerSingleReadSecoSeq <- c(forwardReadSecoSeq, reverseReadSecoSeq)

        # secondarySeq
        forwardReadSecoSeq <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@secondarySeq)
        reverseReadSecoSeq <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@secondarySeq)
        SangerSingleReadSecoSeq <- c(forwardReadSecoSeq, reverseReadSecoSeq)

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
            rownames(basecalls2DF) <- " Second"
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
            rownames(basecalls2DF) <- " Second"
            return(basecalls2DF)
            }
        )
        SangerSingleReadSecoSeqDF <- c(forwardReadSecoSeqDF,
                                       reverseReadSecoSeqDF)


        # traceMatrix
        forwardReadTraceMat <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@traceMatrix)
        reverseReadTraceMat <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@traceMatrix)
        SangerSingleReadTraceMat <- c(forwardReadTraceMat,
                                      reverseReadTraceMat)

        # peakPosMatrix
        forwardReadReadPeakPosMat <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@peakPosMatrix)
        reverseReadReadPeakPosMat <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@peakPosMatrix)
        SangerSingleReadPeakPosMat <- c(forwardReadReadPeakPosMat,
                                        reverseReadReadPeakPosMat)
        # peakAmpMatrix
        forwardReadPeakAmpMat <- sapply(1:forwardReadNum, function(j)
            SangerSingleReadFReadsList[[j]]@peakAmpMatrix)
        reverseReadPeakAmpMat <- sapply(1:reverseReadNum, function(j)
            SangerSingleReadRReadsList[[j]]@peakAmpMatrix)
        SangerSingleReadPeakAmpMat <- c(forwardReadPeakAmpMat,
                                        reverseReadPeakAmpMat)
        return(list(SCName = SCName,
                    SCMinReadsNum = SCMinReadsNum,
                    SCMinReadLength = SCMinReadLength,
                    SCRefAminoAcidSeq = SCRefAminoAcidSeq,
                    SCMinFractionCall = SCMinFractionCall,
                    SCMaxFractionLost = SCMaxFractionLost,
                    SCGeneticCode = SCGeneticCode,
                    SCAcceptStopCodons = SCAcceptStopCodons,
                    SCReadingFrame = SCReadingFrame,
                    SCConsensusRead = SCConsensusRead,
                    SCAlignment = SCAlignment,
                    SCDifferencesDF = SCDifferencesDF,
                    SCDistanceMatrix = SCDistanceMatrix,
                    SCDendrogram = SCDendrogram,
                    SCIndelsDF = SCIndelsDF,
                    SCStopCodonsDF = SCStopCodonsDF,
                    SCSecondaryPeakDF = SCSecondaryPeakDF,
                    SCconsenesusReadName = SCconsenesusReadName,
                    SangerConsensusForRegExp = SangerConsensusForRegExp,
                    SangerConsensusRevRegExp = SangerConsensusRevRegExp,
                    SangerSingleReadFReadsList = SangerSingleReadFReadsList,
                    SangerSingleReadRReadsList = SangerSingleReadRReadsList,
                    forwardReadNum = forwardReadNum,
                    reverseReadNum = reverseReadNum,
                    SangerSingleReadNum = SangerSingleReadNum,
                    SangerConsensusFRReadsList = SangerConsensusFRReadsList,
                    SangerSingleReadBFN = SangerSingleReadBFN,
                    SangerSingleReadAFN = SangerSingleReadAFN,
                    SangerSingleReadFeature = SangerSingleReadFeature,
                    SangerSingleReadAbifRawData = SangerSingleReadAbifRawData,
                    SangerSingleReadQualReport = SangerSingleReadQualReport,
                    SangerSingleReadPrimSeqID = SangerSingleReadPrimSeqID,
                    SangerSingleReadPrimSeq = SangerSingleReadPrimSeq,
                    forwardReadPrimSeqDF = forwardReadPrimSeqDF,
                    reverseReadPrimSeqDF = reverseReadPrimSeqDF,
                    SangerSingleReadPrimSeqDF = SangerSingleReadPrimSeqDF,
                    SangerSingleReadSecoSeqID = SangerSingleReadSecoSeqID,
                    SangerSingleReadSecoSeq = SangerSingleReadSecoSeq,
                    forwardReadSecoSeqDF = forwardReadSecoSeqDF,
                    reverseReadSecoSeqDF = reverseReadSecoSeqDF,
                    SangerSingleReadSecoSeqDF = SangerSingleReadSecoSeqDF,
                    SangerSingleReadTraceMat = SangerSingleReadTraceMat,
                    SangerSingleReadPeakPosMat = SangerSingleReadPeakPosMat,
                    SangerSingleReadPeakAmpMat = SangerSingleReadPeakAmpMat))
    })

    trimmedRV <- reactiveValues(rawSeqLength = 0,
                                rawMeanQualityScore = 0,
                                rawMinQualityScore = 0,
                                trimmedStartPos = 0,
                                trimmedFinishPos = 0,
                                trimmedSeqLength = 0,
                                trimmedMeanQualityScore = 0,
                                trimmedMinQualityScore = 0,
                                remainingRatio = 0)

    ############################################################################
    ### output$ID
    ############################################################################
    dynamicMenuSideBarSCSet(input, output, session, SangerCSetParam)

    output$aligned_consensusRead_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
            h1(paste("You've selected:", input$sidebar_menu))
            box(title = tags$p("Input Parameters: ",
                               style = "font-size: 26px;
                                       font-weight: bold;"),
                solidHeader = TRUE, collapsible = TRUE,
                status = "success", width = 12,
                tags$hr(style = ("border-top: 0.2px hidden #A9A9A9;")),
                fluidRow(
                    column(12,
                           column(3,
                                  h4("Output Directory: ",
                                     style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(shinyDirectory),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Raw ABI Parent Directory: ",
                                     style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(SangerCSetParentDir),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Forward Suffix RegExp: ",
                                     style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(SangerCSetSuffixForwardRegExp),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Reverse Suffix RegExp: ",
                                     style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(SangerCSetSuffixReverseRegExp),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Consensus Read Number: ",
                                     style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(SangerConsensusSetNum),
                           )
                    ),
                )
            )
        } else {
            if (!is.na(as.numeric(sidebar_menu[[1]]))) {
                if (sidebar_menu[[2]] == "Sanger" &&
                    sidebar_menu[[3]] == "Consensus" &&
                    sidebar_menu[[4]] == "Read" &&
                    sidebar_menu[[5]] == "Overview") {
                    consensusReadIndex <- strtoi(sidebar_menu[[1]])

                    SCName <- SangerCSetParam[[consensusReadIndex]]$SCName
                    SCMinReadsNum <-
                        SangerCSetParam[[consensusReadIndex]]$SCMinReadsNum
                    SCMinReadLength <-
                        SangerCSetParam[[consensusReadIndex]]$SCMinReadLength
                    SCRefAminoAcidSeq <-
                        SangerCSetParam[[consensusReadIndex]]$SCRefAminoAcidSeq
                    SCMinFractionCall <-
                        SangerCSetParam[[consensusReadIndex]]$SCMinFractionCall
                    SCMaxFractionLost <-
                        SangerCSetParam[[consensusReadIndex]]$SCMaxFractionLost
                    SCGeneticCode <-
                        SangerCSetParam[[consensusReadIndex]]$SCGeneticCode
                    SCAcceptStopCodons <-
                        SangerCSetParam[[consensusReadIndex]]$SCAcceptStopCodons
                    SCReadingFrame <-
                        SangerCSetParam[[consensusReadIndex]]$SCReadingFrame
                    SCConsensusRead <-
                        SangerCSetParam[[consensusReadIndex]]$SCConsensusRead
                    SCAlignment <-
                        SangerCSetParam[[consensusReadIndex]]$SCAlignment
                    SCDifferencesDF <-
                        SangerCSetParam[[consensusReadIndex]]$SCDifferencesDF
                    SCDistanceMatrix <-
                        SangerCSetParam[[consensusReadIndex]]$SCDistanceMatrix
                    SCDendrogram <-
                        SangerCSetParam[[consensusReadIndex]]$SCDendrogram
                    SCIndelsDF <-
                        SangerCSetParam[[consensusReadIndex]]$SCIndelsDF
                    SCStopCodonsDF <-
                        SangerCSetParam[[consensusReadIndex]]$SCStopCodonsDF
                    SCSecondaryPeakDF <-
                        SangerCSetParam[[consensusReadIndex]]$SCSecondaryPeakDF
                    SCconsenesusReadName <-
                        SangerCSetParam[[consensusReadIndex]]$
                        SCconsenesusReadName
                    SangerConsensusForRegExp <-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerConsensusForRegExp
                    SangerConsensusRevRegExp <-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerConsensusRevRegExp
                    SangerSingleforwardReadNum <-
                        SangerCSetParam[[consensusReadIndex]]$forwardReadNum
                    SangerSinglereverseReadNum <-
                        SangerCSetParam[[consensusReadIndex]]$reverseReadNum
                    SangerSingleReadNum <-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadNum
                    fluidRow(
                        useShinyjs(),
                        box(title = tags$p("Input Parameters: ",
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(
                                style = ("border-top: 0.2px hidden #A9A9A9;")),
                            fluidRow(
                                column(12,
                                       column(3,
                                              h4("Consenesus Read Number: ",
                                                 style="font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SCName),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              h4("Consenesus Read Name: ",
                                                 style="font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SCconsenesusReadName),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              h4("Forward Suffix RegExp: ",
                                                 style="font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SangerCSetSuffixForwardRegExp),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              h4("Forward Read Number: ",
                                                 style="font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SangerSingleforwardReadNum),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              h4("Reverse Suffix RegExp: ",
                                                 style="font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SangerCSetSuffixReverseRegExp),
                                       )
                                ),
                                column(12,
                                       column(3,
                                              h4("Reverse Read Number: ",
                                                 style="font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SangerSinglereverseReadNum),
                                       )
                                ),
                            ),
                            ################################################
                            #### Add this after having reference sample ####
                            ################################################
                            # If it is null
                            # h1(SCRefAminoAcidSeq),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
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
                                           style = paste("height:100%;",
                                                         "overflow-y: hidden;",
                                                         "overflow-x: scroll;")
                                    ),
                                ),
                            ),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                column(3,
                                       valueBox(
                                           subtitle =
                                               tags$p("MinReadsNum",
                                                      style = "font-size: 15px;
                                                      font-weight: bold;"),
                                           value =
                                               tags$p(strtoi(SCMinReadsNum),
                                                     style ="font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       )
                                       # uiOutput("SCMinReadsNum") ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle =
                                               tags$p("MinReadLength",
                                                      style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value =
                                               tags$p(strtoi(SCMinReadLength),
                                                      style="font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                       # uiOutput("SCMinReadLength")  ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle =
                                               tags$p("MinFractionCall",
                                                      style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value =
                                               tags$p(
                                                  as.numeric(SCMinFractionCall),
                                                  style = "font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                ),
                                column(3,
                                       valueBox(
                                           subtitle =
                                               tags$p("MaxFractionLost",
                                                      style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value =
                                               tags$p(
                                                  as.numeric(SCMaxFractionLost),
                                                  style = "font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                ),
                                column(3,
                                       valueBox(
                                           subtitle =
                                               tags$p("AcceptStopCodons",
                                                      style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value =
                                               tags$p(SCAcceptStopCodons,
                                                      style="font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                       # uiOutput("SCAcceptStopCodons") ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle =
                                               tags$p("ReadingFrame",
                                                      style = "font-size: 15px;
                                                      font-weight: bold;"),
                                           value =
                                               tags$p(strtoi(SCReadingFrame),
                                                      style="font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                ),
                            ),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                        ),

                        box(title = tags$p("Results: ",
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(
                                style = ("border-top: 4px hidden #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Alignment",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           htmlOutput("consensusAlignmentHTML"),
                                    ),
                                ),
                            ),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Differences Dataframe",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           dataTableOutput("differencesDF") ,
                                           style = paste0("height:100%;",
                                                          "overflow-y:scroll;",
                                                          "overflow-x: scroll;")
                                    )
                                ),
                            ),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Dendrogram",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           # plot()
                                           plotOutput("dendrogramPlot"),
                                           style = paste0("height:100%;",
                                                          "overflow-y:scroll;",
                                                          "overflow-x: scroll;")
                                    ),
                                    column(width = 12,
                                           dataTableOutput("dendrogramDF"),
                                           style = paste0("height:100%;",
                                                          "overflow-y:scroll;",
                                                          "overflow-x: scroll;")
                                    )
                                ),
                            ),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Secondary Peak Dataframe",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           uiOutput("secondaryPeakDFUI"),
                                           style = paste0("height:100%;",
                                                          "overflow-y:scroll;",
                                                          "overflow-x: scroll;")
                                    )
                                ),
                            ),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Indels Dataframe",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           uiOutput("SCIndelsDFUI"),
                                           style = paste0("height:100%;",
                                                          "overflow-y:scroll;",
                                                          "overflow-x: scroll;")
                                    )
                                ),
                            ),
                            tags$hr(
                                style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Stop Codons Dataframe",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           uiOutput("SCStopCodonsDFUI"),
                                           style = paste0("height:100%;",
                                                          "overflow-y:scroll;",
                                                          "overflow-x: scroll;")
                                    )
                                )
                            ),
                        ),
                    )
                } else if (sidebar_menu[[2]] == "Consensus" &&
                           sidebar_menu[[3]] == "Read" &&
                           sidebar_menu[[4]] == "-" &&
                           !is.na(as.numeric(sidebar_menu[[5]])) &&
                           (sidebar_menu[[6]] == "Forward" ||
                            sidebar_menu[[6]] == "Reverse") &&
                           sidebar_menu[[7]] == "Read") {
                    consensusReadIndex <- strtoi(sidebar_menu[[1]])
                    singleReadIndex <- strtoi(sidebar_menu[[5]])
                    SangerConsensusFRReadsList <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerConsensusFRReadsList
                    SangerSingleReadBFN <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadBFN
                    SangerSingleReadAFN <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadAFN
                    SangerSingleReadFeature <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadFeature
                    SangerSingleReadAbifRawData <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadAbifRawData
                    SangerSingleReadQualReport <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadQualReport
                    SangerSingleReadPrimSeqID <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadPrimSeqID
                    SangerSingleReadPrimSeq <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadPrimSeq
                    forwardReadPrimSeqDF <-
                        SangerCSetParam[[consensusReadIndex]]$
                            forwardReadPrimSeqDF
                    reverseReadPrimSeqDF <-
                        SangerCSetParam[[consensusReadIndex]]$
                            reverseReadPrimSeqDF
                    SangerSingleReadPrimSeqDF <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadPrimSeqDF
                    SangerSingleReadSecoSeqID <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadSecoSeqID
                    SangerSingleReadSecoSeq <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadSecoSeq
                    forwardReadSecoSeqDF <-
                        SangerCSetParam[[consensusReadIndex]]$
                            forwardReadPrimSeqDF
                    reverseReadSecoSeqDF <-
                        SangerCSetParam[[consensusReadIndex]]$
                            reverseReadSecoSeqDF
                    SangerSingleReadSecoSeqDF <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadSecoSeqDF
                    SangerSingleReadTraceMat <-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadTraceMat
                    SangerSingleReadPeakPosMat <-
                        SangerCSetParam[[consensusReadIndex]]$
                            SangerSingleReadPeakPosMat
                    SangerSingleReadPeakAmpMat <-
                        SangerCSetParam[[consensusReadIndex]]$
                        SangerSingleReadPeakAmpMat













                    fluidRow(
                        useShinyjs(),
                        box(title = tags$p("Raw File: ",
                                           style = "font-size: 26px;
                                         font-weight: bold;"),
                            solidHeader = TRUE,
                            status = "success", width = 12,
                            h1(paste0(
                                SangerSingleReadBFN[[singleReadIndex]])),
                            tags$h5(paste("( full path:",
                                          SangerSingleReadAFN[[
                                              singleReadIndex]],
                                          ")"), style = "font-style:italic")),
                        box(title = tags$p("Primary & Secondary Peaks: ",
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                            column(width = 1,
                                   tags$p("Primary",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                                   tags$br(),
                                   tags$br(),
                                   tags$p("Second",
                                          style = "font-size: 15px;
                                       font-weight: bold;"),
                            ),
                            column(width = 11,
                                   excelOutput("primarySeqDF",
                                               width = "100%", height = "50"),
                                   excelOutput("secondSeqDF",
                                               width = "100%", height = "50"),
                                   style = paste("overflow-y: hidden;",
                                                 "overflow-x: scroll;")
                            )
                        ),
                        box(title = tags$p("Quality Report: ",
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                            box(title = tags$p(tagList(shiny::icon("arrow-circle-right"),
                                                       "Trimming Parameters Input"),
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,

                                fluidRow(
                                    column(width = 5,
                                           column(width = 1,
                                           ),
                                           column(width = 11,
                                                  selectInput("TrimmingMethodSelection", label = h4("Select Your Trimming Method"),
                                                              choices = list("Method 1" = "M1", "Method 2" = "M2"),
                                                              selected = SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@TrimmingMethod,
                                                              width = "100%"),
                                                  column(width = 12,
                                                         column(width = 2,
                                                                shiny::icon("mouse-pointer"),
                                                         ),
                                                         column(width = 10,
                                                                textOutput("TrimmingMethodSelectionOutput"),
                                                         ),
                                                  ),
                                           ),
                                    ),
                                    column(width = 7,
                                           uiOutput("TrimmingMethodUI") ,
                                    ),
                                ),
                            ),
                            box(title = tags$p(tagList(shiny::icon("arrow-circle-left"),
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
                                                      uiOutput("rawMeanQualityScore") ,
                                               ),
                                               column(4,
                                                      uiOutput("rawMinQualityScore") ,
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
                                                      uiOutput("trimmedSeqLength") ,
                                               ),
                                               column(4,
                                                      uiOutput("trimmedMeanQualityScore") ,
                                               ),
                                               column(4,
                                                      uiOutput("trimmedMinQualityScore") ,
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
                                tags$hr(style = ("border-top: 6px double #A9A9A9;")),
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
                        box(title = tags$p("Chromatogram: ",
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),

                            box(title = tags$p(tagList(shiny::icon("arrow-circle-right"),
                                                       "Chromatogram Input"),
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(12,
                                       column(3,
                                              sliderInput("ChromatogramBasePerRow",
                                                          label = h3("Slider"), min = 5,
                                                          max = 200, value = 100),
                                       ),
                                       column(3,
                                              uiOutput("ChromatogramtrimmedStartPos"),
                                       ),
                                       column(3,
                                              uiOutput("ChromatogramtrimmedFinishPos"),
                                       ),
                                       column(3,
                                              numericInput(
                                                  "ChromatogramSignalRatioCutoff",
                                                  h3("Signal Ratio Cutoff"),
                                                  value = 0.33),
                                              checkboxInput(
                                                  "ChromatogramCheckShowTrimmed",
                                                  "Whether show trimmed region",
                                                  value = TRUE),)
                                ),
                            ),
                            box(title = tags$p(tagList(shiny::icon("arrow-circle-left"),
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












                    # fluidRow(
                    #     useShinyjs(),
                    #     box(title = tags$p("Raw File: ",
                    #                        style = "font-size: 26px;
                    #                      font-weight: bold;"),
                    #         solidHeader = TRUE,
                    #         status = "success", width = 12,
                    #         h1(paste0(
                    #             SangerSingleReadBFN[[singleReadIndex]])),
                    #         tags$h5(paste("( full path:",
                    #                       SangerSingleReadAFN[[
                    #                           singleReadIndex]],
                    #                       ")"), style = "font-style:italic")),
                    #     box(title = tags$p("Primary & Secondary Peaks: ",
                    #                        style = "font-size: 26px;
                    #                    font-weight: bold;"),
                    #         solidHeader = TRUE, collapsible = TRUE,
                    #         status = "success", width = 12,
                    #         tags$hr(
                    #             style = ("border-top: 4px hidden #A9A9A9;")),
                    #         column(width = 1,
                    #                tags$p("Primary",
                    #                       style = "font-size: 15px;
                    #                    font-weight: bold;"),
                    #                tags$br(),
                    #                tags$br(),
                    #                tags$p("Second",
                    #                       style = "font-size: 15px;
                    #                    font-weight: bold;"),
                    #         ),
                    #         column(width = 11,
                    #                excelOutput("primarySeqDF",
                    #                            width = "100%", height = "50"),
                    #                excelOutput("secondSeqDF",
                    #                            width = "100%", height = "50"),
                    #                style = paste("overflow-y: hidden;",
                    #                              "overflow-x: scroll;")
                    #         )
                    #     ),
                    #     box(title = tags$p("Quality Report: ",
                    #                        style = "font-size: 26px;
                    #                    font-weight: bold;"),
                    #         solidHeader = TRUE, collapsible = TRUE,
                    #         status = "success", width = 12,
                    #         tags$hr(
                    #             style = ("border-top: 4px hidden #A9A9A9;")),
                    #         fluidRow(
                    #             column(3,
                    #                    uiOutput("M2CutoffQualityScore") ,
                    #                    tags$ul(
                    #                        textInput(
                    #                            "M2CutoffQualityScoreText",
                    #                            label = p("Change Value"),
                    #                            value = toString(
                    #                                SangerSingleReadQualReport
                    #                                [[singleReadIndex]]@
                    #                                    M2CutoffQualityScore),
                    #                            width = '90%')
                    #                    ),
                    #             ),
                    #             column(3,
                    #                    uiOutput("M2SlidingWindowSize") ,
                    #                    tags$ul(
                    #                        textInput(
                    #                            "M2SlidingWindowSizeText",
                    #                            label = p("Change Value"),
                    #                            value = toString(
                    #                                SangerSingleReadQualReport
                    #                                [[singleReadIndex]]@
                    #                                    M2SlidingWindowSize),
                    #                            width = '90%')
                    #                    ),
                    #             ),
                    #             column(3,
                    #                    uiOutput("trimmedStartPos") ,
                    #             ),
                    #             column(3,
                    #                    uiOutput("trimmedFinishPos") ,
                    #             )
                    #         ),
                    #         tags$hr(
                    #             style = ("border-top: 6px double #A9A9A9;")),
                    #         fluidRow(
                    #             column(6,
                    #                    uiOutput("trimmedRatio")
                    #             ),
                    #             column(6,
                    #                    uiOutput("remainingBP")
                    #             )
                    #         ),
                    #         box(title = tags$p("Cumulative Ratio Plot",
                    #                            style = "font-size: 24px;
                    #                    font-weight: bold;"),
                    #             collapsible = TRUE,
                    #             status = "success", width = 6,
                    #             plotlyOutput("qualityTrimmingRatioPlot") %>%
                    #                 withSpinner()),
                    #         box(title = tags$p("Cumulative Ratio Plot",
                    #                            style = "font-size: 24px;
                    #                    font-weight: bold;"),
                    #             collapsible = TRUE,
                    #             status = "success", width = 6,
                    #             plotlyOutput("qualityQualityBasePlot") %>%
                    #                 withSpinner()),
                    #     ),
                    #     box(title = tags$p("Chromatogram: ",
                    #                        style = "font-size: 26px;
                    #                    font-weight: bold;"),
                    #         solidHeader = TRUE, collapsible = TRUE,
                    #         status = "success", width = 12,
                    #         tags$hr(
                    #             style = ("border-top: 4px hidden #A9A9A9;")),
                    #         column(12,
                    #                column(3,
                    #                       sliderInput(
                    #                           "ChromatogramBasePerRow",
                    #                           label = h3("Slider"), min = 5,
                    #                           max = 200, value = 100),
                    #                ),
                    #                column(3,
                    #                       uiOutput(
                    #                           "ChromatogramtrimmedStartPos"),
                    #                ),
                    #                column(3,
                    #                       uiOutput(
                    #                           "ChromatogramtrimmedFinishPos"),
                    #                ),
                    #                column(3,
                    #                       numericInput(
                    #                           "ChromatogramSignalRatioCutoff",
                    #                           h3("Signal Ratio Cutoff"),
                    #                           value = 0.33),
                    #                       checkboxInput(
                    #                           "ChromatogramCheckShowTrimmed",
                    #                           "Whether show trimmed region",
                    #                           value = TRUE),)
                    #         ),
                    #         tags$hr(
                    #             style = ("border-top: 6px double #A9A9A9;")),
                    #         column(width = 12,
                    #                tags$hr(
                    #                    style=("border-top:6px double #A9A9A9;")
                    #                ),
                    #                uiOutput("chromatogramUIOutput"),
                    #         )
                    #     )
                    # )
                }
            }
        }
    })






    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        # message("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
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
                message("save SangerConsensus quality S4 object Forward")
                }
            )
            sapply(1:reverseReadNum, function(j) {
                SangerConsensusSet@consensusReadsList[[i]]@
                    reverseReadsList[[j]]@QualityReport <<-
                    SangerCSetParam[[i]]$
                    SangerSingleReadQualReport[[forwardReadNum + j]]
                message("save SangerConsensus quality S4 object Reverse")
            }
            )
        })
        saveRDS(SangerConsensusSet, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        NEW_SANGER_CONSENSUS_READ <<- readRDS(file=newS4Object)
        # shinyOptions(NewSangerConsensusSet = newS4)
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        btn <- input$closeUI
        stopApp()
    })


    output$geneticCodeDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCGeneticCode <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCGeneticCode
        excelTable(data =  t(data.frame(SCGeneticCode)),
                   defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })

    output$differencesDF = renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDifferencesDF <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDifferencesDF
        SCDifferencesDF
    })

    output$dendrogramDF <- renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDendrogram <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDendrogram
        SCDendrogram[[1]]
    })

    output$dendrogramPlot <- renderPlot({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDendrogram <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDendrogram
        plot(SCDendrogram[[2]])
    })

    output$SCIndelsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCIndelsDF <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCIndelsDF
        if (all(dim(SCIndelsDF) == c(0,0))) {
            h4("*** 'Indels' data frame is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCIndelsDF")
        }
    })

    output$SCStopCodonsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCStopCodonsDF <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCStopCodonsDF
        if (all(dim(SCStopCodonsDF) == c(0,0))) {
            h4("*** 'Stop Codons' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCStopCodonsDF")
        }
    })

    output$secondaryPeakDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDistanceMatrix <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDistanceMatrix
        if (all(dim(SCDistanceMatrix) == c(0,0))) {
            h4("*** 'Distance Matrix' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("secondaryPeakDF")
        }
    })

    output$SCIndelsDF <- renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCIndelsDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCIndelsDF
        SCIndelsDF
    })

    output$SCStopCodonsDF <- renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCStopCodonsDF <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCStopCodonsDF
        SCStopCodonsDF
    })

    output$secondaryPeakDF = renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDistanceMatrix <-
            SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDistanceMatrix
        SCDistanceMatrix
    })

    output$primarySeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        excelTable(data =SangerCSetParam[[consensusReadIndex]]$
                       SangerSingleReadPrimSeqDF[[singleReadIndex]],
                   defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)

    })

    output$secondSeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        excelTable(data =
                       SangerCSetParam[[consensusReadIndex]]$
                       SangerSingleReadSecoSeqDF[[singleReadIndex]],
                   defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })





    ##### Need Fix Here
    # valueBoxM1TrimmingCutoff(input, output, session, SangerSingleReadQualReport)
    output$M1TrimmingCutoff <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        message("input$M1TrimmingCutoffText: ", input$M1TrimmingCutoffText)
        if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
            as.numeric(input$M1TrimmingCutoffText) > 0 &&
            as.numeric(input$M1TrimmingCutoffText) <= 1) {
            inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
        } else {
            inputM1TrimmingCutoffText <- 0.0001
        }
        if (SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M1") {
            # message("&&&& Dynamic M1")
            trimmingPos <-
                M1inside_calculate_trimming(
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityPhredScores,
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityBaseScore,
                    as.numeric(inputM1TrimmingCutoffText))
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

                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)

                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    rawSeqLength <<- rawSeqLength
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    rawMeanQualityScore <<- rawMeanQualityScore
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    rawMinQualityScore <<- rawMinQualityScore
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedStartPos <<- trimmedStartPos
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedFinishPos <<- trimmedFinishPos
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedSeqLength <<- trimmedSeqLength
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMeanQualityScore <<- trimmedMeanQualityScore
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMinQualityScore <<- trimmedMinQualityScore
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                    remainingRatio <<- remainingRatio

                trimmedRV[["rawSeqLength"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@remainingRatio * 100, 2)
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

    ##### Need Fix Here
    # valueBoxM2CutoffQualityScore (input, output, session, SangerSingleReadQualReport)
    output$M2CutoffQualityScore <- renderUI({
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
        if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
            strtoi(input$M2SlidingWindowSizeText) > 0 &&
            strtoi(input$M2SlidingWindowSizeText) <= 20 &&
            strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
            inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
        } else {
            inputM2SlidingWindowSizeText <- 5
        }



        if (SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M2") {
            # message("&&&& Dynamic M2")
            if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
                strtoi(input$M2CutoffQualityScoreText) > 0 &&
                strtoi(input$M2CutoffQualityScoreText) <= 60 &&
                strtoi(input$M2CutoffQualityScoreText) %% 1 ==0 &&
                !is.na(strtoi(input$M2SlidingWindowSizeText)) &&
                strtoi(input$M2SlidingWindowSizeText) > 0 &&
                strtoi(input$M2SlidingWindowSizeText) <= 20 &&
                strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
                trimmingPos <- M2inside_calculate_trimming(
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                            qualityPhredScores,
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityBaseScore,
                        strtoi(inputM2CutoffQualityScoreText),
                        strtoi(inputM2SlidingWindowSizeText))
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

                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        M2CutoffQualityScore <<- as.numeric(inputM2CutoffQualityScoreText)
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)

                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        rawSeqLength <<- rawSeqLength
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        rawMeanQualityScore <<- rawMeanQualityScore
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        rawMinQualityScore <<- rawMinQualityScore
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        trimmedStartPos <<- trimmedStartPos
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        trimmedFinishPos <<- trimmedFinishPos
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        trimmedSeqLength <<- trimmedSeqLength
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        trimmedMeanQualityScore <<- trimmedMeanQualityScore
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        trimmedMinQualityScore <<- trimmedMinQualityScore
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                        remainingRatio <<- remainingRatio

                    trimmedRV[["rawSeqLength"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
                    trimmedRV[["rawMeanQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@rawMeanQualityScore
                    trimmedRV[["rawMinQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@rawMinQualityScore
                    trimmedRV[["trimmedStartPos"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedStartPos
                    trimmedRV[["trimmedFinishPos"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedFinishPos
                    trimmedRV[["trimmedSeqLength"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedSeqLength
                    trimmedRV[["trimmedMeanQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedMeanQualityScore
                    trimmedRV[["trimmedMinQualityScore"]] <<-
                        SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmedMinQualityScore
                    trimmedRV[["remainingRatio"]] <<-
                        round(SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@remainingRatio * 100, 2)
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
            inputM1TrimmingCutoffText <- 0.0003
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

                SangerCSetParam[[consensusReadIndex]]$
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    M1TrimmingCutoff <<- strtoi(inputM1TrimmingCutoffText)

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
        }
    })


    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
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
    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
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
                                 round((trimmedFinishPos - trimmedStartPos+1) /
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
            )
    })

    valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)

    # chromatogram
    output$chromatogramUIOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(as.numeric(sidebar_menu[[1]])) &&
            sidebar_menu[[2]] == "Consensus" &&
            sidebar_menu[[3]] == "Read" &&
            sidebar_menu[[4]] == "-" &&
            !is.na(as.numeric(sidebar_menu[[5]])) &&
            (sidebar_menu[[6]] == "Forward" ||
             sidebar_menu[[6]] == "Reverse") &&
            sidebar_menu[[7]] == "Read") {
            chromatogramRowNumAns <-
                chromatogramRowNum(
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerConsensusFRReadsList[[singleReadIndex]],
                    strtoi(input$ChromatogramBasePerRow)) * 200
            message("chromatogramRowNumAns: ", chromatogramRowNumAns)
            plotOutput("chromatogram", height = chromatogramRowNumAns) %>%
                withSpinner()
        }
    })

    output$chromatogram <- renderPlot({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(as.numeric(sidebar_menu[[1]])) &&
            sidebar_menu[[2]] == "Consensus" &&
            sidebar_menu[[3]] == "Read" &&
            sidebar_menu[[4]] == "-" &&
            !is.na(as.numeric(sidebar_menu[[5]])) &&
            (sidebar_menu[[6]] == "Forward" || sidebar_menu[[6]] == "Reverse") &&
            sidebar_menu[[7]] == "Read") {
            rawSeqLength =
                SangerCSetParam[[consensusReadIndex]]$
                SangerSingleReadQualReport[[singleReadIndex]]@rawSeqLength
            hetcalls <-
                makeBaseCalls(
                    SangerCSetParam[[consensusReadIndex]]$
                        SangerConsensusFRReadsList[[singleReadIndex]],
                    ratio = as.numeric(
                        input$ChromatogramSignalRatioCutoff))
            chromatogram(hetcalls,
                         width = strtoi(input$ChromatogramBasePerRow),
                         height = 2, trim5 = trimmedRV[["trimmedStartPos"]],
                         trim3 = rawSeqLength - trimmedRV[["trimmedFinishPos"]],
                         showtrim = (input$ChromatogramCheckShowTrimmed),
                         showcalls = "both")
        }
    })


    output$consensusAlignmentHTML<-renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        SCconsenesusReadName <-
            SangerCSetParam[[consensusReadIndex]]$SCconsenesusReadName
        SCAlignment <-
            SangerCSetParam[[consensusReadIndex]]$SCAlignment
        browseSeqHTML <-
            file.path(shinyDirectory, "BrowseSeqs_html",
                      paste0(sidebar_menu[[1]], "_",
                             SCconsenesusReadName,
                             "_Alignment_BrowseSeqs.html"))
        if (!dir.exists(file.path(shinyDirectory, "BrowseSeqs_html"))) {
            dir.create(file.path(shinyDirectory, "BrowseSeqs_html"))
        }
        if (!file.exists(browseSeqHTML)) {
            BrowseSeqs(SCAlignment,
                       openURL=FALSE, htmlFile=browseSeqHTML)
        }
        includeHTML(browseSeqHTML)
    })





    output$TrimmingMethodUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.null(SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]])) {
            if (input$TrimmingMethodSelection == "M1") {
                trimmingMethodLocal ="Method 1"
                message("Inside Method 1!!")
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod <<- "M1"
                if (is.null(SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff)) {
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff <<-  0.0001
                }


                fluidRow(
                    column(6,
                           uiOutput("M1TrimmingCutoff") ,
                           tags$ul(
                               textInput("M1TrimmingCutoffText",
                                         label = p("Change Value"),
                                         value = toString(
                                             SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                                                 M1TrimmingCutoff),
                                         width = '70%')
                           ),
                    ),
                )






            } else if (input$TrimmingMethodSelection == "M2") {
                trimmingMethodLocal ="Method 2"
                message("Inside Method 2!!")
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod <<- "M2"
                if (is.null(SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore)) {
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore <<-  20
                }
                if (is.null(SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize )) {
                    SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize <<-  5
                }

                fluidRow(
                    column(6,
                           uiOutput("M2CutoffQualityScore") ,
                           tags$ul(
                               textInput("M2CutoffQualityScoreText",
                                         label = p("Change Value"),
                                         value = toString(
                                             SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                                                 M2CutoffQualityScore),
                                         width = '70%')
                           ),
                    ),
                    column(6,
                           uiOutput("M2SlidingWindowSize") ,
                           tags$ul(
                               textInput("M2SlidingWindowSizeText",
                                         label = p("Change Value"),
                                         value = toString(
                                             SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                                                 M2SlidingWindowSize),
                                         width = '70%')
                           ),
                    ),
                )


            }
        }
    })
    output$TrimmingMethodSelectionOutput <- renderText({
        # tags$p(paste("You currently select '", input$TrimmingMethodSelection, "'"),
        #        style = "font-size: 15px; font-weight: bold;")
        if (input$TrimmingMethodSelection == "M1") {
            "Logarithmic Scale Trimming"
        } else if (input$TrimmingMethodSelection == "M2") {
            "Logarithmic Scale Sliding Window Trimming"
        }
    })
}

