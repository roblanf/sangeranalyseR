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
        SangerSingleReadPrimSeqDF <- c(forwardReadPrimSeqDF, reverseReadPrimSeqDF)


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


    trimmedRV <- reactiveValues(trimmedStart = 0, trimmedEnd = 0,
                                remainingBP = 0, trimmedRatio = 0)

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

                    # BrowseSeqs_html
                    if (!dir.exists(file.path(shinyDirectory, "BrowseSeqs_html"))) {
                        dir.create(file.path(shinyDirectory, "BrowseSeqs_html"))
                    }
                    browseSeqHTML <-
                        file.path(shinyDirectory, "BrowseSeqs_html",
                                  paste0(sidebar_menu[[1]], "_",
                                         SCconsenesusReadName,
                                         "_Alignment_BrowseSeqs.html"))
                    if (!file.exists(browseSeqHTML)) {
                        BrowseSeqs(SCAlignment,
                                   openURL=FALSE, htmlFile=browseSeqHTML)
                    }

                    fluidRow(
                        useShinyjs(),
                        # h1(SCconsensusRead),
                        h1(paste("You've selected:", input$sidebar_menu)),
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
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Alignment",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           includeHTML(browseSeqHTML)
                                    ),
                                ),
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
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
                           (sidebar_menu[[6]] == "Forward" || sidebar_menu[[6]] == "Reverse") &&
                           sidebar_menu[[7]] == "Read") {
                    consensusReadIndex <- strtoi(sidebar_menu[[1]])
                    singleReadIndex <- strtoi(sidebar_menu[[5]])
                    SangerConsensusFRReadsList <- SangerCSetParam[[consensusReadIndex]]$SangerConsensusFRReadsList
                    SangerSingleReadBFN <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadBFN
                    SangerSingleReadAFN <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadAFN
                    SangerSingleReadFeature <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadFeature
                    SangerSingleReadAbifRawData <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadAbifRawData
                    SangerSingleReadQualReport <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport
                    SangerSingleReadPrimSeqID <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadPrimSeqID
                    SangerSingleReadPrimSeq <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadPrimSeq
                    forwardReadPrimSeqDF <- SangerCSetParam[[consensusReadIndex]]$forwardReadPrimSeqDF
                    reverseReadPrimSeqDF <- SangerCSetParam[[consensusReadIndex]]$reverseReadPrimSeqDF
                    SangerSingleReadPrimSeqDF <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadPrimSeqDF
                    SangerSingleReadSecoSeqID <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadSecoSeqID
                    SangerSingleReadSecoSeq <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadSecoSeq
                    forwardReadSecoSeqDF <- SangerCSetParam[[consensusReadIndex]]$forwardReadPrimSeqDF
                    reverseReadSecoSeqDF <- SangerCSetParam[[consensusReadIndex]]$reverseReadSecoSeqDF
                    SangerSingleReadSecoSeqDF <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadSecoSeqDF
                    SangerSingleReadTraceMat <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadTraceMat
                    SangerSingleReadPeakPosMat <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadPeakPosMat
                    SangerSingleReadPeakAmpMat <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadPeakAmpMat


                    fluidRow(
                        useShinyjs(),
                        h1("consensusReadIndex: ", consensusReadIndex, "  singleReadIndex: ", singleReadIndex),
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
                            fluidRow(
                                column(3,
                                       uiOutput("cutoffQualityScore") ,
                                       tags$ul(
                                           textInput("cutoffQualityScoreText",
                                                     label = p("Change Value"),
                                                     value = toString(
                                                         SangerSingleReadQualReport
                                                         [[singleReadIndex]]@
                                                             cutoffQualityScore),
                                                     width = '90%')
                                       ),
                                ),
                                column(3,
                                       uiOutput("slidingWindowSize") ,
                                       tags$ul(
                                           textInput("slidingWindowSizeText",
                                                     label = p("Change Value"),
                                                     value = toString(
                                                         SangerSingleReadQualReport
                                                         [[singleReadIndex]]@
                                                             slidingWindowSize),
                                                     width = '90%')
                                       ),
                                ),
                                column(3,
                                       uiOutput("trimmingStartPos") ,
                                ),
                                column(3,
                                       uiOutput("trimmingFinishPos") ,
                                )
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                column(6,
                                       uiOutput("trimmedRatio")
                                ),
                                column(6,
                                       uiOutput("remainingBP")
                                )
                            ),
                            box(title = tags$p("Cumulative Ratio Plot",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 6,
                                plotlyOutput("qualityTrimmingRatioPlot") %>%
                                    withSpinner()),
                            box(title = tags$p("Cumulative Ratio Plot",
                                               style = "font-size: 24px;
                                       font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 6,
                                plotlyOutput("qualityQualityBasePlot") %>%
                                    withSpinner()),
                        ),
                        box(title = tags$p("Chromatogram: ",
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                            column(12,
                                   column(3,
                                          sliderInput("ChromatogramBasePerRow",
                                                      label = h3("Slider"), min = 5,
                                                      max = 200, value = 100),
                                   ),
                                   column(3,
                                          uiOutput("ChromatogramTrimmingStartPos"),
                                   ),
                                   column(3,
                                          uiOutput("ChromatogramTrimmingFinishPos"),
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
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                            column(width = 12,
                                   tags$hr(
                                       style = ("border-top: 6px double #A9A9A9;")
                                   ),
                                   uiOutput("chromatogramUIOutput"),
                            )
                        )
                    )
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
            #     trimmedRV[["trimmedStart"]] <-
            #         SangerSingleReadQualReport[[
            #             strtoi(sidebar_menu[[1]])]]@trimmingStartPos
            #     trimmedRV[["trimmedEnd"]] <-
            #         SangerSingleReadQualReport[[
            #             strtoi(sidebar_menu[[1]])]]@trimmingFinishPos
            #     qualityPhredScores = SangerSingleReadQualReport[[
            #         strtoi(sidebar_menu[[1]])]]@qualityPhredScores
            #
            #     readLen = length(qualityPhredScores)
            #     trimmedRV[["remainingBP"]] <- trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1
            #     trimmedRV[["trimmedRatio"]] <- round(((trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1) / readLen) * 100, digits = 2)
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        btn <- input$saveS4
        id <- paste0('txt', btn)
        newS4Object <- file.path(tempdir(), "SangerAlignedConsensusSet.Rda")
        saveRDS(SangerCSetParam, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        showNotification(paste("New S4 object is store as:", newS4Object),
                         type = "message", duration = 10)
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
        SCGeneticCode <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCGeneticCode
        excelTable(data =  t(data.frame(SCGeneticCode)),
                   defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })

    output$differencesDF = renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDifferencesDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDifferencesDF
        SCDifferencesDF
    })

    output$dendrogramDF <- renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDendrogram <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDendrogram
        SCDendrogram[[1]]
    })

    output$dendrogramPlot <- renderPlot({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDendrogram <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDendrogram
        plot(SCDendrogram[[2]])
    })

    output$SCIndelsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCIndelsDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCIndelsDF
        if (all(dim(SCIndelsDF) == c(0,0))) {
            h4("*** 'Indels' data frame is empty. ***", style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCIndelsDF")
        }
    })

    output$SCStopCodonsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCStopCodonsDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCStopCodonsDF
        if (all(dim(SCStopCodonsDF) == c(0,0))) {
            h4("*** 'Stop Codons' dataframe is empty. ***", style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCStopCodonsDF")
        }
    })

    output$secondaryPeakDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDistanceMatrix <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDistanceMatrix
        if (all(dim(SCDistanceMatrix) == c(0,0))) {
            h4("*** 'Distance Matrix' dataframe is empty. ***", style="font-weight: bold; font-style: italic;")
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
        SCStopCodonsDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCStopCodonsDF
        SCStopCodonsDF
    })

    output$secondaryPeakDF = renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        SCDistanceMatrix <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDistanceMatrix
        SCDistanceMatrix
    })

    output$primarySeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        excelTable(data =
                       SangerCSetParam[[consensusReadIndex]]$SangerSingleReadPrimSeqDF[[singleReadIndex]],
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
                       SangerCSetParam[[consensusReadIndex]]$SangerSingleReadSecoSeqDF[[singleReadIndex]],
                   defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })

    valueBoxCutoffQualityScore (input, output, session)
    valueBoxSlidingWindowSize (input, output, session)
    valueBoxTrimmingStartPos (input, output, session, trimmedRV)
    valueBoxTrimmingFinishPos (input, output, session, trimmedRV)

    valueBoxRemainingBP (input, output, session, trimmedRV)
    valueBoxTrimmedRatio (input, output, session, trimmedRV)

    observeEvent(input$cutoffQualityScoreText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        message("@@@singleReadIndex: ", singleReadIndex)
        if (!is.na(strtoi(input$cutoffQualityScoreText)) &&
            strtoi(input$cutoffQualityScoreText) > 0 &&
            strtoi(input$cutoffQualityScoreText) <= 60 &&
            strtoi(input$cutoffQualityScoreText) %% 1 ==0) {
            inputCutoffQualityScoreText <- input$cutoffQualityScoreText
        } else {
            inputCutoffQualityScoreText <- 20
        }
        trimmingPos <-
            inside_calculate_trimming(
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityBaseScore, strtoi(inputCutoffQualityScoreText),
                SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@slidingWindowSize)
        if (!is.null(trimmingPos[1]) && !is.null(trimmingPos[2])) {
            SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@
                cutoffQualityScore <<- strtoi(inputCutoffQualityScoreText)
            SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingStartPos <<- trimmingPos[1]
            SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingFinishPos <<- trimmingPos[2]
            trimmedRV[["trimmedStart"]] <<- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingStartPos
            trimmedRV[["trimmedEnd"]] <<- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingFinishPos
            qualityPhredScores <-  SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
            readLen = length(qualityPhredScores)
            trimmedRV[["remainingBP"]] <<-
                trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1
            trimmedRV[["trimmedRatio"]] <<-
                round(((trimmedRV[["trimmedEnd"]] -
                            trimmedRV[["trimmedStart"]] + 1) / readLen)*100,
                      digits = 2)
        }
    })

    observeEvent(input$slidingWindowSizeText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        if (!is.na(strtoi(input$slidingWindowSizeText)) &&
            strtoi(input$slidingWindowSizeText) > 0 &&
            strtoi(input$slidingWindowSizeText) <= 20 &&
            strtoi(input$slidingWindowSizeText) %% 1 ==0) {
            inputSlidingWindowSizeText <- input$slidingWindowSizeText
        } else {
            inputSlidingWindowSizeText <- 5
        }
        trimmingPos <-
            inside_calculate_trimming(SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityBaseScore,
                                      SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@cutoffQualityScore,
                strtoi(inputSlidingWindowSizeText))
        if (!is.null(trimmingPos[1]) && !is.null(trimmingPos[2])) {
            SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@slidingWindowSize <<- strtoi(inputSlidingWindowSizeText)
            SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingStartPos <<- trimmingPos[1]
            SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingFinishPos <<- trimmingPos[2]
            trimmedRV[["trimmedStart"]] <<-SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingStartPos
            trimmedRV[["trimmedEnd"]] <<-SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@trimmingFinishPos
            qualityPhredScores = SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
            readLen = length(qualityPhredScores)
            trimmedRV[["remainingBP"]] <<-trimmedRV[["trimmedEnd"]] - trimmedRV[["trimmedStart"]] + 1
            trimmedRV[["trimmedRatio"]] <<-
                round(((trimmedRV[["trimmedEnd"]] -
                            trimmedRV[["trimmedStart"]] + 1) / readLen)*100,
                      digits = 2)
        }
    })
    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        readFeature <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadFeature[[singleReadIndex]]
        trimmingStartPos = trimmedRV[["trimmedStart"]]
        trimmingFinishPos = trimmedRV[["trimmedEnd"]]
        qualityPhredScores <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
        readLen = length(qualityPhredScores)

        stepRatio = 1 / readLen
        trimmingStartPos / readLen
        trimmingFinishPos / readLen

        trimmedPer <- c()
        remainingPer <- c()

        for (i in 1:trimmingStartPos) {
            if (i != trimmingStartPos) {
                trimmedPer <- c(trimmedPer, stepRatio)
            }
        }

        for (i in trimmingStartPos:trimmingFinishPos) {
            trimmedPer <- c(trimmedPer, 0)
        }


        for (i in trimmingFinishPos:readLen) {
            if (i != trimmingFinishPos) {
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
            layout(xaxis = x, yaxis = y, legend = list(orientation = 'h',
                                                       xanchor = "center",  # use center of legend as anchor
                                                       x = 0.5, y = 1.1)) %>%
            add_annotations(
                text = "Trimmed Ratio (Each BP)",
                x = (trimmingStartPos + trimmingFinishPos) / 2,
                y = ((trimmedPer[1] + trimmedPer[length(trimmedPer)]) / 2)
                + 0.06,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Remaining Ratio (Each BP)",
                x = (trimmingStartPos+trimmingFinishPos) / 2,
                y = ((remainingPer[1] + remainingPer[length(remainingPer)]) / 2)
                - 0.06,
                showarrow=FALSE
            )
    })
    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[5]])
        readFeature <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadFeature[[singleReadIndex]]
        trimmingStartPos = trimmedRV[["trimmedStart"]]
        trimmingFinishPos = trimmedRV[["trimmedEnd"]]
        qualityPhredScores <- SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
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
            add_trace(x=seq(trimmingStartPos,
                            trimmingFinishPos,
                            len=trimmingFinishPos-trimmingStartPos+1),
                      y=rep(70, trimmingFinishPos-trimmingStartPos+1),
                      mode="lines", hoverinfo="text",
                      text=paste("Trimmed Reads BP length:",
                                 trimmingFinishPos-trimmingStartPos+1,
                                 "BPs <br>",
                                 "Trimmed Reads BP ratio:",
                                 round((trimmingFinishPos - trimmingStartPos+1)/
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
                   shapes = list(vline(trimmingStartPos),
                                 vline(trimmingFinishPos)),
                   legend = list(orientation = 'h',
                                 xanchor = "center",  # use center of legend as anchor
                                 x = 0.5, y = 1.1)) %>%
            # add_segments(x = trimmingStartPos, xend = trimmingFinishPos, y = 70, yend = 70, inherit = TRUE, width = 10, line = list(width = 8)) %>%
            # add_segments(x = 0, xend = readLen, y = 75, yend = 75, inherit = TRUE, width = 4, line = list(width = 8)) %>%
            add_annotations(
                text = "Trimming Strat <br> BP Index",
                x = trimmingStartPos + 40,
                y = 15,
                showarrow=FALSE
            ) %>%
            add_annotations(
                text = "Trimming End <br> BP Index",
                x = trimmingFinishPos - 40,
                y = 15,
                showarrow=FALSE
            )
        # add_markers(qualityPlotDf, x=~Index, y=~Score)
        # add_segments(x = trimmingStartPos, xend = trimmingFinishPos, y = 70, yend = 70, inherit = TRUE)
    })

    valueBoxChromTrimmingStartPos (input, output, session, trimmedRV)
    valueBoxChromTrimmingFinishPos (input, output, session, trimmedRV)

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
            (sidebar_menu[[6]] == "Forward" || sidebar_menu[[6]] == "Reverse") &&
            sidebar_menu[[7]] == "Read") {
            chromatogramRowNumAns <-
                chromatogramRowNum(
                    SangerCSetParam[[consensusReadIndex]]$SangerConsensusFRReadsList
                    [[singleReadIndex]],
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
            qualityPhredScores = SangerCSetParam[[consensusReadIndex]]$SangerSingleReadQualReport[[singleReadIndex]]@qualityPhredScores
            readLen = length(qualityPhredScores)
            hetcalls <-
                makeBaseCalls(SangerCSetParam[[consensusReadIndex]]$SangerConsensusFRReadsList
                              [[singleReadIndex]],
                    ratio = as.numeric(
                        input$ChromatogramSignalRatioCutoff))
            chromatogram(hetcalls,
                         width = strtoi(input$ChromatogramBasePerRow),
                         height = 2, trim5 = trimmedRV[["trimmedStart"]],
                         trim3 = readLen - trimmedRV[["trimmedEnd"]],
                         showtrim = (input$ChromatogramCheckShowTrimmed),
                         showcalls = "both")
        }
    })
}

