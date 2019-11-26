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
        SangerSingleReadRReadsList <- SangerCSetList[[1]]@reverseReadsList
        forwardReadNum <- length(SangerSingleReadFReadsList)
        reverseReadNum <- length(SangerSingleReadRReadsList)
        SangerSingleReadNum <- forwardReadNum + reverseReadNum

        # Forward + reverse reads list
        SangerConsensusFRReadsList <- c(SangerSingleReadFReadsList,
                                        SangerSingleReadRReadsList)

        # readFileName (basename)
        forwardReadBFN <- sapply(1:forwardReadNum, function(i)
            basename(SangerSingleReadFReadsList[[i]]@readFileName))
        reverseReadBFN <- sapply(1:reverseReadNum, function(i)
            basename(SangerSingleReadRReadsList[[i]]@readFileName))
        SangerSingleReadBFN <- c(forwardReadBFN, reverseReadBFN)

        # readFileName (absolute)
        forwardReadAFN <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@readFileName)
        reverseReadAFN <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@readFileName)
        SangerSingleReadAFN <- c(forwardReadAFN, reverseReadAFN)

        # readFeature
        forwardReadFeature <- sapply(1:forwardReadNum, function(i)
            paste0(i, " ",
                   SangerSingleReadFReadsList[[i]]@readFeature))
        reverseReadFeature <- sapply(1:reverseReadNum, function(i)
            paste0(i+forwardReadNum, " ",
                   SangerSingleReadRReadsList[[i]]@readFeature))
        SangerSingleReadFeature <- c(forwardReadFeature, reverseReadFeature)

        # abifRawData
        forwardReadAbifRawData <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@abifRawData)
        reverseReadAbifRawData <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@abifRawData)
        SangerSingleReadAbifRawData <- c(forwardReadAbifRawData,
                                         reverseReadAbifRawData)

        # QualityReport
        forwardReadQualReport <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@QualityReport)
        reverseReadQualReport <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@QualityReport)
        SangerSingleReadQualReport <- c(forwardReadQualReport,
                                        reverseReadQualReport)

        # primarySeqID
        forwardReadPrimSeqID <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@primarySeqID)
        reverseReadPrimSeqID <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@primarySeqID)
        SangerSingleReadPrimSeqID <- c(forwardReadPrimSeqID, reverseReadPrimSeqID)

        # primarySeq
        forwardReadPrimSeq <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@primarySeq)
        reverseReadPrimSeq <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@primarySeq)
        SangerSingleReadPrimSeq <- c(forwardReadPrimSeq, reverseReadPrimSeq)

        # secondarySeqID
        forwardReadSecoSeqID <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@secondarySeqID)
        reverseReadSecoSeqID <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@secondarySeqID)
        SangerSingleReadSecoSeqID <- c(forwardReadSecoSeqID, reverseReadSecoSeqID)

        # secondarySeq
        forwardReadSecoSeq <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@secondarySeq)
        reverseReadSecoSeq <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@secondarySeq)
        SangerSingleReadSecoSeq <- c(forwardReadSecoSeq, reverseReadSecoSeq)

        # traceMatrix
        forwardReadTraceMat <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@traceMatrix)
        reverseReadTraceMat <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@traceMatrix)
        SangerSingleReadTraceMat <- c(forwardReadTraceMat, reverseReadTraceMat)

        # peakPosMatrix
        forwardReadReadPeakPosMat <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@peakPosMatrix)
        reverseReadReadPeakPosMat <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@peakPosMatrix)
        SangerSingleReadPeakPosMat <- c(forwardReadReadPeakPosMat,
                                        reverseReadReadPeakPosMat)
        # peakAmpMatrix
        forwardReadPeakAmpMat <- sapply(1:forwardReadNum, function(i)
            SangerSingleReadFReadsList[[i]]@peakAmpMatrix)
        reverseReadPeakAmpMat <- sapply(1:reverseReadNum, function(i)
            SangerSingleReadRReadsList[[i]]@peakAmpMatrix)
        SangerSingleReadPeakAmpMat <- c(forwardReadPeakAmpMat, reverseReadPeakAmpMat)
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
                    SangerSingleReadSecoSeqID = SangerSingleReadSecoSeqID,
                    SangerSingleReadSecoSeq = SangerSingleReadSecoSeq,
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
                                  h4("Output Directory: ", style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(shinyDirectory),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Parent Directory: ", style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(SangerCSetParentDir),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Forward Suffix RegExp: ", style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(SangerCSetSuffixForwardRegExp),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Reverse Suffix RegExp: ", style="font-weight: bold;"),
                           ),
                           column(9,
                                  h4(SangerCSetSuffixReverseRegExp),
                           )
                    ),
                    column(12,
                           column(3,
                                  h4("Consensus Read Number: ", style="font-weight: bold;"),
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
                    consensusReadIndex <- sidebar_menu[[1]]
                    h1(paste("You've selected:", input$sidebar_menu))

                    SCName <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCName
                    SCMinReadsNum <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCMinReadsNum
                    SCMinReadLength <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCMinReadLength
                    SCRefAminoAcidSeq <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCRefAminoAcidSeq
                    SCMinFractionCall <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCMinFractionCall
                    SCMaxFractionLost <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCMaxFractionLost
                    SCGeneticCode <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCGeneticCode
                    SCAcceptStopCodons <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCAcceptStopCodons
                    SCReadingFrame <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCReadingFrame
                    SCConsensusRead <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCConsensusRead
                    SCAlignment <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCAlignment
                    SCDifferencesDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDifferencesDF
                    SCDistanceMatrix <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDistanceMatrix
                    SCDendrogram <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCDendrogram
                    SCIndelsDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCIndelsDF
                    SCStopCodonsDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCStopCodonsDF
                    SCSecondaryPeakDF <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCSecondaryPeakDF
                    SCconsenesusReadName <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SCconsenesusReadName
                    SangerConsensusForRegExp <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerConsensusForRegExp
                    SangerConsensusRevRegExp <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerConsensusRevRegExp
                    SangerSingleReadNum <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadNum
                    SangerConsensusFRReadsList <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerConsensusFRReadsList
                    SangerSingleReadBFN <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadBFN
                    SangerSingleReadAFN <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadAFN
                    SangerSingleReadFeature <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadFeature
                    SangerSingleReadAbifRawData <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadAbifRawData
                    SangerSingleReadQualReport <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadQualReport
                    SangerSingleReadPrimSeqID <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadPrimSeqID
                    SangerSingleReadPrimSeq <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadPrimSeq
                    SangerSingleReadSecoSeqID <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadSecoSeqID
                    SangerSingleReadSecoSeq <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadSecoSeq
                    SangerSingleReadTraceMat <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadTraceMat
                    SangerSingleReadPeakPosMat <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadPeakPosMat
                    SangerSingleReadPeakAmpMat <- SangerCSetParam[[strtoi(sidebar_menu[[1]])]]$SangerSingleReadPeakAmpMat

                    # BrowseSeqs_html
                    if (!dir.exists(file.path(shinyDirectory, "BrowseSeqs_html"))) {
                        dir.create(file.path(shinyDirectory, "BrowseSeqs_html"))
                    }
                    browseSeqHTML <- file.path(shinyDirectory, "BrowseSeqs_html", paste0(sidebar_menu[[1]], "_", SCconsenesusReadName, "_Alignment_BrowseSeqs.html"))
                    if (!file.exists(browseSeqHTML)) {
                        BrowseSeqs(SCAlignment, openURL=FALSE, htmlFile=browseSeqHTML)
                    }

                    fluidRow(
                        useShinyjs(),
                        # h1(SCconsensusRead),
                        box(title = tags$p("Input Parameters: ",
                                           style = "font-size: 26px;
                                       font-weight: bold;"),
                            solidHeader = TRUE, collapsible = TRUE,
                            status = "success", width = 12,
                            tags$hr(style = ("border-top: 0.2px hidden #A9A9A9;")),
                            fluidRow(
                                column(12,
                                       column(3,
                                              h4("Consenesus Read Name: ", style="font-weight: bold;"),
                                       ),
                                       column(9,
                                              h4(SCName),
                                       )
                                ),
                            ),
                            ################################################
                            #### Add this after having reference sample ####
                            ################################################
                            # If it is null
                            # h1(SCRefAminoAcidSeq),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
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
                                           excelOutput("geneticCodeDF", width = "100%", height = "50"),
                                           style = "height:60px; overflow-y: hidden;overflow-x: scroll;"
                                    ),
                                ),
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                column(3,
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
                                       # uiOutput("SCMinReadsNum") ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle = tags$p("MinReadLength",
                                                             style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value = tags$p(strtoi(SCMinReadLength),
                                                          style = "font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                       # uiOutput("SCMinReadLength")  ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle = tags$p("MinFractionCall",
                                                             style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value = tags$p(as.numeric(SCMinFractionCall),
                                                          style = "font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                       # uiOutput("SCMinFractionCall") ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle = tags$p("MaxFractionLost",
                                                             style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value = tags$p(as.numeric(SCMaxFractionLost),
                                                          style = "font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                       # uiOutput("SCMaxFractionLost") ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle = tags$p("AcceptStopCodons",
                                                             style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value = tags$p(SCAcceptStopCodons,
                                                          style = "font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                       # uiOutput("SCAcceptStopCodons") ,
                                ),
                                column(3,
                                       valueBox(
                                           subtitle = tags$p("ReadingFrame",
                                                             style = "font-size: 15px;
                                            font-weight: bold;"),
                                           value = tags$p(strtoi(SCReadingFrame),
                                                          style = "font-size: 29px;"),
                                           icon = icon("cut", "fa-sm"),
                                           color = "olive",
                                           width = 12,
                                       ),
                                       # uiOutput("SCReadingFrame") ,
                                ),
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
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
                                           # htmlOutput("consensusAlignmentHTML"),
                                           # sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
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
                                           dataTableOutput("differencesDF") , style = "height:100%; overflow-y: scroll;overflow-x: scroll;"
                                    )
                                ),
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Dendrogram",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           # plot()
                                           plotOutput("dendrogramPlot"),
                                           style = "height:100%; overflow-y: scroll;overflow-x: scroll;"
                                    ),
                                    column(width = 12,
                                           dataTableOutput("dendrogramDF"),
                                           style = "height:100%; overflow-y: scroll;overflow-x: scroll;"
                                    )
                                ),
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Secondary Peak Dataframe",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           uiOutput("secondaryPeakDFUI"),
                                           style = "height:100%; overflow-y: scroll;overflow-x: scroll;"
                                    )
                                ),
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Indels Dataframe",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           uiOutput("SCIndelsDFUI"),
                                           style = "height:100%; overflow-y: scroll;overflow-x: scroll;"
                                    )
                                ),
                            ),
                            tags$hr(style = ("border-top: 6px double #A9A9A9;")),
                            fluidRow(
                                box(title = tags$p("Stop Codons Dataframe",
                                                   style = "font-size: 24px;
                                       font-weight: bold;"),
                                    collapsible = TRUE,
                                    status = "success", width = 12,
                                    column(width = 12,
                                           uiOutput("SCStopCodonsDFUI"),
                                           style = "height:100%; overflow-y: scroll;overflow-x: scroll;"
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
                    consensusReadIndex <- sidebar_menu[[1]]
                    singleReadIndex <- sidebar_menu[[5]]

                    h1(paste("You've selected:", input$sidebar_menu))
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
    observeEventButtonSaveSCSet(input, output, session, SangerCSetParam)
    observeEventButtonClose(input, output, session)


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

    # output$consensusAlignmentHTML<-renderUI({
    #     sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
    #     if (!dir.exists(file.path(shinyDirectory, "BrowseSeqs_html"))) {
    #         dir.create(file.path(shinyDirectory, "BrowseSeqs_html"))
    #     }
    #     localAlignment <- SangerCSetParam[[paste0(sidebar_menu[[1]]$SCAlignment
    #     browseSeqHTML <- file.path(shinyDirectory, "BrowseSeqs_html", paste0(sidebar_menu[[1]], "_", localAlignment, "_Alignment_BrowseSeqs.html"))
    #     if (!file.exists(browseSeqHTML)) {
    #         BrowseSeqs(localAlignment, openURL=FALSE, htmlFile=browseSeqHTML)
    #     }
    #     includeHTML(browseSeqHTML)
    # })
}

