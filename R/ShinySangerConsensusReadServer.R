### ============================================================================
### R shiny consensusRead server function
### ============================================================================
consensusReadServer <- function(input, output, session) {
    # Suppress Warning
    options(warn = -1)

    ### ------------------------------------------------------------------------
    ### SangerConsensusRead parameters initialization.
    ### ------------------------------------------------------------------------
    SangerConsensusRead <- getShinyOption("SangerConsensusRead")
    shinyDirectory <- getShinyOption("shinyDirectory")
    SangerConsensus <- SangerConsensusRead[[1]]

    ### ------------------------------------------------------------------------
    ### ConsensusRead-related parameters initialization.
    ### ------------------------------------------------------------------------
    parentDirectory <- SangerConsensus@parentDirectory
    consenesusReadName <- SangerConsensus@consenesusReadName
    suffixForwardRegExp <- SangerConsensus@suffixForwardRegExp
    suffixReverseRegExp <- SangerConsensus@suffixReverseRegExp
    SCMinReadsNum <- SangerConsensus@minReadsNum
    SCMinReadLength <- SangerConsensus@minReadLength
    SCRefAminoAcidSeq <- SangerConsensus@refAminoAcidSeq
    SCMinFractionCall <- SangerConsensus@minFractionCall
    SCMaxFractionLost <- SangerConsensus@maxFractionLost
    SCGeneticCode <- SangerConsensus@geneticCode
    SCAcceptStopCodons <- SangerConsensus@acceptStopCodons
    SCReadingFrame <- SangerConsensus@readingFrame
    SCAlignment<- SangerConsensus@alignment
    SCDifferencesDF<- SangerConsensus@differencesDF
    SCDistanceMatrix <- SangerConsensus@distanceMatrix
    SCDendrogram <- SangerConsensus@dendrogram
    SCIndelsDF <- SangerConsensus@indelsDF
    SCStopCodonsDF <- SangerConsensus@stopCodonsDF
    SCSecondaryPeakDF <- SangerConsensus@secondaryPeakDF
    SCSecondaryPeakDF <- SangerConsensus@secondaryPeakDF
    SCTrimmingMethod <- SangerConsensus@forwardReadsList[[1]]@
        QualityReport@TrimmingMethod
    if (SCTrimmingMethod == "M1") {
        SCTrimmingMethodName = "Method 1: 'Logarithmic Scale Trimming'"
    } else if (SCTrimmingMethod == "M2") {
        SCTrimmingMethodName = "Method 2: 'Logarithmic Scale Sliding Window Trimming'"
    }

    ### ------------------------------------------------------------------------
    ### Reads-related parameters initialization.
    ### ------------------------------------------------------------------------
    forwardReadNum <- length((SangerConsensus)@forwardReadsList)
    reverseReadNum <- length((SangerConsensus)@reverseReadsList)
    SangerSingleReadNum <- forwardReadNum + reverseReadNum

    # Forward & reverse reads list
    SangerConsensusFRReadsList <- c(SangerConsensus@forwardReadsList,
                                    SangerConsensus@reverseReadsList)

    # readFileName (basename)
    forwardReadBFN <- sapply(1:forwardReadNum, function(i)
        basename(SangerConsensus@forwardReadsList[[i]]@readFileName))
    reverseReadBFN <- sapply(1:reverseReadNum, function(i)
        basename(SangerConsensus@reverseReadsList[[i]]@readFileName))
    SangerSingleReadBFN <- c(forwardReadBFN, reverseReadBFN)

    # readFileName (absolute)
    forwardReadAFN <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@readFileName)
    reverseReadAFN <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@readFileName)
    SangerSingleReadAFN <- c(forwardReadAFN, reverseReadAFN)

    # readFeature
    forwardReadFeature <- sapply(1:forwardReadNum, function(i)
        paste0(i, " ",
               SangerConsensus@forwardReadsList[[i]]@readFeature))
    reverseReadFeature <- sapply(1:reverseReadNum, function(i)
        paste0(i+forwardReadNum, " ",
               SangerConsensus@reverseReadsList[[i]]@readFeature))
    SangerSingleReadFeature <- c(forwardReadFeature, reverseReadFeature)

    # abifRawData
    forwardReadAbifRawData <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@abifRawData)
    reverseReadAbifRawData <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@abifRawData)
    SangerSingleReadAbifRawData <- c(forwardReadAbifRawData,
                                     reverseReadAbifRawData)

    # QualityReport
    forwardReadQualReport <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@QualityReport)
    reverseReadQualReport <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@QualityReport)
    SangerSingleReadQualReport <- c(forwardReadQualReport,
                                    reverseReadQualReport)

    # Quality Score Dataframe
    forwardQualityScoreDF <- lapply(1:forwardReadNum, function(i) {
        PhredScore <- forwardReadQualReport[[i]]@qualityPhredScores
        PhredScoreDF <- data.frame(
            t(data.frame(PhredScore)), stringsAsFactors = FALSE)
        colnames(PhredScoreDF) <- substr(colnames(PhredScoreDF), 2, 100)
        rownames(PhredScoreDF) <- NULL
        return(PhredScoreDF)
        }
    )

    reverseQualityScoreDF <- lapply(1:reverseReadNum, function(i) {
        PhredScore <- reverseReadQualReport[[i]]@qualityPhredScores
        PhredScoreDF <- data.frame(
            t(data.frame(PhredScore)), stringsAsFactors = FALSE)
        colnames(PhredScoreDF) <- substr(colnames(PhredScoreDF), 2, 100)
        rownames(PhredScoreDF) <- NULL
        return(PhredScoreDF)
        }
    )
    SangerSingleReadQSDF <- c(forwardQualityScoreDF, reverseQualityScoreDF)

    # primarySeqID
    forwardReadPrimSeqID <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@primarySeqID)
    reverseReadPrimSeqID <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@primarySeqID)
    SangerSingleReadPrimSeqID <- c(forwardReadPrimSeqID, reverseReadPrimSeqID)

    # primarySeq
    forwardReadPrimSeq <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@primarySeq)
    reverseReadPrimSeq <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@primarySeq)
    SangerSingleReadPrimSeq <- c(forwardReadPrimSeq, reverseReadPrimSeq)
    # SangerSingleReadPrimSeqChara <- as.character(SangerSingleReadPrimSeq)

    # primarySeqDF
    forwardReadPrimSeqDF <- lapply(1:forwardReadNum, function(i) {
        basecalls1 <- unlist(strsplit(
            toString(SangerConsensus@forwardReadsList[[i]]@primarySeq), ""))
        aveposition <- rowMeans(
            SangerConsensus@forwardReadsList[[i]]@peakPosMatrix, na.rm=TRUE)
        basecalls1 <- basecalls1[1:length(aveposition)]
        basecalls1DF <- data.frame(
            t(data.frame(basecalls1)), stringsAsFactors = FALSE)
        colnames(basecalls1DF) <- substr(colnames(basecalls1DF), 2, 100)
        rownames(basecalls1DF) <- NULL
        return(basecalls1DF)
        }
    )
    reverseReadPrimSeqDF <- lapply(1:reverseReadNum, function(i) {
        basecalls1 <- unlist(strsplit(
            toString(SangerConsensus@reverseReadsList[[i]]@primarySeq), ""))
        aveposition <- rowMeans(
            SangerConsensus@reverseReadsList[[i]]@peakPosMatrix, na.rm=TRUE)
        basecalls1 <- basecalls1[1:length(aveposition)]
        basecalls1DF <- data.frame(
            t(data.frame(basecalls1)), stringsAsFactors = FALSE)
        colnames(basecalls1DF) <- substr(colnames(basecalls1DF), 2, 100)
        rownames(basecalls1DF) <- NULL
        return(basecalls1DF)
        }
    )
    SangerSingleReadPrimSeqDF <- c(forwardReadPrimSeqDF, reverseReadPrimSeqDF)

    # secondarySeqID
    forwardReadSecoSeqID <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@secondarySeqID)
    reverseReadSecoSeqID <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@secondarySeqID)
    SangerSingleReadSecoSeqID <- c(forwardReadSecoSeqID, reverseReadSecoSeqID)

    # secondarySeq
    forwardReadSecoSeq <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@secondarySeq)
    reverseReadSecoSeq <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@secondarySeq)
    SangerSingleReadSecoSeq <- c(forwardReadSecoSeq, reverseReadSecoSeq)

    # secondarySeqDF
    forwardReadSecoSeqDF <- lapply(1:forwardReadNum, function(i) {
        basecalls2 <- unlist(strsplit(
            toString(SangerConsensus@forwardReadsList[[i]]@secondarySeq), ""))
        aveposition <- rowMeans(
            SangerConsensus@forwardReadsList[[i]]@peakPosMatrix, na.rm=TRUE)
        basecalls2 <- basecalls2[1:length(aveposition)]
        basecalls2DF <- data.frame(
            t(data.frame(basecalls2)), stringsAsFactors = FALSE)
        colnames(basecalls2DF) <- substr(colnames(basecalls2DF), 2, 100)
        rownames(basecalls2DF) <- NULL
        return(basecalls2DF)
        }
    )
    reverseReadSecoSeqDF <- lapply(1:reverseReadNum, function(i) {
        basecalls2 <- unlist(strsplit(toString(
            SangerConsensus@reverseReadsList[[i]]@secondarySeq), ""))
        aveposition <- rowMeans(
            SangerConsensus@forwardReadsList[[i]]@peakPosMatrix, na.rm=TRUE)
        basecalls2 <- basecalls2[1:length(aveposition)]
        basecalls2DF <- data.frame(
            t(data.frame(basecalls2)), stringsAsFactors = FALSE)
        colnames(basecalls2DF) <- substr(colnames(basecalls2DF), 2, 100)
        rownames(basecalls2DF) <- NULL
        return(basecalls2DF)
        }
    )
    SangerSingleReadSecoSeqDF <- c(forwardReadSecoSeqDF, reverseReadSecoSeqDF)

    # primaryAASeq
    forwardReadPrimAASeq <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@primaryAASeq)
    reverseReadPrimAASeq <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@primaryAASeq)
    SangerSingleReadPrimAASeq <- c(forwardReadPrimAASeq, reverseReadPrimAASeq)

    # primaryAASeqDF
    forwardReadPrimAASeqDF <- lapply(1:forwardReadNum, function(i) {
        AAString <- SangerConsensus@forwardReadsList[[i]]@primaryAASeq
        AAStringDF <- data.frame(
            t(data.frame(AAString)), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
        }
    )
    reverseReadPrimAASeqDF <- lapply(1:reverseReadNum, function(i) {
        AAString <- SangerConsensus@reverseReadsList[[i]]@primaryAASeq
        AAStringDF <- data.frame(
            t(data.frame(AAString)), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
        }
    )
    SangerSingleReadPrimAASeqDF <- c(forwardReadPrimAASeqDF,
                                     reverseReadPrimAASeqDF)

    # traceMatrix
    forwardReadTraceMat <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@traceMatrix)
    reverseReadTraceMat <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@traceMatrix)
    SangerSingleReadTraceMat <- c(forwardReadTraceMat, reverseReadTraceMat)

    # peakPosMatrix
    forwardReadReadPeakPosMat <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@peakPosMatrix)
    reverseReadReadPeakPosMat <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@peakPosMatrix)
    SangerSingleReadPeakPosMat <- c(forwardReadReadPeakPosMat,
                                    reverseReadReadPeakPosMat)
    # peakAmpMatrix
    forwardReadPeakAmpMat <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@peakAmpMatrix)
    reverseReadPeakAmpMat <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@peakAmpMatrix)
    SangerSingleReadPeakAmpMat <- c(forwardReadPeakAmpMat,
                                    reverseReadPeakAmpMat)
    # trimmedQS <- reactiveValues(cuffOffQuality = 0, M2SlidingWindowSize = 0)

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
    ### Functions for all UI page
    ############################################################################
    ### ------------------------------------------------------------------------
    ### dynamic side menu bar
    ### ------------------------------------------------------------------------
    dynamicMenuSideBarSC(input, output, session,
                         SangerSingleReadNum, SangerSingleReadFeature)


    output$consensusRead_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (input$sidebar_menu == "Sanger Consensus Read Overview") {
            ### ----------------------------------------------------------------
            ### Dynamic page navigation: consensusRead content overview
            ### ----------------------------------------------------------------
            # refAminoAcidSeq (string)
            fluidRow(
                useShinyjs(),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Basic Information: "),
                                   style = "font-size: 26px;
                                       font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 0.2px hidden #A9A9A9;")),
                    fluidRow(
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
                                      tags$p(tagList(icon("caret-right"),
                                                     "Raw ABI Parent Directory: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(parentDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Consenesus Read Name: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(consenesusReadName),
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
                                      h4(suffixForwardRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Forward Read Number: "),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(forwardReadNum),
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
                                      h4(suffixReverseRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Reverse Read Number: "),
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
                    # If it is null
                    # h1(SCRefAminoAcidSeq),
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
                                           width = "100%", height = "50"),
                               style = paste("height:100%; ",
                                             "overflow-y: hidden;",
                                             "overflow-x: scroll;")
                        ),
                    ),
                ),

                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Consensus Read Results: "),
                                   style = "font-size: 26px;
                                       font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Alignment",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               htmlOutput("consensusAlignmentHTML"),
                        ),
                    ),
                    box(title = tags$p("Differences Dataframe",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SCDifferencesDFUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Dendrogram",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               # plot()
                               plotOutput("dendrogramPlot"),
                               # dataTableOutput("differencesDF")
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        ),
                        column(width = 12,
                               dataTableOutput("dendrogramDF"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Secondary Peak Dataframe",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SCDistanceMatrixUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Indels Dataframe",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SCIndelsDFUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Stop Codons Dataframe",
                                       style = "font-size: 24px;
                                       font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SCStopCodonsDFUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    )
                ),
            )
        } else {
            if (!is.na(as.numeric(sidebar_menu[[1]]))) {
                fluidRow(
                    useShinyjs(),
                    box(title = tags$p(tagList(icon("dot-circle"),
                                               "Raw File: "),
                                       style = "font-size: 26px;
                                         font-weight: bold;"),
                        solidHeader = TRUE,
                        status = "success", width = 12,
                        h1(paste0(
                            SangerSingleReadBFN[[strtoi(sidebar_menu[[1]])]])),
                        tags$h5(paste("( full path:",
                                      SangerSingleReadAFN[[
                                          strtoi(sidebar_menu[[1]])]],
                                      ")"), style = "font-style:italic")),
                    box(title = tags$p(tagList(icon("dot-circle"),
                                               "Primary & Secondary Sequences: "),
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
                        tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                        box(title = tags$p(tagList(icon("arrow-circle-right"),
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
    })


    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic rightHeader text
    ### ------------------------------------------------------------------------
    observeEventDynamicHeaderSC(input, output, session, trimmedRV,
                              SangerSingleReadQualReport)

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        newS4Object <- file.path(shinyDirectory, "SangerConsensus.Rda")
        showNotification(paste("New S4 object is store as:", newS4Object),
                         type = "message", duration = 10)
        # Read number
        forwardReadNum <- length((SangerConsensus)@forwardReadsList)
        reverseReadNum <- length((SangerConsensus)@reverseReadsList)
        SangerSingleReadNum <- forwardReadNum + reverseReadNum

        sapply(1:forwardReadNum, function(i) {
            SangerConsensus@forwardReadsList[[i]]@QualityReport <<-
                SangerSingleReadQualReport[[i]]
            message("save SangerConsensus quality S4 object Forward")
            }
        )
        sapply(1:reverseReadNum, function(i) {
            SangerConsensus@reverseReadsList[[i]]@QualityReport <<-
                SangerSingleReadQualReport[[forwardReadNum + i]]
            message("save SangerConsensus quality S4 object Reverse")
            }
        )
        saveRDS(SangerConsensus, file=newS4Object)
        message("New S4 object is store as: ", newS4Object)
        NEW_SANGER_CONSENSUS_READ <<- readRDS(file=newS4Object)
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        btn <- input$closeUI
        stopApp()
    })


    ############################################################################
    ### ConsensusRead (Function for Sanger Consensus Read Overview)
    ############################################################################
    output$geneticCodeDF <- renderExcel({
        excelTable(data =  t(data.frame(SCGeneticCode)),
                   defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })

    output$consensusAlignmentHTML<-renderUI({
        browseSeqHTML <-
            file.path(shinyDirectory,
                      paste0(A_chloroticConsensusReads@consenesusReadName,
                             "_Alignment_BrowseSeqs.html"))
        if (!file.exists(browseSeqHTML)) {
            BrowseSeqs(A_chloroticConsensusReads@alignment,
                       openURL=FALSE, htmlFile=browseSeqHTML)
        }
        includeHTML(
            file.path(shinyDirectory,
                      paste0(A_chloroticConsensusReads@consenesusReadName,
                             "_Alignment_BrowseSeqs.html")))
    })

    output$dendrogramDF <- renderDataTable({
        SCDendrogram[[1]]
    })

    output$dendrogramPlot <- renderPlot({
        plot(SCDendrogram[[2]])
    })

    output$SCDifferencesDFUI <- renderUI({
        if (all(dim(SCDifferencesDF) == c(0,0))) {
            h4("*** 'Differences' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDifferencesDF")

        }
    })

    output$SCDistanceMatrixUI <- renderUI({
        if (all(dim(SCDistanceMatrix) == c(0,0))) {
            h4("*** 'Distance' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDistanceMatrix")

        }
    })

    output$SCIndelsDFUI <- renderUI({
        if (all(dim(SCIndelsDF) == c(0,0))) {
            h4("*** 'Indels' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCIndelsDF")

        }
    })

    output$SCStopCodonsDFUI <- renderUI({
        if (all(dim(SCStopCodonsDF) == c(0,0))) {
            h4("*** 'Stop Codons' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCStopCodonsDF")
        }
    })

    output$SCDifferencesDF = renderDataTable({
        SCDifferencesDF
    })

    output$SCDistanceMatrix = renderDataTable({
        SCDistanceMatrix
    })

    output$SCIndelsDF <- renderDataTable({
        SCIndelsDF
    })

    output$SCStopCodonsDF <- renderDataTable({
        SCStopCodonsDF
    })


    ############################################################################
    ### SangerSingleRead (Function for singel read in consensusRead)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### primarySeq & secondSeq related
    ### ------------------------------------------------------------------------
    output$primarySeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        excelTable(data =
                       SangerSingleReadPrimSeqDF[[singleReadIndex]],
                   defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)

    })

    output$secondSeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        excelTable(data =
                       SangerSingleReadSecoSeqDF[[singleReadIndex]],
                   defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })

    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        excelTable(data =
                       SangerSingleReadQSDF[[singleReadIndex]],
                   defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })

    output$PrimAASeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        excelTable(data =
                       SangerSingleReadPrimAASeqDF[[singleReadIndex]],
                   defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
    })

    ### ------------------------------------------------------------------------
    ### Quality trimming related (value box)
    ### ------------------------------------------------------------------------
    valueBoxSCMinReadsNum(input, output, SCMinReadsNum, session)
    valueBoxSCMinReadLength(input, output, SCMinReadLength, session)
    valueBoxSCMinFractionCall(input, output, SCMinFractionCall, session)
    valueBoxSCMaxFractionLost(input, output, SCMaxFractionLost, session)
    valueBoxSCAcceptStopCodons(input, output, SCAcceptStopCodons, session)
    valueBoxSCReadingFrame(input, output, SCReadingFrame, session)

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


    qualityTrimmingRatioPlot (input, output, session, trimmedRV,
                              SangerSingleReadQualReport,
                              SangerSingleReadFeature)
    qualityQualityBasePlot (input, output, session, trimmedRV,
                            SangerSingleReadQualReport, SangerSingleReadFeature)


    ### ------------------------------------------------------------------------
    ### observeEvent: Method 1 trimming parameter (M1TrimmingCutoff)
    ### ------------------------------------------------------------------------
    observeEvent(input$M1TrimmingCutoffText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
            as.numeric(input$M1TrimmingCutoffText) > 0 &&
            as.numeric(input$M1TrimmingCutoffText) <= 1) {
            inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
        } else {
            inputM1TrimmingCutoffText <- 0.0003
        }
        if (SangerSingleReadQualReport[[
            strtoi(sidebar_menu[[1]])]]@TrimmingMethod == "M1") {
            trimmingPos <-
                M1inside_calculate_trimming(
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        qualityPhredScores,
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@qualityBaseScore,
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
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Method 2 trimming parameter (M2CutoffQualityScore)
    ### ------------------------------------------------------------------------
    observeEvent(input$M2CutoffQualityScoreText, {
        message("Inside M2CutoffQualityScoreText:")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
            strtoi(input$M2CutoffQualityScoreText) > 0 &&
            strtoi(input$M2CutoffQualityScoreText) <= 60 &&
            strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
            inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
        } else {
            inputM2CutoffQualityScoreText <- 20
        }
        if (SangerSingleReadQualReport[[
            strtoi(sidebar_menu[[1]])]]@TrimmingMethod == "M2") {
            trimmingPos <-
                M2inside_calculate_trimming(
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        qualityPhredScores,
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        qualityBaseScore,
                    strtoi(inputM2CutoffQualityScoreText),
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
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
                SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                    M2CutoffQualityScore <<- strtoi(inputM2CutoffQualityScoreText)

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
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Method 2 trimming parameter (M2SlidingWindowSize)
    ### ------------------------------------------------------------------------
    observeEvent(input$M2SlidingWindowSizeText, {
        message("Inside M2SlidingWindowSizeText:")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
            strtoi(input$M2SlidingWindowSizeText) > 0 &&
            strtoi(input$M2SlidingWindowSizeText) <= 20 &&
            strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
            inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
        } else {
            inputM2SlidingWindowSizeText <- 5
        }
        if (SangerSingleReadQualReport[[
            strtoi(sidebar_menu[[1]])]]@TrimmingMethod == "M2") {
            trimmingPos <-
                M2inside_calculate_trimming(
                    SangerSingleReadQualReport[[strtoi(sidebar_menu[[1]])]]@
                        qualityPhredScores,
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@qualityBaseScore,
                    SangerSingleReadQualReport[[
                        strtoi(sidebar_menu[[1]])]]@M2CutoffQualityScore,
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
    })

    # chromatogram
    ### ------------------------------------------------------------------------
    ### chromatogram related feature
    ### ------------------------------------------------------------------------
    output$chromatogramUIOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (input$sidebar_menu != "Sanger Consensus Read Overview") {
            if (!is.na(as.numeric(sidebar_menu[[1]]))) {
                chromatogramRowNumAns <-
                    chromatogramRowNum(SangerConsensusFRReadsList[[singleReadIndex]],
                        strtoi(input$ChromatogramBasePerRow)) * 200
                message("chromatogramRowNumAns: ", chromatogramRowNumAns)
                plotOutput("chromatogram", height = chromatogramRowNumAns) %>%
                    withSpinner()
            }
        }
    })

    output$chromatogram <- renderPlot({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (input$sidebar_menu != "Sanger Consensus Read Overview") {
            if (!is.na(as.numeric(sidebar_menu[[1]]))) {
                rawSeqLength <-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawSeqLength

                hetcalls <-
                    makeBaseCalls(SangerConsensusFRReadsList[[
                        singleReadIndex]],
                                  ratio = as.numeric(
                                      input$ChromatogramSignalRatioCutoff))
                chromatogram(hetcalls,
                             width = strtoi(input$ChromatogramBasePerRow),
                             height = 2, trim5 = trimmedRV[["trimmedStartPos"]],
                             trim3 = rawSeqLength -
                                 trimmedRV[["trimmedFinishPos"]],
                             showtrim = (input$ChromatogramCheckShowTrimmed),
                             showcalls = "both")
            }
        }
    })


    ############################################################################
    ### Switch trimming method related function
    ############################################################################
    output$TrimmingMethodUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.null(SangerSingleReadQualReport[[singleReadIndex]])) {
            if (SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M1") {
                if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff)) {
                    SangerSingleReadQualReport[[singleReadIndex]]@M1TrimmingCutoff <<-  0.0001
                }
                fluidRow(
                    column(6,
                           uiOutput("M1TrimmingCutoff") ,
                           tags$ul(
                               textInput("M1TrimmingCutoffText",
                                         label = p("Input Value"),
                                         value = toString(
                                             SangerSingleReadQualReport
                                             [[singleReadIndex]]@
                                                 M1TrimmingCutoff),
                                         width = '70%')
                           ),
                    ),
                )
            } else if (SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M2") {
                message("Inside Method 2!!")
                if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore)) {
                    SangerSingleReadQualReport[[singleReadIndex]]@M2CutoffQualityScore <<-  20
                }
                if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize )) {
                    SangerSingleReadQualReport[[singleReadIndex]]@M2SlidingWindowSize <<-  5
                }
                fluidRow(
                    column(6,
                           uiOutput("M2CutoffQualityScore") ,
                           tags$ul(
                               textInput("M2CutoffQualityScoreText",
                                         label = p("Input Value"),
                                         value = toString(
                                             SangerSingleReadQualReport
                                             [[singleReadIndex]]@
                                                 M2CutoffQualityScore),
                                         width = '70%')
                           ),
                    ),
                    column(6,
                           uiOutput("M2SlidingWindowSize") ,
                           tags$ul(
                               textInput("M2SlidingWindowSizeText",
                                         label = p("Input Value"),
                                         value = toString(
                                             SangerSingleReadQualReport
                                             [[singleReadIndex]]@
                                                 M2SlidingWindowSize),
                                         width = '70%')
                           ),
                    ),
                )
            }
        }
    })

    output$TrimmingMethodSelectionOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.null(SangerSingleReadQualReport[[singleReadIndex]])) {
            if (SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection : 'Logarithmic Scale Trimming'")
            } else if (SangerSingleReadQualReport[[singleReadIndex]]@TrimmingMethod == "M2") {
                tagList(icon("check-circle"),
                        "Your trimming method selection : 'Logarithmic Scale Sliding Window Trimming'")
            }
        }
    })
}

