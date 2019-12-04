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
    SCTrimmingMethod <- SangerConsensus@forwardReadsList[[1]]@
        QualityReport@TrimmingMethod
    if (SCTrimmingMethod == "M1") {
        SCTrimmingMethodName = "Method 1:
                                'Logarithmic Scale Trimming'"
    } else if (SCTrimmingMethod == "M2") {
        SCTrimmingMethodName = "Method 2:
                                'Logarithmic Scale Sliding Window Trimming'"
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
    })
    SangerSingleReadQSDF <- c(forwardQualityScoreDF, reverseQualityScoreDF)

    # ChromatogramParam
    forwardReadChromatogramParam <- sapply(1:forwardReadNum, function(i)
        SangerConsensus@forwardReadsList[[i]]@ChromatogramParam)
    reverseReadChromatogramParam <- sapply(1:reverseReadNum, function(i)
        SangerConsensus@reverseReadsList[[i]]@ChromatogramParam)
    SangerSingleReadChromatogramParam <- c(forwardReadChromatogramParam,
                                           reverseReadChromatogramParam)

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
    })
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
    })
    SangerSingleReadPrimSeqDF <- c(forwardReadPrimSeqDF, reverseReadPrimSeqDF)

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
    })
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
    })
    SangerSingleReadSecoSeqDF <- c(forwardReadSecoSeqDF, reverseReadSecoSeqDF)

    # primaryAASeqDF
    forwardReadPrimAASeqDF <- lapply(1:forwardReadNum, function(i) {
        AAString <- SangerConsensus@forwardReadsList[[i]]@primaryAASeq
        AAStringDF <- data.frame(
            t(data.frame(AAString)), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    reverseReadPrimAASeqDF <- lapply(1:reverseReadNum, function(i) {
        AAString <- SangerConsensus@reverseReadsList[[i]]@primaryAASeq
        AAStringDF <- data.frame(
            t(data.frame(AAString)), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    SangerSingleReadPrimAASeqDF <- c(forwardReadPrimAASeqDF,
                                     reverseReadPrimAASeqDF)

    trimmedRV <- reactiveValues(rawSeqLength            = 0,
                                rawMeanQualityScore     = 0,
                                rawMinQualityScore      = 0,
                                trimmedStartPos         = 0,
                                trimmedFinishPos        = 0,
                                trimmedSeqLength        = 0,
                                trimmedMeanQualityScore = 0,
                                trimmedMinQualityScore  = 0,
                                remainingRatio          = 0)

    consensusParam <- reactiveValues(consensusRead   = NULL,
                                     differencesDF   = NULL,
                                     alignment       = NULL,
                                     distanceMatrix  = NULL,
                                     dendrogram      = NULL,
                                     indelsDF        = NULL,
                                     stopCodonsDF    = NULL,
                                     secondaryPeakDF = NULL)

    trimmedParam <- reactiveValues(M1TrimmingCutoff     = 0,
                                   M2CutoffQualityScore = 0,
                                   M2SlidingWindowSize  = 0)

    ChromatogramParam <- reactiveValues(baseNumPerRow     = 0,
                                        heightPerRow      = 0,
                                        signalRatioCutoff = 0,
                                        showTrimmed       = TRUE)



    ############################################################################
    ### Functions for all UI page
    ############################################################################
    ### ------------------------------------------------------------------------
    ### dynamic side menu bar
    ### ------------------------------------------------------------------------
    dynamicMenuSideBarSC(input, output, session,
                         SangerSingleReadNum, SangerSingleReadFeature)

    ############################################################################
    ### Main page switch
    ############################################################################
    output$consensusRead_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (input$sidebar_menu == "Sanger Consensus Read Overview") {
            message(">>>>>>>> Inside '", input$sidebar_menu, "'")
            ### ----------------------------------------------------------------
            ### Dynamic page navigation: consensusRead content overview
            ### ----------------------------------------------------------------
            ### ----------------------------------------------------------------
            ### 1. Re-calculate consensus read & Update global variable
            ### ----------------------------------------------------------------
            consensusParam[["consensusRead"]] <<- SangerConsensus@consensusRead
            consensusParam[["differencesDF"]] <<- SangerConsensus@differencesDF
            consensusParam[["alignment"]] <<- SangerConsensus@alignment
            consensusParam[["distanceMatrix"]] <<-SangerConsensus@distanceMatrix
            consensusParam[["dendrogram"]] <<- SangerConsensus@dendrogram
            consensusParam[["indelsDF"]] <<- SangerConsensus@indelsDF
            consensusParam[["stopCodonsDF"]] <<- SangerConsensus@stopCodonsDF
            consensusParam[["secondaryPeakDF"]] <<-
                SangerConsensus@secondaryPeakDF
            fluidRow(
                useShinyjs(),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Basic Information: "),
                                   style = "font-size: 26px;
                                   font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 2px hidden #A9A9A9;")),
                    fluidRow(
                        column(width = 12,
                               actionBttn("recalculateButton",
                                          "Re-calculate consensus read",
                                          icon = icon("calculator"),
                                          style = "simple", color = "danger",
                                          block = TRUE, size = "lg")
                        ),
                        column(12,
                               tags$hr(
                                   style = ("border-top: 2px hidden #A9A9A9;")),
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
                                      tags$p(tagList(icon("caret-right"),
                                                "Raw ABI Parent Directory:"),
                                             style = "font-size: 20px;
                                       font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerConsensus@parentDirectory),
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
                                      h4(SangerConsensus@consensusReadName),
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
                                      h4(SangerConsensus@suffixForwardRegExp),
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
                                      h4(SangerConsensus@suffixReverseRegExp),
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
                    uiOutput("SCrefAminoAcidSeq") ,
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
                    box(title = tags$p("Differences Data frame",
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
                               plotOutput("dendrogramPlot"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        ),
                        column(width = 12,
                               tags$hr(
                                   style = ("border-top: 4px hidden #A9A9A9;")),
                        ),
                        column(width = 12,
                               dataTableOutput("dendrogramDF"),
                               style = paste("height:100%; overflow-y:",
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
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        ),
                        column(width = 12,
                               tags$hr(
                                   style = ("border-top: 4px hidden #A9A9A9;")),
                        ),
                        column(width = 12,
                               uiOutput("SCDistanceMatrixUI"),
                               style = paste("height:100%; overflow-y:",
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
                               style = paste("height:100%; overflow-y:",
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
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    )
                ),
            )
        } else {
            if (!is.na(strtoi(singleReadIndex))) {
                ### ------------------------------------------------------------
                ### First assign the ChromatogramParam parameter
                ### ------------------------------------------------------------
                message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                singleReadIndex <- strtoi(sidebar_menu[[1]])
                ChromatogramParam[["baseNumPerRow"]] <<-
                    SangerSingleReadChromatogramParam[[singleReadIndex]]@
                    baseNumPerRow
                ChromatogramParam[["heightPerRow"]] <<-
                    SangerSingleReadChromatogramParam[[singleReadIndex]]@
                    heightPerRow
                ChromatogramParam[["signalRatioCutoff"]] <<-
                    SangerSingleReadChromatogramParam[[singleReadIndex]]@
                    signalRatioCutoff
                ChromatogramParam[["showTrimmed"]] <<-
                    SangerSingleReadChromatogramParam[[singleReadIndex]]@
                    showTrimmed

                trimmedParam[["M1TrimmingCutoff"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    M1TrimmingCutoff
                trimmedParam[["M2CutoffQualityScore"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    M2CutoffQualityScore
                trimmedParam[["M2SlidingWindowSize"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    M2SlidingWindowSize


                trimmedRV[["rawSeqLength"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerSingleReadQualReport[[singleReadIndex]]@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerSingleReadQualReport[[singleReadIndex]]@
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
                    box(title = tags$p(tagList(icon("dot-circle"),
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
                            column(3,
                                   sliderInput("ChromatogramBasePerRow",
                                               label =h4("Base Number Per Row"),
                                               min = 5,
                                               max = 200,
                                               value = ChromatogramParam[["baseNumPerRow"]]),
                                   sliderInput("ChromatogramHeightPerRow",
                                               label = h4("Height Per Row"),
                                               min = 50,
                                               max = 600,
                                               value = ChromatogramParam[["heightPerRow"]]),
                            ),
                            column(3,
                                   tags$hr(
                                       style =
                                           ("border-top: 4px hidden #A9A9A9;")),
                                   numericInput(
                                       "ChromatogramSignalRatioCutoff",
                                       h3("Signal Ratio Cutoff"),
                                       value = ChromatogramParam[["signalRatioCutoff"]]),
                                   checkboxInput(
                                       "ChromatogramCheckShowTrimmed",
                                       "Whether show trimmed region",
                                       value =
                                           ChromatogramParam[["showTrimmed"]])
                            ),
                            column(3,
                                   tags$hr(
                                      style=("border-top: 4px hidden #A9A9A9;")),
                                   uiOutput("ChromatogramtrimmedStartPos"),
                            ),
                            column(3,
                                   tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
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
        ### --------------------------------------------------------------------
        ### Save SangerConsensus quality S4 object
        ### --------------------------------------------------------------------
        forwardReadNum <- length((SangerConsensus)@forwardReadsList)
        reverseReadNum <- length((SangerConsensus)@reverseReadsList)
        SangerSingleReadNum <- forwardReadNum + reverseReadNum
        sapply(1:forwardReadNum, function(i) {
            SangerConsensus@forwardReadsList[[i]]@QualityReport <<-
                SangerSingleReadQualReport[[i]]
            SangerConsensus@forwardReadsList[[i]]@ChromatogramParam <<-
                SangerSingleReadChromatogramParam[[i]]
            message("save SangerConsensus quality S4 object Forward")
        })
        sapply(1:reverseReadNum, function(i) {
            SangerConsensus@reverseReadsList[[i]]@QualityReport <<-
                SangerSingleReadQualReport[[forwardReadNum + i]]
            SangerConsensus@reverseReadsList[[i]]@ChromatogramParam <<-
                SangerSingleReadChromatogramParam[[forwardReadNum + i]]
            message("save SangerConsensus quality S4 object Reverse")
        })
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

    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus read re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButton, {
        message("######## Reactive button clicked !!!")
        message("######## Start recalculating consensus read (SC")
        CSResult<-
            calculateConsensusRead (SangerConsensus@forwardReadsList,
                                    SangerConsensus@reverseReadsList,
                                    SangerConsensus@refAminoAcidSeq,
                                    SangerConsensus@minFractionCall,
                                    SangerConsensus@maxFractionLost,
                                    SangerConsensus@geneticCode,
                                    SangerConsensus@acceptStopCodons,
                                    SangerConsensus@readingFrame)

        SangerConsensus@consensusRead <<- CSResult$consensusGapfree
        SangerConsensus@differencesDF <<- CSResult$diffsDf
        SangerConsensus@alignment <<- CSResult$aln2
        SangerConsensus@distanceMatrix <<- CSResult$dist
        SangerConsensus@dendrogram <<- CSResult$dend
        SangerConsensus@indelsDF <<- CSResult$indels
        SangerConsensus@stopCodonsDF <<- CSResult$stopsDf
        SangerConsensus@secondaryPeakDF <<- CSResult$spDf

        consensusParam[["consensusRead"]] <<- SangerConsensus@consensusRead
        consensusParam[["differencesDF"]] <<- SangerConsensus@differencesDF
        consensusParam[["alignment"]] <<- SangerConsensus@alignment
        consensusParam[["distanceMatrix"]] <<-SangerConsensus@distanceMatrix
        consensusParam[["dendrogram"]] <<- SangerConsensus@dendrogram
        consensusParam[["indelsDF"]] <<- SangerConsensus@indelsDF
        consensusParam[["stopCodonsDF"]] <<- SangerConsensus@stopCodonsDF
        consensusParam[["secondaryPeakDF"]] <<-
            SangerConsensus@secondaryPeakDF
        message("######## Finish recalculation")
    })


    ### ------------------------------------------------------------------------
    ### observeEvent: Button Consensus chromatogram parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$saveChromatogramParam, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        ### ------------------------------------------------------------
        ### Update ChromatogramBasePerRow
        ### ------------------------------------------------------------
        SangerSingleReadChromatogramParam[[singleReadIndex]]@
            baseNumPerRow <<- input$ChromatogramBasePerRow
        SangerSingleReadChromatogramParam[[singleReadIndex]]@
            heightPerRow <<- input$ChromatogramHeightPerRow
        SangerSingleReadChromatogramParam[[singleReadIndex]]@
            signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
        SangerSingleReadChromatogramParam[[singleReadIndex]]@
            showTrimmed <<- input$ChromatogramCheckShowTrimmed

        ### ------------------------------------------------------------
        ### Save SangerConsensus quality S4 object
        ### ------------------------------------------------------------
        forwardReadNum <- length((SangerConsensus)@forwardReadsList)
        reverseReadNum <- length((SangerConsensus)@reverseReadsList)
        SangerSingleReadNum <- forwardReadNum + reverseReadNum

        if (singleReadIndex <= forwardReadNum) {
            # This is forward list
            SangerConsensus@
                forwardReadsList[[singleReadIndex]]@ChromatogramParam <<-
                SangerSingleReadChromatogramParam[[singleReadIndex]]
        } else {
            # This is reverse list
            SangerConsensus@reverseReadsList[[singleReadIndex-forwardReadNum]]@
                ChromatogramParam <<-
                SangerSingleReadChromatogramParam[[singleReadIndex]]
        }
        ChromatogramParam[["baseNumPerRow"]] <<-
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            baseNumPerRow
        ChromatogramParam[["heightPerRow"]] <<-
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            heightPerRow
        ChromatogramParam[["signalRatioCutoff"]] <<-
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            signalRatioCutoff
        ChromatogramParam[["showTrimmed"]] <<-
            SangerSingleReadChromatogramParam[[singleReadIndex]]@
            showTrimmed
    })


    ############################################################################
    ### ConsensusRead (Function for Sanger Consensus Read Overview)
    ############################################################################
    output$geneticCodeDF <- renderExcel({
        suppressMessages(
            excelTable(data =  t(data.frame(SangerConsensus@geneticCode)),
                       defaultColWidth = 50, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })

    output$SCrefAminoAcidSeqDF <- renderExcel({
        refAminoAcidSeqVec <- strsplit(SangerConsensus@refAminoAcidSeq, "")[[1]]
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
    ### refAminoAcidSeq
    ### ------------------------------------------------------------------------
    output$SCrefAminoAcidSeq <- renderUI({
        if (SangerConsensus@refAminoAcidSeq == "") {
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
    ### ------------------------------------------------------------------------
    ### Alignment
    ### ------------------------------------------------------------------------
    output$consensusAlignmentHTML <- renderUI({
        browseSeqHTML <-
            file.path(shinyDirectory,
                      paste0(SangerConsensus@consensusReadName,
                             "_Alignment_BrowseSeqs.html"))
        BrowseSeqs(consensusParam[["alignment"]],
                   openURL=FALSE, htmlFile=browseSeqHTML)
        includeHTML(
            file.path(shinyDirectory,
                      paste0(SangerConsensus@consensusReadName,
                             "_Alignment_BrowseSeqs.html")))
    })

    ### ------------------------------------------------------------------------
    ### Difference
    ### ------------------------------------------------------------------------
    output$SCDifferencesDFUI <- renderUI({
        if (all(dim(consensusParam[["differencesDF"]]) == c(0,0))) {
            h4("*** 'Differences' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDifferencesDF")

        }
    })
    output$SCDifferencesDF = renderDataTable({
        consensusParam[["differencesDF"]]
    })

    ### ------------------------------------------------------------------------
    ### dendrogram
    ### ------------------------------------------------------------------------
    output$dendrogramPlot <- renderPlot({
        # plot(consensusParam[["dendrogram"]][[2]])
        ggdendrogram(consensusParam[["dendrogram"]][[2]], rotate = TRUE)
    })
    output$dendrogramDF <- renderDataTable({
        consensusParam[["dendrogram"]][[1]]
    })

    ### ------------------------------------------------------------------------
    ### Distance
    ### ------------------------------------------------------------------------
    output$SCDistanceMatrixPlotUI <- renderUI({
        if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
            h4("*** 'Distance' dataframe is empty. (Cannot plot)***",
               style="font-weight: bold; font-style: italic;")
        } else {
            plotlyOutput("SCDistanceMatrixPlot")
        }
    })
    output$SCDistanceMatrixPlot <- renderPlotly({
        plot_ly(x = SangerSingleReadBFN,
                y = SangerSingleReadBFN,
                z = consensusParam[["distanceMatrix"]],
                colors = colorRamp(c("white", "#32a852")),
                type = "heatmap")
    })
    output$SCDistanceMatrixUI <- renderUI({
        if (all(dim(consensusParam[["distanceMatrix"]]) == c(0,0))) {
            h4("*** 'Distance' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCDistanceMatrix")
        }
    })
    output$SCDistanceMatrix = renderDataTable({
        consensusParam[["distanceMatrix"]]
    })


    ### ------------------------------------------------------------------------
    ### Indels
    ### ------------------------------------------------------------------------
    output$SCIndelsDFUI <- renderUI({
        if (all(dim(consensusParam[["indelsDF"]]) == c(0,0))) {
            h4("*** 'Indels' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCIndelsDF")

        }
    })
    output$SCIndelsDF <- renderDataTable({
        consensusParam[["indelsDF"]]
    })

    ### ------------------------------------------------------------------------
    ### StopCodons
    ### ------------------------------------------------------------------------
    output$SCStopCodonsDFUI <- renderUI({
        if (all(dim(consensusParam[["stopCodonsDF"]]) == c(0,0))) {
            h4("*** 'Stop Codons' dataframe is empty. ***",
               style="font-weight: bold; font-style: italic;")
        } else {
            dataTableOutput("SCStopCodonsDF")
        }
    })
    output$SCStopCodonsDF <- renderDataTable({
        consensusParam[["stopCodonsDF"]]
    })

    ### ------------------------------------------------------------------------
    ### Valuebox for basic information
    ### ------------------------------------------------------------------------
    valueBoxSCMinReadsNum(input, output,
                          SangerConsensus@minReadsNum, session)
    valueBoxSCMinReadLength(input, output,
                            SangerConsensus@minReadLength, session)
    valueBoxSCMinFractionCall(input, output,
                              SangerConsensus@minFractionCall, session)
    valueBoxSCMaxFractionLost(input, output,
                              SangerConsensus@maxFractionLost, session)
    valueBoxSCAcceptStopCodons(input, output,
                               SangerConsensus@acceptStopCodons, session)
    valueBoxSCReadingFrame(input, output,
                           SangerConsensus@readingFrame, session)

    ############################################################################
    ### SangerSingleRead (Function for singel read in consensusRead)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### primarySeq & secondSeq related
    ### ------------------------------------------------------------------------
    output$primarySeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        suppressMessages(
            excelTable(data =
                           SangerSingleReadPrimSeqDF[[singleReadIndex]],
                       defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })

    output$secondSeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        suppressMessages(
            excelTable(data =
                           SangerSingleReadSecoSeqDF[[singleReadIndex]],
                       defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })

    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        suppressMessages(
            excelTable(data =
                           SangerSingleReadQSDF[[singleReadIndex]],
                       defaultColWidth = 30, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })

    output$PrimAASeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        suppressMessages(
            excelTable(data =
                           SangerSingleReadPrimAASeqDF[[singleReadIndex]],
                       defaultColWidth = 90, editable = TRUE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })

    ### ------------------------------------------------------------------------
    ### Quality trimming related (value box)
    ### ------------------------------------------------------------------------
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
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
            as.numeric(input$M1TrimmingCutoffText) > 0 &&
            as.numeric(input$M1TrimmingCutoffText) <= 1) {
            inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
        } else {
            inputM1TrimmingCutoffText <- 0.0001
        }
        if (SangerSingleReadQualReport[[singleReadIndex]]@
            TrimmingMethod == "M1") {
            trimmingPos <-
                M1inside_calculate_trimming(
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityPhredScores,
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

            SangerSingleReadQualReport[[singleReadIndex]]@
                M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)

            SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength <<- rawSeqLength
            SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore <<- rawMeanQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore <<- rawMinQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos <<- trimmedStartPos
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos <<- trimmedFinishPos
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength <<- trimmedSeqLength
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore <<- trimmedMeanQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore <<- trimmedMinQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                remainingRatio <<- remainingRatio

            trimmedParam[["M1TrimmingCutoff"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                M1TrimmingCutoff

            trimmedRV[["rawSeqLength"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength
            trimmedRV[["rawMeanQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore
            trimmedRV[["rawMinQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore
            trimmedRV[["trimmedStartPos"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos
            trimmedRV[["trimmedFinishPos"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos
            trimmedRV[["trimmedSeqLength"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength
            trimmedRV[["trimmedMeanQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore
            trimmedRV[["trimmedMinQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore
            trimmedRV[["remainingRatio"]] <<-
                round(SangerSingleReadQualReport[[singleReadIndex]]@
                          remainingRatio * 100, 2)
            ### ------------------------------------------------------------
            ### Save SangerConsensus quality S4 object
            ### ------------------------------------------------------------
            forwardReadNum <- length((SangerConsensus)@forwardReadsList)
            reverseReadNum <- length((SangerConsensus)@reverseReadsList)
            SangerSingleReadNum <- forwardReadNum + reverseReadNum
            if (singleReadIndex <= forwardReadNum) {
                # This is forward list
                SangerConsensus@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport <<-
                    SangerSingleReadQualReport[[singleReadIndex]]
            } else {
                # This is reverse list
                SangerConsensus@
                    reverseReadsList[[singleReadIndex-forwardReadNum]]@
                    QualityReport <<-
                    SangerSingleReadQualReport[[singleReadIndex]]
            }
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Method 2 trimming parameter (M2CutoffQualityScore)
    ### ------------------------------------------------------------------------
    observeEvent(input$M2CutoffQualityScoreText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
            strtoi(input$M2CutoffQualityScoreText) > 0 &&
            strtoi(input$M2CutoffQualityScoreText) <= 60 &&
            strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
            inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
        } else {
            inputM2CutoffQualityScoreText <- 20
        }
        if (SangerSingleReadQualReport[[singleReadIndex]]@
            TrimmingMethod == "M2") {
            trimmingPos <-
                M2inside_calculate_trimming(
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityPhredScores,
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityBaseScore,
                    strtoi(inputM2CutoffQualityScoreText),
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
            SangerSingleReadQualReport[[singleReadIndex]]@
                M2CutoffQualityScore <<-
                strtoi(inputM2CutoffQualityScoreText)

            SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength <<- rawSeqLength
            SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore <<- rawMeanQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore <<- rawMinQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos <<- trimmedStartPos
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos <<- trimmedFinishPos
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength <<- trimmedSeqLength
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore <<- trimmedMeanQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore <<- trimmedMinQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                remainingRatio <<- remainingRatio

            trimmedParam[["M2CutoffQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                M2CutoffQualityScore
            trimmedParam[["M2SlidingWindowSize"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                M2SlidingWindowSize

            trimmedRV[["rawSeqLength"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength
            trimmedRV[["rawMeanQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore
            trimmedRV[["rawMinQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore
            trimmedRV[["trimmedStartPos"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos
            trimmedRV[["trimmedFinishPos"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos
            trimmedRV[["trimmedSeqLength"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength
            trimmedRV[["trimmedMeanQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore
            trimmedRV[["trimmedMinQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore
            trimmedRV[["remainingRatio"]] <<-
                round(SangerSingleReadQualReport[[singleReadIndex]]@
                          remainingRatio * 100, 2)
            ### ------------------------------------------------------------
            ### Save SangerConsensus quality S4 object
            ### ------------------------------------------------------------
            forwardReadNum <- length((SangerConsensus)@forwardReadsList)
            reverseReadNum <- length((SangerConsensus)@reverseReadsList)
            SangerSingleReadNum <- forwardReadNum + reverseReadNum
            if (singleReadIndex <= forwardReadNum) {
                # This is forward list
                SangerConsensus@
                    forwardReadsList[[singleReadIndex]]@
                    QualityReport <<-
                    SangerSingleReadQualReport[[singleReadIndex]]
            } else {
                # This is reverse list
                SangerConsensus@
                    reverseReadsList[[singleReadIndex-forwardReadNum]]@
                    QualityReport <<-
                    SangerSingleReadQualReport[[singleReadIndex]]
            }
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Method 2 trimming parameter (M2SlidingWindowSize)
    ### ------------------------------------------------------------------------
    observeEvent(input$M2SlidingWindowSizeText, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
            strtoi(input$M2SlidingWindowSizeText) > 0 &&
            strtoi(input$M2SlidingWindowSizeText) <= 20 &&
            strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
            inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
        } else {
            inputM2SlidingWindowSizeText <- 5
        }
        if (SangerSingleReadQualReport[[singleReadIndex]]@
            TrimmingMethod == "M2") {
            trimmingPos <-
                M2inside_calculate_trimming(
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityPhredScores,
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        qualityBaseScore,
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

            SangerSingleReadQualReport[[singleReadIndex]]@
                M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)

            SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength <<- rawSeqLength
            SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore <<- rawMeanQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore <<- rawMinQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos <<- trimmedStartPos
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos <<- trimmedFinishPos
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength <<- trimmedSeqLength
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore <<- trimmedMeanQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore <<- trimmedMinQualityScore
            SangerSingleReadQualReport[[singleReadIndex]]@
                remainingRatio <<- remainingRatio

            trimmedParam[["M2CutoffQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                M2CutoffQualityScore
            trimmedParam[["M2SlidingWindowSize"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                M2SlidingWindowSize

            trimmedRV[["rawSeqLength"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawSeqLength
            trimmedRV[["rawMeanQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMeanQualityScore
            trimmedRV[["rawMinQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                rawMinQualityScore
            trimmedRV[["trimmedStartPos"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedStartPos
            trimmedRV[["trimmedFinishPos"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedFinishPos
            trimmedRV[["trimmedSeqLength"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedSeqLength
            trimmedRV[["trimmedMeanQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMeanQualityScore
            trimmedRV[["trimmedMinQualityScore"]] <<-
                SangerSingleReadQualReport[[singleReadIndex]]@
                trimmedMinQualityScore
            trimmedRV[["remainingRatio"]] <<-
                round(SangerSingleReadQualReport[[singleReadIndex]]@
                          remainingRatio * 100, 2)
            ### ------------------------------------------------------------
            ### Save SangerConsensus quality S4 object
            ### ------------------------------------------------------------
            forwardReadNum <- length((SangerConsensus)@forwardReadsList)
            reverseReadNum <- length((SangerConsensus)@reverseReadsList)
            SangerSingleReadNum <- forwardReadNum + reverseReadNum
            if (singleReadIndex <= forwardReadNum) {
                # This is forward list
                SangerConsensus@forwardReadsList[[singleReadIndex]]@QualityReport <<-
                    SangerSingleReadQualReport[[singleReadIndex]]
            } else {
                # This is reverse list
                SangerConsensus@
                    reverseReadsList[[singleReadIndex-forwardReadNum]]@
                    QualityReport <<-
                    SangerSingleReadQualReport[[singleReadIndex]]
            }
        }
    })

    ### ------------------------------------------------------------------------
    ### chromatogram related feature
    ### ------------------------------------------------------------------------
    output$chromatogramUIOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (input$sidebar_menu != "Sanger Consensus Read Overview") {
            if (!is.na(as.numeric(sidebar_menu[[1]]))) {
                trimmedRV[["trimmedSeqLength"]]
                chromatogramRowNumAns <-
                    chromatogramRowNum (
                        strtoi(ChromatogramParam[["baseNumPerRow"]]),
                        SangerSingleReadQualReport[[singleReadIndex]]@
                            rawSeqLength,
                        SangerSingleReadQualReport[[singleReadIndex]]@
                            trimmedSeqLength,
                        ChromatogramParam[["showTrimmed"]]) *
                    strtoi(ChromatogramParam[["heightPerRow"]])
                    plotOutput("chromatogram", height = chromatogramRowNumAns)
                        # addSpinner(spin = "circle", color = "#eb4034")
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
                            ChromatogramParam[["signalRatioCutoff"]]))
                    chromatogram(hetcalls,
                                 width = strtoi(
                                     ChromatogramParam[["baseNumPerRow"]]),
                                 height = 2,
                                 trim5 = trimmedRV[["trimmedStartPos"]],
                                 trim3 = rawSeqLength -
                                     trimmedRV[["trimmedFinishPos"]],
                                 showtrim = (ChromatogramParam[["showTrimmed"]]),
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
            if (SangerSingleReadQualReport[[singleReadIndex]]@
                TrimmingMethod == "M1") {
                if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@
                            M1TrimmingCutoff)) {
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        M1TrimmingCutoff <<-  0.0001
                }
                fluidRow(
                    column(6,
                           uiOutput("M1TrimmingCutoff") ,
                           tags$ul(
                               textInput("M1TrimmingCutoffText",
                                         label = p("Input Value"),
                                         value = toString(
                                             trimmedParam[["M1TrimmingCutoff"]]),
                                         width = '70%')
                           ),
                    ),
                )
            } else if (SangerSingleReadQualReport[[singleReadIndex]]@
                       TrimmingMethod == "M2") {
                message("Inside Method 2!!")
                if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@
                            M2CutoffQualityScore)) {
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        M2CutoffQualityScore <<-  20
                }
                if (is.null(SangerSingleReadQualReport[[singleReadIndex]]@
                            M2SlidingWindowSize )) {
                    SangerSingleReadQualReport[[singleReadIndex]]@
                        M2SlidingWindowSize <<-  5
                }

                fluidRow(
                    column(6,
                           uiOutput("M2CutoffQualityScore") ,
                           tags$ul(
                               textInput("M2CutoffQualityScoreText",
                                         label = p("Input Value"),
                                         value = toString(
                                             trimmedParam[["M2CutoffQualityScore"]]),
                                         width = '70%')
                           ),
                    ),
                    column(6,
                           uiOutput("M2SlidingWindowSize") ,
                           tags$ul(
                               textInput("M2SlidingWindowSizeText",
                                         label = p("Input Value"),
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
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.null(SangerSingleReadQualReport[[singleReadIndex]])) {
            if (SangerSingleReadQualReport[[singleReadIndex]]@
                TrimmingMethod == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Logarithmic Scale Trimming'")
            } else if (SangerSingleReadQualReport[[singleReadIndex]]@
                       TrimmingMethod == "M2") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                        'Logarithmic Scale Sliding Window Trimming'")
            }
        }
    })
}

