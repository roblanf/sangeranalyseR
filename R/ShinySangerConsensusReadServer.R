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
        paste0(i, " ",
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
        basecalls1DF <- data.frame(
            t(data.frame(basecalls1)), stringsAsFactors = FALSE)
        colnames(basecalls1DF) <- substr(colnames(basecalls1DF), 2, 100)
        rownames(basecalls1DF) <- NULL
        return(basecalls1DF)
    })
    reverseReadPrimSeqDF <- lapply(1:reverseReadNum, function(i) {
        basecalls1 <- unlist(strsplit(
            toString(SangerConsensus@reverseReadsList[[i]]@primarySeq), ""))
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
        basecalls2DF <- data.frame(
            t(data.frame(basecalls2)), stringsAsFactors = FALSE)
        colnames(basecalls2DF) <- substr(colnames(basecalls2DF), 2, 100)
        rownames(basecalls2DF) <- NULL
        return(basecalls2DF)
    })
    reverseReadSecoSeqDF <- lapply(1:reverseReadNum, function(i) {
        basecalls2 <- unlist(strsplit(toString(
            SangerConsensus@reverseReadsList[[i]]@secondarySeq), ""))
        basecalls2DF <- data.frame(
            t(data.frame(basecalls2)), stringsAsFactors = FALSE)
        colnames(basecalls2DF) <- substr(colnames(basecalls2DF), 2, 100)
        rownames(basecalls2DF) <- NULL
        return(basecalls2DF)
    })
    SangerSingleReadSecoSeqDF <- c(forwardReadSecoSeqDF, reverseReadSecoSeqDF)

    # primaryAASeqS1DF
    forwardReadPrimAASeqS1DF <- lapply(1:forwardReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@forwardReadsList[[i]]@primaryAASeqS1)
        AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    reverseReadPrimAASeqS1DF <- lapply(1:reverseReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@reverseReadsList[[i]]@primaryAASeqS1)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    SangerSingleReadPrimAASeqS1DF <- c(forwardReadPrimAASeqS1DF,
                                       reverseReadPrimAASeqS1DF)

    # primaryAASeqS2DF
    forwardReadPrimAASeqS2DF <- lapply(1:forwardReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@forwardReadsList[[i]]@primaryAASeqS2)
        AAString <- rbind(NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    reverseReadPrimAASeqS2DF <- lapply(1:reverseReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@reverseReadsList[[i]]@primaryAASeqS2)
        AAString <- rbind(NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    SangerSingleReadPrimAASeqS2DF <- c(forwardReadPrimAASeqS2DF,
                                       reverseReadPrimAASeqS2DF)

    # primaryAASeqS3DF
    forwardReadPrimAASeqS3DF <- lapply(1:forwardReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@forwardReadsList[[i]]@primaryAASeqS3)
        AAString <- rbind(NA, NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    reverseReadPrimAASeqS3DF <- lapply(1:reverseReadNum, function(i) {
        AAString <- data.frame(SangerConsensus@reverseReadsList[[i]]@primaryAASeqS3)
        AAString <- rbind(NA, NA, AAString)
        AAStringDF <- data.frame(
            t(AAString), stringsAsFactors = FALSE)
        colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
        rownames(AAStringDF) <- NULL
        return(AAStringDF)
    })
    SangerSingleReadPrimAASeqS3DF <- c(forwardReadPrimAASeqS3DF,
                                       reverseReadPrimAASeqS3DF)


    ### ------------------------------------------------------------------------
    ### ConsensusRead reactiveValue
    ### ------------------------------------------------------------------------
    consensusParam <-
        reactiveValues(
            consensusRead   = SangerConsensus@consensusRead,
            differencesDF   = SangerConsensus@differencesDF,
            alignment       = SangerConsensus@alignment,
            distanceMatrix  = SangerConsensus@distanceMatrix,
            dendrogram      = SangerConsensus@dendrogram,
            indelsDF        = SangerConsensus@indelsDF,
            stopCodonsDF    = SangerConsensus@stopCodonsDF,
            secondaryPeakDF = SangerConsensus@secondaryPeakDF)

    ### ------------------------------------------------------------------------
    ### SingleRead reactiveValue
    ### ------------------------------------------------------------------------
    sequenceParam <- reactiveValues(primarySeq = 0,
                                    secondarySeq = 0,
                                    primaryAASeqS1 = 0,
                                    primaryAASeqS2 = 0,
                                    primaryAASeqS3 = 0)

























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
    ### Functions for all UI page
    ############################################################################
    ### ------------------------------------------------------------------------
    ### dynamic side menu bar
    ### ------------------------------------------------------------------------
    dynamicMenuSideBarSC(input, output, session,
                         forwardReadNum, reverseReadNum,
                         forwardReadFeature, reverseReadFeature)

    ############################################################################
    ### Main page switch
    ############################################################################
    output$consensusRead_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        directionParam <- sidebar_menu[[2]]
        if (input$sidebar_menu == "Sanger Consensus Read Overview") {
            message(">>>>>>>> Inside '", input$sidebar_menu, "'")
            h1(input$sidebar_menu)

        } else if (!is.na(strtoi(singleReadIndex)) &&
                   directionParam == "Forward") {

            message(">>>>>>>> Inside '", input$sidebar_menu, "'")
            h1(input$sidebar_menu)

        } else if (!is.na(strtoi(singleReadIndex)) &&
                   directionParam == "Reverse") {
            message(">>>>>>>> Inside '", input$sidebar_menu, "'")
            h1(input$sidebar_menu)
        }
    })


    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic rightHeader text
    ### ------------------------------------------------------------------------
    #!!!!!!!! Fix
    observeEventDynamicHeaderSC(input, output, session, trimmedRV,
                              SangerSingleReadQualReport)





















}

