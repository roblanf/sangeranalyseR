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


    consensusParamSet <-
        reactiveValues(consensusReadSCSet = SangerConsensusSet@consensusReadSCSet,
                       alignmentSCSet     = SangerConsensusSet@alignmentSCSet,
                       alignmentTreeSCSet = SangerConsensusSet@alignmentTreeSCSet)

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
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        singleReadIndex <- strtoi(sidebar_menu[[4]])
        message(input$sidebar_menu)
        if (input$sidebar_menu == "Sanger Aligned Consensus Set Overview") {
            h1(input$sidebar_menu)

        } else {
            if (!is.na(consensusReadIndex)) {
                if (sidebar_menu[[2]] == "Sanger" &&
                    sidebar_menu[[3]] == "Consensus" &&
                    sidebar_menu[[4]] == "Read" &&
                    sidebar_menu[[5]] == "Overview") {
                    message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                    h1(input$sidebar_menu)


                } else if (sidebar_menu[[2]] == "CR" &&
                           sidebar_menu[[3]] == "-" &&
                           !is.na(singleReadIndex) &&
                           sidebar_menu[[5]] == "Forward" &&
                           sidebar_menu[[6]] == "Read") {
                    message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                    h1(input$sidebar_menu)


                } else if (sidebar_menu[[2]] == "CR" &&
                          sidebar_menu[[3]] == "-" &&
                          !is.na(singleReadIndex) &&
                          sidebar_menu[[5]] == "Reverse" &&
                          sidebar_menu[[6]] == "Read") {
                    message(">>>>>>>> Inside '", input$sidebar_menu, "'")
                    h1(input$sidebar_menu)

                }
            }
        }
    })


}

