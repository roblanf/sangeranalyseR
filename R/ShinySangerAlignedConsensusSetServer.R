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
        SCName <- paste0(i, "_Consensus_Read")
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
        SangerConsensusForRegExp <- SangerCSetList[[i]]@consenesusReadName
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
            paste0(i, "_",
                   SangerSingleReadFReadsList[[i]]@readFeature))
        reverseReadFeature <- sapply(1:reverseReadNum, function(i)
            paste0(i+forwardReadNum, "_",
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
                    SCMaxFractionLost = SCMinFractionCall,
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
                    SangerConsensusForRegExp = SangerConsensusForRegExp,
                    SangerConsensusForRegExp = SangerConsensusForRegExp,
                    SangerConsensusRevRegExp = SangerConsensusRevRegExp,
                    SangerSingleReadFReadsList = SangerSingleReadFReadsList,
                    SangerSingleReadRReadsList = SangerSingleReadRReadsList,
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
    observeEventDynamicHeaderSCSet(input, output, session, trimmedRV,
                                               SangerCSetParam)
    observeEventButtonSaveSCSet(input, output, session, SangerCSetParam)
    observeEventButtonClose(input, output, session)
}
