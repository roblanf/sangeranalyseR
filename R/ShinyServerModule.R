M1inside_calculate_trimming <- function(qualityPhredScores,
                                        qualityBaseScore,
                                        M1TrimmingCutoff) {
    rawSeqLength <- length(qualityBaseScore)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    start = FALSE
    trimmedStartPos = 0
    qualityBaseScoreCutOff = M1TrimmingCutoff - qualityBaseScore
    ### ------------------------------------------------------------------------
    ### calculate cummulative score
    ### if cumulative value < 0, set it to 0
    ### the BioPython implementation always trims the first base,
    ### this implementation does not.
    ### ------------------------------------------------------------------------
    score = qualityBaseScoreCutOff[1]
    if(score < 0){
        score = 0
    }else{
        trimmedStartPos = 1
        start = TRUE
    }
    cummul_score = c(score)
    ### ------------------------------------------------------------------------
    ### trimmedStartPos = value when cummulative score is first > 0
    ### ------------------------------------------------------------------------
    ### ------------------------------------------------------------------------
    ### trimmedFinishPos = index of highest cummulative score,
    ### marking the end of sequence segment with highest cummulative score
    ### ------------------------------------------------------------------------
    for(i in 2:length(qualityBaseScoreCutOff)){
        score = cummul_score[length(cummul_score)] + qualityBaseScoreCutOff[i]
        if (score <= 0) {
            cummul_score = c(cummul_score, 0)
        }else{
            cummul_score = c(cummul_score, score)
            if(start == FALSE){
                trimmedStartPos = i
                start = TRUE
            }
        }
        trimmedFinishPos = which.max(cummul_score)
    }
    ### ------------------------------------------------------------------------
    ### fix an edge case, where all scores are worse than the cutoff
    ### in this case you wouldn't want to keep any bases at all
    ### ------------------------------------------------------------------------
    if(sum(cummul_score)==0){trimmedFinishPos = 0}
    trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
    trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
    trimmedMinQualityScore <- min(trimmedQualityPhredScore)
    trimmedSeqLength = trimmedFinishPos - trimmedStartPos + 1
    remainingRatio = trimmedSeqLength / rawSeqLength

    return(c("rawSeqLength" = rawSeqLength,
             "rawMeanQualityScore" = rawMeanQualityScore,
             "rawMinQualityScore" = rawMinQualityScore,
             "trimmedStartPos" = trimmedStartPos,
             "trimmedFinishPos" = trimmedFinishPos,
             "trimmedSeqLength" = trimmedSeqLength,
             "trimmedMeanQualityScore" = trimmedMeanQualityScore,
             "trimmedMinQualityScore" = trimmedMinQualityScore,
             "remainingRatio" = remainingRatio))
}

M2inside_calculate_trimming <- function(qualityPhredScores,
                                      qualityBaseScore,
                                      M2CutoffQualityScore,
                                      M2SlidingWindowSize) {
    rawSeqLength <- length(qualityBaseScore)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    qualityPbCutoff <- 10** (M2CutoffQualityScore / (-10.0))
    remainingIndex <- c()
    if (M2SlidingWindowSize > 20 || M2SlidingWindowSize < 0 ||
        M2SlidingWindowSize%%1!=0 ||
        M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
        M2CutoffQualityScore%%1!=0) {
        trimmedStartPos = NULL
        trimmedFinishPos = NULL
    } else {
        for (i in 1:(rawSeqLength-M2SlidingWindowSize+1)) {
            meanSLidingWindow <-
                mean(qualityBaseScore[i:(i+M2SlidingWindowSize-1)])
            if (meanSLidingWindow < qualityPbCutoff) {
                remainingIndex <- c(remainingIndex, i)
                # or ==> i + floor(M2SlidingWindowSize/3)
            }
        }
        trimmedStartPos = remainingIndex[1]
        trimmedFinishPos = remainingIndex[length(remainingIndex)]
        if (is.null(trimmedStartPos) || is.null(trimmedFinishPos)) {
            trimmedStartPos <- 0
            trimmedFinishPos <- 0
            trimmedQualityPhredScore <- 0
            trimmedMeanQualityScore <- 0
            trimmedMinQualityScore <- 0
            trimmedSeqLength = 0
            remainingRatio = 0
        } else {
            trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
            trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
            trimmedMinQualityScore <- min(trimmedQualityPhredScore)
            trimmedSeqLength = trimmedFinishPos - trimmedStartPos + 1
            remainingRatio = trimmedSeqLength / rawSeqLength
        }
    }
    return(list("rawSeqLength" = rawSeqLength,
                "rawMeanQualityScore" = rawMeanQualityScore,
                "rawMinQualityScore" = rawMinQualityScore,
                "trimmedStartPos" = trimmedStartPos,
                "trimmedFinishPos" = trimmedFinishPos,
                "trimmedSeqLength" = trimmedSeqLength,
                "trimmedMeanQualityScore" = trimmedMeanQualityScore,
                "trimmedMinQualityScore" = trimmedMinQualityScore,
                "remainingRatio" = remainingRatio))
}

### ============================================================================
### Calculating consensus read for one read set.
### ============================================================================
calculateConsensusRead <- function(forwardReadsList, reverseReadsList,
                                   refAminoAcidSeq, minFractionCall,
                                   maxFractionLost, geneticCode,
                                   acceptStopCodons, readingFrame,
                                   processorsNum) {
    ### --------------------------------------------------------------------
    ### forward & reverse character reads list string creation
    ### --------------------------------------------------------------------

    trimmedStartPos <- forwardReadsList[[1]]@QualityReport@trimmedStartPos
    trimmedFinishPos <- forwardReadsList[[1]]@QualityReport@trimmedFinishPos
    primaryDNA <- as.character(primarySeq(forwardReadsList[[1]]))[1]
    substr(primaryDNA, trimmedStartPos, trimmedFinishPos)


    fRDNAStringSet <- sapply(forwardReadsList, function(forwardRead) {
        trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
        trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
        primaryDNA <- as.character(primarySeq(forwardRead))
        substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
        # Without trimming
        # as.character(primarySeq(forwardRead))
    })
    rRDNAStringSet <- sapply(reverseReadsList, function(reverseRead) {
        trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
        trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
        primaryDNA <- as.character(primarySeq(reverseRead))
        substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
        # Without trimming
        # as.character(primarySeq(reverseRead))
    })

    ### --------------------------------------------------------------------
    ### DNAStringSet storing forward & reverse reads ! (Origin)
    ### --------------------------------------------------------------------
    frReadSet <- DNAStringSet(c(unlist(fRDNAStringSet),
                                unlist(rRDNAStringSet)))
    frReadFeatureList <- c(rep("Forward Reads", length(fRDNAStringSet)),
                           rep("Reverse Reads", length(rRDNAStringSet)))

    if(length(frReadSet) < 2) {
        error <- paste("\n'Valid abif files should be more than 2.\n",
                       sep = "")
        stop(error)
    }
    processorsNum <- getProcessors(processorsNum)

    ### --------------------------------------------------------------------
    ### Amino acid reference sequence CorrectFrameshifts correction
    ### --------------------------------------------------------------------
    if (refAminoAcidSeq != "") {
        message("Correcting frameshifts in reads using amino acid",
                "reference sequence")
        # My test refAminoAcidSeq data
        # no_N_string <- str_replace_all(frReadSet[1], "N", "T")
        # example.dna <- DNAStringSet(c(`IGHV1-18*01`=no_N_string))
        # refAminoAcidSeq <- translate(example.dna)
        # Verbose should be FALSE, but I get error when calling it
        corrected =
            CorrectFrameshifts(myXStringSet = frReadSet,
                               myAAStringSet = AAStringSet(refAminoAcidSeq),
                               geneticCode = geneticCode,
                               type = 'both',
                               processors = processorsNum)
        frReadSet = corrected$sequences
        indels = getIndelDf(corrected$indels)
        stops = as.numeric(unlist(mclapply(frReadSet, countStopSodons,
                                           readingFrame, geneticCode,
                                           mc.cores = processorsNum)))
        stopsDf = data.frame("read" = names(frReadSet),
                             "stop.codons" = stops)
        frReadSetLen = unlist(lapply(frReadSet, function(x) length(x)))
        frReadSet = frReadSet[which(frReadSetLen>0)]
    } else {
        indels = data.frame()
        stopsDf = data.frame()
    }
    if(length(frReadSet) < 2) {
        error <- paste("\n'After running 'CorrectFrameshifts' function, ",
                       "forward and reverse reads should be more than 2.\n",
                       sep = "")
        stop(error)
    }

    ### --------------------------------------------------------------------
    ### Reads with stop codons elimination
    ### --------------------------------------------------------------------
    ### ----------------------------------------------------------------
    ### Remove reads with stop codons
    ### ----------------------------------------------------------------
    if (!acceptStopCodons) {
        print("Removing reads with stop codons")
        if(refAminoAcidSeq == ""){ # otherwise we already did it above
            stops =
                as.numeric(unlist(mclapply(frReadSet,
                                           countStopSodons,
                                           readingFrame, geneticCode,
                                           mc.cores = processorsNum)))
            stopsDf = data.frame("read" = names(frReadSet),
                                 "stopCodons" = stops)
        }
        old_length = length(frReadSet)
        frReadSet = frReadSet[which(stops==0)]
        # Modify
        message(old_length - length(frReadSet),
                "reads with stop codons removed")
    }

    if(length(frReadSet) < 2) {
        error <- paste("\n'After removing reads with stop codons, ",
                       "forward and reverse reads should be more than 2.\n",
                       sep = "")
        stop(error)
    }

    ### --------------------------------------------------------------------
    ### Start aligning reads
    ### --------------------------------------------------------------------
    if (refAminoAcidSeq != "") {
        aln = AlignTranslation(frReadSet, geneticCode = geneticCode,
                               processors = processorsNum, verbose = FALSE)
    } else {
        aln = AlignSeqs(frReadSet,
                        processors = processorsNum, verbose = FALSE)
    }
    names(aln) = paste(1:length(aln), "Read",
                       basename(names(aln)), sep="_")
    consensus = ConsensusSequence(aln,
                                  minInformation = minFractionCall,
                                  includeTerminalGaps = TRUE,
                                  ignoreNonBases = TRUE,
                                  threshold = maxFractionLost,
                                  noConsensusChar = "-",
                                  ambiguity = TRUE
    )[[1]]

    diffs = mclapply(aln, nPairwiseDiffs,
                     subject = consensus, mc.cores = processorsNum)
    diffs = do.call(rbind, diffs)
    diffsDf = data.frame("name" = names(aln),
                         "pairwise.diffs.to.consensus" = diffs[,1],
                         "unused.chars" = diffs[,2])
    rownames(diffsDf) = NULL

    # get a dendrogram
    dist = DistanceMatrix(aln, correction = "Jukes-Cantor",
                          penalizeGapLetterMatches = FALSE,
                          processors = processorsNum, verbose = FALSE)
    dend = IdClusters(dist, type = "both",
                      showPlot = FALSE,
                      processors = processorsNum, verbose = FALSE)

    # add consensus to alignment
    aln2 = c(aln, DNAStringSet(consensus))
    names(aln2)[length(aln2)] = "Consensus"
    # strip gaps from consensus (must be an easier way!!)
    consensusGapfree = RemoveGaps(DNAStringSet(consensus))[[1]]

    # count columns in the alignment with >1 coincident secondary peaks
    spDf = countCoincidentSp(aln, processors = processorsNum)
    if (is.null(spDf)) {
        spDf = data.frame()
    }
    return(list("consensusGapfree" = consensusGapfree,
                "diffsDf"          = diffsDf,
                "aln2"             = aln2,
                "dist"             = dist,
                "dend"             = dend,
                "indels"           = indels,
                "stopsDf"          = stopsDf,
                "spDf"             = spDf))
}

### ============================================================================
### Aligning consensus reads into a new consensus read for all reads
### ============================================================================
alignConsensusReads <- function(SangerConsensusReadList,
                                geneticCode, refAminoAcidSeq,
                                minFractionCallSCSet, maxFractionLostSCSet,
                                processorsNum) {
    ### --------------------------------------------------------------------
    ### Creating SangerConsensusReadList DNAStringSet
    ### --------------------------------------------------------------------
    SangerConsensusReadDNAList <-
        sapply(SangerConsensusReadList, function(SangerConsensusRead) {
            as.character(SangerConsensusRead@consensusRead)
        })

    SangerConsensusReadDNASet <- DNAStringSet(SangerConsensusReadDNAList)

    ### --------------------------------------------------------------------
    ### Aligning consensus reads
    ### --------------------------------------------------------------------
    if(length(SangerConsensusReadDNASet) > 1) {
        message("Aligning consensus reads ... ")
        if(!is.null(refAminoAcidSeq)){
            aln = AlignTranslation(SangerConsensusReadDNASet,
                                   geneticCode = geneticCode,
                                   processors = processorsNum,
                                   verbose = FALSE)
        }else{
            aln = AlignSeqs(SangerConsensusReadDNASet,
                            processors = processorsNum,
                            verbose = FALSE)
        }

        # Making a rough NJ tree. Labels are rows in the summary df
        neat.labels = match(names(aln),
                            as.character(names(SangerConsensusReadDNASet))
        )
        aln2 = aln
        names(aln2) = neat.labels


        aln.bin = as.DNAbin(aln2)

        aln.dist = dist.dna(aln.bin, pairwise.deletion = TRUE)

        # Making a rough NJ tree. Labels are rows in the summary df
        #    (If tree cannot be created ==> NULL)
        aln.tree = NULL
        try({
            aln.tree = bionjs(aln.dist)
            aln.tree$tip.label <- names(aln)
            # deal with -ve branches
            # This is not necessarily accurate, but it is good enough to judge seuqences using the tree
            aln.tree$edge.length[which(aln.tree$edge.length<0)] = abs(aln.tree$edge.length[which(aln.tree$edge.length<0)])            },
            silent = TRUE
        )

        # Get consensus read and add to alignment result
        consensus = ConsensusSequence(aln,
                                      minInformation = minFractionCallSCSet,
                                      includeTerminalGaps = TRUE,
                                      ignoreNonBases = TRUE,
                                      threshold = maxFractionLostSCSet,
                                      noConsensusChar = "-",
                                      ambiguity = TRUE
        )[[1]]
    } else {
        aln = NULL
        aln.tree = NULL
    }
    return(list("consensus" = consensus,
                "aln" = aln,
                "aln.tree" = aln.tree))
}
### ============================================================================
### Adding dynamic menu to sidebar.
### ============================================================================
dynamicMenuSideBarSC <- function(input, output, session,
                               SangerSingleReadNum, SangerSingleReadFeature) {
    output$singleReadMenu <- renderMenu({
        menu_list <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(text = SangerSingleReadFeature[i],
                          tabName = SangerSingleReadFeature[i],
                          selected = TRUE, icon = icon("minus")))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Sanger Consensus Read Overview")})
}

dynamicMenuSideBarSCSet <- function(input, output, session, SangerCSetParam) {
    output$singleReadMenu <- renderMenu({
        SangerCSNum <- length(SangerCSetParam)
        menu_list <- sapply(1:SangerCSNum, function(i) {
            SangerSingleReadNum <- SangerCSetParam[[i]]$SangerSingleReadNum
            SangerCSMenuSubItem <- sapply(1:SangerSingleReadNum, function(j) {
                list(menuSubItem(text = SangerCSetParam[[i]]$SangerSingleReadFeature[[j]],
                                 tabName = paste0(i, " Consensus Read - ", SangerCSetParam[[i]]$SangerSingleReadFeature[[j]])))
            })
            SangerCSMenuSubItem <- c(list(menuSubItem(text = paste(SangerCSetParam[[i]]$SCName, "Overview"),
                                                      tabName = paste0(i, " Sanger Consensus Read Overview"), icon = icon("align-left"))),
                                     SangerCSMenuSubItem)
            SangerCSetParam[[i]]$SCName
            list(menuItem(text = SangerCSetParam[[i]]$SCName,
                          tabName = SangerCSetParam[[i]]$SCName,
                          icon = icon("minus"), SangerCSMenuSubItem))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Sanger Aligned Consensus Set Overview")})
}

### ============================================================================
### observeEvent: Adding dynamic rightHeader text
### ============================================================================
observeEventDynamicHeaderSC <- function(input, output, session, trimmedRV,
                                           SangerSingleReadQualReport) {
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
            trimmedRV[["trimmedStartPos"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmedStartPos
            trimmedRV[["trimmedFinishPos"]] <-
                SangerSingleReadQualReport[[
                    strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
            qualityPhredScores = SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityPhredScores

            readLen = length(qualityPhredScores)
            trimmedRV[["trimmedSeqLength"]] <- trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1
            trimmedRV[["remainingRatio"]] <- round(((trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1) / readLen) * 100, digits = 2)
        }
    })
}






















# observeEventDynamicHeaderSCSet <- function(input, output, session, trimmedRV,
#                                            SangerCSetParam) {
#
#     output$res <- renderText({
#         paste("You've selected:", input$sidebar_menu)
#     })
#
#     observeEvent(input$sidebar_menu, {
#         menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
#         html("rightHeader", menuItem)
#         sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
#         # message("strtoi(sidebar_menu[[1]]): ", strtoi(sidebar_menu[[1]]))
#         if (!is.na(suppressWarnings(as.numeric(sidebar_menu[[1]])))) {
#         #     trimmedRV[["trimmedStartPos"]] <-
#         #         SangerSingleReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedStartPos
#         #     trimmedRV[["trimmedFinishPos"]] <-
#         #         SangerSingleReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
#         #     qualityPhredScores = SangerSingleReadQualReport[[
#         #         strtoi(sidebar_menu[[1]])]]@qualityPhredScores
#         #
#         #     readLen = length(qualityPhredScores)
#         #     trimmedRV[["trimmedSeqLength"]] <- trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1
#         #     trimmedRV[["remainingRatio"]] <- round(((trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1) / readLen) * 100, digits = 2)
#         }
#     })
# }

### ============================================================================
### valueBox: SCMinReadsNum
### ============================================================================
valueBoxSCMinReadsNum <- function(input, output, SCMinReadsNum, session) {
    output$SCMinReadsNum <- renderUI({
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
    })
}

### ============================================================================
### valueBox: SCMinReadLength
### ============================================================================
valueBoxSCMinReadLength <- function(input, output, SCMinReadLength, session) {
    output$SCMinReadLength <- renderUI({
        valueBox(
            subtitle = tags$p("MinReadLength",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCMinReadLength),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMinFractionCall
### ============================================================================
valueBoxSCMinFractionCall <- function(input, output, SCMinFractionCall, session) {
    output$SCMinFractionCall <- renderUI({
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SCMinFractionCall),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMaxFractionLost
### ============================================================================
valueBoxSCMaxFractionLost <- function(input, output, SCMaxFractionLost, session) {
    output$SCMaxFractionLost <- renderUI({
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SCMaxFractionLost),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCAcceptStopCodons
### ============================================================================
valueBoxSCAcceptStopCodons <- function(input, output, SCAcceptStopCodons, session) {
    output$SCAcceptStopCodons <- renderUI({
        valueBox(
            subtitle = tags$p("AcceptStopCodons",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(SCAcceptStopCodons,
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCReadingFrame
### ============================================================================
valueBoxSCReadingFrame <- function(input, output, SCReadingFrame, session) {
    output$SCReadingFrame <- renderUI({
        valueBox(
            subtitle = tags$p("ReadingFrame",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCReadingFrame),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMinReadsNum
### ============================================================================
valueBoxSCMinReadsNumCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SCMinReadsNum <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])

        SCMinReadsNum <- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@minReadsNum
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
    })
}

### ============================================================================
### valueBox: SCMinReadLength
### ============================================================================
valueBoxSCMinReadLengthCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SCMinReadLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        SCMinReadLength <- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@minReadLength
        valueBox(
            subtitle = tags$p("MinReadLength",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCMinReadLength),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMinFractionCall
### ============================================================================
valueBoxSCMinFractionCallCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SCMinFractionCall <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        SCMinFractionCall <- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@minFractionCall
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SCMinFractionCall),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCMaxFractionLost
### ============================================================================
valueBoxSCMaxFractionLostCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SCMaxFractionLost <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        SCMaxFractionLost <- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@maxFractionLost
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SCMaxFractionLost),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCAcceptStopCodons
### ============================================================================
valueBoxSCAcceptStopCodonsCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SCAcceptStopCodons <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        SCAcceptStopCodons <- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@acceptStopCodons
        valueBox(
            subtitle = tags$p("AcceptStopCodons",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(SCAcceptStopCodons,
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SCReadingFrame
### ============================================================================
valueBoxSCReadingFrameCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SCReadingFrame <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        consensusReadIndex <- strtoi(sidebar_menu[[1]])
        SCReadingFrame <- SangerConsensusSet@consensusReadsList[[consensusReadIndex]]@readingFrame
        valueBox(
            subtitle = tags$p("ReadingFrame",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SCReadingFrame),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

################################################################################
### Each Read
################################################################################
### ============================================================================
### valueBox: Change M1TrimmingCutoff
### ============================================================================
valueBoxM1TrimmingCutoff <- function(input, output, session) {
    output$M1TrimmingCutoff <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        singleReadIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
            as.numeric(input$M1TrimmingCutoffText) > 0 &&
            as.numeric(input$M1TrimmingCutoffText) <= 1) {
            inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
        } else {
            inputM1TrimmingCutoffText <- 0.0001
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
}

### ============================================================================
### valueBox: Change M2CutoffQualityScore
### ============================================================================
valueBoxM2CutoffQualityScore <- function(input, output, session) {
    output$M2CutoffQualityScore <- renderUI({
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
}

### ============================================================================
### valueBox: Change M2SlidingWindowSize
### ============================================================================
valueBoxM2SlidingWindowSize <- function(input, output, session) {
    output$M2SlidingWindowSize <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
            strtoi(input$M2SlidingWindowSizeText) > 0 &&
            strtoi(input$M2SlidingWindowSizeText) <= 20 &&
            strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
            inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
        } else {
            inputM2SlidingWindowSizeText <- 5
        }
        valueBox(
            # strtoi(input$M2CutoffQualityScoreText
            subtitle = tags$p("Sliding Window Size ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(inputM2SlidingWindowSizeText),
                           style = "font-size: 29px;"),
            icon = icon("expand", "fa-sm"),
            color = "olive", width = 10,
        )
    })
}

### ============================================================================
### valueBox: Change rawSeqLength
### ============================================================================
valueBoxRawSeqLength <- function(input, output, session, trimmedRV) {
    output$rawSeqLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Raw Seqence Len ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["rawSeqLength"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change rawMeanQualityScore
### ============================================================================
valueBoxRawMeanQualityScore <- function(input, output, session, trimmedRV) {
    output$rawMeanQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Raw Mean Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(round(trimmedRV[["rawMeanQualityScore"]], 2)),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change rawMinQualityScore
### ============================================================================
valueBoxRawMinQualityScore <- function(input, output, session, trimmedRV) {
    output$rawMinQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Raw Min Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["rawMinQualityScore"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedStartPos
### ============================================================================
valueBoxTrimmedStartPos <- function(input, output, session, trimmedRV) {
    output$trimmedStartPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimming Start Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedStartPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedFinishPos
### ============================================================================
valueBoxTrimmedFinishPos <- function(input, output, session, trimmedRV) {
    output$trimmedFinishPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimming End Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedFinishPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("times-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedSeqLength
### ============================================================================
valueBoxTrimmedSeqLength <- function(input, output, session, trimmedRV) {
    output$trimmedSeqLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimmed Seqence Length ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedSeqLength"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedMeanQualityScore
### ============================================================================
valueBoxTrimmedMeanQualityScore <- function(input, output, session, trimmedRV) {
    output$trimmedMeanQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimmed Mean Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(round(trimmedRV[["trimmedMeanQualityScore"]], 2)),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedMinQualityScore
### ============================================================================
valueBoxTrimmedMinQualityScore <- function(input, output, session, trimmedRV) {
    output$trimmedMinQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimmed Min Quality Score ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedMinQualityScore"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change remainingRatio
### ============================================================================
valueBoxRemainingRatio <- function(input, output, session, trimmedRV) {
    output$remainingRatio <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Remaining Ratio",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(paste(trimmedRV[["remainingRatio"]], "%"),
                           style = "font-size: 32px;"),
            icon = icon("divide", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedStartPos
### ============================================================================
valueBoxChromTrimmedStartPos <- function(input, output, session, trimmedRV) {
    output$ChromatogramtrimmedStartPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimming Start Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedStartPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("check-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### valueBox: Change trimmedFinishPos
### ============================================================================
valueBoxChromTrimmedFinishPos <- function(input, output, session, trimmedRV) {
    output$ChromatogramtrimmedFinishPos <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        valueBox(
            subtitle = tags$p("Trimming End Pos ",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(toString(trimmedRV[["trimmedFinishPos"]]),
                           style = "font-size: 29px;"),
            icon = icon("times-circle", "fa-sm"),
            color = "olive", width = 12,
        )
    })
}

### ============================================================================
### qualityTrimmingRatioPlot
### ============================================================================
qualityTrimmingRatioPlot <- function(input, output, session, trimmedRV,
                                     SangerSingleReadQualReport,
                                     SangerSingleReadFeature) {
    output$qualityTrimmingRatioPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(strtoi(sidebar_menu[[1]])) &&
            !is.null(SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]])) {
            readFeature <- SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]]
            trimmedStartPos = trimmedRV[["trimmedStartPos"]]
            trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
            qualityPhredScores = SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityPhredScores
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
            # Change font setting
            # f <- list(
            #     family = "Courier New, monospace",
            #     size = 18,
            #     color = "#7f7f7f"
            # )
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
        }
    })
}

### ============================================================================
### qualityQualityBasePlot
### ============================================================================
qualityQualityBasePlot <- function(input, output, session, trimmedRV,
                                   SangerSingleReadQualReport,
                                   SangerSingleReadFeature) {
    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (!is.na(strtoi(sidebar_menu[[1]])) &&
            !is.null(SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]])) {
            readFeature <- SangerSingleReadFeature[[strtoi(sidebar_menu[[1]])]]
            trimmedStartPos = trimmedRV[["trimmedStartPos"]]
            trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
            qualityPhredScores = SangerSingleReadQualReport[[
                strtoi(sidebar_menu[[1]])]]@qualityPhredScores
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
                                     round((trimmedFinishPos - trimmedStartPos+1)/
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
                                     xanchor = "center",  # use center of legend as anchor
                                     x = 0.5, y = 1.1)) %>%
                # add_segments(x = trimmedStartPos, xend = trimmedFinishPos, y = 70, yend = 70, inherit = TRUE, width = 10, line = list(width = 8)) %>%
                # add_segments(x = 0, xend = readLen, y = 75, yend = 75, inherit = TRUE, width = 4, line = list(width = 8)) %>%
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
        }
    })
}

### ============================================================================
### chromatogram row number counting
### ============================================================================
# chromatogramRowNum <- function(obj, width) {
chromatogramRowNum <- function(width, rawLength, trimmedLength, showTrimmed) {
    if (showTrimmed) {
        numplots = ceiling(rawLength / width)
    } else {
        numplots = ceiling(trimmedLength / width)
    }
}


vline <- function(x = 0, color = "red") {
    list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = x,
        x1 = x,
        line = list(color = color)
    )
}

