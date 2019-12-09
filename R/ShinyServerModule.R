M1inside_calculate_trimming <- function(qualityPhredScores,
                                        qualityBaseScores,
                                        M1TrimmingCutoff) {
    rawSeqLength <- length(qualityBaseScores)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    start = FALSE
    trimmedStartPos = 0
    qualityBaseScoresCutOff = M1TrimmingCutoff - qualityBaseScores
    ### ------------------------------------------------------------------------
    ### calculate cummulative score
    ### if cumulative value < 0, set it to 0
    ### the BioPython implementation always trims the first base,
    ### this implementation does not.
    ### ------------------------------------------------------------------------
    score = qualityBaseScoresCutOff[1]
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
    for(i in 2:length(qualityBaseScoresCutOff)){
        score = cummul_score[length(cummul_score)] + qualityBaseScoresCutOff[i]
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
    if (trimmedFinishPos - trimmedStartPos == 0) {
        trimmedStartPos = 1
        trimmedFinishPos = 2
    }
    trimmedSeqLength = trimmedFinishPos - trimmedStartPos
    trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
    trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
    trimmedMinQualityScore <- min(trimmedQualityPhredScore)
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
                                      qualityBaseScores,
                                      M2CutoffQualityScore,
                                      M2SlidingWindowSize) {
    rawSeqLength <- length(qualityBaseScores)
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
                mean(qualityBaseScores[i:(i+M2SlidingWindowSize-1)])
            if (meanSLidingWindow < qualityPbCutoff) {
                remainingIndex <- c(remainingIndex, i)
                # or ==> i + floor(M2SlidingWindowSize/3)
            }
        }
        trimmedStartPos = remainingIndex[1]
        trimmedFinishPos = remainingIndex[length(remainingIndex)]
        if (is.null(trimmedStartPos) || is.null(trimmedFinishPos)) {
            trimmedStartPos <- 1
            trimmedFinishPos <- 2
        }
        trimmedQualityPhredScore <- qualityPhredScores[trimmedStartPos:trimmedFinishPos]
        trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
        trimmedMinQualityScore <- min(trimmedQualityPhredScore)
        trimmedSeqLength = trimmedFinishPos - trimmedStartPos
        remainingRatio = trimmedSeqLength / rawSeqLength
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
calculateContigSeq <- function(forwardReadsList, forwardReadList,
                                   refAminoAcidSeq, minFractionCall,
                                   maxFractionLost, geneticCode,
                                   acceptStopCodons, readingFrame,
                                   processorsNum) {
    ### ------------------------------------------------------------------------
    ### forward & reverse character reads list string creation
    ### ------------------------------------------------------------------------
    fRDNAStringSet <- sapply(forwardReadsList, function(forwardRead) {
        trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
        trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
        primaryDNA <- as.character(forwardRead@primarySeq)
        substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
    })
    rRDNAStringSet <- sapply(forwardReadList, function(reverseRead) {
        trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
        trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
        primaryDNA <- as.character(reverseRead@primarySeq)
        substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
    })

    ### --------------------------------------------------------------------
    ### DNAStringSet storing forward & reverse reads ! (Origin)
    # ### --------------------------------------------------------------------
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
alignContigs <- function(SangerContigList,
                                geneticCode, refAminoAcidSeq,
                                minFractionCallSA, maxFractionLostSA,
                                processorsNum) {
    ### --------------------------------------------------------------------
    ### Creating SangerContigList DNAStringSet
    ### --------------------------------------------------------------------
    SangerContigDNAList <-
        sapply(SangerContigList, function(SangerContig) {
            as.character(SangerContig@contigSeq)
        })

    SangerContigDNASet <- DNAStringSet(SangerContigDNAList)

    ### --------------------------------------------------------------------
    ### Aligning consensus reads
    ### --------------------------------------------------------------------
    if(length(SangerContigDNASet) > 1) {
        message("Aligning consensus reads ... ")
        if(refAminoAcidSeq != ""){
            aln = AlignTranslation(SangerContigDNASet,
                                   geneticCode = geneticCode,
                                   processors = processorsNum,
                                   verbose = FALSE)
        }else{
            aln = AlignSeqs(SangerContigDNASet,
                            processors = processorsNum,
                            verbose = FALSE)
        }

        # Making a rough NJ tree. Labels are rows in the summary df
        neat.labels = match(names(aln),
                            as.character(names(SangerContigDNASet))
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
                                      minInformation = minFractionCallSA,
                                      includeTerminalGaps = TRUE,
                                      ignoreNonBases = TRUE,
                                      threshold = maxFractionLostSA,
                                      noConsensusChar = "-",
                                      ambiguity = TRUE
        )[[1]]
    } else {
        aln = NULL
        aln.tree = NULL
    }
    return(list("consensus" = consensus,
                "aln"       = aln,
                "aln.tree"  = aln.tree))
}
### ============================================================================
### Adding dynamic menu to sidebar.
### ============================================================================
dynamicMenuSideBarSA <- function(input, output, session,
                                  forwardReadNum, reverseReadNum,
                                  forwardReadFeature, reverseReadFeature) {
    output$singleReadMenu <- renderMenu({

        fmenuSub_list <- sapply(1:forwardReadNum, function(i) {
            list(menuSubItem(text = forwardReadFeature[i],
                          tabName = forwardReadFeature[i],
                          icon = icon("minus")))
        })
        rmenuSub_list <- sapply(1:reverseReadNum, function(i) {
            list(menuSubItem(text = reverseReadFeature[i],
                          tabName = reverseReadFeature[i],
                          icon = icon("minus")))
        })
        fmenu_list <- menuItem(text = "Forward Reads",
                               tabName = "forwardReads",
                               icon = icon("minus"), fmenuSub_list)

        rmenu_list <- menuItem(text = "Reverse Reads",
                               tabName = "reverseReads",
                               icon = icon("minus"), rmenuSub_list)
        sidebarMenu(.list = list(fmenu_list, rmenu_list))
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu",
                            "Sanger Contig Overview")})
}

dynamicMenuSideBarSA <- function(input, output, session, SangerAlignmentParam) {
    output$singleReadMenu <- renderMenu({
        SangerCSNum <- length(SangerAlignmentParam)
        menu_list <- sapply(1:SangerCSNum, function(i) {
            forwardReadNum <- SangerAlignmentParam[[i]]$forwardReadNum
            reverseReadNum <- SangerAlignmentParam[[i]]$reverseReadNum

            forwardReadFeature <- SangerAlignmentParam[[i]]$forwardReadFeature
            reverseReadFeature <- SangerAlignmentParam[[i]]$reverseReadFeature
            fmenuSub_list <- sapply(1:forwardReadNum, function(j) {
                list(menuSubItem(text = paste0(i, " CR - ", forwardReadFeature[j]),
                                 tabName = paste0(i, " CR - ", forwardReadFeature[j]),
                                 icon = icon("minus")))
            })
            rmenuSub_list <- sapply(1:reverseReadNum, function(j) {
                list(menuSubItem(text = paste0(i, " CR - ", reverseReadFeature[j]),
                                 tabName = paste0(i, " CR - ", reverseReadFeature[j]),
                                 icon = icon("minus")))
            })
            fmenu_list <- menuItem(text = "Forward Reads",
                                   tabName = "forwardReads",
                                   icon = icon("minus"), fmenuSub_list)

            rmenu_list <- menuItem(text = "Reverse Reads",
                                   tabName = "reverseReads",
                                   icon = icon("minus"), rmenuSub_list)
            SangerCSMenuSubItem <- list(fmenu_list, rmenu_list)

            SangerCSMenuSubItem <- c(list(menuSubItem(text = paste(SangerAlignmentParam[[i]]$SCName, "Overview"),
                                                      tabName = paste0(i, " Sanger Consensus Read Overview"),
                                                      icon = icon("align-left"))),
                                     SangerCSMenuSubItem)

            SangerAlignmentParam[[i]]$SCName
            list(menuItem(text = SangerAlignmentParam[[i]]$SCName,
                          tabName = SangerAlignmentParam[[i]]$SCName,
                          icon = icon("minus"), SangerCSMenuSubItem))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Sanger Alignment Overview")})
}

### ============================================================================
### observeEvent: Adding dynamic rightHeader text
### ============================================================================
observeEventDynamicHeaderSA <- function(input, output, session, trimmedRV) {
    observeEvent(input$sidebar_menu, {
        menuItem <- switch(input$sidebar_menu, input$sidebar_menu)
        message("menuItem: ", menuItem)
        html("rightHeader", menuItem)
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
    })
}






















# observeEventDynamicHeaderSA <- function(input, output, session, trimmedRV,
#                                            SangerAlignmentParam) {
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
#         #         SangerReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedStartPos
#         #     trimmedRV[["trimmedFinishPos"]] <-
#         #         SangerReadQualReport[[
#         #             strtoi(sidebar_menu[[1]])]]@trimmedFinishPos
#         #     qualityPhredScores = SangerReadQualReport[[
#         #         strtoi(sidebar_menu[[1]])]]@qualityPhredScores
#         #
#         #     readLen = length(qualityPhredScores)
#         #     trimmedRV[["trimmedSeqLength"]] <- trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1
#         #     trimmedRV[["remainingRatio"]] <- round(((trimmedRV[["trimmedFinishPos"]] - trimmedRV[["trimmedStartPos"]] + 1) / readLen) * 100, digits = 2)
#         }
#     })
# }

### ============================================================================
### valueBox: SAMinReadsNum
### ============================================================================
valueBoxSAMinReadsNum <- function(input, output, SAMinReadsNum, session) {
    output$SAMinReadsNum <- renderUI({
        valueBox(
            subtitle = tags$p("MinReadsNum",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SAMinReadsNum),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SAMinReadLength
### ============================================================================
valueBoxSAMinReadLength <- function(input, output, SAMinReadLength, session) {
    output$SAMinReadLength <- renderUI({
        valueBox(
            subtitle = tags$p("MinReadLength",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SAMinReadLength),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SAMinFractionCall
### ============================================================================
valueBoxSAMinFractionCall <- function(input, output, SAMinFractionCall, session) {
    output$SAMinFractionCall <- renderUI({
        valueBox(
            subtitle = tags$p("MinFractionCall",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SAMinFractionCall),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SAMaxFractionLost
### ============================================================================
valueBoxSAMaxFractionLost <- function(input, output, SAMaxFractionLost, session) {
    output$SAMaxFractionLost <- renderUI({
        valueBox(
            subtitle = tags$p("MaxFractionLost",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(as.numeric(SAMaxFractionLost),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SAAcceptStopCodons
### ============================================================================
valueBoxSAAcceptStopCodons <- function(input, output, SAAcceptStopCodons, session) {
    output$SAAcceptStopCodons <- renderUI({
        valueBox(
            subtitle = tags$p("AcceptStopCodons",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(SAAcceptStopCodons,
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SAReadingFrame
### ============================================================================
valueBoxSAReadingFrame <- function(input, output, SAReadingFrame, session) {
    output$SAReadingFrame <- renderUI({
        valueBox(
            subtitle = tags$p("ReadingFrame",
                              style = "font-size: 15px;
                                            font-weight: bold;"),
            value = tags$p(strtoi(SAReadingFrame),
                           style = "font-size: 29px;"),
            icon = icon("cut", "fa-sm"),
            color = "olive",
            width = 12,
        )
    })
}

### ============================================================================
### valueBox: SAMinReadsNum
### ============================================================================
valueBoxSAMinReadsNumCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SAMinReadsNum <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMinReadsNum <- SangerConsensusSet@contigList[[contigIndex]]@minReadsNum
            valueBox(
                subtitle = tags$p("MinReadsNum",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(SAMinReadsNum),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}

### ============================================================================
### valueBox: SAMinReadLength
### ============================================================================
valueBoxSAMinReadLengthCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SAMinReadLength <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMinReadLength <- SangerConsensusSet@contigList[[contigIndex]]@minReadLength
            valueBox(
                subtitle = tags$p("MinReadLength",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(SAMinReadLength),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}

### ============================================================================
### valueBox: SAMinFractionCall
### ============================================================================
valueBoxSAMinFractionCallCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SAMinFractionCall <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMinFractionCall <- SangerConsensusSet@contigList[[contigIndex]]@minFractionCall
            valueBox(
                subtitle = tags$p("MinFractionCall",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(as.numeric(SAMinFractionCall),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}

### ============================================================================
### valueBox: SAMaxFractionLost
### ============================================================================
valueBoxSAMaxFractionLostCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SAMaxFractionLost <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAMaxFractionLost <- SangerConsensusSet@contigList[[contigIndex]]@maxFractionLost
            valueBox(
                subtitle = tags$p("MaxFractionLost",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(as.numeric(SAMaxFractionLost),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}

### ============================================================================
### valueBox: SAAcceptStopCodons
### ============================================================================
valueBoxSAAcceptStopCodonsCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SAAcceptStopCodons <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAAcceptStopCodons <- SangerConsensusSet@contigList[[contigIndex]]@acceptStopCodons
            valueBox(
                subtitle = tags$p("AcceptStopCodons",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(SAAcceptStopCodons,
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
    })
}

### ============================================================================
### valueBox: SAReadingFrame
### ============================================================================
valueBoxSAReadingFrameCSSet <- function(input, output, SangerConsensusSet, session) {
    output$SAReadingFrame <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SAReadingFrame <-
                SangerConsensusSet@
                contigList[[contigIndex]]@readingFrame
            valueBox(
                subtitle = tags$p("ReadingFrame",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(SAReadingFrame),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 12,
            )
        }
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
        readIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(readIndex))) {
            if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
                as.numeric(input$M1TrimmingCutoffText) > 0 &&
                as.numeric(input$M1TrimmingCutoffText) <= 1) {
                inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
            } else {
                inputM1TrimmingCutoffText <- 0.0001
            }
            valueBox(
                subtitle = tags$p("Cutoff Score",
                                  style = "font-size: 15px;
                                       font-weight: bold;"),
                value = tags$p(as.numeric(inputM1TrimmingCutoffText),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 10,
            )
        }
    })
}

### ============================================================================
### valueBox: Change M2CutoffQualityScore
### ============================================================================
valueBoxM2CutoffQualityScore <- function(input, output, session) {
    output$M2CutoffQualityScore <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(readIndex))) {
            if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
                strtoi(input$M2CutoffQualityScoreText) > 0 &&
                strtoi(input$M2CutoffQualityScoreText) <= 60 &&
                strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
                inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
            } else {
                inputM2CutoffQualityScoreText <- 20
            }
            valueBox(
                subtitle = tags$p("Cutoff Quality Score",
                                  style = "font-size: 15px;
                                            font-weight: bold;"),
                value = tags$p(strtoi(inputM2CutoffQualityScoreText),
                               style = "font-size: 29px;"),
                icon = icon("cut", "fa-sm"),
                color = "olive",
                width = 10,
            )
        }
    })
}

### ============================================================================
### valueBox: Change M2SlidingWindowSize
### ============================================================================
valueBoxM2SlidingWindowSize <- function(input, output, session) {
    output$M2SlidingWindowSize <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        readIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(strtoi(readIndex))) {
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
        }
    })
}

### ============================================================================
### valueBox: Change rawSeqLength
### ============================================================================
valueBoxRawSeqLength <- function(input, output, session, trimmedRV) {
    output$rawSeqLength <- renderUI({
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
qualityTrimmingRatioPlotDisplay <- function(input, output, session,
                                            trimmedRV, qualityPhredScores) {
    trimmedStartPos = trimmedRV[["trimmedStartPos"]]
    trimmedFinishPos = trimmedRV[["trimmedFinishPos"]]
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
    )
    PerDataPlot <- melt(PerData, id.vars = c("Base"))
    suppressPlotlyMessage(
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
            ))
}


### ============================================================================
### qualityQualityBasePlot
### ============================================================================
qualityQualityBasePlotDisplay <- function(input, output, session,
                                          trimmedRV, qualityPhredScores) {
    trimmedStartPos <- trimmedRV[["trimmedStartPos"]]
    trimmedFinishPos <- trimmedRV[["trimmedFinishPos"]]
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
    suppressPlotlyMessage(
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
            ))
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




### ============================================================================
### SangerRead Sequence Render Function (DNA / AA) dynamic
### ============================================================================
primarySeqDisplay <- function(sequenceParam) {
    primarySeq <- unlist(strsplit(
        sequenceParam[["primarySeq"]], ""))
    primarySeqDF <- data.frame(
        t(data.frame(primarySeq)), stringsAsFactors = FALSE)
    colnames(primarySeqDF) <- substr(colnames(primarySeqDF), 2, 100)
    rownames(primarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(primarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(primarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(primarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(primarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = primarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
secondarySeqDisplay <- function(sequenceParam) {
    secondarySeq <-unlist(strsplit(
        sequenceParam[["secondarySeq"]], ""))
    secondarySeqDF <- data.frame(
        t(data.frame(secondarySeq)), stringsAsFactors = FALSE)
    colnames(secondarySeqDF) <- substr(colnames(secondarySeqDF), 2, 100)
    rownames(secondarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(secondarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(secondarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(secondarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(secondarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = secondarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
qualityScoreDisplay <- function(PhredScore) {
    PhredScoreDF <- data.frame(
        t(data.frame(PhredScore)), stringsAsFactors = FALSE)
    colnames(PhredScoreDF) <- substr(colnames(PhredScoreDF), 2, 100)
    rownames(PhredScoreDF) <- NULL
    styleList <- SetAllStyleList(PhredScoreDF, "#ecffd9")
    suppressMessages(
        excelTable(data =
                       PhredScoreDF, defaultColWidth = 30, editable = FALSE,
                   rowResize = FALSE, columnResize = FALSE,
                   allowInsertRow = FALSE, allowInsertColumn = FALSE,
                   allowDeleteRow = FALSE, allowDeleteColumn = FALSE,
                   style = styleList, allowRenameColumn = FALSE,
                   loadingSpin = TRUE)
    )
}
PrimAASeqS1Display <- function(sequenceParam) {
    AAString <- data.frame(AAString(sequenceParam[["primaryAASeqS1"]]))
    AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
    colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
    rownames(AAStringDF) <- NULL
    width <- rep(90, length(AAStringDF))
    styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
    styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
    styleList <- c(styleList1, styleList2)
    suppressMessages(
        excelTable(data = AAStringDF, columns = data.frame(width = width),
                   defaultColWidth = 90, editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
PrimAASeqS2Display <- function(sequenceParam) {
    AAString <- data.frame(AAString(sequenceParam[["primaryAASeqS2"]]))
    AAString <- rbind(NA, AAString)
    AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
    colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
    rownames(AAStringDF) <- NULL
    width <- rep(90, length(AAStringDF) - 1)
    width <- c(30, width)
    styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
    styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
    styleList <- c(styleList1, styleList2)
    styleList[['A1']] <- 'background-color: black;'
    suppressMessages(
        excelTable(data = AAStringDF, columns = data.frame(width = width),
                   defaultColWidth = 90, editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
PrimAASeqS3Display <- function(sequenceParam) {
    AAString <- data.frame(AAString(sequenceParam[["primaryAASeqS3"]]))
    AAString <- rbind(NA, NA, AAString)
    AAStringDF <- data.frame(t(AAString), stringsAsFactors = FALSE)
    colnames(AAStringDF) <- substr(colnames(AAStringDF), 2, 100)
    rownames(AAStringDF) <- NULL
    width <- rep(90, length(AAStringDF) - 2)
    width <- c(30, 30, width)
    styleList1 <- SetAllStyleList(AAStringDF, "#ecffd9")
    styleList2 <- SetCharStyleList (AAStringDF, "*", "#cf0000")
    styleList <- c(styleList1, styleList2)
    styleList[['A1']] <- 'background-color: black;'
    styleList[['B1']] <- 'background-color: black;'
    suppressMessages(
        excelTable(data = AAStringDF, columns = data.frame(width = width),
                   defaultColWidth = 90, editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}
primarySeqTrimmedDisplay <- function(input, output, session,
                                     sequenceParam, trimmedRV) {
    primarySeq <- unlist(strsplit(
        substr(sequenceParam[["primarySeq"]],
               trimmedRV[["trimmedStartPos"]] + 1,
               trimmedRV[["trimmedFinishPos"]])
        , ""))
    primarySeqDF <- data.frame(
        t(data.frame(primarySeq)), stringsAsFactors = FALSE)
    if ((trimmedRV[["trimmedFinishPos"]]-trimmedRV[["trimmedStartPos"]]) == 1) {
        colnames(primarySeqDF) <- "1"
    } else {
        colnames(primarySeqDF) <- substr(colnames(primarySeqDF), 2, 100)
    }
    rownames(primarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(primarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(primarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(primarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(primarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = primarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}

secondSeqTrimmedDisplay <- function(input, output, session,
                                    sequenceParam, trimmedRV) {
    secondarySeq <- unlist(strsplit(
        substr(sequenceParam[["secondarySeq"]],
               trimmedRV[["trimmedStartPos"]] + 1,
               trimmedRV[["trimmedFinishPos"]])
        , ""))
    secondarySeqDF <- data.frame(
        t(data.frame(secondarySeq)), stringsAsFactors = FALSE)
    if ((trimmedRV[["trimmedFinishPos"]]-trimmedRV[["trimmedStartPos"]]) == 1) {
        colnames(secondarySeqDF) <- "1"
    } else {
        colnames(secondarySeqDF) <- substr(colnames(secondarySeqDF), 2, 100)
    }
    rownames(secondarySeqDF) <- NULL
    AstyleList <- SetCharStyleList(secondarySeqDF, "A", "#1eff00")
    TstyleList <- SetCharStyleList(secondarySeqDF, "T", "#ff7a7a")
    CstyleList <- SetCharStyleList(secondarySeqDF, "C", "#7ac3ff")
    GstyleList <- SetCharStyleList(secondarySeqDF, "G", "#c9c9c9")
    styleList <- c(AstyleList, TstyleList, CstyleList, GstyleList)
    suppressMessages(
        excelTable(data = secondarySeqDF, defaultColWidth = 30,
                   editable = FALSE, rowResize = FALSE,
                   columnResize = FALSE, allowInsertRow = FALSE,
                   allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                   allowDeleteColumn = FALSE, allowRenameColumn = FALSE,
                   style = styleList, loadingSpin = TRUE)
    )
}

