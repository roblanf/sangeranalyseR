### ============================================================================
### Global helper functions
### ============================================================================

### Phase 6: BPPARAM resolution.
###
### Translates a (processorsNum, BPPARAM) input pair into a single BPPARAM
### object suitable for `bplapply` and `bpnworkers`. Honours the historical
### contract that `processorsNum = 1` means "run serially":
###   * BPPARAM supplied -> use it (caller is responsible).
###   * processorsNum == 1 (or NULL on Windows) -> SerialParam().
###   * processorsNum >= 2 -> MulticoreParam(workers = N) on Unix,
###                          SnowParam(workers = N) on Windows.
###   * processorsNum == NULL on Unix -> bpparam() (registered default,
###     usually MulticoreParam with all cores).
.resolveBPPARAM <- function(processorsNum = NULL, BPPARAM = NULL) {
    if (!is.null(BPPARAM)) return(BPPARAM)
    on_windows <- identical(.Platform$OS.type, "windows")
    if (is.null(processorsNum)) {
        if (on_windows) return(BiocParallel::SerialParam())
        return(BiocParallel::bpparam())
    }
    if (!is.numeric(processorsNum) || length(processorsNum) != 1L ||
        processorsNum < 1L) {
        return(BiocParallel::SerialParam())
    }
    n <- as.integer(processorsNum)
    if (n == 1L) return(BiocParallel::SerialParam())
    if (on_windows) return(BiocParallel::SnowParam(workers = n))
    BiocParallel::MulticoreParam(workers = n)
}

### Back-compat shim: keep `getProcessors()` returning an integer count,
### now derived from a BPPARAM. External callers (e.g. DECIPHER's
### `processors=` argument) consume the integer; internal call sites can
### either keep using the integer or migrate to BPPARAM directly.
getProcessors <- function(processors = NULL) {
    if (!is.null(processors) && is.numeric(processors) &&
        length(processors) == 1L && processors >= 1L) {
        return(as.integer(processors))
    }
    BiocParallel::bpnworkers(.resolveBPPARAM(processors))
}

suppressPlotlyMessage <- function(p) {
    suppressMessages(plotly_build(p))
}

# <---------------------------------------------------------------------------->

### ============================================================================
### SangerAlignment related helper functions
### ============================================================================
### ----------------------------------------------------------------------------
### Aligning SangerContigs (SangerAlignment)
### ----------------------------------------------------------------------------
alignContigs <- function(SangerContigList, geneticCode, refAminoAcidSeq,
                         minFractionCallSA, maxFractionLostSA, processorsNum,
                         BPPARAM = NULL) {
    if (!is.null(BPPARAM)) processorsNum <- BiocParallel::bpnworkers(BPPARAM)
    ### ------------------------------------------------------------------------
    ### Creating SangerContigList DNAStringSet
    ### ------------------------------------------------------------------------
    SangerContigDNAList <-
        vapply(SangerContigList, function(SangerContig) {
            as.character(SangerContig@contigSeq)
        }, FUN.VALUE = character(1))
    SangerContigDNASet <- DNAStringSet(SangerContigDNAList)
    ### ------------------------------------------------------------------------
    ### Aligning consensus reads
    ### ------------------------------------------------------------------------
    if(length(SangerContigDNASet) > 1) {
        log_info("Aligning consensus reads ... ")
        if(refAminoAcidSeq != ""){
            aln = AlignTranslation(SangerContigDNASet,
                                   geneticCode = geneticCode,
                                   processors = processorsNum,
                                   verbose = FALSE)
        }else{
            log_info('Before building!!')
            aln = AlignSeqs(SangerContigDNASet, processors = processorsNum,
                            verbose = FALSE)
            log_info('After building!!')
        }
        # Making a rough NJ tree. Labels are rows in the summary df

        if (length(aln) > 2) {
            neat.labels = match(names(aln),
                                as.character(names(SangerContigDNASet)))
            aln2 = aln
            names(aln2) = neat.labels

            aln.bin = as.DNAbin(aln2)

            aln.dist = dist.dna(aln.bin, pairwise.deletion = TRUE)
            # Making a rough NJ tree. Labels are rows in the summary df
            #    (If tree cannot be created ==> NULL)
            aln.tree = read.tree(text="();")
            tryCatch({
                aln.tree = bionjs(aln.dist)
                aln.tree$tip.label <- names(aln)
                # deal with -ve branches
                # This is not necessarily accurate, but it is good enough to
                # judge seuquences using the tree
                aln.tree$edge.length[which(aln.tree$edge.length<0)] =
                    abs(aln.tree$edge.length[which(aln.tree$edge.length<0)])
            }, warning = function(warning_condition) {
                log_info("The number of contigs is less than 3 or quality of reads ",
                        "are too low. 'Contigs Tree' cannot be created.")
            }, error = function(error_condition) {
                log_info("The number of contigs is less than 3 or quality of reads ",
                        "are too low. 'Contigs Tree' cannot be created.")
                aln.tree = read.tree(text="();")
            })
        } else {
            log_info("The number of contigs is less than 3 or quality of reads ",
                    "are too low. 'Contigs Tree' cannot be created.")
            aln.tree = read.tree(text="();")
        }

        # Get consensus read and add to alignment result
        consensus = ConsensusSequence(aln,
                                      minInformation = minFractionCallSA,
                                      includeTerminalGaps = TRUE,
                                      threshold = maxFractionLostSA,
                                      noConsensusChar = "-",
                                      ambiguity = TRUE)[[1]]
    } else {
        consensus = NULL
        aln = NULL
        phyloT <- as.phylo(rtree(n = 2))
        phyloT$edge <- matrix(c(0,0,0,0), 2, 2)
        phyloT$tip.label <- NULL
        phyloT$edge.length <- NULL
        phyloT$Nnode <- NULL
        aln.tree = phyloT
    }
    return(list("consensus" = consensus,
                "aln"       = aln,
                "aln.tree"  = aln.tree))
}

# <---------------------------------------------------------------------------->

### ============================================================================
### SangerContig related helper functions
### ============================================================================
getIndelDf <- function(indelList){
    r = lapply(indelList, indelRow)
    indelDf = data.frame(matrix(unlist(r), byrow = TRUE, nrow = length(indelList)))
    indelDf = cbind(names(indelList), indelDf)
    names(indelDf) = c('read', 'insertions', 'deletions', 'distance')
    return(indelDf)
}
indelRow <- function(row){
    nIns = length(row$insertions)
    nDel = length(row$deletions)
    dist = row$distance
    return(c(nIns, nDel, dist))
}
nPairwiseDiffs <- function(pattern, subject){
    # pairwise differences assuming pattern and subject are aligned
    comp = compareStrings(pattern, subject)
    qs = str_count(comp, '\\?')
    ps = str_count(comp, '\\+')
    return(c(qs, ps))
}
countCoincidentSp <- function(aln, processorsNum = NULL){
    # Phase 6: each oneAmbiguousColumn call is microsecond-fast -- Phase-2 audit
    # showed mclapply fork overhead dominated the actual work. Plain serial
    # lapply is faster on realistic alignment widths.
    is = seq_len(aln@ranges@width[1])
    r = lapply(is, oneAmbiguousColumn, aln=aln)
    r = Filter(Negate(is.null), r)

    if(length(r)>0){
        r = as.data.frame(matrix(unlist(r), nrow=length(r), byrow=TRUE))
        names(r) = c('column.number', 'ambiguities', 'column')
        return(r)
    }else{
        return(NULL)
    }
}
oneAmbiguousColumn <- function(i, aln){
    ambiguous = names(IUPAC_CODE_MAP)[5:length(names(IUPAC_CODE_MAP))]
    col = as.character(subseq(aln, i, i))
    str = paste(col, sep="", collapse="")
    amb = sum(col %in% ambiguous)
    if(amb>1){
        return(c(i, amb, str))
    }
}
### ----------------------------------------------------------------------------
### Phase 17: alternative consensus base-callers (Sprint 3)
###
### `.computeConsensusMajority(aln, weights = NULL)` — at each alignment
### column, picks the base with the highest count (or highest summed
### weight if `weights` is non-NULL). Returns a list with:
###   * `consensus`: a single DNAString with one character per column
###     ("-" if every read has a gap there).
###   * `qualityScores`: an integer vector matching the consensus length,
###     reporting either the synthetic agreement-Phred (no weights) or
###     the mean Phred of agreeing reads (with weights).
###
### Used by Issue #87 (majority rule) and Issue #48 (Phred-aware
### resolution). The synthetic / averaged scores power Issue #33
### (consensus quality scores reported on `@contigSeq` via attr()).
### ----------------------------------------------------------------------------
.computeConsensusMajority <- function(aln, weights = NULL) {
    mat   <- as.matrix(aln)
    nrows <- nrow(mat)
    ncols <- ncol(mat)
    cons_chars <- character(ncols)
    qscores    <- integer(ncols)
    for (j in seq_len(ncols)) {
        col   <- mat[, j]
        keep  <- col != "-"
        bases <- col[keep]
        if (length(bases) == 0L) {
            cons_chars[j] <- "-"
            qscores[j]    <- 0L
            next
        }
        if (is.null(weights)) {
            tab    <- table(bases)
            winner <- names(tab)[which.max(tab)]
            qscores[j] <- as.integer(round(40 * max(tab) / length(bases)))
        } else {
            w    <- weights[keep, j]
            tab  <- tapply(w, bases, sum)
            winner <- names(tab)[which.max(tab)]
            agree_w <- weights[col == winner, j]
            qscores[j] <- if (length(agree_w) > 0L)
                              as.integer(round(mean(agree_w)))
                          else 0L
        }
        cons_chars[j] <- winner
    }
    list(
        consensus     = DNAString(paste(cons_chars, collapse = "")),
        qualityScores = qscores
    )
}

### `.buildQualityMatrix(aln, qualityPhredScoresList)` — map each read's
### per-base Phred score onto its aligned columns. Returns a numeric
### matrix the same shape as `as.matrix(aln)` where entry (i, j) is the
### Phred score of read i at alignment column j (0 in gap columns).
.buildQualityMatrix <- function(aln, qualityPhredScoresList) {
    mat    <- as.matrix(aln)
    nrows  <- nrow(mat)
    ncols  <- ncol(mat)
    out    <- matrix(0L, nrow = nrows, ncol = ncols)
    rnms   <- rownames(mat)
    for (i in seq_len(nrows)) {
        ## DECIPHER prefixes alignment names like "1_Read_<file>";
        ## try the prefixed form first, then fall back to the basename.
        key       <- rnms[i]
        clean_key <- sub("^[0-9]+_Read_", "", key)
        q <- qualityPhredScoresList[[key]]
        if (is.null(q)) q <- qualityPhredScoresList[[clean_key]]
        if (is.null(q)) {
            ## No quality known — fall back to flat Phred 30 across the
            ## non-gap columns of this read.
            n_nongap <- sum(mat[i, ] != "-")
            q <- rep(30L, n_nongap)
        }
        pos <- 0L
        for (j in seq_len(ncols)) {
            if (mat[i, j] != "-") {
                pos <- pos + 1L
                if (pos <= length(q)) out[i, j] <- as.integer(q[pos])
            }
        }
    }
    out
}

### ----------------------------------------------------------------------------
### Calculating SangerContig
### ----------------------------------------------------------------------------
calculateContigSeq <- function(inputSource, forwardReadList, reverseReadList,
                               refAminoAcidSeq, minFractionCall,
                               maxFractionLost, geneticCode,
                               acceptStopCodons, readingFrame,
                               processorsNum = NULL, printLevel="",
                               BPPARAM = NULL,
                               minOverlapFraction     = 0.0,
                               minOverlapBases        = 0L,
                               alignSeqsParams        = list(),
                               consensusMethod        = "strict",
                               qualityAware           = FALSE,
                               qualityPhredScoresList = NULL) {
    BPPARAM <- .resolveBPPARAM(processorsNum, BPPARAM)
    processorsNum <- BiocParallel::bpnworkers(BPPARAM)
    ### ------------------------------------------------------------------------
    ### forward & reverse character reads list string creation
    ### ------------------------------------------------------------------------
    fRDNAStringSet <- lapply(forwardReadList, function(forwardRead) {
        primaryDNA <- as.character(forwardRead@primarySeq)
        if (inputSource == "ABIF") {
            trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
            primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
        }
        return(primaryDNA)
    })
    rRDNAStringSet <- lapply(reverseReadList, function(reverseRead) {
        DNALen <- length(reverseRead@primarySeq)
        primaryDNA <- as.character(reverseComplement(reverseRead@primarySeq))
        if (inputSource == "ABIF") {
            trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
            ## Trimming on the original reads (not reverseComplement).
            primaryDNA <- substr(primaryDNA, DNALen - trimmedFinishPos + 1,
                                 DNALen - trimmedStartPos)
        }
        return(primaryDNA)
    })
    ### ------------------------------------------------------------------------
    ### DNAStringSet storing forward & reverse reads ! (Origin)
    ### ------------------------------------------------------------------------
    frReadSet <- DNAStringSet(c(unlist(fRDNAStringSet),
                                unlist(rRDNAStringSet)))
    frReadFeatureList <- c(rep("Forward Reads", length(fRDNAStringSet)),
                           rep("Reverse Reads", length(rRDNAStringSet)))
    # Read number in each contig can be 1!
    # if(length(frReadSet) < 2) {
    #     error <- paste("\n'Valid abif files should be more than 2.\n",
    #                    sep = "")
    #     log_error(error)
    # }
    processorsNum <- getProcessors(processorsNum)

    ### ------------------------------------------------------------------------
    ### Amino acid reference sequence CorrectFrameshifts correction
    ### ------------------------------------------------------------------------
    if (refAminoAcidSeq != "") {
        if (printLevel == "SangerContig") {
            log_info("Correcting frameshifts in reads using amino acid",
                     "reference sequence")   
        }
        # verbose_print <- printLevel == "SangerContig"
        # Verbose should be FALSE, but I get error when calling it
        corrected =
            CorrectFrameshifts(myXStringSet = frReadSet,
                               myAAStringSet = AAStringSet(refAminoAcidSeq),
                               geneticCode = geneticCode,
                               type = 'both',
                               processors = processorsNum)
        
        
        
        
        
        
        
        
        
        
        
        
        frReadSet = corrected$sequences
        indels = getIndelDf(corrected$indels)
        stops = as.numeric(unlist(BiocParallel::bplapply(
            frReadSet, countStopSodons,
            readingFrame, geneticCode,
            BPPARAM = BPPARAM)))
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
        log_error(error)
    }
    ### ------------------------------------------------------------------------
    ### Reads with stop codons elimination
    ### ------------------------------------------------------------------------
    ### ------------------------------------------------------------------------
    ### Remove reads with stop codons
    ### ------------------------------------------------------------------------
    if (!acceptStopCodons) {
        log_info("Removing reads with stop codons")
        if(refAminoAcidSeq == ""){ # otherwise we already did it above
            stops =
                as.numeric(unlist(BiocParallel::bplapply(
                    frReadSet, countStopSodons,
                    readingFrame, geneticCode,
                    BPPARAM = BPPARAM)))
            stopsDf = data.frame("read" = names(frReadSet),
                                 "stopCodons" = stops)
        }
        old_length = length(frReadSet)
        frReadSet = frReadSet[which(stops==0)]
        # Modify
        log_info(old_length - length(frReadSet),
                "reads with stop codons removed")
    }

    if(length(frReadSet) < 2) {
        error <- paste("\n'After removing reads with stop codons, ",
                       "forward and reverse reads should be more than 2.\n",
                       sep = "")
        log_error(error)
    }

    ### ------------------------------------------------------------------------
    ### Issue #42: defensive pre-alignment filter — drop reads whose
    ### trimmed primary sequence is shorter than 2 bp. The Phase-3 / -4
    ### `minReadLength` filter at the SangerContig level handles the
    ### typical case, but on aggressively-trimmed degenerate inputs (M1
    ### produces trimmedFinishPos = 0 for some windows; FASTA reads can
    ### be length 1 if the user supplied a very short fragment) a length-
    ### 1 entry can survive into `frReadSet` and silently break the
    ### downstream alignment / consensus.
    ### ------------------------------------------------------------------------
    too_short <- BiocGenerics::width(frReadSet) < 2L
    if (any(too_short)) {
        log_warn(">> Dropping ", sum(too_short),
                 " read(s) with trimmed length < 2 bp ",
                 "(MIN_READ_LENGTH_DEFENSIVE_DROP).")
        frReadSet <- frReadSet[!too_short]
    }
    if (length(frReadSet) < 2L) {
        ## Issue #42: when the defensive filter (or upstream filtering)
        ## drops us below the AlignSeqs minimum, we must NOT enter
        ## DECIPHER — it errors on length-1 inputs. Return a degenerate
        ## but well-formed result so the SangerContig.initialize caller
        ## can fold this into a controlled READ_NUMBER_ERROR instead of
        ## crashing the whole alignment.
        log_warn(">> Fewer than 2 usable reads after defensive filter; ",
                 "returning empty consensus.")
        if (length(frReadSet) == 1L) {
            consensusGapfree <- frReadSet[[1L]]
        } else {
            consensusGapfree <- DNAString()
        }
        return(list("consensusGapfree" = consensusGapfree,
                    "diffsDf"          = data.frame(),
                    "aln2"             = DNAStringSet(),
                    "dist"             = matrix(),
                    "dend"             = list(),
                    "indels"           = indels,
                    "stopsDf"          = stopsDf,
                    "spDf"             = data.frame()))
    }

    ### ------------------------------------------------------------------------
    ### Start aligning reads
    ###
    ### Issue #94: users with low-overlap F+R reads (~50 bp overlap on
    ### 800 bp reads, common in 16S barcoding) need fine control over the
    ### DECIPHER alignment parameters. Pass `alignSeqsParams` through to
    ### `AlignSeqs` / `AlignTranslation` so callers can tune
    ### `iterations`, `gapOpening`, `refinements`, etc.
    ### ------------------------------------------------------------------------
    base_align_args <- list(myXStringSet = frReadSet,
                            processors    = processorsNum,
                            verbose       = FALSE)
    extra_align_args <- alignSeqsParams[
        !names(alignSeqsParams) %in% names(base_align_args)]
    if (refAminoAcidSeq != "") {
        aln = do.call(AlignTranslation,
                      c(base_align_args,
                        list(geneticCode = geneticCode),
                        extra_align_args))
    } else {
        aln = do.call(AlignSeqs,
                      c(base_align_args, extra_align_args))
    }
    names(aln) = paste(seq_len(length(aln)), "Read",
                       basename(names(aln)), sep="_")

    ### ------------------------------------------------------------------------
    ### Issue #94 / #66: post-alignment overlap-quality check.
    ###
    ### Compute, for each pair of aligned reads, the number of alignment
    ### columns where BOTH reads have a non-gap base. That count is the
    ### "shared overlap length". If the minimum shared overlap across all
    ### pairs is below `minOverlapFraction * shorter_read_length` AND
    ### below `minOverlapBases` (when supplied), log a `LOW_OVERLAP_WARN`.
    ###
    ### The default thresholds (0.0 / 0L) preserve pre-Phase-16 behaviour
    ### unless the caller opts in. Callers that opt in get protection
    ### against the silent IUPAC-ambiguity-soup outputs reported in #66.
    ### ------------------------------------------------------------------------
    overlap_min_obs <- NA_integer_
    overlap_min_pair <- NA_character_
    if (length(aln) >= 2L &&
        (minOverlapFraction > 0 || minOverlapBases > 0L)) {
        nongap <- as.matrix(aln) != "-"
        ## Pairwise overlap counts via crossprod.
        ## (rows = reads, cols = alignment columns)
        ovl <- crossprod(t(nongap))     # nreads x nreads
        diag(ovl) <- NA_integer_
        ## Each read's effective length (non-gap column count).
        read_len <- rowSums(nongap)
        overlap_min_obs <- min(ovl, na.rm = TRUE)
        ix <- which(ovl == overlap_min_obs, arr.ind = TRUE)
        if (nrow(ix) > 0L) {
            overlap_min_pair <- paste(rownames(ovl)[ix[1L, 1L]],
                                       colnames(ovl)[ix[1L, 2L]],
                                       sep = " <-> ")
        }
        shorter_pair_len <- min(read_len[ix[1L, ]])
        threshold_frac <- minOverlapFraction * shorter_pair_len
        threshold      <- max(threshold_frac, as.numeric(minOverlapBases))
        if (overlap_min_obs < threshold) {
            log_warn(">> LOW_OVERLAP_WARN: smallest pairwise overlap is ",
                     overlap_min_obs, " bp (between ", overlap_min_pair,
                     "); required threshold is ", round(threshold, 1L),
                     " bp. Consensus may contain spurious IUPAC ",
                     "ambiguity codes; review carefully or tighten ",
                     "trimming parameters before merging.")
        }
    }

    ### ------------------------------------------------------------------------
    ### Issues #87 / #48 / #33: pluggable consensus base-callers.
    ###
    ### `consensusMethod`:
    ###   * "strict"           — pre-Phase-17 default; uses DECIPHER's
    ###                          ConsensusSequence with IUPAC ambiguity codes
    ###                          for disagreeing bases (issue #87 reporter
    ###                          calls this "ambiguous bases in consensus").
    ###   * "majority"         — at each column pick the most-frequent base
    ###                          (plurality vote); ties break by alphabetical
    ###                          order of the base. Synthesises per-position
    ###                          Phred = 40 * (winner_count / total_count).
    ###   * "quality_weighted" — same as majority but votes are weighted by
    ###                          source-read Phred scores. Per-position
    ###                          consensus Phred is the mean of agreeing
    ###                          reads' scores at that column. (issue #48)
    ###
    ### `qualityAware = TRUE` is shorthand for `consensusMethod =
    ### "quality_weighted"`.
    ###
    ### Consensus quality scores (issue #33) are returned both as a
    ### top-level list element AND attached to the gap-free consensus via
    ### `attr(..., "qualityScores")` so they're discoverable from
    ### `attributes(sc@contigSeq)$qualityScores`.
    ### ------------------------------------------------------------------------
    if (isTRUE(qualityAware)) consensusMethod <- "quality_weighted"
    if (!consensusMethod %in% c("strict", "majority", "quality_weighted")) {
        stop("`consensusMethod` must be one of 'strict', 'majority', ",
             "or 'quality_weighted'.")
    }

    if (consensusMethod == "strict") {
        consensus <- ConsensusSequence(aln,
                                        minInformation = minFractionCall,
                                        includeTerminalGaps = TRUE,
                                        threshold = maxFractionLost,
                                        noConsensusChar = "-",
                                        ambiguity = TRUE)[[1L]]
        consensusQualityScores <- integer(0)
    } else {
        weights_mat <- NULL
        if (consensusMethod == "quality_weighted") {
            if (is.null(qualityPhredScoresList)) {
                log_warn(">> consensusMethod = 'quality_weighted' requested ",
                         "but no qualityPhredScoresList supplied; falling ",
                         "back to flat Phred-30 weights.")
            }
            qpls <- if (is.null(qualityPhredScoresList)) list()
                    else qualityPhredScoresList
            weights_mat <- .buildQualityMatrix(aln, qpls)
        }
        majority_res <- .computeConsensusMajority(aln, weights = weights_mat)
        consensus              <- majority_res$consensus
        consensusQualityScores <- majority_res$qualityScores
    }

    # Phase 6: nPairwiseDiffs is microsecond-fast per read; serial lapply
    # avoids fork overhead that previously dominated this call site.
    diffs = lapply(aln, nPairwiseDiffs, subject = consensus)
    diffs = do.call(rbind, diffs)
    diffsDf = data.frame("name" = names(aln),
                         "pairwise.diffs.to.consensus" = diffs[,1],
                         "unused.chars" = diffs[,2])
    rownames(diffsDf) = NULL

    # get a dendrogram
    dist = DistanceMatrix(aln, correction = "JC69",
                          penalizeGapLetterMatches = FALSE,
                          processors = processorsNum, verbose = FALSE)
    dend = Treeline(myDistMatrix=dist, method="UPGMA",
                    type = "both", showPlot = FALSE, processors = processorsNum, verbose = FALSE)

    # add consensus to alignment
    aln2 = c(aln, DNAStringSet(consensus))
    names(aln2)[length(aln2)] = "Consensus"
    # strip gaps from consensus (must be an easier way!!)
    consensusGapfree = RemoveGaps(DNAStringSet(consensus))[[1]]

    ### ------------------------------------------------------------------------
    ### Issue #33: align the consensus quality vector to the gap-stripped
    ### consensus. `consensusQualityScores` initially has one entry per
    ### alignment column including the gap columns we just stripped; we
    ### subset to the non-gap positions so length(qualityScores) ==
    ### length(consensusGapfree).
    ### ------------------------------------------------------------------------
    if (length(consensusQualityScores) > 0L) {
        cons_str <- as.character(consensus)
        cons_chars <- strsplit(cons_str, "", fixed = TRUE)[[1L]]
        keep <- cons_chars != "-"
        cons_qs_gapfree <- consensusQualityScores[keep]
    } else {
        cons_qs_gapfree <- integer(0)
    }
    attr(consensusGapfree, "qualityScores") <- cons_qs_gapfree

    # count columns in the alignment with >1 coincident secondary peaks
    spDf = countCoincidentSp(aln, processorsNum = processorsNum)
    if (is.null(spDf)) {
        spDf = data.frame()
    }
    return(list("consensusGapfree"        = consensusGapfree,
                "diffsDf"                 = diffsDf,
                "aln2"                    = aln2,
                "dist"                    = dist,
                "dend"                    = dend,
                "indels"                  = indels,
                "stopsDf"                 = stopsDf,
                "spDf"                    = spDf,
                "consensusQualityScores"  = cons_qs_gapfree))
}
### ----------------------------------------------------------------------------
### MakeBaseCalls related function
### ----------------------------------------------------------------------------
MakeBaseCallsInside <- function(traceMatrix, peakPosMatrixRaw,
                                qualityPhredScoresRaw,
                                signalRatioCutoff, readFeature, printLevel) {
    if (printLevel == "SangerRead") {
        log_info("          * Making basecall !!")
    }
    #get peaks for each base
    Apeaks <- getpeaks(traceMatrix[,1])
    Cpeaks <- getpeaks(traceMatrix[,2])
    Gpeaks <- getpeaks(traceMatrix[,3])
    Tpeaks <- getpeaks(traceMatrix[,4])
    
    #get window around primary basecall peaks
    primarypeaks <- peakPosMatrixRaw[,1]
    diffs <- diff(c(0,primarypeaks))
    starts <- primarypeaks - 0.5*diffs
    stops <- c(primarypeaks[seq_len((length(primarypeaks)-1))] +
                   0.5*diffs[2:length(diffs)],
               primarypeaks[length(diffs)] + 0.5*diffs[length(diffs)]
    )
    #hack for last peak. Just uses distance preceding peak
    #as distance after peak

    #Now get max peak value for each channel in each peak window.
    #If no peak return 0
    primary <- NULL
    secondary <- NULL
    tempPosMatrix <- matrix(nrow=length(starts), ncol=4)
    tempAmpMatrix <- matrix(nrow=length(starts), ncol=4)
    indexBaseCall <- c()

    # Phase 7: batch the peak lookups. peakvalues_batch_cpp processes all
    # peak windows for one channel in a single .Call, eliminating the
    # per-window R-to-C++ marshalling overhead that dominated when the
    # function was called once per (channel, window) pair.
    AbatchOut <- peakvalues_batch_cpp(Apeaks, starts, stops)
    CbatchOut <- peakvalues_batch_cpp(Cpeaks, starts, stops)
    GbatchOut <- peakvalues_batch_cpp(Gpeaks, starts, stops)
    TbatchOut <- peakvalues_batch_cpp(Tpeaks, starts, stops)

    for(i in seq_len(length(starts))) {
        Apeak <- AbatchOut[, i]
        Cpeak <- CbatchOut[, i]
        Gpeak <- GbatchOut[, i]
        Tpeak <- TbatchOut[, i]
        if(is.na(Apeak[2]) &
           is.na(Cpeak[2]) &
           is.na(Gpeak[2]) &
           is.na(Tpeak[2])) {
            next #rare case where no peak found
        }
        ### --------------------------------------------------------------------
        ### My modification here: Tracking BaseCall index
        ###     Add qualtiy score when making basecall
        ### --------------------------------------------------------------------
        indexBaseCall <- c(indexBaseCall, i)
        signals <- c(Apeak[1], Cpeak[1], Gpeak[1], Tpeak[1])

        tempAmpMatrix[i,] <- signals
        positions <- c(Apeak[2], Cpeak[2], Gpeak[2], Tpeak[2])
        tempPosMatrix[i,] <- positions
        signalratios <- signals/max(signals, na.rm=TRUE)
        Bases <- c("A", "C", "G", "T")
        Bases[signalratios < signalRatioCutoff] <- NA
        #sort by decreasing signal strength
        Bases <- Bases[order(signals, decreasing=TRUE)]
        positions <- positions[order(signals, decreasing=TRUE)]
        if(length(Bases[!is.na(Bases)]) == 4
           | length(Bases[!is.na(Bases)]) == 0) {
            primary <- c(primary, "N")
            secondary <- c(secondary, "N")
        } else if (length(Bases[!is.na(Bases)]) > 1) {
            primary <- c(primary, Bases[1])
            Bases2 <- Bases[2:4]
            sortedBase2<- sort(Bases2[!is.na(Bases2)])
            dicValue <- paste(sortedBase2, collapse="")
            secondaryLetter <- names(IUPAC_CODE_MAP[IUPAC_CODE_MAP == dicValue])
            secondary <- c(secondary, secondaryLetter)
        }
        else {
            primary <- c(primary, Bases[1])
            secondary <- c(secondary, Bases[1])
        }
    }
    if (readFeature == "Forward Read") {
        qualityPhredScores <- qualityPhredScoresRaw[indexBaseCall]
        primarySeq <- DNAString(paste(primary, collapse=""))
        secondarySeq <- DNAString(paste(secondary, collapse=""))
    } else if (readFeature == "Reverse Read") {
        qualityPhredScores <- qualityPhredScoresRaw[indexBaseCall]
        primarySeq <- DNAString(paste(primary, collapse=""))
        secondarySeq <- DNAString(paste(secondary, collapse=""))
    }
    peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    if (printLevel == "SangerRead") {
        log_info("          * Updating slots in 'SangerRead' instance !!")
    }
    return(list("qualityPhredScores" = qualityPhredScores,
                "peakPosMatrix" = peakPosMatrix,
                "peakAmpMatrix" = peakAmpMatrix,
                "primarySeq" = primarySeq,
                "secondarySeq" = secondarySeq))
}
getpeaks <- function(trace) {
    r <- rle(trace)
    indexes <- which(rep(diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2,
                         times = r$lengths))
    cbind(indexes, trace[indexes])
}
### Phase 7: kept as a private helper for the equivalence test in
### tests/testthat/test-Rcpp-peakvalues.R. Production paths use the C++
### implementation via peakvalues_cpp() in src/peakvalues.cpp (~30x faster
### per call on typical Sanger reads).
.peakvalues_r <- function(x, pstart, pstop) {
    region <- x[x[,1] > pstart & x[,1] < pstop, ,drop=FALSE]
    if (length(region[,1]) == 0) return(c(0, NA))
    else return(c(max(region[,2], na.rm=TRUE), region[which.max(region[,2]),1]))
}

### ----------------------------------------------------------------------------
### MakeBasecall secondary peak finding helper function
### ----------------------------------------------------------------------------
IUPAC_CODE_MAP <- c(
    A="A",
    C="C",
    G="G",
    T="T",
    M="AC",
    R="AG",
    W="AT",
    S="CG",
    Y="CT",
    K="GT",
    V="ACG",
    H="ACT",
    D="AGT",
    B="CGT",
    N="ACGT"
)

### ----------------------------------------------------------------------------
### chromatogram related function
### ----------------------------------------------------------------------------
chromatogramRowNum <- function(width, rawLength, trimmedLength, showTrimmed) {
    if (showTrimmed) {
        numplots = ceiling(rawLength / width)
    } else {
        numplots = ceiling(trimmedLength / width)
    }
}

# <---------------------------------------------------------------------------->

### ============================================================================
### SangerRead related helper functions
### ============================================================================
SetCharStyleList <- function(AASeqDF, selectChar, colorCode) {
    stopIndex <- AASeqDF %>% `==` (selectChar) %>% which()
    stopExcelIndex <- int2col(stopIndex)
    stopExcelIndexName <- paste0(stopExcelIndex, "1")
    styleList <-
        as.list(rep(paste('background-color:', colorCode,
                          "; font-weight: bold;"), length(stopExcelIndex)))
    if (length(stopIndex) != 0) {
        names(styleList) <- stopExcelIndexName
    }
    return(styleList)
}
SetAllStyleList <- function(AASeqDF, colorCode) {
    Index <- strtoi(names(AASeqDF))
    ExcelIndex <- int2col(Index)
    ExcelIndexName <- paste0(ExcelIndex, "1")
    styleList <-
        as.list(rep(paste('background-color:', colorCode,
                          "; font-weight: bold;"), length(ExcelIndex)))
    names(styleList) <- ExcelIndexName
    return(styleList)
}
countStopSodons <- function(sequence,
                            readingFrame = 1, geneticCode = GENETIC_CODE){
    l = length(sequence) + 1 - readingFrame
    if(l < 3){
        sprintf("Cannot calculate stop codons on sequence of length %d",
                " in reading frame %d",length(sequence), readingFrame)
        # return(NULL)
        error <- paste("\nCannot calculate stop codons on sequence of length ",
                       length(sequence), " in reading frame ", readingFrame,
                       ".\n", sep = "")
        log_error(error)
    }
    # this comes almost straight from the BioStrings manual
    tri = trinucleotideFrequency(sequence[readingFrame:length(sequence)],step=3)
    names(tri) <- geneticCode[names(tri)]
    freqs = lapply(split(tri, names(tri)), sum)
    stops = freqs["*"]
    return(as.numeric(stops))
}
calculateAASeq <- function(primarySeq, trimmedStartPos,
                           trimmedFinishPos, geneticCode) {
    DNASeqshift0 <- DNAString(substr(as.character(primarySeq),
                                     1+trimmedStartPos, trimmedFinishPos))
    primaryAASeqS1 <-
        suppressWarnings(Biostrings::translate(DNASeqshift0,
                                   genetic.code = geneticCode,
                                   no.init.codon=TRUE,
                                   if.fuzzy.codon="solve"))

    DNASeqshift1 <- DNAString(substr(as.character(primarySeq),
                                     2+trimmedStartPos, trimmedFinishPos))
    primaryAASeqS2 <-
        suppressWarnings(Biostrings::translate(DNASeqshift1,
                                   genetic.code = geneticCode,
                                   no.init.codon=TRUE,
                                   if.fuzzy.codon="solve"))

    DNASeqshift2 <- DNAString(substr(as.character(primarySeq),
                                     3+trimmedStartPos, trimmedFinishPos))
    primaryAASeqS3 <-
        suppressWarnings(Biostrings::translate(DNASeqshift2,
                                   genetic.code = geneticCode,
                                   no.init.codon=TRUE,
                                   if.fuzzy.codon="solve"))
    return(list("primaryAASeqS1" = primaryAASeqS1,
                "primaryAASeqS2" = primaryAASeqS2,
                "primaryAASeqS3" = primaryAASeqS3))
}
### ----------------------------------------------------------------------------
### Quality trimming related parameter
### ----------------------------------------------------------------------------
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
    trimmedQualityPhredScore <-
        qualityPhredScores[(trimmedStartPos+1):trimmedFinishPos]
    trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
    trimmedMinQualityScore <- min(trimmedQualityPhredScore)
    remainingRatio = trimmedSeqLength / rawSeqLength

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
M2inside_calculate_trimming <-function(qualityPhredScores,
                                       M2CutoffQualityScore,
                                       M2SlidingWindowSize) {
    rawSeqLength <- length(qualityPhredScores)
    rawMeanQualityScore <- mean(qualityPhredScores)
    rawMinQualityScore <- min(qualityPhredScores)
    if (M2SlidingWindowSize > 40 || M2SlidingWindowSize < 0 ||
        M2SlidingWindowSize%%1!=0 ||
        M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
        M2CutoffQualityScore%%1!=0) {
        trimmedStartPos = NULL
        trimmedFinishPos = NULL
    } else {
        ### ------------------------------------------------------------------------
        ### Find the trimming start point
        ###  First window that the average score is bigger than threshold score.
        ###  (Whole window will be kept)
        ### ------------------------------------------------------------------------
        totalThresholdScore <- M2CutoffQualityScore * M2SlidingWindowSize
        trimmedStartPos = 1
        trimmedFinishPos = 2
        for (i in seq_len((rawSeqLength-M2SlidingWindowSize+1))) {
            totalScore <-
                sum(qualityPhredScores[i:(i+M2SlidingWindowSize-1)])
            if (totalScore > totalThresholdScore) {
                trimmedStartPos = i
                break
            }
        }
        qualityPhredScoresRev <- rev(qualityPhredScores)
        for (i in seq_len((rawSeqLength-M2SlidingWindowSize+1))) {
            totalScore <-
                sum(qualityPhredScoresRev[i:(i+M2SlidingWindowSize-1)])
            if (totalScore > totalThresholdScore) {
                trimmedFinishPos = i
                break
            }
        }
        trimmedFinishPos <- length(qualityPhredScoresRev) - trimmedFinishPos + 1
        # for (i in (trimmedStartPos+M2SlidingWindowSize-1):(rawSeqLength-M2SlidingWindowSize+1)) {
        #     totalScore <-
        #         sum(qualityPhredScores[i:(i+M2SlidingWindowSize-1)])
        #     if (totalScore < totalThresholdScore) {
        #         # Keep all base pairs in the previous window.
        #         trimmedFinishPos = i + M2SlidingWindowSize - 2
        #         break
        #     }
        # }
        if (trimmedStartPos == (rawSeqLength-M2SlidingWindowSize+1) ||
            trimmedStartPos == trimmedFinishPos) {
            trimmedStartPos = 1
            trimmedFinishPos = 2
        }
        trimmedSeqLength = trimmedFinishPos - trimmedStartPos
        trimmedQualityPhredScore <-
            qualityPhredScores[(trimmedStartPos+1):trimmedFinishPos]
        trimmedMeanQualityScore <- mean(trimmedQualityPhredScore)
        trimmedMinQualityScore <- min(trimmedQualityPhredScore)
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

### ----------------------------------------------------------------------------
### Quality score base pair plot functions
### ----------------------------------------------------------------------------
QualityBasePlotly <- function(trimmedStartPos, trimmedFinishPos,
                              readLen, qualityPlotDf, x,  y) {
    p <- suppressPlotlyMessage(
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
                                 xanchor = "center",
                                 x = 0.5, y = 1.1)) %>%
            add_annotations(
                text = "Trimming Start <br> BP Index",
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
    return(p)
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
### Phase 8: Plotly + WebGL chromatogram renderer.
###
### A full Sanger trace can have ~10^4 points per channel. The legacy
### `chromatogram_overwrite` uses base-R graphics via `polygon()` which is
### fine for static images but freezes the browser when wrapped in
### Shiny/htmlwidgets at full resolution.
###
### `chromatogram_plotly()` returns a single Plotly htmlwidget that:
###   * uses `scattergl` (WebGL) traces -- keeps the browser responsive
###     even at >50k points per channel,
###   * downsamples to `max_points` per channel by uniform-stride
###     subsampling when the trace is longer (preserves peak silhouettes
###     well; for production use one would prefer LTTB, but stride is
###     deterministic and zero-dep),
###   * supports the same `colors` argument as `chromatogram_overwrite`
###     ("default" / "cb_friendly" / a 5-vector of hex colours).
###
### Returned object is a `plotly::plotly` htmlwidget that the Shiny app
### can render with `plotly::renderPlotly`.
### ============================================================================

#' Render a Sanger chromatogram as an interactive Plotly widget
#'
#' Wraps the four trace channels (A/C/G/T) of a sangerseq / SangerRead
#' object into a single \code{plotly} htmlwidget that renders via WebGL
#' (\code{scattergl}). Intended for embedding in Shiny dashboards where
#' the static \code{\link{chromatogram_overwrite}} would be too heavy.
#'
#' @param obj A sangerseq or SangerRead instance with a populated
#'   \code{traceMatrix}.
#' @param trim5 Integer; if \code{showtrim} is TRUE, shade the first
#'   \code{trim5} positions to indicate the 5' trim region.
#' @param trim3 Integer; if \code{showtrim} is TRUE, shade the last
#'   \code{trim3} positions to indicate the 3' trim region.
#' @param max_points Integer cap on the number of points rendered per
#'   channel. When the trace exceeds \code{max_points} it is downsampled
#'   by uniform stride.
#' @param showtrim Logical; whether to overlay shaded trim regions.
#' @param colors Either \code{"default"}, \code{"cb_friendly"}, or a
#'   length-5 character vector of hex colours for (A, T, C, G, other).
#'
#' @return A \code{plotly} htmlwidget. The returned object carries a
#'   \code{downsample_info} attribute reporting the original and rendered
#'   point counts plus the stride.
#'
#' @examples
#' data(sangerReadFData)
#' \donttest{
#' chromatogram_plotly(sangerReadFData)
#' }
#' @export
chromatogram_plotly <- function(obj,
                                 trim5      = 0,
                                 trim3      = 0,
                                 max_points = 8000L,
                                 showtrim   = FALSE,
                                 colors     = "default") {
    if (!is(obj, "sangerseq")) {
        stop("'obj' must be a sangerseq (or SangerRead) S4 object.")
    }

    palette <- if (identical(colors, "default")) {
        c(A = "#2ca02c", T = "#1f77b4", C = "#000000", G = "#d62728",
          other = "#9467bd")
    } else if (identical(colors, "cb_friendly")) {
        c(A = "#000000", T = "#c7c7c7", C = "#0072b2", G = "#d55e00",
          other = "#cc79a7")
    } else if (is.character(colors) && length(colors) == 5L) {
        setNames(colors, c("A", "T", "C", "G", "other"))
    } else {
        stop("'colors' must be \"default\", \"cb_friendly\", or a length-5 character vector")
    }

    trace_mat <- obj@traceMatrix
    if (is.null(trace_mat) || nrow(trace_mat) == 0L) {
        stop("`obj@traceMatrix` is empty -- no chromatogram to render.")
    }
    n_total <- nrow(trace_mat)

    # Uniform-stride downsample if the trace is longer than max_points.
    if (n_total > max_points) {
        stride <- ceiling(n_total / max_points)
        idx    <- seq.int(1L, n_total, by = stride)
    } else {
        idx <- seq.int(1L, n_total)
    }
    x_axis <- idx

    # Order matches sangerseqR convention: traceMatrix columns are A, C, G, T.
    p <- plotly::plot_ly()
    p <- plotly::add_trace(p,
        x = x_axis, y = trace_mat[idx, 1L],
        type = "scattergl", mode = "lines",
        line = list(color = palette[["A"]], width = 1),
        name = "A")
    p <- plotly::add_trace(p,
        x = x_axis, y = trace_mat[idx, 2L],
        type = "scattergl", mode = "lines",
        line = list(color = palette[["C"]], width = 1),
        name = "C")
    p <- plotly::add_trace(p,
        x = x_axis, y = trace_mat[idx, 3L],
        type = "scattergl", mode = "lines",
        line = list(color = palette[["G"]], width = 1),
        name = "G")
    p <- plotly::add_trace(p,
        x = x_axis, y = trace_mat[idx, 4L],
        type = "scattergl", mode = "lines",
        line = list(color = palette[["T"]], width = 1),
        name = "T")

    # Optional shaded trim region.
    if (showtrim && (trim5 > 0 || trim3 > 0)) {
        if (trim5 > 0) {
            p <- plotly::add_trace(p,
                x = c(0, trim5), y = c(0, 0),
                type = "scattergl", mode = "lines",
                line = list(color = "rgba(200, 200, 200, 0.5)", width = 30),
                name = "5' trimmed", hoverinfo = "skip")
        }
        if (trim3 > 0) {
            p <- plotly::add_trace(p,
                x = c(n_total - trim3, n_total), y = c(0, 0),
                type = "scattergl", mode = "lines",
                line = list(color = "rgba(200, 200, 200, 0.5)", width = 30),
                name = "3' trimmed", hoverinfo = "skip")
        }
    }

    p <- plotly::layout(p,
        xaxis  = list(title = "Trace position",
                      range = c(1L, n_total)),
        yaxis  = list(title = "Signal"),
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.1),
        hovermode = "x unified")

    attr(p, "downsample_info") <- list(
        original_points     = n_total,
        rendered_points     = length(idx),
        downsample_stride   = if (exists("stride", inherits = FALSE)) stride else 1L
    )
    p
}

SangerReadInnerTrimming <- function(SangerReadInst, inputSource) {
    primaryDNA <- as.character(SangerReadInst@primarySeq)
    if (inputSource == "ABIF") {
        trimmedStartPos <- 
            SangerReadInst@QualityReport@trimmedStartPos
        trimmedFinishPos <- 
            SangerReadInst@QualityReport@trimmedFinishPos
        primaryDNA <- 
            substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
    }
    return(primaryDNA)
}

#' Static base-R chromatogram renderer with a corrected color palette.
#'
#' Reimplementation of \code{sangerseqR::chromatogram} with a fix for
#' base color rendering. Intended for static (PDF / PNG) export. For
#' interactive embedding in Shiny see \code{\link{chromatogram_plotly}}.
#'
#' @param obj A sangerseq or SangerRead instance.
#' @param trim5 Integer; number of bases to mark as 5' trimmed.
#' @param trim3 Integer; number of bases to mark as 3' trimmed.
#' @param showcalls One of \code{"primary"}, \code{"secondary"},
#'   \code{"both"}, or \code{"none"}.
#' @param width Bases per row.
#' @param height Plot height per row (relative units).
#' @param cex.mtext Text size for marginal annotations.
#' @param cex.base Text size for base-call labels.
#' @param ylim Maximum y-axis multiplier (relative to robust mean).
#' @param filename Optional path to write a PDF.
#' @param showtrim Logical; if TRUE, shade the trim regions.
#' @param showhets Logical; if TRUE, mark heterozygous positions.
#' @param colors Either \code{"default"}, \code{"cb_friendly"}, or a
#'   length-5 character vector of hex colours.
#'
#' @return Invisibly returns NULL; called for its side effect of
#'   plotting (or writing to \code{filename}).
#'
#' @examples
#' data(sangerReadFData)
#' \donttest{
#' chromatogram_overwrite(sangerReadFData)
#' }
#' @export
chromatogram_overwrite <- function(obj, trim5=0, trim3=0,
                                   showcalls=c("primary", "secondary", "both", "none"),
                                   width=100, height=2, cex.mtext=1, cex.base=1, ylim=3,
                                   filename=NULL, showtrim=FALSE, showhets=TRUE, colors="default") {
    if (colors == "default") {
        A_color = "green"
        T_color = "blue"
        C_color = "black"
        G_color = "red"
        unknown_color = "purple"
    } else if (colors == "cb_friendly") {
        A_color = rgb(0, 0, 0, maxColorValue = 255)
        T_color = rgb(199, 199, 199, maxColorValue = 255)
        C_color = rgb(0, 114, 178, maxColorValue = 255)
        G_color = rgb(213, 94, 0, maxColorValue = 255)
        unknown_color = rgb(204, 121, 167, maxColorValue = 255)
    } else {
        A_color = colors[1]
        T_color = colors[2]
        C_color = colors[3]
        G_color = colors[4]
        unknown_color = colors[5]
    }
    originalpar <- par(no.readonly=TRUE)
    showcalls <- showcalls[1]
    traces <- obj@traceMatrix
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    basecalls1 <- basecalls1[1:length(aveposition)] 
    basecalls2 <- basecalls2[1:length(aveposition)] 
    if(showtrim == FALSE) {
        if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
        else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
        if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
        else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
        aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
    }
    indexes <- 1:length(basecalls1)
    trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all 
    #false if not trimmed
    if (!is.null(trim3)) {
        traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, 
                                nrow(traces))), ]
    }
    if (!is.null(trim5)) {
        offset <- max(c(1, aveposition[1] - 10))
        traces <- traces[offset:nrow(traces),]
        aveposition <- aveposition - (offset-1)
    }
    maxsignal <- apply(traces, 1, max)
    ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))           
    p <- c(0, aveposition, nrow(traces))
    midp <- diff(p)/2
    starts <- aveposition - midp[1:(length(midp)-1)]
    starthets <- starts
    starthets[basecalls1 == basecalls2] <- NA
    ends <- aveposition + midp[2:(length(midp))]
    endhets <- ends
    endhets[basecalls1 == basecalls2] <- NA
    starttrims <- starts
    starttrims[!trimmed] <- NA
    endtrims <- ends
    endtrims[!trimmed] <- NA
    
    colortranslate <- c(A=A_color, C=C_color, G=G_color, T=T_color)
    colorvector1 <- unname(colortranslate[basecalls1])
    colorvector1[is.na(colorvector1)] <- unknown_color
    colorvector2 <- unname(colortranslate[basecalls2])
    colorvector2[is.na(colorvector2)] <- unknown_color
    
    valuesperbase <- nrow(traces)/length(basecalls1)
    tracewidth <- width*valuesperbase
    breaks <- seq(1,nrow(traces), by=tracewidth)
    numplots <- length(breaks)
    if(!is.null(filename)) {
        if (! dir.exists(dirname(filename))) {
            dir.create(dirname(filename))
        }
        pdf(filename, width=8.5, height=height*numplots) 
    }
    par(mar=c(2,2,2,1), mfrow=c(numplots, 1))
    basecallwarning1 = 0
    basecallwarning2 = 0
    j = 1
    
    for(i in breaks) {
        range <- aveposition >= i & aveposition < (i+tracewidth)
        starthet <- starthets[range] - tracewidth*(j-1)
        starthet[starthet < 0] <- 0
        endhet <- endhets[range] - tracewidth*(j-1)
        endhet[endhet > tracewidth] <- tracewidth
        lab1 <- basecalls1[range]
        lab2 <- basecalls2[range]
        pos <- aveposition[range] - tracewidth*(j-1)
        colors1 <- colorvector1[range]
        colors2 <- colorvector2[range]
        starttrim <- starttrims[range] - tracewidth*(j-1)
        endtrim <- endtrims[range] - tracewidth*(j-1)
        plotrange <- i:min(i+tracewidth, nrow(traces))
        plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n", 
             bty="n", xlab="", yaxt="n", , xlim=c(1,tracewidth))
        if (showhets==TRUE) {
            rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')
        }
        if (showtrim==TRUE) {
            rect(starttrim, 0, endtrim, ylims[2], col='red', border='transparent', 
                 density=15)
        }
        lines(traces[plotrange,1], col=A_color)
        lines(traces[plotrange,2], col=C_color )
        lines(traces[plotrange,3], col=G_color)
        lines(traces[plotrange,4], col=T_color)
        mtext(as.character(which(range)[1]), side=2, line=0, cex=cex.mtext)
        
        for(k in 1:length(lab1)) {
            if (showcalls=="primary" | showcalls=="both") {
                if (is.na(basecalls1[1]) & basecallwarning1==0) {
                    warning("Primary basecalls missing")
                    basecallwarning1 = 1
                } 
                else if (length(lab1) > 0) {   
                    axis(side=3, at=pos[k], labels=lab1[k], col.axis=colors1[k], 
                         family="mono", cex=cex.base, line=ifelse(showcalls=="both", 0, 
                                                                  -1), tick=FALSE)
                }
            }
            if (showcalls=="secondary" | showcalls=="both") {
                if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
                    warning("Secondary basecalls missing")
                    basecallwarning2 = 1
                } 
                else if (length(lab2) > 0) { 
                    axis(side=3, at=pos[k], labels=lab2[k], col.axis=colors2[k], 
                         family="mono", cex=cex.base, line=-1, tick=FALSE)
                }
            }
        }
        j = j + 1
    }
    if(!is.null(filename)) {
        dev.off()
        log_info("Chromatogram saved to ", filename,
                 " in the current working directory")
    }
    else par(originalpar)
}
