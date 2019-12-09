### ============================================================================
### Plotting trimmed and remaining ratio for "SangerContig" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaConsensusRead.RDdata")
#' trimmingRatioPlot(A_chloroticaConsensusRead)
setMethod("trimmingRatioPlot",  "SangerContig", function(object){
    ### ------------------------------------------------------------------------
    ### Trimmed ratio plot for forward read
    ### ------------------------------------------------------------------------
    fdQualityReportObject = object@forwardReadSangerseq@QualityReport
    fdPlotting <- preTrimmingRatioPlot(fdQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Trimmed ratio plot for reverse read
    ### ------------------------------------------------------------------------
    rvQualityReportObject = object@reverseReadSangerseq@QualityReport
    rvPlotting <- preTrimmingRatioPlot(rvQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Grid plotting for forward and reverse reads
    ### ------------------------------------------------------------------------
    grid.arrange(fdPlotting, rvPlotting, ncol=2)
})



### ============================================================================
### Plotting quality for each base for "SangerContig" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaConsensusRead.RDdata")
#' qualityBasePlot(A_chloroticaConsensusRead)
setMethod("qualityBasePlot",  "SangerContig", function(object){
    ### ------------------------------------------------------------------------
    ### Quality base plot for forward read
    ### ------------------------------------------------------------------------
    fdQualityReportObject = object@forwardReadSangerseq@QualityReport
    fdPlotting <- preQualityBasePlot(fdQualityReportObject)

    ### ------------------------------------------------------------------------
    ### Grid plotting for forward and reverse reads
    ### ------------------------------------------------------------------------
    rvQualityReportObject = object@reverseReadSangerseq@QualityReport
    rvPlotting <- preQualityBasePlot(rvQualityReportObject)

    grid.arrange(fdPlotting, rvPlotting, ncol=2)
})


## ============================================================================
## Updating quality parameters for SangerContig object.
## ============================================================================
#' @example
#' load("data/A_chloroticaConsensusRead.RDdata")
#' trimmingRatioPlot(A_chloroticaConsensusRead)
#' qualityBasePlot(A_chloroticaConsensusRead)
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@TrimmingMethod
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M1TrimmingCutoff
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2CutoffQualityScore
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2SlidingWindowSize
#'
#' A_chloroticaConsensusRead <-
#'                 updateQualityParam(A_chloroticaConsensusRead,
#'                                    "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(A_chloroticaConsensusRead)
#' qualityBasePlot(A_chloroticaConsensusRead)
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@TrimmingMethod
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M1TrimmingCutoff
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2CutoffQualityScore
#' A_chloroticaConsensusRead@forwardReadSangerseq@QualityReport@M2SlidingWindowSize
setMethod("updateQualityParam",  "SangerContig",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
    ### ------------------------------------------------------------------------
    ### Updating forward read quality parameters
    ### ------------------------------------------------------------------------
    # qualityBaseScoresFD <-
    #     object@forwardReadSangerseq@QualityReport@qualityBaseScores
    object@forwardReadSangerseq <-
        updateQualityParam(object@forwardReadSangerseq,
                           TrimmingMethod,
                           M1TrimmingCutoff,
                           M2CutoffQualityScore,
                           M2SlidingWindowSize)

    ### ------------------------------------------------------------------------
    ### Updating reverse read quality parameters
    ### ------------------------------------------------------------------------
    object@reverseReadSangerseq <-
        updateQualityParam(object@reverseReadSangerseq,
                           TrimmingMethod,
                           M1TrimmingCutoff,
                           M2CutoffQualityScore,
                           M2SlidingWindowSize)
    return(object)
})

setMethod("writeFASTA", "SangerContig", function(obj, outputDir, compress,
                                                        compression_level,
                                                        selection = "all") {
    ### ------------------------------------------------------------------------
    ### selection can be 'all', 'alignment' 'allReads' 'contig'
    ### ------------------------------------------------------------------------
    consensusReadName = obj@consensusReadName
    message("Start to write '", consensusReadName,
            "' to FASTA format ...")
    ### ------------------------------------------------------------------------
    ### Writing alignment result to FASTA file (Exclude consensus read)
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "alignment") {
        message("\n    >> Writing alignment to FASTA ...")
        alignmentObject = obj@alignment
        alignmentObject$Consensus <- NULL
        writeXStringSet(alignmentObject,
                        file.path(outputDir,
                                  paste0(consensusReadName, "_alignment.fa")))
    }


    ### ------------------------------------------------------------------------
    ### Writing all single read into FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "allReads") {
        message("\n    >> Writing all single reads to FASTA ...")
        fRDNAStringSet <- sapply(obj@forwardReadsList, function(forwardRead) {
            trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
            primaryDNA <- as.character(forwardRead@primarySeq)
            substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
        })
        names(fRDNAStringSet) <- basename(names(fRDNAStringSet))
        rRDNAStringSet <- sapply(obj@reverseReadsList, function(reverseRead) {
            trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
            primaryDNA <- as.character(reverseRead@primarySeq)
            substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
        })
        names(rRDNAStringSet) <- basename(names(rRDNAStringSet))
        frReadSet <- DNAStringSet(c(unlist(fRDNAStringSet),
                                    unlist(rRDNAStringSet)))
        writeXStringSet(frReadSet,
                        file.path(outputDir,
                                  paste0(consensusReadName,
                                         "_all_trimmed_reads.fa")))
    }


    ### ------------------------------------------------------------------------
    ### Writing consensus read into FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "contig") {
        message("\n    >> Writing consensus read to FASTA ...")

        writeTarget <- DNAStringSet(obj@consensusRead)
        names(writeTarget) <- paste0(consensusReadName, "_contig")
        writeXStringSet(writeTarget,
                        file.path(outputDir,
                                  paste0(consensusReadName, "_contig.fa")))
    }
    message("\nFinish writing '", consensusReadName, "' to FASTA format")
})
