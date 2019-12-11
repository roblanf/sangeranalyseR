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
    if (selection != "all" && selection != "alignment" &&
        selection != "allReads" && selection != "contig") {
        stop("\nSelection must be 'all', 'alignment', 'allReads' or 'contig'.")
    }
    contigName = obj@contigName
    message("Start to write '", contigName, "' to FASTA format ...")
    ### ------------------------------------------------------------------------
    ### Writing alignment result to FASTA file (Exclude consensus read)
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "alignment") {
        message("\n    >> Writing alignment to FASTA ...")
        alignmentObject = obj@alignment
        alignmentObject$Consensus <- NULL
        writeXStringSet(alignmentObject,
                        file.path(outputDir,
                                  paste0(contigName, "_reads_alignment.fa")),
                        compress = compress,
                        compression_level = compression_level)
    }

    ### ------------------------------------------------------------------------
    ### Writing all single read into FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "allReads") {
        message("\n    >> Writing all single reads to FASTA ...")
        fRDNAStringSet <- sapply(obj@forwardReadList, function(forwardRead) {
            trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
            trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
            primaryDNA <- as.character(forwardRead@primarySeq)
            substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
        })
        names(fRDNAStringSet) <- basename(names(fRDNAStringSet))
        rRDNAStringSet <- sapply(obj@reverseReadList, function(reverseRead) {
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
                                  paste0(contigName,
                                         "_all_trimmed_reads.fa")),
                        compress = compress,
                        compression_level = compression_level)
    }

    ### ------------------------------------------------------------------------
    ### Writing consensus read into FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "contig") {
        message("\n    >> Writing consensus read to FASTA ...")

        writeTarget <- DNAStringSet(obj@contigSeq)
        names(writeTarget) <- paste0(contigName, "_contig")
        writeXStringSet(writeTarget,
                        file.path(outputDir,paste0(contigName, "_contig.fa")),
                        compress = compress,
                        compression_level = compression_level)
    }
    message("\nFinish writing '", contigName, "' to FASTA format")
})

#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' contigName <- "RBNII395-13[C_LepFolF,C_LepFolR]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' A_chloroticContig <- SangerContig(
#'                                  parentDirectory       = inputFilesParentDir,
#'                                  contigName            = contigName,
#'                                  suffixForwardRegExp   = suffixForwardRegExp,
#'                                  suffixReverseRegExp   = suffixReverseRegExp,
#'                                  TrimmingMethod        = "M1",
#'                                  M1TrimmingCutoff      = 0.0001,
#'                                  M2CutoffQualityScore  = NULL,
#'                                  M2SlidingWindowSize   = NULL,
#'                                  baseNumPerRow         = 100,
#'                                  heightPerRow          = 200,
#'                                  signalRatioCutoff     = 0.33,
#'                                  showTrimmed           = TRUE)
#' RShinyCS <- launchAppSangerContig(A_chloroticContig)
setMethod("launchAppSangerContig", "SangerContig",
          function(obj, outputDir = NULL) {
              ### ------------------------------------------------------------------------
              ### Checking SangerContig input parameter is a list containing
              ### one S4 object.
              ### ------------------------------------------------------------------------
              if (is.null(outputDir)) {
                  outputDir <- tempdir()
                  suppressWarnings(dir.create(outputDir))
              }
              if (dir.exists(outputDir)) {
                  shinyOptions(sangerContig = list(obj))
                  shinyOptions(shinyDirectory = outputDir)
                  newSangerContig <- shinyApp(SangerContigUI, SangerContigServer)
                  return(newSangerContig)
              } else {
                  stop("'", outputDir, "' is not valid. Please check again")
              }
})

setMethod("generateReport", "SangerContig", function(obj, outputDir) {
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir))
    }

    # obj <- A_chloroticContig
    outputDirSR <- file.path(outputDir, "SangerRead_Report")
    if (!dir.exists(outputDirSR)) {
        suppressWarnings(dir.create(outputDirSR, recursive = TRUE))
    }

    allSangerRead <- c(obj@forwardReadList, obj@reverseReadList)

    outputFN <- sapply(allSangerRead, generateReport, outputDir)

    outputDirSC <- file.path(outputDir, "SangerRead_Report")

    rootDir <- system.file(package = "sangeranalyseR")
    originRmd <- file.path(rootDir, "vignettes", "SangerContig_Report.Rmd")
    outputRmd <- file.path(outputDirSC, "Sanger_Report.Rmd")

    res <- render(input = "/Users/chaokuan-hao/Documents/ANU_2019_Semester_2/Lanfear_Lab/sangeranalyseR/vignettes/SangerContig_Report.Rmd",
                  output_dir = outputDirSC,
                  params = list(SangerContig = obj,
                                outputDir = outputDir,
                                outputFN = outputFN))
})
