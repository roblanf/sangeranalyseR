## =============================================================================
## Updating quality parameters for SangerContig object.
## =============================================================================
#' @examples
#' load("data/sangerContig.RData")
#' updateQualityParam(sangerContig,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
setMethod("updateQualityParam",  "SangerContig",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
    ### ------------------------------------------------------------------------
    ### Updating forward read quality parameters
    ### Quality parameters is checked in 'QualityReport' method
    ### ------------------------------------------------------------------------
    errors <- character()
    errors <- checkTrimParam(TrimmingMethod,
                             M1TrimmingCutoff,
                             M2CutoffQualityScore,
                             M2SlidingWindowSize,
                             errors)
    if (length(errors) == 0) {
        newForwardReadList <- sapply(object@forwardReadList,
                                     function(forwardRead) {
            forwardRead <-
                updateQualityParam(forwardRead,
                                   TrimmingMethod         = TrimmingMethod,
                                   M1TrimmingCutoff       = M1TrimmingCutoff,
                                   M2CutoffQualityScore   = M2CutoffQualityScore,
                                   M2SlidingWindowSize    = M2SlidingWindowSize)
        })
        object@forwardReadList <- newForwardReadList
        newReverseReadList <- sapply(object@reverseReadList,
                                     function(reverseRead) {
            updforwardRead <-
                updateQualityParam(reverseRead,
                                   TrimmingMethod         = TrimmingMethod,
                                   M1TrimmingCutoff       = M1TrimmingCutoff,
                                   M2CutoffQualityScore   = M2CutoffQualityScore,
                                   M2SlidingWindowSize    = M2SlidingWindowSize)
        })
        object@reverseReadList <- newReverseReadList
        return(object)
    } else {
        stop(errors)
    }
})

#' @examples
#' load("data/sangerContig.RData")
#' RShinySC <- launchAppSC(sangerContig)
setMethod("launchAppSC", "SangerContig", function(obj, outputDir = NULL) {
    ### --------------------------------------------------------------
    ### Checking SangerContig input parameter is a list containing
    ### one S4 object.
    ### --------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    message(">>> outputDir : ", outputDir)
    if (dir.exists(outputDir)) {
        shinyOptions(sangerContig = list(obj))
        shinyOptions(shinyDirectory = outputDir)
        newSangerContig <-
            shinyApp(SangerContigUI, SangerContigServer)
        return(newSangerContig)
    } else {
        stop("'", outputDir, "' is not valid. Please check again")
    }
})

## =============================================================================
## Writing primary sequence into FASTA format
## =============================================================================
#' @examples
#' load("data/sangerContig.RData")
#' writeFastaSC(sangerContig)
setMethod("writeFastaSC", "SangerContig", function(obj, outputDir, compress,
                                                   compression_level,
                                                   selection = "all") {
    ### ------------------------------------------------------------------------
    ### selection can be 'all', 'alignment' 'allReads' 'contig'
    ### ------------------------------------------------------------------------
    if (selection != "all" && selection != "alignment" &&
        selection != "allReads" && selection != "contig") {
        stop("\nSelection must be 'all', 'alignment', 'allReads' or 'contig'.")
    }
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    message(">>> outputDir : ", outputDir)
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

## =============================================================================
## Generating report for SangerContig
## =============================================================================
#' @examples
#' load("data/sangerContig.RData")
#' generateReportSC(sangerContig)
setMethod("generateReportSC", "SangerContig",
          function(obj, outputDir, includeSangerRead = TRUE,
                   navigationAlignmentFN = NULL) {
    ### ------------------------------------------------------------------------
    ### Make sure the input directory is not NULL
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    message(">>> outputDir : ", outputDir)
    ### ------------------------------------------------------------------------
    ### Make sure the directory is exist (SangerContig level)
    ###  => SangerRead level directory will be created recursively
    ### ------------------------------------------------------------------------
    outputDirSC <- file.path(outputDir, obj@contigName)
    if (!dir.exists(outputDirSC)) {
        suppressWarnings(dir.create(outputDirSC, recursive = TRUE))
    }

    rootDir <- system.file(package = "sangeranalyseR")
    originRmd <- file.path(rootDir, "rmd", "SangerContig_Report.Rmd")
    outputHtml <- file.path(outputDirSC, "SangerContig_Report.html")

    forwardReads <- obj@forwardReadList
    reverseReads <- obj@reverseReadList

    if(includeSangerRead) {
        forwardReadFN <- sapply(forwardReads, generateReportSR,
                                outputDir = outputDirSC,
                                navigationContigFN = outputHtml,
                                navigationAlignmentFN = navigationAlignmentFN)
        reverseReadFN <- sapply(reverseReads, generateReportSR,
                                outputDir = outputDirSC,
                                navigationContigFN = outputHtml,
                                navigationAlignmentFN = navigationAlignmentFN)
    } else {
        forwardReadFN = NULL
        reverseReadFN = NULL
    }
    res <- render(input = originRmd,
                  output_dir = outputDirSC,
                  params = list(SangerContig = obj,
                                outputDir = outputDirSC,
                                forwardReadFN = forwardReadFN,
                                reverseReadFN = reverseReadFN,
                                navigationAlignmentFN = navigationAlignmentFN))
    return(outputHtml)
})
