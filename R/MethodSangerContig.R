## =============================================================================
## Updating quality parameters for SangerContig object.
## =============================================================================
#' @title updateQualityParam
#' @name SangerContig-class-updateQualityParam
#' @rdname SangerContig-class-updateQualityParam
#' @examples
#' load("data/sangerContig.RData")
#' updateQualityParam(sangerContig,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)
setMethod("updateQualityParam",  "SangerContig",function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
    if (object@inputSource == "ABIF") {
        ### --------------------------------------------------------------------
        ### Updating forward read quality parameters
        ### Quality parameters is checked in 'QualityReport' method
        ### --------------------------------------------------------------------
        errors <- character()
        errors <- checkTrimParam(TrimmingMethod,
                                 M1TrimmingCutoff,
                                 M2CutoffQualityScore,
                                 M2SlidingWindowSize,
                                 errors)
        if (length(errors) == 0) {
            newForwardReadList <-
                sapply(object@forwardReadList,
                       function(forwardRead) {
                           forwardRead <-
                               updateQualityParam(
                                   forwardRead,
                                   TrimmingMethod         = TrimmingMethod,
                                   M1TrimmingCutoff       = M1TrimmingCutoff,
                                   M2CutoffQualityScore   = M2CutoffQualityScore,
                                   M2SlidingWindowSize    = M2SlidingWindowSize)
                       })
            object@forwardReadList <- newForwardReadList
            newReverseReadList <-
                sapply(object@reverseReadList,
                       function(reverseRead) {
                           updforwardRead <-
                               updateQualityParam(
                                   reverseRead,
                                   TrimmingMethod         = TrimmingMethod,
                                   M1TrimmingCutoff       = M1TrimmingCutoff,
                                   M2CutoffQualityScore   = M2CutoffQualityScore,
                                   M2SlidingWindowSize    = M2SlidingWindowSize)
                       })
            object@reverseReadList <- newReverseReadList
            object@trimmingMethodSC <- TrimmingMethod
            CSResult <-
                calculateContigSeq (inputSource      = object@inputSource,
                                    forwardReadList  = object@forwardReadList,
                                    reverseReadList  = object@reverseReadList,
                                    refAminoAcidSeq  = object@refAminoAcidSeq,
                                    minFractionCall  = object@minFractionCall,
                                    maxFractionLost  = object@maxFractionLost,
                                    geneticCode      = object@geneticCode,
                                    acceptStopCodons = object@acceptStopCodons,
                                    readingFrame     = object@readingFrame,
                                    processorsNum    = getProcessors(NULL))
            object@contigSeq <- CSResult$consensusGapfree
            object@differencesDF <- CSResult$diffsDf
            object@alignment <- CSResult$aln2
            object@distanceMatrix <- CSResult$dist
            object@dendrogram <- CSResult$dend
            object@indelsDF <- CSResult$indels
            object@stopCodonsDF <- CSResult$stopsDf
            object@secondaryPeakDF <- CSResult$spDf
            return(object)
        } else {
            stop(errors)
        }
    } else if (object@inputSource == "FASTA") {
        message("SangerContig with 'FASTA' inputSource ",
                "cannot update quality parameters")
    }
})

#' @title launchAppSC
#' @name SangerContig-class-launchAppSC
#' @rdname SangerContig-class-launchAppSC
#' @examples
#' load("data/sangerContig.RData")
#' RShinySC <- launchAppSC(sangerContig)
setMethod("launchAppSC", "SangerContig", function(obj, outputDir = NULL) {
    if (obj@inputSource == "ABIF") {
        ### --------------------------------------------------------------------
        ### Checking SangerContig input parameter is a list containing
        ### one S4 object.
        ### --------------------------------------------------------------------
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
    } else if (obj@inputSource == "FASTA") {
        message("SangerContig with 'FASTA' inputSource ",
                "cannot run Shiny app\n (You don't need to ",
                "do trimming or base calling)")
    }
})

## =============================================================================
## Writing primary sequence into FASTA format
## =============================================================================
#' @title writeFastaSC
#' @name SangerContig-class-writeFastaSC
#' @rdname SangerContig-class-writeFastaSC
#' @examples
#' load("data/sangerContig.RData")
#' writeFastaSC(sangerContig, "/Users/chaokuan-hao/Desktop/sangeranalyseR_fasta/SangerContig")
setMethod("writeFastaSC", "SangerContig", function(obj, outputDir, compress,
                                                   compression_level,
                                                   selection = "all") {
    ### ------------------------------------------------------------------------
    ### selection can be 'all', 'reads_alignment' 'reads_unalignment' 'contig'
    ### ------------------------------------------------------------------------
    if (selection != "all" && selection != "reads_alignment" &&
        selection != "reads_unalignment" && selection != "contig") {
        stop("\nSelection must be 'all',
             'reads_alignment', 'reads_unalignment' or 'contig'.")
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
    if (selection == "all" || selection == "reads_alignment") {
        message("\n    >> Writing alignment to FASTA ...")
        alignmentObject = obj@alignment
        alignmentObject$Consensus <- NULL
        writeAlignment <- append(alignmentObject, list(obj@contigSeq))
        names(writeAlignment) <- sub("^[0-9]*_", "", names(writeAlignment))
        names(writeAlignment)[length(writeAlignment)] <-
            paste0(obj@contigName, "_contig")
        writeXStringSet(writeAlignment,
                        file.path(outputDir,
                                  paste0(contigName, "_reads_alignment.fa")),
                        compress = compress,
                        compression_level = compression_level)
    }

    ### ------------------------------------------------------------------------
    ### Writing all single read into FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "reads_unalignment") {
        message("\n    >> Writing all single reads to FASTA ...")
        fRDNAStringSet <- sapply(obj@forwardReadList, function(forwardRead) {
            primaryDNA <- as.character(forwardRead@primarySeq)
            ### ----------------------------------------------------------------
            ### Only read in ABIF file format needs to do trimming
            ### ----------------------------------------------------------------
            if (obj@inputSource == "ABIF") {
                trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
                trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
                primaryDNA <- substr(primaryDNA,
                                     trimmedStartPos+1, trimmedFinishPos)
            }
            return(primaryDNA)
        })
        names(fRDNAStringSet) <- basename(names(fRDNAStringSet))
        rRDNAStringSet <- sapply(obj@reverseReadList, function(reverseRead) {
            ### ----------------------------------------------------------------
            ### Only read in ABIF file format needs to do trimming
            ### ----------------------------------------------------------------
            primaryDNA <- as.character(reverseRead@primarySeq)
            if (obj@inputSource == "ABIF") {
                # Trim first and then reverse complement
                trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
                trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
                primaryDNA <- substr(primaryDNA,
                                     trimmedStartPos+1, trimmedFinishPos)
                primaryDNA <- as.character(reverseComplement(DNAString(primaryDNA)))
            }
            return(primaryDNA)
        })
        names(rRDNAStringSet) <- basename(names(rRDNAStringSet))
        frReadSet <- DNAStringSet(c(unlist(fRDNAStringSet),
                                    unlist(rRDNAStringSet)))
        writeXStringSet(frReadSet,
                        file.path(outputDir,
                                  paste0(contigName,
                                         "_reads_unalignment.fa")),
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
#' @title generateReportSC
#' @name SangerContig-class-generateReportSC
#' @rdname SangerContig-class-generateReportSC
#' @examples
#' load("data/sangerContig.RData")
#' generateReportSC(sangerContig)
setMethod("generateReportSC", "SangerContig",
          function(obj, outputDir, includeSangerRead = TRUE,
                   navigationAlignmentFN = NULL) {





    # if (object@inputSource == "ABIF") {
    #
    # } else if (object@inputSource == "FASTA") {
    #     message("SangerContig with 'FASTA' inputSource ",
    #             "cannot run Shiny app\n (You don't need to ",
    #             "do trimming or base calling)")
    # }




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

    if (obj@inputSource == "ABIF") {
        outputHtml <- file.path(outputDirSC, "SangerContig_Report.html")
    } else if (obj@inputSource == "FASTA") {
        outputHtml <- file.path(outputDirSC, "SangerContig_Report.html")
    }
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
