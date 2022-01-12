## =============================================================================
## Updating quality parameters for SangerContig object.
## =============================================================================
#' A SangerContig method which updates QualityReport parameter for each the SangerRead instance inside SangerContig.
#'
#' @title updateQualityParam
#' @name SangerContig-class-updateQualityParam
#' @aliases updateQualityParam,SangerContig-method
#'
#' @param object A SangerContig S4 instance.
#' @param TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#' @param processorsNum The number of processors to use, or NULL (the default) for all available processors.
#'
#' @return A SangerContig instance.
#'
#' @examples
#' data("sangerContigData")
#' \dontrun{
#' updateQualityParam(sangerContigData,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)}
setMethod("updateQualityParam",  "SangerContig",function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL,
                   processorsNum          = NULL){
    if (object@inputSource == "ABIF") {
        ### --------------------------------------------------------------------
        ### Updating forward read quality parameters
        ### Quality parameters is checked in 'QualityReport' method
        ### --------------------------------------------------------------------
        errors <- list(character(0), character(0))
        errors <- checkTrimParam(TrimmingMethod,
                                 M1TrimmingCutoff,
                                 M2CutoffQualityScore,
                                 M2SlidingWindowSize,
                                 errors[[1]], errors[[2]])
        if (length(errors[[1]]) == 0) {
            newForwardReadList <-
                lapply(object@forwardReadList,
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
                lapply(object@reverseReadList,
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
                                    processorsNum    = getProcessors(processorsNum))
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
            sapply(paste0(errors[[2]], errors[[1]], '\n') , 
                   log_error, simplify = FALSE)
        }
    } else if (object@inputSource == "FASTA") {
        log_info("SangerContig with 'FASTA' inputSource ",
                "cannot update quality parameters")
    }
})

#' A SangerContig method which launches Shiny app for SangerContig instance.
#'
#' @title launchAppSC
#' @name SangerContig-class-launchAppSC
#' @aliases launchAppSC,SangerContig-method
#'
#' @param object A SangerContig S4 instance.
#' @param outputDir The output directory of the saved new SangerContig S4 instance.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#'
#' @return A \code{shiny.appobj} object.
#'
#' @examples 
#' data("sangerContigData")
#' RShinySC <- launchAppSC(sangerContigData)
#' RShinySC <- launchAppSC(sangerContigData, colors="cb_friendly")
setMethod("launchAppSC", "SangerContig", function(object, outputDir = NULL, colors = "default") {
    if (object@inputSource == "ABIF") {
        ### --------------------------------------------------------------------
        ### Checking SangerContig input parameter is a list containing
        ### one S4 object.
        ### --------------------------------------------------------------------
        if (is.null(outputDir)) {
            outputDir <- tempdir()
        }
        if (!dir.exists(outputDir)) {
            suppressWarnings(dir.create(outputDir, recursive = TRUE))
        }
        log_info(">>> outputDir : ", outputDir)
        if (dir.exists(outputDir)) {
            shinyOptions(sangerContig = list(object))
            shinyOptions(shinyDirectory = outputDir)
            shinyOptions(colors = colors)
            newSangerContig <-
                shinyApp(SangerContigUI, SangerContigServer)
            return(newSangerContig)
        } else {
            stop("'", outputDir, "' is not valid. Please check again")
        }
    } else if (object@inputSource == "FASTA") {
        log_info("SangerContig with 'FASTA' inputSource ",
                "cannot run Shiny app\n (You don't need to ",
                "do trimming or base calling)")
    }
})

## =============================================================================
## Writing primary sequence into FASTA format
## =============================================================================
#' A SangerContig method which writes sequences into Fasta files.
#'
#' @title writeFastaSC
#' @name SangerContig-class-writeFastaSC
#' @aliases writeFastaSC,SangerContig-method
#'
#' @param object A SangerContig S4 instance.
#' @param outputDir The output directory of generated FASTA files.
#' @param compress Like for the \code{save} function in base R, must be \code{TRUE} or \code{FALSE} (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip". This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param compression_level This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param selection This value can be \code{all}, \code{reads_alignment}, \code{reads_unalignment} or \code{contig}. It generates reads and the contig FASTA files.
#'
#' @return The output directory of FASTA files.
#'
#' @examples
#' data("sangerContigData")
#' writeFastaSC(sangerContigData)
setMethod("writeFastaSC", "SangerContig", function(object, outputDir, compress,
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
    outputDir <- normalizePath(outputDir)
    log_info(">>> outputDir : ", outputDir)
    contigName = object@contigName
    log_info("Start to write '", contigName, "' to FASTA format ...")
    ### ------------------------------------------------------------------------
    ### Writing alignment result to FASTA file (Exclude consensus read)
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "reads_alignment") {
        log_info("\n    >> Writing alignment to FASTA ...")
        alignmentObject = object@alignment
        alignmentObject$Consensus <- NULL
        writeAlignment <- append(alignmentObject, list(object@contigSeq))
        names(writeAlignment) <- sub("^[0-9]*_", "", names(writeAlignment))
        names(writeAlignment)[length(writeAlignment)] <-
            paste0(object@contigName, "_contig")
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
        log_info("\n    >> Writing all single reads to FASTA ...")
        fRDNAStringSet <- lapply(object@forwardReadList, function(forwardRead) {
            primaryDNA <- as.character(forwardRead@primarySeq)
            ### ----------------------------------------------------------------
            ### Only read in ABIF file format needs to do trimming
            ### ----------------------------------------------------------------
            if (object@inputSource == "ABIF") {
                trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
                trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
                primaryDNA <- substr(primaryDNA,
                                     trimmedStartPos+1, trimmedFinishPos)
            }
            return(primaryDNA)
        })
        names(fRDNAStringSet) <- basename(names(fRDNAStringSet))
        rRDNAStringSet <- lapply(object@reverseReadList, function(reverseRead) {
            ### ----------------------------------------------------------------
            ### Only read in ABIF file format needs to do trimming
            ### ----------------------------------------------------------------
            primaryDNA <- as.character(reverseRead@primarySeq)
            if (object@inputSource == "ABIF") {
                # Trim first and then reverse complement
                trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
                trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
                primaryDNA <- substr(primaryDNA,
                                     trimmedStartPos+1, trimmedFinishPos)
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
        log_info("\n    >> Writing consensus read to FASTA ...")

        writeTarget <- DNAStringSet(object@contigSeq)
        names(writeTarget) <- paste0(contigName, "_contig")
        writeXStringSet(writeTarget,
                        file.path(outputDir,paste0(contigName, "_contig.fa")),
                        compress = compress,
                        compression_level = compression_level)
    }
    log_info("\nFinish writing '", contigName, "' to FASTA format")
})

## =============================================================================
## Generating report for SangerContig
## =============================================================================
#' A SangerContig method which generates final reports of the SangerContig instance.
#'
#' @title generateReportSC
#' @name SangerContig-class-generateReportSC
#' @aliases generateReportSC,SangerContig-method
#'
#' @param object A SangerContig S4 instance.
#' @param outputDir The output directory of the generated HTML report.
#' @param includeSangerRead The parameter that decides whether to include SangerRead level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param navigationAlignmentFN The internal parameter passed to HTML report. Users should not modify this parameter on their own.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#'
#' @return The output absolute path to the SangerContig's HTML file.
#'
#' @examples
#' data("sangerContigData")
#' \dontrun{
#' generateReportSC(sangerContigData)
#' generateReportSC(sangerContigData, colors="cb_friendly")}
setMethod("generateReportSC", "SangerContig",
          function(object, outputDir, includeSangerRead = TRUE, 
                   colors, navigationAlignmentFN = NULL) {
    # if (object@inputSource == "ABIF") {
    #
    # } else if (object@inputSource == "FASTA") {
    #     log_info("SangerContig with 'FASTA' inputSource ",
    #             "cannot run Shiny app\n (You don't need to ",
    #             "do trimming or base calling)")
    # }
    ### ------------------------------------------------------------------------
    ### Make sure the input directory is not NULL
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
    }
    if (!dir.exists(outputDir)) {
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    outputDir <- normalizePath(outputDir)
    log_info(">>> outputDir : ", outputDir)
    ### ------------------------------------------------------------------------
    ### Make sure the directory is exist (SangerContig level)
    ###  => SangerRead level directory will be created recursively
    ### ------------------------------------------------------------------------
    outputDirSC <- file.path(outputDir, object@contigName)
    if (!dir.exists(outputDirSC)) {
        suppressWarnings(dir.create(outputDirSC, recursive = TRUE))
    }

    rootDir <- system.file(package = "sangeranalyseR")
    originRmd <- file.path(rootDir, "rmd", "SangerContig_Report.Rmd")

    if (object@inputSource == "ABIF") {
        outputHtml <- file.path(outputDirSC, "SangerContig_Report.html")
    } else if (object@inputSource == "FASTA") {
        outputHtml <- file.path(outputDirSC, "SangerContig_Report.html")
    }
    forwardReads <- object@forwardReadList
    reverseReads <- object@reverseReadList

    if(includeSangerRead) {
        forwardReadFN <- lapply(forwardReads, generateReportSR,
                                outputDir = outputDirSC,
                                colors = colors,
                                navigationContigFN = outputHtml,
                                navigationAlignmentFN = navigationAlignmentFN)
        reverseReadFN <- lapply(reverseReads, generateReportSR,
                                outputDir = outputDirSC,
                                colors = colors,
                                navigationContigFN = outputHtml,
                                navigationAlignmentFN = navigationAlignmentFN)
    } else {
        forwardReadFN = NULL
        reverseReadFN = NULL
    }
    res <- render(input = originRmd,
                  output_dir = outputDirSC,
                  params = list(SangerContig = object,
                                outputDir = outputDirSC,
                                forwardReadFN = forwardReadFN,
                                reverseReadFN = reverseReadFN,
                                navigationAlignmentFN = navigationAlignmentFN,
                                colors = colors))
    return(outputHtml)
})







## =============================================================================
## Generating summary table for SangerContig instance
## =============================================================================
#' A SangerContig method which generates summary table for SangerContig instance
#'
#' @title readTable
#' @name SangerContig-class-readTable
#' @aliases readTable,SangerContig-method
#'
#' @param object A SangerContig S4 instance.
#' @param indentation The indentation for different level printing.
#'
#' @return None
#'
#' @examples
#' data(sangerReadFData)
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' readTable(sangerReadFData)
#' readTable(sangerContigData)
#' readTable(sangerAlignmentData)
#' }
setMethod("readTable", "SangerContig", function(object, indentation = 0) {
    space <- paste(rep(' ', indentation), collapse = "")
    inputSource <- object@inputSource
    CSV_NamesConversion <- object@CSV_NamesConversion
    contigName <- object@contigName
    REGEX_SuffixForward <- object@REGEX_SuffixForward
    REGEX_SuffixReverse <- object@REGEX_SuffixReverse
    minReadsNum <- object@minReadsNum
    minReadLength <- object@minReadLength
    minFractionCall <- object@minFractionCall
    refAminoAcidSeq <- object@refAminoAcidSeq
    maxFractionLost <- object@maxFractionLost
    acceptStopCodons <- object@acceptStopCodons
    readingFrame <- object@readingFrame
    contigSeq <- object@contigSeq
    
    forwardReadList <- object@forwardReadList
    reverseReadList <- object@reverseReadList
    
    differencesDF <- object@differencesDF
    distanceMatrix <- object@distanceMatrix
    indelsDF <- object@indelsDF
    stopCodonsDF <- object@stopCodonsDF
    secondaryPeakDF <- object@secondaryPeakDF
    
    if (indentation == 0) {
        log_info("**********************************************")
        log_info("******** SangerContig readTable print ********")
        log_info("**********************************************")
    }
    cat(" ================================\n", 
        "=== SangerContig S4 instance ===\n", 
        "================================\n")
    cat(space, "                 Input source : ", inputSource, "\n",
        space, "    Names conversion CSV file : ", CSV_NamesConversion, "\n",
        space, "                  Contig name : ", contigName, "\n",
        space, "         Suffix forward regex : ", REGEX_SuffixForward, "\n",
        space, "         Suffix reverse regex : ", REGEX_SuffixReverse, "\n",
        space, "         Minimum reads number : ", minReadsNum, "\n",
        space, "          Minimum read length : ", minReadLength, "\n",
        space, "        Minimum fraction call : ", minFractionCall, "\n",
        space, "Reference amino acid sequence : ", refAminoAcidSeq, "\n",
        space, "        Maximum fraction lost : ", maxFractionLost, "\n",
        space, "           Accept stop codons : ", acceptStopCodons, "\n",
        space, "                Reading frame : ", readingFrame, "\n",
        space, "              Contig sequence : ", as.character(contigSeq), "\n\n"
    )
    lapply(forwardReadList, readTable, indentation = 1)
    lapply(reverseReadList, readTable, indentation = 1)
})
