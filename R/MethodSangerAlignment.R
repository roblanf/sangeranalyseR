## =============================================================================
## Updating quality parameters for SangerAlignment object.
## =============================================================================
#' A SangerAlignment method which updates QualityReport parameter for each the SangerRead instance inside SangerAlignment.
#'
#' @title updateQualityParam
#' @name SangerAlignment-class-updateQualityParam
#' @aliases updateQualityParam,SangerAlignment-method
#'
#' @param object A SangerAlignment S4 instance.
#' @param TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#' @param processorsNum The number of processors to use, or NULL (the default) for all available processors.
#'
#' @return A SangerAlignment instance.
#'
#' @examples
#' data("sangerAlignmentData")
#' \dontrun{
#' updateQualityParam(sangerAlignmentData,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)}
setMethod("updateQualityParam",  "SangerAlignment",
          function(object,
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
        errors <- character()
        errors <- checkTrimParam(TrimmingMethod,
                                 M1TrimmingCutoff,
                                 M2CutoffQualityScore,
                                 M2SlidingWindowSize,
                                 errors)
        if (length(errors) == 0) {
            newContigList <-
                lapply(object@contigList,
                       function(contig) {
                           contig <-
                               updateQualityParam(
                                   contig,
                                   TrimmingMethod         = TrimmingMethod,
                                   M1TrimmingCutoff       = M1TrimmingCutoff,
                                   M2CutoffQualityScore   = M2CutoffQualityScore,
                                   M2SlidingWindowSize    = M2SlidingWindowSize,
                                   processorsNum          = processorsNum)
                       })
            object@contigList <- newContigList
            object@trimmingMethodSA <- TrimmingMethod
            acResult <- alignContigs(object@contigList, object@geneticCode,
                                     object@refAminoAcidSeq,
                                     object@minFractionCallSA,
                                     object@maxFractionLostSA,
                                     getProcessors(processorsNum))
            object@contigsConsensus <- acResult[["consensus"]]
            object@contigsAlignment <- acResult[["aln"]]
            object@contigsTree <- acResult[["aln.tree"]]
            return(object)
        } else {
            log_error(paste(errors, collapse = ""))
        }
    } else if (object@inputSource == "FASTA") {
        log_info("SangerAlignment with 'FASTA' inputSource ",
                "cannot update quality parameters")
    }
})

#' A SangerAlignment method which launches Shiny app for SangerAlignment instance.
#'
#' @title launchAppSA
#' @name SangerAlignment-class-launchAppSA
#' @aliases launchAppSA,SangerAlignment-method
#'
#' @param object A SangerAlignment S4 instance.
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
#' data("sangerAlignmentData")
#' RShinySA <- launchAppSA(sangerAlignmentData)
#' RShinySA <- launchAppSA(sangerAlignmentData, colors="cb_friendly")
setMethod("launchAppSA", "SangerAlignment", function(object, outputDir = NULL, 
                                                     colors = "default") {
    if (object@inputSource == "ABIF") {
        ### --------------------------------------------------------------------
        ### Checking SangerAlignment input parameter is a list containing
        ### one S4 object.
        ### --------------------------------------------------------------------
        if (is.null(outputDir)) {
            outputDir <- tempdir()
            suppressWarnings(dir.create(outputDir, recursive = TRUE))
        }
        log_info(">>> outputDir : ", outputDir)

        if (dir.exists(outputDir)) {
            shinyOptions(sangerAlignment = list(object))
            shinyOptions(shinyDirectory = outputDir)
            shinyOptions(colors = colors)
            newSangerAlignment <- shinyApp(SangerAlignmentUI, SangerAlignmentServer)
            return(newSangerAlignment)
        } else {
            log_error("'", outputDir, "' is not valid. Please check again")
        }
    } else if (object@inputSource == "FASTA") {
        log_info("SangerAlignment with 'FASTA' inputSource ",
                "cannot run Shiny app\n (You don't need to ",
                "do trimming or base calling)")
    }
})

## =============================================================================
## Writing primary sequence into FASTA format
## =============================================================================
#' A SangerAlignment method which writes sequences into Fasta files.
#'
#' @title writeFastaSA
#' @name SangerAlignment-class-writeFastaSA
#' @aliases writeFastaSA,SangerAlignment-method
#'
#' @param object A SangerAlignment S4 instance.
#' @param outputDir The output directory of generated FASTA files.
#' @param compress Like for the \code{save} function in base R, must be \code{TRUE} or \code{FALSE} (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip". This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param compression_level This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param selection This value can be \code{all}, \code{contigs_alignment}, \code{contigs_unalignment} or \code{all_reads}. It generates reads and contigs FASTA files.
#'
#' @return The output directory of FASTA files.
#'
#' @examples
#' data("sangerAlignmentData")
#' writeFastaSA(sangerAlignmentData)
setMethod("writeFastaSA", "SangerAlignment", function(object, outputDir, compress,
                                                 compression_level,
                                                 selection = "all") {
    ### ------------------------------------------------------------------------
    ### selection can be 'all', 'alignment' 'all_reads' 'contig'
    ### ------------------------------------------------------------------------
    if (selection != "all" && selection != "contigs_alignment" &&
        selection != "contigs_unalignment" && selection != "all_reads") {
        log_error(paste0("\nSelection must be 'all', 'contigs_alignment',",
                    " 'contigs_unalignment' or 'all_reads'."))
    }
    ### ------------------------------------------------------------------------
    ### Make sure the input directory is not NULL
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    log_info(">>> outputDir : ", outputDir)
    log_info("Start to write 'SangerAlignment' to FASTA format ...")
    ### ------------------------------------------------------------------------
    ### Writing 'contigs alignment' result to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "contigs_alignment") {
        log_info("\n    >> Writing 'alignment' to FASTA ...")
        alignmentObject = object@contigsAlignment
        writeXStringSet(alignmentObject,
                        file.path(outputDir, "Sanger_contigs_alignment.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    ### ------------------------------------------------------------------------
    ### Writing 'all contigs' to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "contigs_unalignment") {
        log_info("\n    >> Writing 'contigs' to FASTA ...")
        contigsList <- lapply(object@contigList, function(contig) {
            contig@contigSeq
        })
        contigsListDNASet<- DNAStringSet(contigsList)
        writeXStringSet(contigsListDNASet,
                        file.path(outputDir, "Sanger_contigs_unalignment.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    # ### ----------------------------------------------------------------------
    # ### Writing 'contigs consensus read' to FASTA file
    # ### ----------------------------------------------------------------------
    # if (selection == "all" || selection == "consensusRead") {
    #     log_info("\n    >> Writing 'consensusRead' to FASTA ...")
    #     contigsConsensusDNASet<- DNAStringSet(object@contigsConsensus)
    #     names(contigsConsensusDNASet) <- "Sanger Consensus Read"
    #     writeXStringSet(contigsConsensusDNASet,
    #                     file.path(outputDir, "Sanger_consensus_read.fa"),
    #                     compress = compress,
    #                     compression_level = compression_level)
    # }

    ### ------------------------------------------------------------------------
    ### Writing 'all reads' to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "all_reads") {
        log_info("\n    >> Writing all single reads to FASTA ...")
        fRDNASet <- vapply(object@contigList, function(contig) {
            fRDNAStringSet <- vapply(contig@forwardReadList, function(forwardRead) {
                primaryDNA <- as.character(forwardRead@primarySeq)
                if (object@inputSource == "ABIF") {
                    trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
                    trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
                    primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
                }
                return(primaryDNA)
            }, character(1))
            names(fRDNAStringSet) <- basename(names(fRDNAStringSet))
            fRDNAStringSet
        }, character(1))
        rRDNASet <- vapply(object@contigList, function(contig) {
            rRDNAStringSet <- vapply(contig@reverseReadList, function(reverseRead) {
                primaryDNA <- as.character(reverseRead@primarySeq)
                if (object@inputSource == "ABIF") {
                    # Trim first and then reverse complement
                    trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
                    trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
                    primaryDNA <- substr(primaryDNA,
                                         trimmedStartPos+1, trimmedFinishPos)
                }
                return(primaryDNA)
            }, character(1))
            names(rRDNAStringSet) <- basename(names(rRDNAStringSet))
            rRDNAStringSet
        }, character(1))
        allDNASet <- DNAStringSet(c(fRDNASet, rRDNASet))
        writeXStringSet(allDNASet,
                        file.path(outputDir, "Sanger_all_trimmed_reads.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    log_info("\nFinish writing 'SangerAlignment' to FASTA format")
})

## =============================================================================
## Generating report for SangerContig
## =============================================================================
#' A SangerAlignment method which generates final reports of the SangerContig instance.
#'
#' @title generateReportSA
#' @name SangerAlignment-class-generateReportSA
#' @aliases generateReportSA,SangerAlignment-method
#'
#' @param object A SangerAlignment S4 instance.
#' @param outputDir The output directory of the generated HTML report.
#' @param includeSangerContig The parameter that decides whether to include SangerContig level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param includeSangerRead The parameter that decides whether to include SangerRead level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#'
#' @return The output absolute path to the SangerAlignment's HTML file.
#'
#' @examples
#' data("sangerAlignmentData")
#' \dontrun{
#' generateReportSA(sangerAlignmentData)
#' generateReportSA(sangerAlignmentData, colors="cb_friendly")}
setMethod("generateReportSA", "SangerAlignment",
          function(object, outputDir,
                   includeSangerContig = TRUE, includeSangerRead = TRUE,
                   colors) {
    ### ------------------------------------------------------------------------
    ### Make sure the input directory is not NULL
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    log_info(">>> outputDir : ", outputDir)
    ### ------------------------------------------------------------------------
    ### Make sure the directory is exist (SangerAlignment level)
    ###  => SangerContig, SangerRead level directory will be created recursively
    ### ------------------------------------------------------------------------
    outputDirSA <- file.path(outputDir, "SangerAlignment")
    if (!dir.exists(outputDirSA)) {
        suppressWarnings(dir.create(outputDirSA, recursive = TRUE))
    }
    rootDir <- system.file(package = "sangeranalyseR")
    originRmd <- file.path(rootDir, "rmd", "SangerAlignment_Report.Rmd")
    outputHtml <- file.path(outputDirSA, "SangerAlignment_Report.html")

    # Start for loop
    if (includeSangerContig) {
        contigsFN <- lapply(object@contigList, function (objContig) {
            log_info("!!! outputHtml: ", outputHtml)
            generateReportSC(objContig, outputDir = outputDirSA,
                             includeSangerRead = includeSangerRead,
                             colors=colors,
                             navigationAlignmentFN = outputHtml)
        })
    } else {
        contigsFN <- NULL
    }
    res <- render(input = originRmd,
                  output_dir = outputDirSA,
                  params = list(SangerAlignment = object,
                                outputDir = outputDirSA,
                                contigsFN = contigsFN, 
                                colors = colors))
    return(outputHtml)
})
