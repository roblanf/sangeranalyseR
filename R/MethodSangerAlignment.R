## =============================================================================
## Updating quality parameters for SangerAlignment object.
## =============================================================================
#' A SangerAlignment method which updates QualityReport parameter for each the SangerRead instance inside SangerAlignment.
#'
#' @title updateQualityParam
#' @name SangerAlignment-class-updateQualityParam
#' @rdname SangerAlignment-Method
#'
#' @docType methods
#' @examples
#' \dontrun{load("data/sangerAlignment.RData")
#' updateQualityParam(sangerAlignment,
#'                    TrimmingMethod         = "M2",
#'                    M1TrimmingCutoff       = NULL,
#'                    M2CutoffQualityScore   = 40,
#'                    M2SlidingWindowSize    = 15)}
setMethod("updateQualityParam",  "SangerAlignment",
          function(object,
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
            newContigList <-
                sapply(object@contigList,
                       function(contig) {
                           contig <-
                               updateQualityParam(
                                   contig,
                                   TrimmingMethod         = TrimmingMethod,
                                   M1TrimmingCutoff       = M1TrimmingCutoff,
                                   M2CutoffQualityScore   = M2CutoffQualityScore,
                                   M2SlidingWindowSize    = M2SlidingWindowSize)
                       })
            object@contigList <- newContigList
            object@trimmingMethodSA <- TrimmingMethod
            acResult <- alignContigs(object@contigList, object@geneticCode,
                                     object@refAminoAcidSeq,
                                     object@minFractionCallSA,
                                     object@maxFractionLostSA,
                                     getProcessors(NULL))
            object@contigsConsensus <- acResult[["consensus"]]
            object@contigsAlignment <- acResult[["aln"]]
            object@contigsTree <- acResult[["aln.tree"]]
            return(object)
        } else {
            stop(errors)
        }
    } else if (object@inputSource == "FASTA") {
        message("SangerAlignment with 'FASTA' inputSource ",
                "cannot update quality parameters")
    }
})

#' A SangerAlignment method which launches Shiny app for SangerAlignment instance.
#'
#' @title launchAppSA
#' @name SangerAlignment-class-launchAppSA
#' @rdname SangerAlignment-Method
#'
#' @docType methods
#' @examples
#' \dontrun{load("data/sangerAlignment.RData")
#' RShinySA <- launchAppSA(sangerAlignment)}
setMethod("launchAppSA", "SangerAlignment", function(obj, outputDir = NULL) {
    if (obj@inputSource == "ABIF") {
        ### ------------------------------------------------------------------------
        ### Checking SangerAlignment input parameter is a list containing
        ### one S4 object.
        ### ------------------------------------------------------------------------
        if (is.null(outputDir)) {
            outputDir <- tempdir()
            suppressWarnings(dir.create(outputDir, recursive = TRUE))
        }
        message(">>> outputDir : ", outputDir)

        if (dir.exists(outputDir)) {
            shinyOptions(sangerAlignment = list(obj))
            shinyOptions(shinyDirectory = outputDir)
            newSangerAlignment <- shinyApp(SangerAlignmentUI, SangerAlignmentServer)
            return(newSangerAlignment)
        } else {
            stop("'", outputDir, "' is not valid. Please check again")
        }
    } else if (obj@inputSource == "FASTA") {
        message("SangerAlignment with 'FASTA' inputSource ",
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
#' @rdname SangerAlignment-Method
#'
#' @docType methods
#' @examples
#' \dontrun{load("data/sangerAlignment.RData")
#' writeFastaSA(sangerAlignment, "/Users/chaokuan-hao/Desktop/sangeranalyseR_fasta/SangerAlignment")}
setMethod("writeFastaSA", "SangerAlignment", function(obj, outputDir, compress,
                                                 compression_level,
                                                 selection = "all") {
    ### ------------------------------------------------------------------------
    ### selection can be 'all', 'alignment' 'all_reads' 'contig'
    ### ------------------------------------------------------------------------
    if (selection != "all" && selection != "contigs_alignment" &&
        selection != "contigs_unalignment" && selection != "all_reads") {
        stop(paste0("\nSelection must be 'all', 'contigs_alignment',",
                    " 'contigs_unalignment' or 'all_reads'."))
    }
    ### ------------------------------------------------------------------------
    ### Make sure the input directory is not NULL
    ### ------------------------------------------------------------------------
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir, recursive = TRUE))
    }
    message(">>> outputDir : ", outputDir)
    message("Start to write 'SangerAlignment' to FASTA format ...")
    ### ------------------------------------------------------------------------
    ### Writing 'contigs alignment' result to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "contigs_alignment") {
        message("\n    >> Writing 'alignment' to FASTA ...")
        alignmentObject = obj@contigsAlignment
        writeXStringSet(alignmentObject,
                        file.path(outputDir, "Sanger_contigs_alignment.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    ### ------------------------------------------------------------------------
    ### Writing 'all contigs' to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "contigs_unalignment") {
        message("\n    >> Writing 'contigs' to FASTA ...")
        contigsList <- sapply(obj@contigList, function(contig) {
            contig@contigSeq
        })
        contigsListDNASet<- DNAStringSet(contigsList)
        writeXStringSet(contigsListDNASet,
                        file.path(outputDir, "Sanger_contigs_unalignment.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    # ### ------------------------------------------------------------------------
    # ### Writing 'contigs consensus read' to FASTA file
    # ### ------------------------------------------------------------------------
    # if (selection == "all" || selection == "consensusRead") {
    #     message("\n    >> Writing 'consensusRead' to FASTA ...")
    #     contigsConsensusDNASet<- DNAStringSet(obj@contigsConsensus)
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
        message("\n    >> Writing all single reads to FASTA ...")
        fRDNASet <- sapply(obj@contigList, function(contig) {
            fRDNAStringSet <- sapply(contig@forwardReadList, function(forwardRead) {
                primaryDNA <- as.character(forwardRead@primarySeq)
                if (obj@inputSource == "ABIF") {
                    trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
                    trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
                    primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
                }
                return(primaryDNA)
            })
            names(fRDNAStringSet) <- basename(names(fRDNAStringSet))
            fRDNAStringSet
        })
        rRDNASet <- sapply(obj@contigList, function(contig) {
            rRDNAStringSet <- sapply(contig@reverseReadList, function(reverseRead) {
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
            rRDNAStringSet
        })
        allDNASet <- DNAStringSet(c(fRDNASet, rRDNASet))
        writeXStringSet(allDNASet,
                        file.path(outputDir, "Sanger_all_trimmed_reads.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    message("\nFinish writing 'SangerAlignment' to FASTA format")
})

## =============================================================================
## Generating report for SangerContig
## =============================================================================
#' A SangerAlignment method which generates final reports of the SangerContig instance.
#'
#' @title generateReportSA
#' @name SangerAlignment-class-generateReportSA
#' @rdname SangerAlignment-Method
#'
#' @docType methods
#' @examples
#' \dontrun{load("data/sangerAlignment.RData")
#' generateReportSA(sangerAlignment)}
setMethod("generateReportSA", "SangerAlignment",
          function(obj, outputDir,
                   includeSangerContig = TRUE,
                   includeSangerRead = TRUE) {

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
        contigsFN <- sapply(obj@contigList, function (objContig) {
            message("!!! outputHtml: ", outputHtml)
            generateReportSC(objContig, outputDir = outputDirSA,
                             navigationAlignmentFN = outputHtml,
                             includeSangerRead = includeSangerRead)
        })
    } else {
        contigsFN <- NULL
    }
    res <- render(input = originRmd,
                  output_dir = outputDirSA,
                  params = list(SangerAlignment = obj,
                                outputDir = outputDirSA,
                                contigsFN = contigsFN))
    return(outputHtml)
})
