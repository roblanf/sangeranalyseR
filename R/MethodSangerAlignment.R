#' @title launchAppSA
#' @name SangerAlignment-class-launchAppSA
#' @rdname SangerAlignment-class-launchAppSA
#' @examples
#' load("data/sangerAlignment.RData")
#' RShinySA <- launchAppSA(sangerAlignment)
setMethod("launchAppSA", "SangerAlignment", function(obj, outputDir = NULL) {
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
})

## =============================================================================
## Writing primary sequence into FASTA format
## =============================================================================
#' @title writeFastaSA
#' @name SangerAlignment-class-writeFastaSA
#' @rdname SangerAlignment-class-writeFastaSA
#' @examples
#' load("data/sangerAlignment.RData")
#' writeFastaSA(sangerAlignment)
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
                trimmedStartPos <- forwardRead@QualityReport@trimmedStartPos
                trimmedFinishPos <- forwardRead@QualityReport@trimmedFinishPos
                primaryDNA <- as.character(forwardRead@primarySeq)
                substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
            })
            names(fRDNAStringSet) <- basename(names(fRDNAStringSet))
            fRDNAStringSet
        })
        rRDNASet <- sapply(obj@contigList, function(contig) {
            rRDNAStringSet <- sapply(contig@reverseReadList, function(reverseRead) {
                trimmedStartPos <- reverseRead@QualityReport@trimmedStartPos
                trimmedFinishPos <- reverseRead@QualityReport@trimmedFinishPos
                primaryDNA <- as.character(reverseRead@primarySeq)
                substr(primaryDNA, trimmedStartPos, trimmedFinishPos)
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
#' @title generateReportSA
#' @name SangerAlignment-class-generateReportSA
#' @rdname SangerAlignment-class-generateReportSA
#' @examples
#' load("data/sangerAlignment.RData")
#' generateReportSA(sangerAlignment)
setMethod("generateReportSA", "SangerAlignment",
          function(obj, outputDir,
                   includeSangerContig = TRUE,
                   includeSangerRead = TRUE) {
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
