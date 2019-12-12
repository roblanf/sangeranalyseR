setMethod("writeFASTA", "SangerAlignment", function(obj, outputDir, compress,
                                                 compression_level,
                                                 selection = "all") {
    ### ------------------------------------------------------------------------
    ### selection can be 'all', 'alignment' 'allReads' 'contig'
    ### ------------------------------------------------------------------------
    if (selection != "all" && selection != "alignment" &&
        selection != "contigs" && selection != "consensusRead" &&
        selection != "allReads") {
        stop(paste0("\nSelection must be 'all', ",
                    "'alignment' 'contigs','consensusRead' or 'allReads'."))
    }
    message("Start to write 'SangerAlignment' to FASTA format ...")
    ### ------------------------------------------------------------------------
    ### Writing 'contigs alignment' result to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "alignment") {
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
    if (selection == "all" || selection == "contigs") {
        message("\n    >> Writing 'contigs' to FASTA ...")
        contigsList <- sapply(obj@contigList, function(contig) {
            contig@contigSeq
        })
        contigsListDNASet<- DNAStringSet(contigsList)
        writeXStringSet(contigsListDNASet,
                        file.path(outputDir, "Sanger_contigs.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    ### ------------------------------------------------------------------------
    ### Writing 'contigs consensus read' to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "consensusRead") {
        message("\n    >> Writing 'consensusRead' to FASTA ...")
        contigsConsensusDNASet<- DNAStringSet(obj@contigsConsensus)
        names(contigsConsensusDNASet) <- "Sanger Consensus Read"
        writeXStringSet(contigsConsensusDNASet,
                        file.path(outputDir, "Sanger_consensus_read.fa"),
                        compress = compress,
                        compression_level = compression_level)
    }

    ### ------------------------------------------------------------------------
    ### Writing 'all reads' to FASTA file
    ### ------------------------------------------------------------------------
    if (selection == "all" || selection == "allReads") {
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

#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' newAlignment <- SangerAlignment(
#'                     parentDirectory       = rawDataDir,
#'                     suffixForwardRegExp   = suffixForwardRegExp,
#'                     suffixReverseRegExp   = suffixReverseRegExp,
#'                     refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                     TrimmingMethod        = "M1",
#'                     M1TrimmingCutoff      = 0.0001,
#'                     M2CutoffQualityScore  = NULL,
#'                     M2SlidingWindowSize   = NULL,
#'                     baseNumPerRow         = 100,
#'                     heightPerRow          = 200,
#'                     signalRatioCutoff     = 0.33,
#'                     showTrimmed           = TRUE)
#' RShinyCSSet <- launchAppSangerAlignment(newAlignment)
setMethod("launchAppSangerAlignment", "SangerAlignment",
          function(obj, outputDir = NULL) {
              ### ------------------------------------------------------------------------
              ### Checking SangerAlignment input parameter is a list containing
              ### one S4 object.
              ### ------------------------------------------------------------------------
              if (is.null(outputDir)) {
                  outputDir <- tempdir()
                  suppressWarnings(dir.create(outputDir))
              }
              if (dir.exists(outputDir)) {
                  shinyOptions(sangerAlignment = list(obj))
                  shinyOptions(shinyDirectory = outputDir)
                  newSangerAlignment <- shinyApp(SangerAlignmentUI, SangerAlignmentServer)
                  return(newSangerAlignment)
              } else {
                  stop("'", outputDir, "' is not valid. Please check again")
              }
})


setMethod("generateReport", "SangerAlignment", function(obj, outputDir) {
    if (is.null(outputDir)) {
        outputDir <- tempdir()
        suppressWarnings(dir.create(outputDir))
    }

    outputDirSA <- file.path(outputDir, "SangerAlignment")
    if (!dir.exists(outputDirSA)) {
        suppressWarnings(dir.create(outputDirSA, recursive = TRUE))
    }

    outputHtml <- file.path(outputDirSA, "SangerAlignment_Report.html")

    # Start for loop
    contigsFN <- sapply(obj@contigList, function (objContig) {
        message("!!! outputHtml: ", outputHtml)
        outputDirSC <- file.path(outputDirSA, objContig@contigName)
        if (!dir.exists(outputDirSC)) {
            suppressWarnings(dir.create(outputDirSC, recursive = TRUE))
        }
        generateReport(objContig, outputDir = outputDirSA,
                       navigationAlignmentFN = outputHtml)
    })

    res <- render(input = "/Users/chaokuan-hao/Documents/ANU_2019_Semester_2/Lanfear_Lab/sangeranalyseR/vignettes/SangerAlignment_Report.Rmd",
                  output_dir = outputDirSA,
                  params = list(SangerAlignment = obj,
                                outputDir = outputDir,
                                contigsFN = contigsFN))


    # rootDir <- system.file(package = "sangeranalyseR")
    # originRmd <- file.path(rootDir, "vignettes", "SangerContig_Report.Rmd")
    # outputHtml <- file.path(outputDirSC, "SangerContig_Report.html")
    #
    # forwardReads <- objContig@forwardReadList
    # reverseReads <- objContig@reverseReadList
    #
    # forwardReadFN <- sapply(forwardReads, generateReport, outputDir, outputHtml)
    # reverseReadFN <- sapply(reverseReads, generateReport, outputDir, outputHtml)
    #
    # res <- render(input = "/Users/chaokuan-hao/Documents/ANU_2019_Semester_2/Lanfear_Lab/sangeranalyseR/vignettes/SangerContig_Report.Rmd",
    #               output_dir = outputDirSC,
    #               params = list(SangerContig = objContig,
    #                             outputDir = outputDir,
    #                             forwardReadFN = forwardReadFN,
    #                             reverseReadFN = reverseReadFN))
    return(contigsFN)
})
