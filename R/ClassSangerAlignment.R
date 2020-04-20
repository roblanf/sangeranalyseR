#' @export
setOldClass("phylo")

setClassUnion("DNAStringORNULL", c("DNAString", "NULL"))
setClassUnion("DNAStringSetORNULL", c("DNAStringSet", "NULL"))

#' @title SangerAlignment
#'
#' @description  An S4 class containing SangerContigs lists and contigs alignment results which corresponds to a final alignment in Sanger sequencing.
#'
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot fastaFileName If \code{inputSource} is \code{"FASTA"}, then this value has to be the name of the FASTA file; if \code{inputSource} is \code{"ABIF"}, then this value is \code{""} by default.
#' @slot namesConversionCSV The file path to the CSV file that provides read names that follow the naming regulation. If \code{inputSource} is \code{"FASTA"}, then users need to prepare the csv file or make sure the original names inside FASTA file are valid; if \code{inputSource} is \code{"ABIF"}, then this value is \code{NULL} by default.
#' @slot parentDirectory If \code{inputSource} is \code{"ABIF"}, then this value is the path of the parent directory storing all reads in ABIF format you wish to analyse and cannot be NULL. In SangerAlignment, all reads in subdirectories will be scanned recursively. If \code{inputSource} is \code{"FASTA"}, then this value is \code{NULL} by default.
#' @slot suffixForwardRegExp The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be \code{"_F.ab1"}.
#' @slot suffixReverseRegExp The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be \code{"_R.ab1"}.
#' @slot trimmingMethodSA The read trimming method for all SangerRead S4 instances in SangerAlignment. The value must be \code{"M1"} (the default) or \code{'M2'}. All SangerReads must have the same trimming method.
#' @slot minFractionCallSA Minimum fraction of the sequences required to call a consensus sequence for SangerAlignment at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
#' @slot maxFractionLostSA Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerAlignment (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
#' @slot geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @slot refAminoAcidSeq An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is \code{""}.
#' @slot contigList A list storing all SangerContigs S4 instances.
#' @slot contigsConsensus The consensus read of all SangerContig S4 instances in DNAString object.
#' @slot contigsAlignment The alignment of all SangerContig S4 instances with the called consensus sequence in DNAStringSet object. Users can use BrowseSeqs() to view the alignment.
#' @slot contigsTree A phylo instance returned by bionj function in ape package. It can be used to draw the tree.
#'
#' @name SangerAlignment-class
#' @exportClass SangerAlignment
#'
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerRead.R
#' @examples
#' ## Input From ABIF file format (Regex)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica', 'ACHLO')
#' suffixForwardRegExp <- "_[0-9]*_F.ab1"
#' suffixReverseRegExp <- "_[0-9]*_R.ab1"
#' sangerAlignment <- new("SangerAlignment",
#'                        inputSource           = "ABIF",
#'                        parentDirectory       = parentDir,
#'                        suffixForwardRegExp   = suffixForwardRegExp,
#'                        suffixReverseRegExp   = suffixReverseRegExp,
#'                        refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                        TrimmingMethod        = "M1",
#'                        M1TrimmingCutoff      = 0.0001,
#'                        M2CutoffQualityScore  = NULL,
#'                        M2SlidingWindowSize   = NULL,
#'                        baseNumPerRow         = 100,
#'                        heightPerRow          = 200,
#'                        signalRatioCutoff     = 0.33,
#'                        showTrimmed           = TRUE,
#'                        processorsNum         = 2)
#'
#' ## Input From ABIF file format (Csv three column)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica', 'ACHLO')
#' namesConversionCSV <- file.path(rawDataDir, "ab1", "SangerAlignment",
#' "names_conversion.csv")
#' sangerAlignment <- new("SangerAlignment",
#'                        inputSource           = "ABIF",
#'                        parentDirectory       = parentDir,
#'                        namesConversionCSV    = namesConversionCSV,
#'                        refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                        TrimmingMethod        = "M1",
#'                        M1TrimmingCutoff      = 0.0001,
#'                        M2CutoffQualityScore  = NULL,
#'                        M2SlidingWindowSize   = NULL,
#'                        baseNumPerRow         = 100,
#'                        heightPerRow          = 200,
#'                        signalRatioCutoff     = 0.33,
#'                        showTrimmed           = TRUE,
#'                        processorsNum         = 2)
#'
#' ## Input From FASTA file format (No Csv - Regex)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerAlignment", "Sanger_all_reads.fa")
#' suffixForwardRegExpFa <- "_[0-9]*_F$"
#' suffixReverseRegExpFa <- "_[0-9]*_R$"
#' sangerAlignmentFa <- new("SangerAlignment",
#'                          inputSource           = "FASTA",
#'                          fastaFileName         = fastaFN,
#'                          suffixForwardRegExp   = suffixForwardRegExpFa,
#'                          suffixReverseRegExp   = suffixReverseRegExpFa,
#'                          refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                          processorsNum         = 2)
#'
#' ## Input From FASTA file format (Csv three column method)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerAlignment", "Sanger_all_reads.fa")
#' namesConversionCSV <- file.path(rawDataDir, "fasta",
#'                                 "SangerAlignment", "names_conversion.csv")
#' sangerAlignmentFa <- new("SangerAlignment",
#'                          inputSource           = "FASTA",
#'                          fastaFileName         = fastaFN,
#'                          namesConversionCSV    = namesConversionCSV,
#'                          refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                          processorsNum         = 2)
setClass("SangerAlignment",
         # Users need to name their ab1 files in a systematic way. Here is the
         # regulation:
         #  1. Naming: XXXXX_F[0-9]*.ab1 / XXXXX_R[0-9]*.ab1
         #        For reads in same contig, XXXXX must be same.
         #  2. Users can set
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerAlignment'
         ### -------------------------------------------------------------------
         representation(inputSource                 = "character",
                        fastaFileName               = "characterORNULL",
                        namesConversionCSV          = "characterORNULL",
                        parentDirectory             = "characterORNULL",
                        suffixForwardRegExp         = "characterORNULL",
                        suffixReverseRegExp         = "characterORNULL",
                        trimmingMethodSA            = "character",
                        minFractionCallSA           = "numeric",
                        maxFractionLostSA           = "numeric",
                        geneticCode                 = "character",
                        refAminoAcidSeq             = "character",
                        contigList                  = "list",
                        contigsConsensus            = "DNAStringORNULL",
                        contigsAlignment            = "DNAStringSetORNULL",
                        contigsTree                 = "phylo"
         ),
)

### ============================================================================
### Overwrite initialize for 'SangerContig' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerAlignment",
          function(.Object,
                   inputSource            = "ABIF",
                   fastaFileName          = NULL,
                   namesConversionCSV     = NULL,
                   parentDirectory        = NULL,
                   suffixForwardRegExp    = NULL,
                   suffixReverseRegExp    = NULL,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL,
                   baseNumPerRow          = 100,
                   heightPerRow           = 200,
                   signalRatioCutoff      = 0.33,
                   showTrimmed            = TRUE,
                   refAminoAcidSeq        = "",
                   minReadsNum            = 2,
                   minReadLength          = 20,
                   minFractionCall        = 0.5,
                   maxFractionLost        = 0.5,
                   geneticCode            = GENETIC_CODE,
                   acceptStopCodons       = TRUE,
                   readingFrame           = 1,
                   minFractionCallSA      = 0.5,
                   maxFractionLostSA      = 0.5,
                   processorsNum          = NULL) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()
    ##### ------------------------------------------------------------------
    ##### Input parameter prechecking for SangerContig parameter
    ##### ------------------------------------------------------------------
    errors <- checkMinReadsNum(minReadsNum, errors)
    errors <- checkMinReadLength(minReadLength, errors)
    errors <- checkMinFractionCall(minFractionCall, errors)
    errors <- checkMaxFractionLost(maxFractionLost, errors)
    errors <- checkGeneticCode(geneticCode, errors)
    errors <- checkAcceptStopCodons(acceptStopCodons, errors)
    errors <- checkReadingFrame(readingFrame, errors)
    ##### ----------------------------------------------------------------------
    ##### Input parameter prechecking for processorsNum
    ##### ----------------------------------------------------------------------
    errors <- checkProcessorsNum(processorsNum, errors)
    log_info('*************************************************')
    log_info('**** Start creating SangerAlignment instance ****')
    log_info('*************************************************')
    if (length(errors) == 0) {
        processorsNum <- getProcessors (processorsNum)
        if (inputSource == "ABIF") {
            ### ----------------------------------------------------------------
            ##### 'parentDirectory' prechecking
            ### ----------------------------------------------------------------
            errors <- checkParentDirectory (parentDirectory, errors)
            ### ----------------------------------------------------------------
            ##### Input parameter prechecking for TrimmingMethod.
            ### ----------------------------------------------------------------
            errors <- checkTrimParam(TrimmingMethod, M1TrimmingCutoff,
                                     M2CutoffQualityScore,
                                     M2SlidingWindowSize, errors)
            ### ----------------------------------------------------------------
            ##### Input parameter prechecking for ChromatogramParam
            ### ----------------------------------------------------------------
            errors <- checkBaseNumPerRow (baseNumPerRow, errors)
            errors <- checkHeightPerRow (baseNumPerRow, errors)
            errors <- checkSignalRatioCutoff (signalRatioCutoff, errors)
            errors <- checkShowTrimmed (showTrimmed, errors)
            if (length(errors) != 0) {
                log_error(errors)
            }
            ab1RegexChecker <- is.null(namesConversionCSV) &&
                !is.null(suffixForwardRegExp) && !is.null(suffixReverseRegExp)
            ab1CSVChecker <- !is.null(namesConversionCSV)
            trimmingMethodSA <- TrimmingMethod
            if (ab1RegexChecker) {
                errors <- 
                    checkNamesConversionCSV(TRUE, parentDirectory, 
                                            fastaFileName, namesConversionCSV, 
                                            inputSource, "ab1Regex", errors)
                if (length(errors) != 0) {
                    log_error(errors)
                }
                ### ------------------------------------------------------------
                ### Automatically finding contig name by forward&reverse suffix
                ### ------------------------------------------------------------
                log_info("**** You are using Regex Method to group AB1 files!")
                parentDirFiles <- list.files(parentDirectory, recursive = TRUE)
                forwardSelectInputFiles <- 
                    parentDirFiles[grepl(suffixForwardRegExp, parentDirFiles)]
                reverseSelectInputFiles <- 
                    parentDirFiles[grepl(suffixReverseRegExp, parentDirFiles)]
                
                # Find possible consensus Name for forward and reverse reads
                forwardContigName <-
                    unlist(str_split(forwardSelectInputFiles, 
                                     suffixForwardRegExp, n = Inf, 
                                     simplify = FALSE))[c(TRUE, FALSE)]
                reverseContigName <-
                    unlist(str_split(reverseSelectInputFiles, 
                                     suffixReverseRegExp, n = Inf, 
                                     simplify = FALSE))[c(TRUE, FALSE)]

                contigName <- union(forwardContigName, reverseContigName)
                contigNumber <- length(contigName)

                # Create contig for all list of contigNumber
                ### ----------------------------------------------------------------
                ##### Creating each SangerContig (store as SangerContigList)
                ### ----------------------------------------------------------------
                SangerContigList <-
                    sapply(contigName,
                           function(eachConsRead) {
                               insideDirName<- dirname(eachConsRead)
                               insideContigName <- basename(eachConsRead)
                               newSangerContig <- 
                                   new("SangerContig",
                                       inputSource          = inputSource,
                                       fastaFileName        = fastaFileName,
                                       parentDirectory      =
                                           file.path(parentDirectory, insideDirName),
                                       contigName           = insideContigName,
                                       suffixForwardRegExp  = suffixForwardRegExp,
                                       suffixReverseRegExp  = suffixReverseRegExp,
                                       TrimmingMethod       = TrimmingMethod,
                                       M1TrimmingCutoff     = M1TrimmingCutoff,
                                       M2CutoffQualityScore = M2CutoffQualityScore,
                                       M2SlidingWindowSize  = M2SlidingWindowSize,
                                       baseNumPerRow        = baseNumPerRow,
                                       heightPerRow         = heightPerRow,
                                       signalRatioCutoff    = signalRatioCutoff,
                                       showTrimmed          = showTrimmed,
                                       refAminoAcidSeq      = refAminoAcidSeq,
                                       minReadsNum          = minReadsNum,
                                       minReadLength        = minReadLength,
                                       minFractionCall      = minFractionCall,
                                       maxFractionLost      = maxFractionLost,
                                       geneticCode          = geneticCode,
                                       acceptStopCodons     = acceptStopCodons,
                                       readingFrame         = readingFrame,
                                       processorsNum        = processorsNum, 
                                       logLevel             = FALSE)
                               forwardNumber <- 
                                   length(newSangerContig@forwardReadList)
                               reverseNumber <- 
                                   length(newSangerContig@reverseReadList)
                               readNumber <- forwardNumber + reverseNumber
                               if (readNumber >= minReadsNum && 
                                   readNumber >= 2) {
                                   newSangerContig
                               } else {
                                   NULL
                               } 
                           })
            } else if (ab1CSVChecker) {
                errors <- 
                    checkNamesConversionCSV(TRUE, parentDirectory, 
                                            fastaFileName, namesConversionCSV, 
                                            inputSource, "ab1CSV", errors)
                if (length(errors) != 0) {
                    log_error(errors)
                }
                log_info("**** You are using CSV Name Conversion Method ",
                        "to group AB1 files!")
                csvFile <- read.csv(namesConversionCSV, header = TRUE)
                log_info("**** Contig number in your Csv file is ", length(unique(csvFile$contig)))
                contigNames <- as.character(unique(csvFile$contig))
                SangerContigList <- sapply(contigNames, function(contigName) {
                    newSangerContig <- 
                        new("SangerContig",
                            inputSource          = inputSource,
                            fastaFileName        = fastaFileName,
                            namesConversionCSV   = namesConversionCSV,
                            parentDirectory      = parentDirectory,
                            contigName           = contigName,
                            refAminoAcidSeq      = refAminoAcidSeq,
                            minReadsNum          = minReadsNum,
                            minReadLength        = minReadLength,
                            minFractionCall      = minFractionCall,
                            maxFractionLost      = maxFractionLost,
                            geneticCode          = geneticCode,
                            acceptStopCodons     = acceptStopCodons,
                            readingFrame         = readingFrame,
                            processorsNum        = processorsNum, 
                            logLevel             = FALSE)
                    forwardNumber <- length(newSangerContig@forwardReadList)
                    reverseNumber <- length(newSangerContig@reverseReadList)
                    if ((forwardNumber + reverseNumber) >= minReadsNum) {
                        newSangerContig
                    } else {
                        NULL
                    } 
                })
            }
        } else if (inputSource == "FASTA") {
            errors <- checkFastaFileName(fastaFileName, errors)
            if(length(errors) != 0) {
                log_error(errors)
            }
            csvRegexChecker <- is.null(namesConversionCSV) &&
                !is.null(suffixForwardRegExp) &&
                !is.null(suffixReverseRegExp)
            csvCSVChecker <- !is.null(namesConversionCSV)
            readFasta <- read.fasta(fastaFileName, as.string = TRUE)
            trimmingMethodSA <- ""
            if (csvRegexChecker) {
                errors <- 
                    checkNamesConversionCSV(TRUE, parentDirectory, 
                                            fastaFileName, namesConversionCSV, 
                                            inputSource, "csvRegex", errors)
                if (length(errors) != 0) {
                    log_error(errors)
                }
                log_info("**** You are using Regex Method ",
                        "to group reads in FASTA file (No CSV file)!")
                readFastaNames <- names(readFasta)
                forwardSelectInputFiles <- readFastaNames[grepl(suffixForwardRegExp,
                                                                readFastaNames)]
                reverseSelectInputFiles <- readFastaNames[grepl(suffixReverseRegExp,
                                                                readFastaNames)]
                # Find possible consensus Name for forward and reverse reads
                forwardContigName <-
                    unlist(str_split(forwardSelectInputFiles, suffixForwardRegExp,
                                     n = Inf, simplify = FALSE))[c(TRUE, FALSE)]
                reverseContigName <-
                    unlist(str_split(reverseSelectInputFiles, suffixReverseRegExp,
                                     n = Inf, simplify = FALSE))[c(TRUE, FALSE)]
                contigNames <- union(forwardContigName, reverseContigName)
                SangerContigList <- sapply(contigNames, function(contigName) {
                    newSangerContig <- 
                        new("SangerContig",
                            inputSource          = inputSource,
                            fastaFileName        = fastaFileName,
                            namesConversionCSV   = namesConversionCSV,
                            parentDirectory      = parentDirectory,
                            contigName           = contigName,
                            suffixForwardRegExp  = suffixForwardRegExp,
                            suffixReverseRegExp  = suffixReverseRegExp,
                            refAminoAcidSeq      = refAminoAcidSeq,
                            minReadsNum          = minReadsNum,
                            minReadLength        = minReadLength,
                            minFractionCall      = minFractionCall,
                            maxFractionLost      = maxFractionLost,
                            geneticCode          = geneticCode,
                            acceptStopCodons     = acceptStopCodons,
                            readingFrame         = readingFrame,
                            processorsNum        = processorsNum, 
                            logLevel             = FALSE)
                    forwardNumber <- length(newSangerContig@forwardReadList)
                    reverseNumber <- length(newSangerContig@reverseReadList)
                    if ((forwardNumber + reverseNumber) >= minReadsNum) {
                        newSangerContig
                    } else {
                        NULL
                    } 
                })
            } else if (csvCSVChecker) {
                errors <- 
                    checkNamesConversionCSV(TRUE, parentDirectory, 
                                            fastaFileName, namesConversionCSV, 
                                            inputSource, "csvCSV", errors)
                if (length(errors) != 0) {
                    log_error(errors)
                }
                log_info("**** You are using CSV Name Conversion Method ",
                        "to group reads in FASTA file (with CSV file)!")
                csvFile <- read.csv(namesConversionCSV, header = TRUE)
                contigNames <- unique(as.character(csvFile$contig))
                SangerContigList <- sapply(contigNames, function(contigName) {
                    newSangerContig <- 
                        new("SangerContig",
                            inputSource          = inputSource,
                            fastaFileName        = fastaFileName,
                            namesConversionCSV   = namesConversionCSV,
                            parentDirectory      = parentDirectory,
                            contigName           = contigName,
                            suffixForwardRegExp  = suffixForwardRegExp,
                            suffixReverseRegExp  = suffixReverseRegExp,
                            refAminoAcidSeq      = refAminoAcidSeq,
                            minReadsNum          = minReadsNum,
                            minReadLength        = minReadLength,
                            minFractionCall      = minFractionCall,
                            maxFractionLost      = maxFractionLost,
                            geneticCode          = geneticCode,
                            acceptStopCodons     = acceptStopCodons,
                            readingFrame         = readingFrame,
                            processorsNum        = processorsNum,
                            logLevel             = FALSE)
                    forwardNumber <- length(newSangerContig@forwardReadList)
                    reverseNumber <- length(newSangerContig@reverseReadList)
                    if ((forwardNumber + reverseNumber) >= minReadsNum) {
                        newSangerContig
                    } else {
                        NULL
                    } 
                })
            }
        }
        # message("SangerContigList length: ", length(SangerContigList))
        SangerContigList <- Filter(Negate(is.null), SangerContigList)
        # message("SangerContigList length: ", length(SangerContigList))
        acResult <- alignContigs(SangerContigList, geneticCode,
                                 refAminoAcidSeq, minFractionCallSA,
                                 maxFractionLostSA, processorsNum)
        consensus <- acResult[["consensus"]]
        aln <- acResult[["aln"]]
        aln.tree <- acResult[["aln.tree"]]
        log_success("  >> 'SangerAlignment' S4 instance is created !!")
    } else {
        log_error(errors)
    }
    callNextMethod(.Object,
                   inputSource           = inputSource,
                   fastaFileName         = fastaFileName,
                   namesConversionCSV    = namesConversionCSV,
                   parentDirectory       = parentDirectory,
                   suffixForwardRegExp   = suffixForwardRegExp,
                   suffixReverseRegExp   = suffixReverseRegExp,
                   trimmingMethodSA      = trimmingMethodSA,
                   contigList            = SangerContigList,
                   minFractionCallSA     = minFractionCallSA,
                   maxFractionLostSA     = maxFractionLostSA,
                   geneticCode           = geneticCode,
                   contigsConsensus      = consensus,
                   refAminoAcidSeq       = refAminoAcidSeq,
                   contigsAlignment      = aln,
                   contigsTree           = aln.tree)
})
