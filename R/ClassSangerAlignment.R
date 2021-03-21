#' @export
setOldClass("phylo")

setClassUnion("DNAStringORNULL", c("DNAString", "NULL"))
setClassUnion("DNAStringSetORNULL", c("DNAStringSet", "NULL"))

#' @title SangerAlignment
#'
#' @description  An S4 class containing SangerContigs lists and contigs alignment results which corresponds to a final alignment in Sanger sequencing.
#'
#' @slot objectResults
#' 
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot processMethod The method to create a contig from reads. The value is \code{"REGEX"}, \code{"CSV"}
#' 
#' 
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
#'                        processMethod         = "REGEX",
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
#'                        processMethod         = "CSV",
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
#'                          processMethod         = "REGEX"
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
#'                          processMethod         = "CSV"
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
         representation(objectResults               = "ObjectResults",
                        inputSource                 = "character",
                        processMethod               = "character",
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
                   printLevel             = "SangerAlignment",
                   inputSource            = "ABIF",
                   processMethod          = "REGEX",
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
    creationResult <- TRUE
    errors <- list(character(0), character(0))
    warnings <- c(character(0))
    readResultTableName <- c("readName","creationResult", "errorType", 
                             "errorMessage", "inputSource", "direction")
    readResultTable <- data.frame()
    ############################################################################
    ### First layer of pre-checking: SangerAlignment input parameter prechecking
    ############################################################################
    if (printLevel == "SangerAlignment") {
        errors <- checkInputSource(inputSource, errors[[1]], errors[[2]])
        errors <- checkProcessMethod(inputSource, processMethod, errors[[1]], errors[[2]])
        errors <- checkFastaFileName(inputSource, fastaFileName,
                                     errors[[1]], errors[[2]])
        errors <- checkNamesConversionCSV(processMethod, namesConversionCSV, 
                                          "SangerAlignment", suffixForwardRegExp,
                                          suffixReverseRegExp, inputSource, 
                                          errors[[1]], errors[[2]])
        errors <- checkRefAAS(refAminoAcidSeq, errors[[1]], errors[[2]])
        errors <- checkMinReadsNum(minReadsNum, errors[[1]], errors[[2]])
        errors <- checkMinReadLength(minReadLength, errors[[1]], errors[[2]])
        errors <- checkMinFractionCall(minFractionCall, errors[[1]], errors[[2]])
        errors <- checkMaxFractionLost(maxFractionLost, errors[[1]], errors[[2]])
        errors <- checkMinFractionCall(minFractionCallSA, errors[[1]], errors[[2]])
        errors <- checkMaxFractionLost(maxFractionLostSA, errors[[1]], errors[[2]])
        errors <- checkGeneticCode(geneticCode, errors[[1]], errors[[2]])
        errors <- checkAcceptStopCodons(acceptStopCodons, errors[[1]], errors[[2]])
        errors <- checkReadingFrame(readingFrame, errors[[1]], errors[[2]])
        errors <- checkProcessorsNum(processorsNum, errors[[1]], errors[[2]])
        if (inputSource == "ABIF") {
            ### --------------------------------------------------------------------
            ### 'ABIF' condition checking!
            ### --------------------------------------------------------------------
            ########################################################################
            ### Second layer of pre-checking: 'ABIF' condition checking!
            ########################################################################
            errors <- checkParentDirectory (parentDirectory, errors[[1]], errors[[2]])
            errors <- checkTrimParam(TrimmingMethod,
                                     M1TrimmingCutoff,
                                     M2CutoffQualityScore,
                                     M2SlidingWindowSize,
                                     errors[[1]], errors[[2]])
            errors <- checkBaseNumPerRow (baseNumPerRow, errors[[1]], errors[[2]])
            errors <- checkHeightPerRow (heightPerRow, errors[[1]], errors[[2]])
            errors <- checkSignalRatioCutoff (signalRatioCutoff, errors[[1]], errors[[2]])
            errors <- checkShowTrimmed (showTrimmed, errors[[1]], errors[[2]])
            trimmingMethodSA <- TrimmingMethod    
        } else if (inputSource == "FASTA") {
            ### --------------------------------------------------------------------
            ### 'FASTA' condition checking!
            ### --------------------------------------------------------------------
            ########################################################################
            ### Second layer of pre-checking: 'FASTA' condition checking!
            ########################################################################
            readFasta <- read.fasta(fastaFileName, as.string = TRUE)
            fastaNames <- names(readFasta)
        }
        if (processMethod=="CSV") {
            errors <- checkAb1FastaCsv(parentDirectory, fastaFileName, 
                                       namesConversionCSV, inputSource, 
                                       errors[[1]], errors[[2]])
        }
    }
    
    if (length(errors[[1]]) == 0 ) {
        log_info('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        log_info('%%%% Start creating SangerAlignment instance %%%%')
        log_info('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        processorsNum <- getProcessors (processorsNum)
        if (inputSource == "ABIF" && processMethod == "REGEX") {
            log_info("  >> You are using Regular Expression Method",
                     " to group AB1 files!")
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
            contigNames <- union(forwardContigName, reverseContigName)
            contigNumber <- length(contigNames)
            
            # Create contig for all list of contigNumber
            ### ----------------------------------------------------------------
            ##### Creating each SangerContig (store as SangerContigList)
            ### ----------------------------------------------------------------
            SangerContigList <-
                lapply(contigNames,
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
                               readNumber >= 1) {
                               newSangerContig
                           } else {
                               NULL
                           }
                       })
            names(SangerContigList) <- contigNames
        } else if (inputSource == "ABIF" && processMethod == "CSV") {
            log_info("**** You are using CSV Name Conversion Method ",
                     "to group AB1 files!")
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            log_info("**** Contig number in your Csv file is ", 
                     length(unique(csvFile$contig)))
            contigNames <- as.character(unique(csvFile$contig))
            SangerContigList <- lapply(contigNames, function(contigName) {
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
            names(SangerContigList) <- contigNames
        }
        
        if (inputSource == "FASTA" && processMethod == "REGEX") {
            log_info("**** You are using Regular Expression Method ",
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
            SangerContigList <- lapply(contigNames, function(contigName) {
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
            names(SangerContigList) <- contigNames
        } else if (inputSource == "FASTA" && processMethod == "CSV") {
            log_info("**** You are using CSV Name Conversion Method ",
                     "to group reads in FASTA file (with Csv file)!")
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            contigNames <- unique(as.character(csvFile$contig))
            SangerContigList <- lapply(contigNames, function(contigName) {
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
            names(SangerContigList) <- contigNames
        }
    }
    
    if (length(errors[[1]]) == 0 ) {
       
        # message("SangerContigList length: ", length(SangerContigList))
        SangerContigList <- Filter(Negate(is.null), SangerContigList)
        
        invisible(
            lapply(SangerContigList, function(contig) {
            lapply(contig@forwardReadList, function(read) {
                if (read@objectResults@creationResult) {
                    trimmedLen <- read@QualityReport@trimmedFinishPos -
                        read@QualityReport@trimmedStartPos
                    if (trimmedLen >= minReadLength) {
                        ### Success: readResultTable (SangerContig Level)
                        row <- data.frame(basename(read@readFileName),
                                          read@objectResults@creationResult,
                                          "None", "None",
                                          read@inputSource, read@readFeature)
                        names(row) <- readResultTableName
                        readResultTable <<- rbind(readResultTable, row)
                    } else {
                        msg <- paste0("  * >> ", read@readFileName,
                                      " is shorter than 'minReadLength' ",
                                      minReadLength,". This read is created but skipped!\n")
                        warnings <<- c(warnings, msg)
                        log_warn(msg)
                        ### Failed: readResultTable (SangerContig Level)
                        row <- data.frame(basename(read@readFileName),
                                          FALSE, "TRIMMED_READ_ERROR", msg,
                                          read@inputSource,  read@readFeature)
                        names(row) <- readResultTableName
                        readResultTable <<- rbind(readResultTable, row)
                    }
                } else {
                    ### Failed: readResultTable (SangerRead Level)
                    row <- data.frame(basename(read@readFileName),
                                      read@objectResults@creationResult,
                                      read@objectResults@errorTypes,
                                      read@objectResults@errorMessages,
                                      read@inputSource, read@readFeature)
                    names(row) <- readResultTableName
                    readResultTable <<- rbind(readResultTable, row)
                }
            })
            lapply(contig@reverseReadList, function(read) {
                if (read@objectResults@creationResult) {
                    trimmedLen <- read@QualityReport@trimmedFinishPos -
                        read@QualityReport@trimmedStartPos
                    if (trimmedLen >= minReadLength) {
                        ### Success: readResultTable (SangerContig Level)
                        row <- data.frame(basename(read@readFileName),
                                          read@objectResults@creationResult,
                                          "None", "None",
                                          read@inputSource, read@readFeature)
                        names(row) <- readResultTableName
                        readResultTable <<- rbind(readResultTable, row)
                    } else {
                        msg <- paste0("  * >> ", read@readFileName,
                                      " is shorter than 'minReadLength' ",
                                      minReadLength,". This read is created but skipped!\n")
                        warnings <<- c(warnings, msg)
                        log_warn(msg)
                        ### Failed: readResultTable (SangerContig Level)
                        row <- data.frame(basename(read@readFileName),
                                          FALSE, "TRIMMED_READ_ERROR", msg,
                                          read@inputSource,  read@readFeature)
                        names(row) <- readResultTableName
                        readResultTable <<- rbind(readResultTable, row)
                    }
                } else {
                    ### Failed: readResultTable (SangerRead Level)
                    row <- data.frame(basename(read@readFileName),
                                      read@objectResults@creationResult,
                                      read@objectResults@errorTypes,
                                      read@objectResults@errorMessages,
                                      read@inputSource, read@readFeature)
                    names(row) <- readResultTableName
                    readResultTable <<- rbind(readResultTable, row)
                }
            })
        })
        )


        
        
        
        
        
        
        
        # message("SangerContigList length: ", length(SangerContigList))
        
        # alignContigs <- function(SangerContigList, geneticCode, refAminoAcidSeq,
        #                          minFractionCallSA, maxFractionLostSA, processorsNum) {
        acResult <- alignContigs(SangerContigList, geneticCode,
                                 refAminoAcidSeq, minFractionCallSA,
                                 maxFractionLostSA, processorsNum)
        consensus <- acResult[["consensus"]]
        aln <- acResult[["aln"]]
        aln.tree <- acResult[["aln.tree"]]
        
        
        contigNum <- length(SangerContigList)
        
        # 100 reads detected
        # 12 contigs detected from [regular expression / csv file]
        # 25 forward reads assigned to 12 contigs according to [regular expression / csv file]
        # 3 reverse reads assigned to 2 contigs according to [regular expression / csv file]
        # for more information see [object]
        log_success("#############################################################")
        log_success("######## 'SangerAlignment' S4 instance is created !! ########")
        log_success("#############################################################")
        if (is.null(namesConversionCSV)) {
            log_success("  * >> ", contigNum, " contigs detected from 'regular expression'.")
        } else {
            log_success("  * >> ", contigNum, " contigs detected from 'csv file'.")
        }
        readNum <- 0
        for (contig in SangerContigList) {
            log_success("      * >> Contig '", contig@contigName, "':")
            log_success("          * >> ", length(contig@forwardReadList), " forward reads.")
            log_success("          * >> ", length(contig@reverseReadList), " reverse reads.")
            readNum <- readNum + length(contig@forwardReadList) + length(contig@reverseReadList)
        }
        log_success("  * >> ", readNum, " reads created from ", inputSource, " file.")
        if (TrimmingMethod == "M1") {
            log_success("  * >> Read is trimmed by 'M1 - Mott’s trimming algorithm'.")
        } else if (TrimmingMethod == "M2") {
            log_success("  * >> Read is trimmed by 'M2 - sliding window method'.")
        }
        # log_success("  * >> For more information, please run 'readTable(object)'.")
        
        # Add reads checking.

    }
    # else {
    #     log_error(paste(errors, collapse = ""))
    # }
    
    
    
    
    
    
    
    
    if (length(errors[[1]]) != 0) {
        creationResult <- FALSE
        sapply(paste0(errors[[2]], errors[[1]], '\n') , 
               log_error, simplify = FALSE)
        inputSource            = ""
        processMethod          = ""
        fastaFileName          = NULL
        namesConversionCSV     = NULL
        parentDirectory        = NULL
        contigName             = NULL
        suffixForwardRegExp    = NULL
        suffixReverseRegExp    = NULL
        forwardReadListFilter  = list()
        reverseReadListFilter  = list()
        trimmingMethodSC       = ""
        minReadsNum            = 0
        minReadLength          = 0
        refAminoAcidSeq        = ""
        minFractionCall        = 0
        maxFractionLost        = 0
        geneticCode            = ""
        acceptStopCodons       = FALSE
        readingFrame           = 0
    }
    if (nrow(readResultTable) != 0 && ncol(readResultTable) != 0) {
        names(readResultTable) <- readResultTableName
        log_debug("   >> For more information, please run 'object' or 'readTable(object)'.")
        log_debug("   >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads")
    }
    objectResults <- new("ObjectResults", creationResult = creationResult,
                         errorMessages = errors[[1]], errorTypes = errors[[2]],
                         warningMessages = character(0), warningTypes = character(0),
                         printLevel = printLevel, readResultTable = readResultTable)
    
    
    callNextMethod(.Object,
                   objectResults          = objectResults,
                   inputSource            = inputSource,
                   processMethod          = processMethod,
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
