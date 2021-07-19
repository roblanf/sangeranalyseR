#' @export
setOldClass("phylo")

setClassUnion("DNAStringORNULL", c("DNAString", "NULL"))
setClassUnion("DNAStringSetORNULL", c("DNAStringSet", "NULL"))

#' @title SangerAlignment
#'
#' @description  An S4 class containing SangerContigs lists and contigs alignment results which corresponds to a final alignment in Sanger sequencing.
#'
#' @slot objectResults This is the object that stores all information of the creation result.
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot processMethod The method to create a contig from reads. The value is \code{"REGEX"} or \code{"CSV"}. The default value is \code{"REGEX"}.
#' @slot ABIF_Directory If \code{inputSource} is \code{"ABIF"}, then this value is the path of a parent directory storing all reads in ABIF format you want to analyse. If \code{inputSource} is \code{"FASTA"}, then this value has to be \code{NULL} by default.
#' @slot FASTA_File If \code{inputSource} is \code{"FASTA"}, then this value has to be the path to a valid FASTA file ; if \code{inputSource} is \code{"ABIF"}, then this value has to be \code{NULL} by default.
#' @slot REGEX_SuffixForward The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be \code{"_F.ab1"}.
#' @slot REGEX_SuffixReverse The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be \code{"_R.ab1"}.
#' @slot CSV_NamesConversion The file path to the CSV file that provides read names, directions, and their contig groups. If \code{processMethod} is \code{"CSV"}, then this value has to be the path to a valid CSV file; if \code{processMethod} is \code{"REGEX"}, then this value has to be \code{NULL} by default.
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
#' ## Simple example
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica', 'ACHLO')
#' my_aligned_contigs <- new("SangerAlignment",
#'                           ABIF_Directory     = parentDir,
#'                           REGEX_SuffixForward = "_[0-9]*_F.ab1$",
#'                           REGEX_SuffixReverse = "_[0-9]*_R.ab1$")
#'                           
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica', 'ACHLO')
#' CSV_NamesConversion <- file.path(rawDataDir, "ab1", "SangerAlignment", "names_conversion.csv")
#' sangerAlignment <- new("SangerAlignment",
#'                        processMethod          = "CSV",
#'                        ABIF_Directory         = parentDir,
#'                        CSV_NamesConversion    = CSV_NamesConversion)
#' 
#' ## Input From ABIF file format (Regex)
#' REGEX_SuffixForward <- "_[0-9]*_F.ab1$"
#' REGEX_SuffixReverse <- "_[0-9]*_R.ab1$"
#' sangerAlignment <- new("SangerAlignment",
#'                        printLevel            = "SangerAlignment",
#'                        inputSource           = "ABIF",
#'                        processMethod         = "REGEX",
#'                        FASTA_File            = NULL,
#'                        CSV_NamesConversion   = NULL,
#'                        ABIF_Directory        = parentDir,
#'                        REGEX_SuffixForward   = REGEX_SuffixForward,
#'                        REGEX_SuffixReverse   = REGEX_SuffixReverse,
#'                        TrimmingMethod        = "M1",
#'                        M1TrimmingCutoff      = 0.0001,
#'                        M2CutoffQualityScore  = NULL,
#'                        M2SlidingWindowSize   = NULL,
#'                        baseNumPerRow         = 100,
#'                        heightPerRow          = 200,
#'                        signalRatioCutoff     = 0.33,
#'                        showTrimmed           = TRUE,
#'                        refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                        minReadsNum           = 2,
#'                        minReadLength         = 20,
#'                        minFractionCall       = 0.5,
#'                        maxFractionLost       = 0.5,
#'                        geneticCode           = GENETIC_CODE,
#'                        acceptStopCodons      = TRUE,
#'                        readingFrame          = 1,
#'                        processorsNum         = 2)
#'
#' ## Input From ABIF file format (Csv three column)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica', 'ACHLO')
#' CSV_NamesConversion <- file.path(rawDataDir, "ab1", "SangerAlignment", 
#' "names_conversion_all.csv")
#' sangerAlignment <- new("SangerAlignment",
#'                        inputSource           = "ABIF",
#'                        processMethod         = "CSV",
#'                        ABIF_Directory        = parentDir,
#'                        CSV_NamesConversion   = CSV_NamesConversion,
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
#' REGEX_SuffixForwardFa <- "_[0-9]*_F$"
#' REGEX_SuffixReverseFa <- "_[0-9]*_R$"
#' sangerAlignmentFa <- new("SangerAlignment",
#'                          inputSource           = "FASTA",
#'                          processMethod         = "REGEX",
#'                          FASTA_File            = fastaFN,
#'                          REGEX_SuffixForward   = REGEX_SuffixForwardFa,
#'                          REGEX_SuffixReverse   = REGEX_SuffixReverseFa,
#'                          refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                          processorsNum         = 2)
#'
#' ## Input From FASTA file format (Csv three column method)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerAlignment", "Sanger_all_reads.fa")
#' CSV_NamesConversion <- file.path(rawDataDir, "fasta",
#'                                 "SangerAlignment", "names_conversion.csv")
#' sangerAlignmentFa <- new("SangerAlignment",
#'                          inputSource           = "FASTA",
#'                          processMethod         = "CSV",
#'                          FASTA_File            = fastaFN,
#'                          CSV_NamesConversion   = CSV_NamesConversion,
#'                          refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                          processorsNum         = 2)
setClass("SangerAlignment",
         # Users need to name their ab1 files in a systematic way. Here is the
         # regulation:
         #  1. Naming: XXXXX_F_[0-9]*.ab1 / XXXXX_R_[0-9]*.ab1
         #        For reads in same contig, XXXXX must be same.
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerAlignment'
         ### -------------------------------------------------------------------
         representation(objectResults               = "ObjectResults",
                        inputSource                 = "character",
                        processMethod               = "character",
                        ABIF_Directory              = "characterORNULL",
                        FASTA_File                  = "characterORNULL",
                        REGEX_SuffixForward         = "characterORNULL",
                        REGEX_SuffixReverse         = "characterORNULL",
                        CSV_NamesConversion         = "characterORNULL",
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
                   ABIF_Directory         = NULL,
                   FASTA_File             = NULL,
                   REGEX_SuffixForward    = NULL,
                   REGEX_SuffixReverse    = NULL,
                   CSV_NamesConversion    = NULL,
                   geneticCode            = GENETIC_CODE,
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
                   acceptStopCodons       = TRUE,
                   readingFrame           = 1,
                   processorsNum          = 1) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    creationResult <- TRUE
    errors <- list(character(0), character(0))
    warnings <- list(character(0), character(0))
    readResultTableName <- c("readName","creationResult", "errorType", 
                             "errorMessage", "inputSource", "direction")
    readResultTable <- data.frame()
    ############################################################################
    ### First layer of pre-checking: SangerAlignment input parameter prechecking
    ############################################################################
    if (printLevel == "SangerAlignment") {
        errors <- checkInputSource(inputSource, errors[[1]], errors[[2]])
        errors <- checkProcessMethod(inputSource, processMethod, errors[[1]], errors[[2]])
        errors <- checkGeneticCode(geneticCode, errors[[1]], errors[[2]])

        errors <- checkRefAAS(refAminoAcidSeq, errors[[1]], errors[[2]])
        errors <- checkMinReadsNum(minReadsNum, errors[[1]], errors[[2]])
        errors <- checkMinReadLength(minReadLength, errors[[1]], errors[[2]])
        errors <- checkMinFractionCall(minFractionCall, errors[[1]], errors[[2]])
        errors <- checkMaxFractionLost(maxFractionLost, errors[[1]], errors[[2]])
        errors <- checkAcceptStopCodons(acceptStopCodons, errors[[1]], errors[[2]])
        errors <- checkReadingFrame(readingFrame, errors[[1]], errors[[2]])
        errors <- checkProcessorsNum(processorsNum, errors[[1]], errors[[2]])
        if (inputSource == "ABIF") {
            ### ----------------------------------------------------------------
            ### 'ABIF' condition checking!
            ### ----------------------------------------------------------------
            # Check ABIF_Directory and set FASTA_File to NULL
            errors <- checkABIF_Directory(ABIF_Directory, 
                                          errors[[1]], errors[[2]])
            FASTA_File <- NULL
            # Check the trimming method
            errors <- checkTrimParam(TrimmingMethod,
                                     M1TrimmingCutoff,
                                     M2CutoffQualityScore,
                                     M2SlidingWindowSize,
                                     errors[[1]], errors[[2]])
            if (TrimmingMethod == "M1") {
                M2CutoffQualityScore = NULL
                M2SlidingWindowSize = NULL
            } else if (TrimmingMethod == "M2") {
                M1TrimmingCutoff = NULL
            }
            # Check chromatogram parameter
            errors <- checkBaseNumPerRow (baseNumPerRow, errors[[1]], errors[[2]])
            errors <- checkHeightPerRow (heightPerRow, errors[[1]], errors[[2]])
            errors <- checkSignalRatioCutoff (signalRatioCutoff, 
                                              errors[[1]], errors[[2]])
            errors <- checkShowTrimmed (showTrimmed, errors[[1]], errors[[2]])
        } else if (inputSource == "FASTA") {
            ### ----------------------------------------------------------------
            ### 'FASTA' condition checking!
            ### ----------------------------------------------------------------
            # Check FASTA_File and set ABIF_Directory to NULL
            errors <- checkFASTA_File(inputSource, FASTA_File,
                                      errors[[1]], errors[[2]])
            ABIF_Directory <- NULL
            # Set trimming parameters to NULL
            TrimmingMethod = ""
            M1TrimmingCutoff = NULL
            M2CutoffQualityScore = NULL
            M2SlidingWindowSize = NULL
            if (length(errors[[1]]) == 0 ) {
                # Read FASTA file
                readFasta <- read.fasta(FASTA_File, as.string = TRUE)
                fastaNames <- names(readFasta)
            }
        }
        if (processMethod=="REGEX") {
            ### ----------------------------------------------------------------
            ### 'REGEX' condition checking!
            ### ----------------------------------------------------------------
            # Check REGEX_SuffixForward and REGEX_SuffixReverse 
            #  and set CSV_NamesConversion to NULL
            errors <- checkREGEX_SuffixForward(REGEX_SuffixForward, 
                                               errors[[1]], errors[[2]])
            errors <- checkREGEX_SuffixReverse(REGEX_SuffixReverse, 
                                               errors[[1]], errors[[2]])
            CSV_NamesConversion <- NULL
        } else if (processMethod=="CSV") {
            ### ----------------------------------------------------------------
            ### 'CSV' condition checking!
            ### ----------------------------------------------------------------
            # Check CSV_NamesConversion and set 
            #  REGEX_SuffixForward and REGEX_SuffixReverse to NULL
            errors <- checkCSV_NamesConversion(CSV_NamesConversion, 
                                               errors[[1]], errors[[2]])
            REGEX_SuffixForward <- NULL
            REGEX_SuffixReverse <- NULL
            errors <- checkAb1FastaCsv(ABIF_Directory, FASTA_File, 
                                       CSV_NamesConversion, inputSource, 
                                       errors[[1]], errors[[2]])
        }
    }
    if (length(errors[[1]]) == 0 ) {
        log_info('#################################################')
        log_info('#### Start creating SangerAlignment instance ####')
        log_info('#################################################')
        processorsNum <- getProcessors (processorsNum)
        if (inputSource == "ABIF" && processMethod == "REGEX") {
            log_info("  >> You are using Regular Expression Method",
                     " to group AB1 files!")
            parentDirFiles <- list.files(ABIF_Directory, recursive = TRUE)
            forwardSelectInputFiles <-
                parentDirFiles[grepl(REGEX_SuffixForward, parentDirFiles)]
            reverseSelectInputFiles <-
                parentDirFiles[grepl(REGEX_SuffixReverse, parentDirFiles)]
            # cat("*** forwardSelectInputFiles: ", forwardSelectInputFiles)
            warnings <- checkGreplForward(forwardSelectInputFiles, 
                                          warnings[[1]], warnings[[2]])
            warnings <- checkGreplReverse(reverseSelectInputFiles, 
                                          warnings[[1]], warnings[[2]])
            # Find possible consensus Name for forward and reverse reads
            forwardContigName <-
                unlist(str_split(forwardSelectInputFiles,
                                 REGEX_SuffixForward, n = Inf,
                                 simplify = FALSE))[c(TRUE, FALSE)]
            reverseContigName <-
                unlist(str_split(reverseSelectInputFiles,
                                 REGEX_SuffixReverse, n = Inf,
                                 simplify = FALSE))[c(TRUE, FALSE)]
            contigNames <- union(forwardContigName, reverseContigName)
            contigNumber <- length(contigNames)
            # cat("*** contigNames: ", contigNames)
            # Create contig for all list of contigNumber
            ### ----------------------------------------------------------------
            ##### Creating each SangerContig (store as SangerContigList)
            ### ----------------------------------------------------------------
            SangerContigList <-
                lapply(contigNames,
                       function(eachConsRead) {
                           insideDirName<- dirname(eachConsRead)
                           insideContigName <- basename(eachConsRead)
                           # cat("**** eachConsRead: ", eachConsRead)
                           # cat("**** insideDirName: ", insideDirName)
                           # cat("**** insideContigName: ", insideContigName)
                           newSangerContig <-
                               new("SangerContig",
                                   printLevel           = printLevel,
                                   inputSource          = inputSource,
                                   processMethod        = processMethod,
                                   ABIF_Directory       =
                                       file.path(ABIF_Directory, insideDirName),
                                   FASTA_File           = FASTA_File,
                                   REGEX_SuffixForward  = REGEX_SuffixForward,
                                   REGEX_SuffixReverse  = REGEX_SuffixReverse,
                                   CSV_NamesConversion  = CSV_NamesConversion,
                                   contigName           = insideContigName,
                                   geneticCode          = geneticCode,
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
                                   acceptStopCodons     = acceptStopCodons,
                                   readingFrame         = readingFrame,
                                   processorsNum        = processorsNum)
                           readResultTable <<- rbind(readResultTable, newSangerContig@objectResults@readResultTable)
                           if (newSangerContig@objectResults@creationResult) {
                               newSangerContig
                           } else {
                               NULL
                           }
                       })
            names(SangerContigList) <- contigNames
        } else if (inputSource == "ABIF" && processMethod == "CSV") {
            log_info("**** You are using CSV Name Conversion Method ",
                     "to group AB1 files!")
            parentDirFiles <- list.files(ABIF_Directory, recursive = TRUE)
            csvFile <- read.csv(CSV_NamesConversion, header = TRUE)
            log_info("**** Contig number in your Csv file is ", 
                     length(unique(csvFile$contig)))
            contigNames <- as.character(unique(csvFile$contig))
            contigNames <- lapply(contigNames, function(contigName) {
                contigNameSelectInputFiles <-
                    parentDirFiles[grepl(contigName, parentDirFiles)]   
                inside_contigNames <- file.path(dirname(contigNameSelectInputFiles), contigName)
                inside_contigNames
            })
            contigNames <- unique(unlist(contigNames, recursive = TRUE))
            SangerContigList <- lapply(contigNames, function(contigName) {
                insideDirName<- dirname(contigName)
                insideContigName <- basename(contigName)
                newSangerContig <-
                    new("SangerContig",
                        printLevel           = printLevel,
                        inputSource          = inputSource,
                        processMethod        = processMethod,
                        ABIF_Directory       =
                            file.path(ABIF_Directory, insideDirName),
                        FASTA_File           = FASTA_File,
                        REGEX_SuffixForward  = REGEX_SuffixForward,
                        REGEX_SuffixReverse  = REGEX_SuffixReverse,
                        CSV_NamesConversion  = CSV_NamesConversion,
                        contigName           = insideContigName,
                        geneticCode          = geneticCode,
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
                        acceptStopCodons     = acceptStopCodons,
                        readingFrame         = readingFrame,
                        processorsNum        = processorsNum)
                readResultTable <<- rbind(readResultTable, newSangerContig@objectResults@readResultTable)
                if (newSangerContig@objectResults@creationResult) {
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
            fastaNames <- names(readFasta)
            forwardSelectInputFiles <- fastaNames[grepl(REGEX_SuffixForward,
                                                        fastaNames)]
            reverseSelectInputFiles <- fastaNames[grepl(REGEX_SuffixReverse,
                                                        fastaNames)]
            warnings <- checkGreplForward(forwardSelectInputFiles, 
                                          warnings[[1]], warnings[[2]])
            warnings <- checkGreplReverse(reverseSelectInputFiles, 
                                          warnings[[1]], warnings[[2]])
            # Find possible consensus Name for forward and reverse reads
            forwardContigName <-
                unlist(str_split(forwardSelectInputFiles, REGEX_SuffixForward,
                                 n = Inf, simplify = FALSE))[c(TRUE, FALSE)]
            reverseContigName <-
                unlist(str_split(reverseSelectInputFiles, REGEX_SuffixReverse,
                                 n = Inf, simplify = FALSE))[c(TRUE, FALSE)]
            contigNames <- union(forwardContigName, reverseContigName)
            SangerContigList <- lapply(contigNames, function(contigName) {
                newSangerContig <-
                    new("SangerContig",
                        printLevel           = printLevel,
                        inputSource          = inputSource,
                        processMethod        = processMethod,
                        ABIF_Directory       = ABIF_Directory,
                        FASTA_File           = FASTA_File,
                        REGEX_SuffixForward  = REGEX_SuffixForward,
                        REGEX_SuffixReverse  = REGEX_SuffixReverse,
                        CSV_NamesConversion  = CSV_NamesConversion,
                        contigName           = contigName,
                        geneticCode          = geneticCode,
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
                        acceptStopCodons     = acceptStopCodons,
                        readingFrame         = readingFrame,
                        processorsNum        = processorsNum)
                readResultTable <<- rbind(readResultTable, newSangerContig@objectResults@readResultTable)
                if (newSangerContig@objectResults@creationResult) {
                    newSangerContig
                } else {
                    NULL
                }
            })
            names(SangerContigList) <- contigNames
        } else if (inputSource == "FASTA" && processMethod == "CSV") {
            log_info("**** You are using CSV Name Conversion Method ",
                     "to group reads in FASTA file (with Csv file)!")
            csvFile <- read.csv(CSV_NamesConversion, header = TRUE)
            contigNames <- unique(as.character(csvFile$contig))
            SangerContigList <- lapply(contigNames, function(contigName) {
                newSangerContig <-
                    new("SangerContig",
                        printLevel           = printLevel,
                        inputSource          = inputSource,
                        processMethod        = processMethod,
                        ABIF_Directory       = ABIF_Directory,
                        FASTA_File           = FASTA_File,
                        REGEX_SuffixForward  = REGEX_SuffixForward,
                        REGEX_SuffixReverse  = REGEX_SuffixReverse,
                        CSV_NamesConversion  = CSV_NamesConversion,
                        contigName           = contigName,
                        geneticCode          = geneticCode,
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
                        acceptStopCodons     = acceptStopCodons,
                        readingFrame         = readingFrame,
                        processorsNum        = processorsNum)
                readResultTable <<- rbind(readResultTable, newSangerContig@objectResults@readResultTable)
                if (newSangerContig@objectResults@creationResult) {
                    newSangerContig
                } else {
                    NULL
                }
            })
            names(SangerContigList) <- contigNames
        }
    }
    
    if (length(errors[[1]]) == 0 ) {
        SangerContigList <- Filter(Negate(is.null), SangerContigList)
        acResult <- alignContigs(SangerContigList, geneticCode,
                                 refAminoAcidSeq, minFractionCall,
                                 maxFractionLost, processorsNum)
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
        
        if (contigNum > 0) {
            if (processMethod == "REGEX") {
                log_info("  >> ", contigNum, " contigs detected from 'regular expression'.")
            } else if (processMethod == "CSV") {
                log_info("  >> ", contigNum, " contigs detected from 'csv file'.")
            }
            readNum <- 0
            for (contig in SangerContigList) {
                log_info("      >> Contig '", contig@contigName, "':")
                log_info("          >> ", length(contig@forwardReadList), " forward reads.")
                log_info("          >> ", length(contig@reverseReadList), " reverse reads.")
                readNum <- readNum + length(contig@forwardReadList) + length(contig@reverseReadList)
            }
            log_info("  >> ", readNum, " reads created from ", inputSource, " file.")
            if (TrimmingMethod == "M1") {
                log_info("  >> Reads are trimmed by 'M1 - Mott’s trimming algorithm'.")
            } else if (TrimmingMethod == "M2") {
                log_info("  >> Reads are trimmed by 'M2 - sliding window method'.")
            }   
        } else if (contigNum == 0) {
            ### ----------------------------------------------------------------
            ### Failed: SangerAlignment Level
            ### ----------------------------------------------------------------
            msg <- paste0("The number of your total contig is 0.",
                          "\nPlease check your name matching parameters.")
            errors[[1]] <- c(errors[[1]], msg)
            errors[[2]] <- c(errors[[2]], "CONTIG_NUMBER_ZERO_ERROR")
        }
    }
    
    if (length(errors[[1]]) != 0) {
        creationResult <- FALSE
        sapply(paste0(errors[[2]], '\n',errors[[1]], '\n') , 
               log_error, simplify = FALSE)
        inputSource            = ""
        processMethod          = ""
        FASTA_File         = NULL
        CSV_NamesConversion    = NULL
        ABIF_Directory       = NULL
        REGEX_SuffixForward   = NULL
        REGEX_SuffixReverse   = NULL
        TrimmingMethod        = ""
        SangerContigList      = list()
        minFractionCall       = 0
        maxFractionLost       = 0
        geneticCode           = ""
        consensus             = DNAString()
        refAminoAcidSeq       = ""
        aln                   = DNAStringSet()
        aln.tree              = read.tree(text="();")
    }
    if (printLevel == "SangerAlignment") {
        if (nrow(readResultTable) != 0 && ncol(readResultTable) != 0) {
            names(readResultTable) <- readResultTableName
        }
        log_debug("   >> For more information, please run 'object'.")
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
                   FASTA_File         = FASTA_File,
                   CSV_NamesConversion    = CSV_NamesConversion,
                   ABIF_Directory       = ABIF_Directory,
                   REGEX_SuffixForward   = REGEX_SuffixForward,
                   REGEX_SuffixReverse   = REGEX_SuffixReverse,
                   contigList            = SangerContigList,
                   geneticCode           = geneticCode,
                   contigsConsensus      = consensus,
                   refAminoAcidSeq       = refAminoAcidSeq,
                   contigsAlignment      = aln,
                   contigsTree           = aln.tree)
})

