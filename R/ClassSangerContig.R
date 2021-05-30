#' @title SangerContig
#'
#' @description  An S4 class containing forward and reverse SangerRead lists and alignment, consensus read results which corresponds to a contig in Sanger sequencing.
#' @slot objectResults This is the object that stores all information of the creation result.
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot processMethod The method to create a contig from reads. The value is \code{"REGEX"} or \code{"CSV"}. The default value is \code{"REGEX"}.
#' @slot fastaFileName If \code{inputSource} is \code{"FASTA"}, then this value has to be the path to a valid FASTA file ; if \code{inputSource} is \code{"ABIF"}, then this value has to be \code{NULL} by default.
#' @slot parentDirectory If \code{inputSource} is \code{"ABIF"}, then this value is the path of a parent directory storing all reads in ABIF format you want to analyse. If \code{inputSource} is \code{"FASTA"}, then this value has to be \code{NULL} by default.
#' @slot namesConversionCSV The file path to the CSV file that provides read names, directions, and their contig groups. If \code{processMethod} is \code{"CSV"}, then this value has to be the path to a valid CSV file; if \code{processMethod} is \code{"REGEX"}, then this value has to be \code{NULL} by default.
#' @slot contigName The contig name of all the reads in \code{parentDirectory}.
#' @slot suffixForwardRegExp The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented.
#' @slot suffixReverseRegExp The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented.
#' @slot geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @slot forwardReadList The list of SangerRead S4 instances which are all forward reads.
#' @slot reverseReadList The list of SangerRead S4 instances which are all reverse reads.
#' @slot trimmingMethodSC The read trimming method for all SangerRead S4 instances in SangerContig. The value must be \code{"M1"} (the default) or \code{'M2'}. All SangerRead must have the same trimming method.
#' @slot minReadsNum The minimum number of reads required to make a consensus sequence, must be 2 or more. The default value is \code{2}.
#' @slot minReadLength Reads shorter than this will not be included in the readset. The default \code{20} means that all reads with length of 20 or more will be included. Note that this is the length of a read after it has been trimmed.
#' @slot refAminoAcidSeq An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is \code{""}.
#' @slot minFractionCall Minimum fraction of the sequences required to call a consensus sequence for SangerContig at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
#' @slot maxFractionLost Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerContig (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
#' @slot acceptStopCodons The logical value \code{TRUE} or \code{FALSE}. \code{TRUE} (the defualt): keep all reads, regardless of whether they have stop codons; \code{FALSE}: reject reads with stop codons. If \code{FALSE} is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
#' @slot readingFrame \code{1}, \code{2}, or \code{3}. Only used if \code{accept.stop.codons == FALSE}. This specifies the reading frame that is used to determine stop codons. If you use a \code{refAminoAcidSeq}, then the frame should always be \code{1}, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame.
#' @slot contigSeq The consensus read of all SangerRead S4 instances in DNAString object.
#' @slot alignment The alignment of all SangerRead S4 instances with the called consensus sequence in DNAStringSet object. Users can use BrowseSeqs() to view the alignment.
#' @slot differencesDF A data frame of the number of pairwise differences between each read and the consensus sequence, as well as the number of bases in each input read that did not contribute to the consensus sequence. It can assist in detecting incorrect reads, or reads with a lot of errors.
#' @slot distanceMatrix A distance matrix of genetic distances (corrected with the JC model) between all of the input reads.
#' @slot dendrogram A list storing cluster groups in a data frame and a dendrogram object depicting the distance.matrix. Users can use plot() to see the dendrogram.
#' @slot indelsDF If users specified a reference sequence via \code{refAminoAcidSeq}, then this will be a data frame describing the number of indels and deletions that were made to each of the input reads in order to correct frameshift mutations.
#' @slot stopCodonsDF If users specified a reference sequence via \code{refAminoAcidSeq}, then this will be a data frame describing the number of stop codons in each read.
#' @slot secondaryPeakDF A data frame with one row for each column in the alignment that contained more than one secondary peak. The data frame has three columns: the column number of the alignment; the number of secondary peaks in that column; and the bases (with IUPAC ambiguity codes representing secondary peak calls) in that column represented as a string.
#'
#' @name SangerContig-class
#' @exportClass SangerContig
#'
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerRead.R
#' @examples
#' ## Simple example
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")
#' contigName <- "Achl_RBNII384-13"
#' suffixForwardRegExp <- "_[0-9]*_F.ab1"
#' suffixReverseRegExp <- "_[0-9]*_R.ab1"
#' sangerContig <- new("SangerContig",
#'                      parentDirectory       = parentDir,
#'                      contigName            = contigName,
#'                      suffixForwardRegExp   = suffixForwardRegExp,
#'                      suffixReverseRegExp   = suffixReverseRegExp)
#'                      
#' ## forward / reverse reads match error
#' ## Input From ABIF file format (Regex)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")
#' contigName <- "Achl_RBNII384-13"
#' suffixForwardRegExp <- "_[0-9]*_F.ab1"
#' suffixReverseRegExp <- "_[0-9]*_R.ab1"
#' sangerContig <- new("SangerContig",
#'                      inputSource           = "ABIF",
#'                      processMethod         = "REGEX",
#'                      parentDirectory       = parentDir,
#'                      contigName            = contigName,
#'                      suffixForwardRegExp   = suffixForwardRegExp,
#'                      suffixReverseRegExp   = suffixReverseRegExp,
#'                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                      TrimmingMethod        = "M1",
#'                      M1TrimmingCutoff      = 0.0001,
#'                      baseNumPerRow         = 100,
#'                      heightPerRow          = 200,
#'                      signalRatioCutoff     = 0.33,
#'                      showTrimmed           = TRUE,
#'                      minReadsNum           = 1,
#'                      processorsNum         = 2)
#'
#' ## Input From ABIF file format (Csv three column method)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")
#' namesConversionCSV <- file.path(rawDataDir, "ab1", "SangerContig", "names_conversion_2.csv")
#' sangerContig <- new("SangerContig",
#'                      inputSource           = "ABIF",
#'                      processMethod         = "CSV",
#'                      parentDirectory       = parentDir,
#'                      namesConversionCSV    = namesConversionCSV,
#'                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                      TrimmingMethod        = "M1",
#'                      M1TrimmingCutoff      = 0.000001,
#'                      baseNumPerRow         = 100,
#'                      heightPerRow          = 200,
#'                      signalRatioCutoff     = 0.33,
#'                      showTrimmed           = TRUE,
#'                      processorsNum         = 2)
#'
#'
#' ## Input From FASTA file format (Regex)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerContig", "Achl_ACHLO006-09.fa")
#' contigName <- "Achl_ACHLO006-09"
#' suffixForwardRegExpFa <- "_[0-9]*_F$"
#' suffixReverseRegExpFa <- "_[0-9]*_R$"
#' sangerContigFa <- new("SangerContig",
#'                       inputSource           = "FASTA",
#'                       processMethod         = "REGEX",
#'                       fastaFileName         = fastaFN,
#'                       contigName            = contigName,
#'                       suffixForwardRegExp   = suffixForwardRegExpFa,
#'                       suffixReverseRegExp   = suffixReverseRegExpFa,
#'                       refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                       processorsNum         = 2)
#'
#' ## Input From FASTA file format (Csv - Csv three column method)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerContig", "Achl_ACHLO006-09.fa")
#' namesConversionCSV <- file.path(rawDataDir, "fasta", "SangerContig", "names_conversion_1.csv")
#' sangerContigFa <- new("SangerContig",
#'                       inputSource           = "FASTA",
#'                       processMethod         = "CSV",
#'                       fastaFileName         = fastaFN,
#'                       namesConversionCSV    = namesConversionCSV,
#'                       refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                       processorsNum         = 2)
setClass("SangerContig",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerContig'
         ### -------------------------------------------------------------------
         representation(objectResults             = "ObjectResults",
                        inputSource               = "character",
                        processMethod             = "character",
                        fastaFileName             = "characterORNULL",
                        namesConversionCSV        = "characterORNULL",
                        parentDirectory           = "characterORNULL",
                        contigName                = "characterORNULL",
                        suffixForwardRegExp       = "characterORNULL",
                        suffixReverseRegExp       = "characterORNULL",
                        geneticCode               = "character",
                        forwardReadList           = "list",
                        reverseReadList           = "list",
                        trimmingMethodSC          = "character",
                        minReadsNum               = "numeric",
                        minReadLength             = "numeric",
                        refAminoAcidSeq           = "character",
                        minFractionCall           = "numeric",
                        maxFractionLost           = "numeric",
                        acceptStopCodons          = "logical",
                        readingFrame              = "numeric",
                        contigSeq                 = "DNAString",
                        alignment                 = "DNAStringSet",
                        differencesDF             = "data.frame",
                        distanceMatrix            = "matrix",
                        dendrogram                = "list",
                        indelsDF                  = "data.frame",
                        stopCodonsDF              = "data.frame",
                        secondaryPeakDF           = "data.frame"),
)

### ============================================================================
### Overwrite initialize for 'SangerContig' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerContig",
          function(.Object,
                   printLevel             = "SangerContig",
                   inputSource            = "ABIF",
                   processMethod          = "REGEX",
                   fastaFileName          = NULL,
                   parentDirectory        = NULL,
                   namesConversionCSV     = NULL,
                   contigName             = NULL,
                   suffixForwardRegExp    = NULL,
                   suffixReverseRegExp    = NULL,
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
                   processorsNum          = 1,
                   logLevel               = TRUE) {
    creationResult <- TRUE
    errors <- list(character(0), character(0))
    warnings <- c(character(0))
    readResultTableName <- c("readName","creationResult", "errorType", 
                             "errorMessage", "inputSource", "direction")
    readResultTable <- data.frame()
    ############################################################################
    ### First layer of pre-checking: SangerContig input parameter prechecking
    ############################################################################
    if (printLevel == "SangerContig") {
        errors <- checkInputSource(inputSource, errors[[1]], errors[[2]])
        errors <- checkProcessMethod(inputSource, processMethod, 
                                     errors[[1]], errors[[2]])
        errors <- checkFastaFileName(inputSource, fastaFileName,
                                     errors[[1]], errors[[2]])
        errors <- checkNamesConversionCSV(processMethod, namesConversionCSV, 
                                          contigName, suffixForwardRegExp,
                                          suffixReverseRegExp, inputSource, 
                                          errors[[1]], errors[[2]])
        errors <- checkRefAAS(refAminoAcidSeq, errors[[1]], errors[[2]])
        errors <- checkMinReadsNum(minReadsNum, errors[[1]], errors[[2]])
        errors <- checkMinReadLength(minReadLength, errors[[1]], errors[[2]])
        errors <- checkMinFractionCall(minFractionCall, errors[[1]], errors[[2]])
        errors <- checkMaxFractionLost(maxFractionLost, errors[[1]], errors[[2]])
        errors <- checkGeneticCode(geneticCode, errors[[1]], errors[[2]])
        errors <- checkAcceptStopCodons(acceptStopCodons, errors[[1]], errors[[2]])
        errors <- checkReadingFrame(readingFrame, errors[[1]], errors[[2]])
        errors <- checkProcessorsNum(processorsNum, errors[[1]], errors[[2]])
        if (inputSource == "ABIF") {
            ### ----------------------------------------------------------------
            ### 'ABIF' condition checking!
            ### ----------------------------------------------------------------
            ####################################################################
            ### Second layer of pre-checking: 'ABIF' condition checking!
            ####################################################################
            errors <- checkParentDirectory (parentDirectory, 
                                            errors[[1]], errors[[2]])
            errors <- checkTrimParam(TrimmingMethod,
                                     M1TrimmingCutoff,
                                     M2CutoffQualityScore,
                                     M2SlidingWindowSize,
                                     errors[[1]], errors[[2]])
            errors <- checkBaseNumPerRow (baseNumPerRow, errors[[1]], errors[[2]])
            errors <- checkHeightPerRow (heightPerRow, errors[[1]], errors[[2]])
            errors <- checkSignalRatioCutoff (signalRatioCutoff, 
                                              errors[[1]], errors[[2]])
            errors <- checkShowTrimmed (showTrimmed, errors[[1]], errors[[2]])
        } else if (inputSource == "FASTA") {
            ### ----------------------------------------------------------------
            ### 'FASTA' condition checking!
            ### ----------------------------------------------------------------
            ####################################################################
            ### Second layer of pre-checking: 'FASTA' condition checking!
            ####################################################################
        }
        if (length(errors[[1]]) == 0 && processMethod=="CSV") {
            errors <- checkAb1FastaCsv(parentDirectory, fastaFileName, 
                                       namesConversionCSV, inputSource, 
                                       errors[[1]], errors[[2]])
        }   
    }
    if (length(errors[[1]]) == 0 ) {
        log_info("========================================================")
        log_info("================ Creating 'SangerContig' ===============")
        log_info("========================================================")
        log_info("  >> Contig Name: '", contigName, "'")
        processorsNum <- getProcessors (processorsNum)
        if (inputSource == "ABIF") {
        } else if (inputSource == "FASTA") {
            readFasta <- read.fasta(fastaFileName, as.string = TRUE)
            fastaNames <- names(readFasta)
            TrimmingMethod <- ""
            M1TrimmingCutoff <- NULL
            M2CutoffQualityScore <- NULL
            M2SlidingWindowSize <- NULL
        }
        if (inputSource == "ABIF" && processMethod == "REGEX") {
            if (printLevel == "SangerContig") {
                log_info("  >> You are using Regular Expression Method",
                         " to group AB1 files!")
            }
            parentDirFiles <- list.files(parentDirectory)
            forwardSelectInputFiles <-
                parentDirFiles[grepl(paste0(contigName, suffixForwardRegExp),
                                     parentDirFiles)]
            reverseSelectInputFiles <-
                parentDirFiles[grepl(paste0(contigName, suffixReverseRegExp),
                                     parentDirFiles)]
            forwardAllReads <- lapply(parentDirectory, file.path,
                                      forwardSelectInputFiles)
            reverseAllReads <- lapply(parentDirectory, file.path,
                                      reverseSelectInputFiles)
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(forwardAllReads[[1]], function(forwardN){
                newSangerRead<- new("SangerRead",
                                    printLevel           = printLevel,
                                    inputSource          = inputSource,
                                    readFeature          = "Forward Read",
                                    readFileName         = forwardN,
                                    geneticCode          = geneticCode,
                                    TrimmingMethod       = TrimmingMethod,
                                    M1TrimmingCutoff     = M1TrimmingCutoff,
                                    M2CutoffQualityScore = M2CutoffQualityScore,
                                    M2SlidingWindowSize  = M2SlidingWindowSize,
                                    baseNumPerRow        = baseNumPerRow,
                                    heightPerRow         = heightPerRow,
                                    signalRatioCutoff    = signalRatioCutoff,
                                    showTrimmed          = showTrimmed)
            })
            names(forwardReadList) <- forwardAllReads[[1]]
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseReadList <- lapply(reverseAllReads[[1]], function(reverseN){
                newSangerRead <- new("SangerRead",
                                     printLevel           = printLevel,
                                     inputSource          = inputSource,
                                     readFeature          = "Reverse Read",
                                     readFileName         = reverseN,
                                     geneticCode          = geneticCode,
                                     TrimmingMethod       = TrimmingMethod,
                                     M1TrimmingCutoff     = M1TrimmingCutoff,
                                     M2CutoffQualityScore = M2CutoffQualityScore,
                                     M2SlidingWindowSize  = M2SlidingWindowSize,
                                     baseNumPerRow        = baseNumPerRow,
                                     heightPerRow         = heightPerRow,
                                     signalRatioCutoff    = signalRatioCutoff,
                                     showTrimmed          = showTrimmed)
            })
            names(reverseReadList) <- reverseAllReads[[1]]
        } else if (inputSource == "ABIF" && processMethod == "CSV") {
            if (printLevel == "SangerContig") {
                log_info("  >> You are using CSV Name Conversion Method ",
                         "to group AB1 files!")
            }
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            if (is.null(contigName)) {
                if (length(unique(csvFile$contig)) != 1) {
                    log_error("Error! There are ",
                              length(unique(csvFile$contig)) ,
                              " contigName in the CSV file. ",
                              "There should be only one contigName in the CSV file.")
                } else {
                    if (printLevel == "SangerContig") {
                        log_info(">> Your contig name is ", contigName)
                    }
                    contigName <- unique(csvFile$contig)
                    contigNameCSV <- as.character(unique(csvFile$contig))
                    selectedCsvFile <- csvFile[csvFile$contig == contigNameCSV, ]
                    forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F", ]
                    reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R", ]
                }
            } else if (!is.null(contigName)) {
                if (printLevel == "SangerContig") {
                    log_info(">> Your contig name is ", contigName)
                }
                selectedCsvFile <- csvFile[csvFile$contig == contigName, ]
                forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F",]
                reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R",]
            }
            
            forwardOriginal <- as.character(forwardCsv$reads)
            fAbsoluteAB1 <- file.path(parentDirectory, forwardOriginal)
            if (!all(file.exists(fAbsoluteAB1))) {
                log_error("One of your 'reads' in the csv file ",
                          "does not match any ab1 files in '", parentDirectory, "'")
            }
            # lapply to create forward SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(fAbsoluteAB1, function(forwardN){
                newSangerRead <- new("SangerRead",
                                     printLevel           = printLevel,
                                     inputSource          = inputSource,
                                     readFeature          = "Forward Read",
                                     readFileName         = forwardN,
                                     geneticCode          = geneticCode,
                                     TrimmingMethod       = TrimmingMethod,
                                     M1TrimmingCutoff     = M1TrimmingCutoff,
                                     M2CutoffQualityScore = M2CutoffQualityScore,
                                     M2SlidingWindowSize  = M2SlidingWindowSize,
                                     baseNumPerRow        = baseNumPerRow,
                                     heightPerRow         = heightPerRow,
                                     signalRatioCutoff    = signalRatioCutoff,
                                     showTrimmed          = showTrimmed)
            })
            names(forwardReadList) <- fAbsoluteAB1
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseOriginal <- as.character(reverseCsv$reads)
            rAbsoluteAB1 <- file.path(parentDirectory, reverseOriginal)
            if (!all(file.exists(rAbsoluteAB1))) {
                log_error("One of your 'reads' in the csv file ",
                          "does not match any ab1 files in '", parentDirectory, "'")
            }
            # lapply to create reverse SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            reverseReadList <- lapply(rAbsoluteAB1, function(reverseN){
                newSangerRead <-
                    new("SangerRead",
                        printLevel           = printLevel,
                        inputSource          = inputSource,
                        readFeature          = "Reverse Read",
                        readFileName         = reverseN,
                        geneticCode          = geneticCode,
                        TrimmingMethod       = TrimmingMethod,
                        M1TrimmingCutoff     = M1TrimmingCutoff,
                        M2CutoffQualityScore = M2CutoffQualityScore,
                        M2SlidingWindowSize  = M2SlidingWindowSize,
                        baseNumPerRow        = baseNumPerRow,
                        heightPerRow         = heightPerRow,
                        signalRatioCutoff    = signalRatioCutoff,
                        showTrimmed          = showTrimmed)
            })
            names(reverseReadList) <- rAbsoluteAB1
        }
        if (inputSource == "ABIF") {
            readNum <- length(forwardReadList) + length(reverseReadList)
            log_info("   >> The number of reads detected: ", readNum)
            forwardReadListFilter <- lapply(forwardReadList, function(read) {
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
                        read
                    } else {
                        msg <- paste0("  >> ", read@readFileName, 
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
                        NULL
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
                    NULL
                }
            })
            reverseReadListFilter <- lapply(reverseReadList, function(read) {
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
                        read
                    } else {
                        msg <- paste0("  >> ", read@readFileName, 
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
                        NULL
                    }
                } else {
                    ### --------------------------------------------------------
                    ### Failed: readResultTable (SangerRead Level)
                    ### --------------------------------------------------------
                    ### Failed: readResultTable (SangerRead Level)
                    row <- data.frame(basename(read@readFileName), 
                                      read@objectResults@creationResult, 
                                      read@objectResults@errorTypes, 
                                      read@objectResults@errorMessages, 
                                      read@inputSource, read@readFeature)
                    names(row) <- readResultTableName
                    readResultTable <<- rbind(readResultTable, row)
                    NULL
                }                 
            })            
        }
        if (inputSource == "FASTA" && processMethod == "REGEX") {
            if (printLevel == "SangerContig") {
                log_info("  >> You are using Regular Expression Method ",
                         "to group reads in FASTA file (No CSV file)!")
            }
            ### ----------------------------------------------------------------
            ### Find names with the given contigName
            ### ----------------------------------------------------------------
            contigSubGroupNames <-
                fastaNames[grepl(contigName, fastaNames, fixed=TRUE)]
            ### ----------------------------------------------------------------
            ### Among them, find the forward names
            ### ----------------------------------------------------------------
            forwardSelectNames <-
                contigSubGroupNames[grepl(suffixForwardRegExp,
                                          contigSubGroupNames)]
            ### ----------------------------------------------------------------
            ### Among them, find the reverse names
            ### ----------------------------------------------------------------
            reverseSelectNames <-
                contigSubGroupNames[grepl(suffixReverseRegExp,
                                          contigSubGroupNames)]
            # lapply to create SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(forwardSelectNames, function(forwardName){
                newSangerRead <- new("SangerRead",
                                     printLevel         = printLevel,
                                     inputSource        = inputSource,
                                     readFeature        = "Forward Read",
                                     readFileName       = fastaFileName,
                                     fastaReadName      = forwardName,
                                     geneticCode        = geneticCode)
            })
            names(forwardReadList) <- forwardSelectNames
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseReadList <- lapply(reverseSelectNames, function(reverseName){
                newSangerRead <- new("SangerRead",
                                     printLevel         = printLevel,
                                     inputSource        = inputSource,
                                     readFeature        = "Reverse Read",
                                     readFileName       = fastaFileName,
                                     fastaReadName      = reverseName,
                                     geneticCode        = geneticCode)
            })
            names(reverseReadList) <- reverseSelectNames
        } else if (inputSource == "FASTA" && processMethod == "CSV") {
            if (printLevel == "SangerContig") {
                log_info("  >> You are using CSV Name Conversion Method ",
                         "to group reads in FASTA file (with Csv file)!")
            }
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            if (is.null(contigName)) {
                if (length(unique(csvFile$contig)) != 1) {
                    log_error("Error! There are ",
                              length(unique(csvFile$contig)) ,
                              " contigName in the CSV file. ",
                              "There should be only one contigName in the CSV file.")
                } else {
                    if (printLevel == "SangerContig") {
                        log_info(">> Your contig name is ", contigName)
                    }
                    contigName <- unique(csvFile$contig)
                    contigNameCSV <- as.character(unique(csvFile$contig))
                    selectedCsvFile <- csvFile[csvFile$contig == contigNameCSV, ]
                    forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F", ]
                    reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R", ]
                }
            } else if (!is.null(contigName)) {
                if (printLevel == "SangerContig") {
                    log_info(">> Your contig name is ", contigName)
                }
                selectedCsvFile <- csvFile[csvFile$contig == contigName, ]
                forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F", ]
                reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R", ]
            }
            forwardOriginal <- as.character(forwardCsv$reads)
            reverseOriginal <- as.character(reverseCsv$reads)
            if (!(forwardOriginal %in% names(readFasta) &&
                  reverseOriginal %in% names(readFasta))) {
                log_error("The 'reads' in forwardOriginal is ",
                          "different from read names in FASTA file ('",
                          fastaFileName, "')\n  *'reads': ",
                          paste(forwardOriginal, sep = " "), "\n  *read names ",
                          "in FASTA file: ", paste(names(readFasta), sep = " "))
            }
            # lapply to create SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(forwardOriginal, function(forwardName){
                newSangerRead <- new("SangerRead",
                                     printLevel         = printLevel,
                                     inputSource        = inputSource,
                                     readFeature        = "Forward Read",
                                     readFileName       = fastaFileName,
                                     fastaReadName      = forwardName,
                                     geneticCode        = geneticCode)
            })
            names(forwardReadList) <- forwardOriginal
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseReadList <- lapply(reverseOriginal, function(reverseName){
                newSangerRead <- new("SangerRead",
                                     printLevel         = printLevel,
                                     inputSource        = inputSource,
                                     readFeature        = "Reverse Read",
                                     readFileName       = fastaFileName,
                                     fastaReadName      = reverseName,
                                     geneticCode        = geneticCode)
            })
            names(reverseReadList) <- reverseOriginal
        }
        
        if (inputSource == "FASTA") {
            readNum <- length(forwardReadList) + length(reverseReadList)
            log_info("The number of reads detected: ", readNum)
            forwardReadListFilter <- lapply(forwardReadList, function(read) {
                seqLen <- length(read@primarySeq)
                if (seqLen >= minReadLength) {
                    row <- data.frame(basename(read@readFileName), 
                                      read@objectResults@creationResult, 
                                      "None", "None", 
                                      read@inputSource, read@readFeature)
                    names(row) <- readResultTableName
                    readResultTable <<- rbind(readResultTable, row)
                    read
                } else {
                    msg <- paste0("  >> ", read@readFileName, 
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
                    NULL
                }
            })
            reverseReadListFilter <- lapply(reverseReadList, function(read) {
                seqLen <- length(read@primarySeq)
                if (seqLen >= minReadLength) {
                    row <- data.frame(basename(read@readFileName), 
                                      read@objectResults@creationResult, 
                                      "None", "None", 
                                      read@inputSource, read@readFeature)
                    names(row) <- readResultTableName
                    readResultTable <<- rbind(readResultTable, row)
                    read
                } else {
                    msg <- paste0("  >> ", read@readFileName, 
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
                    NULL
                }
            })
        }
        ### ------------------------------------------------------------------------
        ### 'forwardNumber' + 'reverseNumber' number >= minReadsNum && 2
        ### ------------------------------------------------------------------------
        forwardReadListFilter <- Filter(Negate(is.null), forwardReadListFilter)
        reverseReadListFilter <- Filter(Negate(is.null), reverseReadListFilter)
        forwardNumber <- length(forwardReadListFilter)
        reverseNumber <- length(reverseReadListFilter)
        readNumber <- forwardNumber + reverseNumber
        
        if (readNumber >= minReadsNum) {
            msg <- ""
            if (readNumber >= 2) {
                CSResult <- calculateContigSeq (inputSource      = inputSource,
                                                forwardReadList  = forwardReadListFilter,
                                                reverseReadList  = reverseReadListFilter,
                                                refAminoAcidSeq  = refAminoAcidSeq,
                                                minFractionCall  = minFractionCall,
                                                maxFractionLost  = maxFractionLost,
                                                geneticCode      = geneticCode,
                                                acceptStopCodons = acceptStopCodons,
                                                readingFrame     = readingFrame,
                                                processorsNum    = processorsNum,
                                                printLevel       = printLevel)
                contigGapfree <- CSResult$consensusGapfree
                contigLen <- length(contigGapfree)
                ## This is the only part that is correct!
                diffsDf <- CSResult$diffsDf
                aln2 <- CSResult$aln2
                dist <- CSResult$dist
                dend <- CSResult$dend
                indels <- CSResult$indels
                stopsDf <- CSResult$stopsDf
                spDf <- CSResult$spDf
            } else if (readNumber == 1) {
                if (forwardNumber == 1) {
                    forwardReadListFilter[[1]]
                    primaryDNA <- 
                        SangerReadInnerTrimming(forwardReadListFilter[[1]], inputSource)
                    contigGapfree <- DNAString(primaryDNA)
                } else if (reverseNumber == 1) {
                    primaryDNA <- 
                        SangerReadInnerTrimming(reverseReadListFilter[[1]], inputSource)
                    contigGapfree <- DNAString(primaryDNA)
                }
                diffsDf <- data.frame()
                aln2 <- DNAStringSet()
                dist <- matrix()
                dend <- list()
                indels <- data.frame()
                stopsDf <- data.frame()
                spDf <- data.frame()
            }
            log_success("==========================================================")
            log_success("======== 'SangerContig' S4 instance is created !! ========")
            log_success("==========================================================")
            log_info("   >> ", readNumber, " read(s) created from ", inputSource, " file.")
            if (is.null(namesConversionCSV)) {
                log_info("     >> ", forwardNumber, " reads assigned to 'forward reads' according to 'regular expression'.")
                log_info("     >> ", reverseNumber, " reads assigned to 'reverse reads' according to 'regular expression'.")
            } else {
                log_info("     >> ", forwardNumber, " reads assigned to 'forward reads' according to 'csv file'.")
                log_info("     >> ", reverseNumber, " reads assigned to 'reverse reads' according to 'csv file'.")
            }
            # Here are warning messages
            if (forwardNumber == 0) {
                log_warn("  >> No 'forward read' is detected!")
            }
            if (reverseNumber == 0) {
                log_warn("  >> No 'reverse read' is detected!")
            }
            if (readNumber == 1) {
                log_warn("   >> There is only one read in your SangerContig.")
            }
            
            if (printLevel == "SangerContig") {
                if (TrimmingMethod == "M1") {
                    log_info("   >> Trimmed by 'M1 - Mott’s trimming algorithm'.")
                } else if (TrimmingMethod == "M2") {
                    log_info("   >> Trimmed by 'M2 - sliding window method'.")
                }
            }
        } else {
            msg <- paste0("The number of your total reads is ", readNumber, ".",
                          "\nNumber of total reads has to be equal or more than ",
                          minReadsNum, " ('minReadsNum' that you set)")
            errors[[1]] <- c(errors[[1]], msg)
            errors[[2]] <- c(errors[[2]], "READ_NUMBER_ERROR")
            
            lapply(forwardReadListFilter, function(read) {
                row <- data.frame(basename(read@readFileName), FALSE,
                                  "READ_NUMBER_ERROR", msg,
                                  read@inputSource, read@readFeature)
                names(row) <- readResultTableName
                readResultTable[which(readResultTable$readName == basename(read@readFileName)),] <<- row
            })
            lapply(reverseReadListFilter, function(read) {
                row <- data.frame(basename(read@readFileName), FALSE,
                                  "READ_NUMBER_ERROR", msg,
                                  read@inputSource, read@readFeature)
                names(row) <- readResultTableName
                readResultTable[which(readResultTable$readName == basename(read@readFileName)),] <<- row
            })
        }
    }
    
    if (length(errors[[1]]) != 0) {
        creationResult <- FALSE
        sapply(paste0(errors[[2]], '\n', errors[[1]], '\n') , 
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
        TrimmingMethod         = ""
        minReadsNum            = 0
        minReadLength          = 0
        refAminoAcidSeq        = ""
        minFractionCall        = 0
        maxFractionLost        = 0
        geneticCode            = ""
        acceptStopCodons       = FALSE
        readingFrame           = 0
        contigGapfree          = DNAString()
        diffsDf                = data.frame()
        aln2                   = DNAStringSet()
        dist                   = matrix()
        dend                   = list()
        indels                 = data.frame()
        stopsDf                = data.frame()
        spDf                   = data.frame()
    }
    if (printLevel == "SangerContig") {
        if (nrow(readResultTable) != 0 && ncol(readResultTable) != 0) {
            names(readResultTable) <- readResultTableName
            log_debug("   >> For more information, please run 'object'")
            log_debug("   >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads")
        }
    }
    objectResults <- new("ObjectResults", creationResult = creationResult,
                         errorMessages = errors[[1]], errorTypes = errors[[2]],
                         warningMessages = character(0), warningTypes = character(0),
                         printLevel = printLevel, readResultTable = readResultTable)
    callNextMethod(.Object,
                   objectResults          = objectResults,
                   inputSource            = inputSource,
                   processMethod          = processMethod,
                   fastaFileName          = fastaFileName,
                   namesConversionCSV     = namesConversionCSV,
                   parentDirectory        = parentDirectory,
                   contigName             = contigName,
                   suffixForwardRegExp    = suffixForwardRegExp,
                   suffixReverseRegExp    = suffixReverseRegExp,
                   forwardReadList        = forwardReadListFilter,
                   reverseReadList        = reverseReadListFilter,
                   trimmingMethodSC       = TrimmingMethod,
                   minReadsNum            = minReadsNum,
                   minReadLength          = minReadLength,
                   refAminoAcidSeq        = refAminoAcidSeq,
                   minFractionCall        = minFractionCall,
                   maxFractionLost        = maxFractionLost,
                   geneticCode            = geneticCode,
                   acceptStopCodons       = acceptStopCodons,
                   readingFrame           = readingFrame,
                   contigSeq              = contigGapfree,
                   differencesDF          = diffsDf,
                   alignment              = aln2,
                   distanceMatrix         = dist,
                   dendrogram             = dend,
                   indelsDF               = indels,
                   stopCodonsDF           = stopsDf,
                   secondaryPeakDF        = spDf)
})
