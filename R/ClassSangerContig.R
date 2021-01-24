#' @title SangerContig
#'
#' @description  An S4 class containing forward and reverse SangerRead lists and alignment, consensus read results which corresponds to a contig in Sanger sequencing.
#'
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot fastaFileName If \code{inputSource} is \code{"FASTA"}, then this value has to be the FASTA file; if \code{inputSource} is \code{"ABIF"}, then this value is \code{NULL} by default.
#' @slot namesConversionCSV The file path to the CSV file that provides read names that follow the naming regulation. If \code{inputSource} is \code{"FASTA"}, then users need to prepare the csv file or make sure the original names inside FASTA file are valid; if \code{inputSource} is \code{"ABIF"}, then this value is \code{NULL} by default.
#' @slot parentDirectory If \code{inputSource} is \code{"ABIF"}, then this value is the path of the parent directory storing all reads in ABIF format you wish to analyse and cannot be NULL. In SangerContig, all reads must be in the first layer in this directory. If \code{inputSource} is \code{"FASTA"}, then this value is \code{NULL} by default.
#' @slot contigName The contig name of all the reads in \code{parentDirectory}.
#' @slot suffixForwardRegExp The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be \code{"_F.ab1"}.
#' @slot suffixReverseRegExp The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be \code{"_R.ab1"}.
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
#' ## Input From ABIF file format (Regex)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")
#' contigName <- "Achl_RBNII384-13"
#' suffixForwardRegExp <- "_[0-9]*_F.ab1"
#' suffixReverseRegExp <- "_[0-9]*_R.ab1"
#' sangerContig <- new("SangerContig",
#'                      inputSource           = "ABIF",
#'                      parentDirectory       = parentDir,
#'                      contigName            = contigName,
#'                      suffixForwardRegExp   = suffixForwardRegExp,
#'                      suffixReverseRegExp   = suffixReverseRegExp,
#'                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                      TrimmingMethod        = "M1",
#'                      M1TrimmingCutoff      = 0.0001,
#'                      M2CutoffQualityScore  = NULL,
#'                      M2SlidingWindowSize   = NULL,
#'                      baseNumPerRow         = 100,
#'                      heightPerRow          = 200,
#'                      signalRatioCutoff     = 0.33,
#'                      showTrimmed           = TRUE,
#'                      processorsNum         = 2)
#'
#' ## Input From ABIF file format (Csv three column method)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")
#' namesConversionCSV <- file.path(rawDataDir, "ab1", "SangerContig", "names_conversion_2.csv")
#' sangerContig <- new("SangerContig",
#'                      inputSource           = "ABIF",
#'                      parentDirectory       = parentDir,
#'                      namesConversionCSV    = namesConversionCSV,
#'                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                      TrimmingMethod        = "M1",
#'                      M1TrimmingCutoff      = 0.0001,
#'                      M2CutoffQualityScore  = NULL,
#'                      M2SlidingWindowSize   = NULL,
#'                      baseNumPerRow         = 100,
#'                      heightPerRow          = 200,
#'                      signalRatioCutoff     = 0.33,
#'                      showTrimmed           = TRUE,
#'                      processorsNum         = 2)
#'
#'
#' ## Input From FASTA file format (No Csv - Regex)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerContig", "Achl_ACHLO006-09.fa")
#' contigName <- "Achl_ACHLO006-09"
#' suffixForwardRegExpFa <- "_[0-9]*_F$"
#' suffixReverseRegExpFa <- "_[0-9]*_R$"
#' sangerContigFa <- new("SangerContig",
#'                       inputSource           = "FASTA",
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
#'                       fastaFileName         = fastaFN,
#'                       namesConversionCSV    = namesConversionCSV,
#'                       refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                       processorsNum         = 2)
setClass("SangerContig",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerContig'
         ### -------------------------------------------------------------------
         representation(inputSource               = "character",
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
                   inputSource            = "ABIF",
                   fastaFileName          = NULL,
                   namesConversionCSV     = NULL,
                   parentDirectory        = NULL,
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
                   processorsNum          = NULL,
                   logLevel               = TRUE) {
    errors <- character()
    warnings <- character()
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- checkInputSource (inputSource, errors)
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking for contigSeq parameter
    ### ------------------------------------------------------------------------
    errors <- checkMinReadsNum(minReadsNum, errors)
    errors <- checkMinReadLength(minReadLength, errors)
    errors <- checkMinFractionCall(minFractionCall, errors)
    errors <- checkMaxFractionLost(maxFractionLost, errors)
    errors <- checkGeneticCode(geneticCode, errors)
    errors <- checkAcceptStopCodons(acceptStopCodons, errors)
    errors <- checkReadingFrame(readingFrame, errors)
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking for processorsNum
    ### ------------------------------------------------------------------------
    errors <- checkProcessorsNum(processorsNum, errors)
    if (length(errors) != 0 ) {
        log_error(paste(errors, collapse = ""))
    }
    processorsNum <- getProcessors (processorsNum)
    log_info("******** Contig Name: ", contigName)
    if (inputSource == "ABIF") {
        ### --------------------------------------------------------------------
        ### 'ABIF' condition checking!
        ### --------------------------------------------------------------------
        ab1RegexChecker <- is.null(namesConversionCSV) &&
            !is.null(contigName) &&
            !is.null(suffixForwardRegExp) &&
            !is.null(suffixReverseRegExp)
        ab1CSVChecker <- !is.null(namesConversionCSV)
        errors <- checkParentDirectory (parentDirectory, errors)
        if (ab1RegexChecker) {
            errors <- checkNamesConversionCSV(logLevel, parentDirectory,
                                              fastaFileName, namesConversionCSV,
                                              inputSource, "ab1Regex", errors)
            if (length(errors) != 0 ) {
                log_error(paste(errors, collapse = ""))
            }
            # Regex input
            log_info("**** You are using Regular Expression Method",
                     " to group AB1 files!")
            parentDirFiles <- list.files(parentDirectory)

            forwardSelectInputFiles <-
                parentDirFiles[grepl(paste0(contigName, suffixForwardRegExp),
                                     parentDirFiles)]
            reverseSelectInputFiles <-
                parentDirFiles[grepl(paste0(contigName, suffixReverseRegExp),
                                     parentDirFiles)]
            # print("** forwardSelectInputFiles: ")
            # print(forwardSelectInputFiles)
            # print("** reverseSelectInputFiles")
            # print(reverseSelectInputFiles)
            
            forwardAllReads <- lapply(parentDirectory, file.path,
                                      forwardSelectInputFiles)
            reverseAllReads <- lapply(parentDirectory, file.path,
                                      reverseSelectInputFiles)
            ### ----------------------------------------------------------------
            ### 'forwardAllReads'  files prechecking (must exist)
            ### ----------------------------------------------------------------
            forwardAllErrorMsg <-
                lapply(c(forwardAllReads[[1]]),
                       function(filePath) {
                           if (!file.exists(filePath)) {
                               msg <- paste("\n'", filePath,
                                            "' forward read file does ",
                                            "not exist.\n", sep = "")
                               return(msg)
                           }
                           return()
                       })
            errors <- c(errors, unlist(forwardAllErrorMsg), use.names = FALSE)

            ### ----------------------------------------------------------------
            ### 'reverseAllReads'  files prechecking (must exist)
            ### ----------------------------------------------------------------
            reverseAllErrorMsg <-
                lapply(c(reverseAllReads[[1]]),
                       function(filePath) {
                           if (!file.exists(filePath)) {
                               msg <-
                                   paste("\n'", filePath, "'",
                                         " reverse read file does not exist.\n",
                                         sep = "")
                               return(msg)
                           }
                           return()
                       })
            errors <- c(errors, unlist(reverseAllErrorMsg), use.names = FALSE)
        } else if (ab1CSVChecker) {
            errors <- checkNamesConversionCSV(logLevel, parentDirectory,
                                              fastaFileName, namesConversionCSV,
                                              inputSource, "ab1CSV", errors)
            if (length(errors) != 0 ) {
                log_error(paste(errors, collapse = ""))
            }
            # CSV input
            log_info("**** You are using CSV Name Conversion Method ",
                    "to group AB1 files!")
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            if (is.null(contigName)) {
                if (length(unique(csvFile$contig)) != 1) {
                    log_error("Error! There are ",
                         length(unique(csvFile$contig)) ,
                         " contigName in the CSV file. ",
                         "There should be only one contigName in the CSV file.")
                } else {
                    log_info("**** Contig number in your Csv file is ",
                             length(unique(csvFile$contig)))
                    contigNameCSV <- as.character(unique(csvFile$contig))
                    selectedCsvFile <- csvFile[csvFile$contig == contigNameCSV, ]
                    forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F", ]
                    reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R", ]
                }
            } else if (!is.null(contigName)) {
                log_info("**** Your contig name is ", contigName)
                selectedCsvFile <- csvFile[csvFile$contig == contigName, ]
                forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F",]
                reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R",]
            }
        }
        ### ----------------------------------------------------------------
        ### Input parameter prechecking for TrimmingMethod.
        ### ----------------------------------------------------------------
        errors <- checkTrimParam(TrimmingMethod,
                                 M1TrimmingCutoff,
                                 M2CutoffQualityScore,
                                 M2SlidingWindowSize,
                                 errors)
        ### ----------------------------------------------------------------
        ##### Input parameter prechecking for ChromatogramParam
        ### ----------------------------------------------------------------
        errors <- checkBaseNumPerRow (baseNumPerRow, errors)
        errors <- checkHeightPerRow (baseNumPerRow, errors)
        errors <- checkSignalRatioCutoff (signalRatioCutoff, errors)
        errors <- checkShowTrimmed (showTrimmed, errors)
        if (length(errors) != 0) {
            log_error(paste(errors, collapse = ""))
        }
        ### ----------------------------------------------------------------
        ### Prechecking success. Start to create multiple reads.
        ### ----------------------------------------------------------------
        trimmingMethodSC <- TrimmingMethod
        if (ab1RegexChecker) {
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(forwardAllReads[[1]], function(forwardN){
                newSangerRead<- new("SangerRead",
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
        } else if (ab1CSVChecker) {
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
        forwardReadListFilter <- lapply(forwardReadList, function(read) {
            trimmedLen <- read@QualityReport@trimmedFinishPos -
                read@QualityReport@trimmedStartPos
            if (trimmedLen >= minReadLength) {
                read
            } else {
                NULL
            }
        })
        reverseReadListFilter <- lapply(reverseReadList, function(read) {
            trimmedLen <- read@QualityReport@trimmedFinishPos -
                read@QualityReport@trimmedStartPos
            if (trimmedLen >= minReadLength) {
                read
            } else {
                NULL
            }
        })
    } else if (inputSource == "FASTA") {
        ### --------------------------------------------------------------------
        ### 'FASTA' condition checking!
        ### --------------------------------------------------------------------
        csvRegexChecker <- is.null(namesConversionCSV) &&
            !is.null(contigName) &&
            !is.null(suffixForwardRegExp) &&
            !is.null(suffixReverseRegExp)
        csvCSVChecker <- !is.null(namesConversionCSV)
        errors <- checkFastaFileName(fastaFileName, errors)
        if(length(errors) != 0) {
            log_error(paste(errors, collapse = ""))
        }
        readFasta <- read.fasta(fastaFileName, as.string = TRUE)
        fastaNames <- names(readFasta)
        if (csvRegexChecker) {
            errors <- checkNamesConversionCSV(logLevel, parentDirectory,
                                              fastaFileName, namesConversionCSV,
                                              inputSource, "csvRegex", errors)
            if (length(errors) != 0 ) {
                log_error(paste(errors, collapse = ""))
            }
            # Csv-Regex input
            log_info("**** You are using Regular Expression Method ",
                    "to group reads in FASTA file (No CSV file)!")
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
            trimmingMethodSC <- ""
            # lapply to create SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(forwardSelectNames, function(forwardName){
                newSangerRead <- new("SangerRead",
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
                                     inputSource        = inputSource,
                                     readFeature        = "Reverse Read",
                                     readFileName       = fastaFileName,
                                     fastaReadName      = reverseName,
                                     geneticCode        = geneticCode)
            })
            names(reverseReadList) <- reverseSelectNames
        } else if (csvCSVChecker) {
            errors <- checkNamesConversionCSV(logLevel, parentDirectory,
                                              fastaFileName, namesConversionCSV,
                                              inputSource, "csvCSV", errors)
            if (length(errors) != 0 ) {
                log_error(paste(errors, collapse = ""))
            }
            # Csv-CSV input
            log_info("**** You are using CSV Name Conversion Method ",
                    "to group reads in FASTA file (with Csv file)!")
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            if (is.null(contigName)) {
                if (length(unique(csvFile$contig)) != 1) {
                    log_error("Error! There are ",
                         length(unique(csvFile$contig)) ,
                         " contigName in the CSV file. ",
                         "There should be only one contigName in the CSV file.")
                } else {
                    log_info("**** Contig number in your Csv file is ", length(unique(csvFile$contig)))
                    contigNameCSV <- as.character(unique(csvFile$contig))
                    selectedCsvFile <- csvFile[csvFile$contig == contigNameCSV, ]
                    forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F", ]
                    reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R", ]
                }
            } else if (!is.null(contigName)) {
                log_info("**** Your contig name is ", contigName)
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
                                     inputSource        = inputSource,
                                     readFeature        = "Reverse Read",
                                     readFileName       = fastaFileName,
                                     fastaReadName      = reverseName,
                                     geneticCode        = geneticCode)
            })
            names(reverseReadList) <- reverseOriginal
            trimmingMethodSC <- ""
        }

        forwardReadListFilter <- lapply(forwardReadList, function(read) {
            seqLen <- length(read@primarySeq)
            if (seqLen >= minReadLength) {
                read
            } else {
                log_info("   * Read length is shorter than 'minReadLength' ",
                        minReadLength, ".\n     This read is skipped!!")
                NULL
            }
        })
        reverseReadListFilter <- lapply(reverseReadList, function(read) {
            seqLen <- length(read@primarySeq)
            if (seqLen >= minReadLength) {
                read
            } else {
                log_info("   * Read length is shorter than 'minReadLength' ",
                         minReadLength, ".\n     This read is skipped!!")
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
    if (readNumber >= minReadsNum && readNumber >= 2) {
        CSResult <- calculateContigSeq (inputSource      = inputSource,
                                        forwardReadList  = forwardReadListFilter,
                                        reverseReadList  = reverseReadListFilter,
                                        refAminoAcidSeq  = refAminoAcidSeq,
                                        minFractionCall  = minFractionCall,
                                        maxFractionLost  = maxFractionLost,
                                        geneticCode      = geneticCode,
                                        acceptStopCodons = acceptStopCodons,
                                        readingFrame     = readingFrame,
                                        processorsNum    = processorsNum)
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
        log_success("**********************************************************")
        log_success("******** 'SangerContig' S4 instance is created !! ********")
        log_success("**********************************************************")
        log_success("  * >> ", readNumber, " reads are created from ", inputSource, " file.")
        if (is.null(namesConversionCSV)) {
            log_success("  * >> ", forwardNumber, " reads assigned to 'forward reads' according to 'regular expression'.")
            log_success("  * >> ", reverseNumber, " reads assigned to 'reverse reads' according to 'regular expression'.")
        } else {
            log_success("  * >> ", forwardNumber, " reads assigned to 'forward reads' according to 'csv file'.")
            log_success("  * >> ", reverseNumber, " reads assigned to 'reverse reads' according to 'csv file'.")
        }
        if (TrimmingMethod == "M1") {
            log_success("  * >> Read is trimmed by 'M1 - Mott’s trimming algorithm'.")
        } else if (TrimmingMethod == "M2") {
            log_success("  * >> Read is trimmed by 'M2 - sliding window method'.")
        }
        log_success("  * >> For more information, please run 'readTable(object)'.")
    } else if (readNumber >= minReadsNum && readNumber == 1) {
        msg <- paste("There is only one read in your SangerContig.")
        warnings <- c(warnings, msg)
        invisible(lapply(warnings, log_warn))
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
        log_success("**********************************************************")
        log_success("******** 'SangerContig' S4 instance is created !! ********")
        log_success("**********************************************************")
        log_success("  * >> 1 read is created from ", inputSource, " file.")
        if (is.null(namesConversionCSV)) {
            log_success("  * >> ", forwardNumber, " reads assigned to 'forward reads' according to 'regular expression'.")
            log_success("  * >> ", reverseNumber, " reads assigned to 'reverse reads' according to 'regular expression'.")
        } else {
            log_success("  * >> ", forwardNumber, " reads assigned to 'forward reads' according to 'csv file'.")
            log_success("  * >> ", reverseNumber, " reads assigned to 'reverse reads' according to 'csv file'.")
        }
        if (TrimmingMethod == "M1") {
            log_success("  * >> Read is trimmed by 'M1 - Mott’s trimming algorithm'.")
        } else if (TrimmingMethod == "M2") {
            log_success("  * >> Read is trimmed by 'M2 - sliding window method'.")
        }
        log_success("  * >> For more information, please run 'readTable(object)'.")
    } else {
        msg <- paste("The number of your total reads is ", readNumber, ".",
                     "\nNumber of total reads has to be equal or more than ",
                     minReadsNum, " ('minReadsNum' that you set)", sep = "")
        warnings <- c(warnings, msg)
        invisible(lapply(warnings, log_warn))
        contigGapfree <- DNAString()
        diffsDf <- data.frame()
        aln2 <- DNAStringSet()
        dist <- matrix()
        dend <- list()
        indels <- data.frame()
        stopsDf <- data.frame()
        spDf <- data.frame()
    }
    if (length(errors) != 0) {
        log_error(paste(errors, collapse = ""))
    }
    callNextMethod(.Object,
                   inputSource            = inputSource,
                   fastaFileName          = fastaFileName,
                   namesConversionCSV     = namesConversionCSV,
                   parentDirectory        = parentDirectory,
                   contigName             = contigName,
                   suffixForwardRegExp    = suffixForwardRegExp,
                   suffixReverseRegExp    = suffixReverseRegExp,
                   forwardReadList        = forwardReadListFilter,
                   reverseReadList        = reverseReadListFilter,
                   trimmingMethodSC       = trimmingMethodSC,
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
