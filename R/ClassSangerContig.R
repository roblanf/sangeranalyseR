#' @title SangerContig
#'
#' @description  An S4 class containing forward and reverse SangerRead lists and alignment, consensus read results which corresponds to a contig in Sanger sequencing.
#' @slot objectResults This is the object that stores all information of the creation result.
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot processMethod The method to create a contig from reads. The value is \code{"REGEX"} or \code{"CSV"}. The default value is \code{"REGEX"}.
#' @slot ABIF_Directory If \code{inputSource} is \code{"ABIF"}, then this value is the path of a parent directory storing all reads in ABIF format you want to analyse. If \code{inputSource} is \code{"FASTA"}, then this value has to be \code{NULL} by default.
#' @slot FASTA_File If \code{inputSource} is \code{"FASTA"}, then this value has to be the path to a valid FASTA file ; if \code{inputSource} is \code{"ABIF"}, then this value has to be \code{NULL} by default.
#' @slot REGEX_SuffixForward The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented.
#' @slot REGEX_SuffixReverse The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented.
#' @slot CSV_NamesConversion The file path to the CSV file that provides read names, directions, and their contig groups. If \code{processMethod} is \code{"CSV"}, then this value has to be the path to a valid CSV file; if \code{processMethod} is \code{"REGEX"}, then this value has to be \code{NULL} by default.
#' @slot contigName The contig name of all the reads in \code{ABIF_Directory}.
#' @slot geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @slot forwardReadList The list of SangerRead S4 instances which are all forward reads.
#' @slot reverseReadList The list of SangerRead S4 instances which are all reverse reads.
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
#' REGEX_SuffixForward <- "_[0-9]*_F.ab1$"
#' REGEX_SuffixReverse <- "_[0-9]*_R.ab1$"
#' sangerContig <- new("SangerContig",
#'                      ABIF_Directory       = parentDir,
#'                      contigName            = contigName,
#'                      REGEX_SuffixForward   = REGEX_SuffixForward,
#'                      REGEX_SuffixReverse   = REGEX_SuffixReverse)
#'                      
#' ## forward / reverse reads match error
#' ## Input From ABIF file format (Regex)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "ACHLO")
#' contigName <- "Achl_ACHLO006-09"
#' REGEX_SuffixForward <- "_[0-9]*_F.ab1$"
#' REGEX_SuffixReverse <- "_[0-9]*_R.ab1$"
#' sangerContig <- new("SangerContig",
#'                      inputSource           = "ABIF",
#'                      processMethod         = "REGEX",
#'                      ABIF_Directory       = parentDir,
#'                      contigName            = contigName,
#'                      REGEX_SuffixForward   = REGEX_SuffixForward,
#'                      REGEX_SuffixReverse   = REGEX_SuffixReverse,
#'                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                      TrimmingMethod        = "M1",
#'                      M1TrimmingCutoff      = 0.0001,
#'                      baseNumPerRow         = 100,
#'                      heightPerRow          = 200,
#'                      signalRatioCutoff     = 0.33,
#'                      showTrimmed           = TRUE,
#'                      minReadsNum           = 2,
#'                      processorsNum         = 2)
#'
#' ## Input From ABIF file format (Csv three column method)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")
#' CSV_NamesConversion <- file.path(rawDataDir, "ab1", "SangerContig", "names_conversion_2.csv")
#' sangerContig <- new("SangerContig",
#'                      inputSource           = "ABIF",
#'                      processMethod         = "CSV",
#'                      ABIF_Directory        = parentDir,
#'                      CSV_NamesConversion   = CSV_NamesConversion,
#'                      contigName            = "Achl_RBNII384-13",
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
#' REGEX_SuffixForwardFa <- "_[0-9]*_F$"
#' REGEX_SuffixReverseFa <- "_[0-9]*_R$"
#' sangerContigFa <- new("SangerContig",
#'                       inputSource           = "FASTA",
#'                       processMethod         = "REGEX",
#'                       FASTA_File         = fastaFN,
#'                       contigName            = contigName,
#'                       REGEX_SuffixForward   = REGEX_SuffixForwardFa,
#'                       REGEX_SuffixReverse   = REGEX_SuffixReverseFa,
#'                       refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                       processorsNum         = 2)
#'
#' ## Input From FASTA file format (Csv - Csv three column method)
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerContig", "Achl_ACHLO006-09.fa")
#' CSV_NamesConversion <- file.path(rawDataDir, "fasta", "SangerContig", "names_conversion_1.csv")
#' sangerContigFa <- new("SangerContig",
#'                       inputSource           = "FASTA",
#'                       processMethod         = "CSV",
#'                       FASTA_File         = fastaFN,
#'                       CSV_NamesConversion    = CSV_NamesConversion,
#'                       contigName            = "Achl_ACHLO006-09",
#'                       refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                       processorsNum         = 2)
setClass("SangerContig",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerContig'
         ### -------------------------------------------------------------------
         representation(objectResults             = "ObjectResults",
                        inputSource               = "character",
                        processMethod             = "character",
                        ABIF_Directory            = "characterORNULL",
                        FASTA_File                = "characterORNULL",
                        REGEX_SuffixForward       = "characterORNULL",
                        REGEX_SuffixReverse       = "characterORNULL",
                        CSV_NamesConversion       = "characterORNULL",
                        contigName                = "characterORNULL",
                        geneticCode               = "character",
                        forwardReadList           = "list",
                        reverseReadList           = "list",
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
                   ABIF_Directory         = NULL,
                   FASTA_File             = NULL,
                   REGEX_SuffixForward    = NULL,
                   REGEX_SuffixReverse    = NULL,
                   CSV_NamesConversion    = NULL,
                   contigName             = NULL,
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
    creationResult <- TRUE
    errors <- list(character(0), character(0))
    warnings <- list(character(0), character(0))
    readResultTableName <- c("readName","creationResult", "errorType", 
                             "errorMessage", "inputSource", "direction")
    readResultTable <- data.frame()
    if (printLevel == "SangerContig") {
        ########################################################################
        ### SangerContig input parameter pre-checking
        ########################################################################
        errors <- checkInputSource(inputSource, errors[[1]], errors[[2]])
        errors <- checkProcessMethod(inputSource, processMethod, 
                                     errors[[1]], errors[[2]])
        errors <- checkContigName(contigName, errors[[1]], errors[[2]])
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
        ########################################################################
        ### SangerContig object creation
        ########################################################################
        log_info("========================================================")
        log_info("================ Creating 'SangerContig' ===============")
        log_info("========================================================")
        log_info("  >> Contig Name: '", contigName, "'")
        processorsNum <- getProcessors (processorsNum)
        if (inputSource == "ABIF" && processMethod == "REGEX") {
            if (printLevel == "SangerContig") {
                log_info("  >> You are using Regular Expression Method",
                         " to group AB1 files!")
                log_info(">> Your contig name is ", contigName)
            }
            parentDirFiles <- list.files(ABIF_Directory)
            forwardSelectInputFiles <-
                parentDirFiles[grepl(paste0(contigName, REGEX_SuffixForward),
                                     parentDirFiles)]
            reverseSelectInputFiles <-
                parentDirFiles[grepl(paste0(contigName, REGEX_SuffixReverse),
                                     parentDirFiles)]
            warnings <- checkGreplForward(forwardSelectInputFiles, 
                                          warnings[[1]], warnings[[2]])
            warnings <- checkGreplReverse(reverseSelectInputFiles, 
                                          warnings[[1]], warnings[[2]])
            forwardAllReads <- lapply(ABIF_Directory, file.path,
                                      forwardSelectInputFiles)
            reverseAllReads <- lapply(ABIF_Directory, file.path,
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
                                    fastaReadName        = NULL,
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
                                     fastaReadName        = NULL,
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
                log_info(">> Your contig name is ", contigName)
            }
            parentDirFiles <- list.files(ABIF_Directory)
            csvFile <- read.csv(CSV_NamesConversion, header = TRUE)
            selectedCsvFile <- csvFile[csvFile$contig == contigName, ]
            # forward reads CSV matching
            forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F",]
            forwardCsvReads <- as.character(forwardCsv$reads)
            forwardReads <- intersect(parentDirFiles, forwardCsvReads)
            warnings <- checkCSVConvForward(forwardReads, 
                                            warnings[[1]], warnings[[2]])
            fAbsoluteAB1 <- file.path(ABIF_Directory, forwardReads)
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
                                     fastaReadName        = NULL,
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
            # reverse reads CSV matching
            reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R",]
            reverseCsvReads <- as.character(reverseCsv$reads)
            reverseReads <- intersect(parentDirFiles, reverseCsvReads)
            warnings <- checkCSVConvReverse(reverseReads, 
                                            warnings[[1]], warnings[[2]])
            rAbsoluteAB1 <- file.path(ABIF_Directory, reverseReads)
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseReadList <- lapply(rAbsoluteAB1, function(reverseN){
                newSangerRead <- new("SangerRead",
                                     printLevel           = printLevel,
                                     inputSource          = inputSource,
                                     readFeature          = "Reverse Read",
                                     readFileName         = reverseN,
                                     fastaReadName        = NULL,
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
        if (inputSource == "FASTA" && processMethod == "REGEX") {
            if (printLevel == "SangerContig") {
                log_info("  >> You are using Regular Expression Method ",
                         "to group reads in FASTA file (No CSV file)!")
                log_info(">> Your contig name is ", contigName)
            }
            # Read FASTA file
            readFasta <- read.fasta(FASTA_File, as.string = TRUE)
            fastaNames <- names(readFasta)
            ### ----------------------------------------------------------------
            ### Find names with the given contigName
            ### ----------------------------------------------------------------
            contigSubGroupNames <-
                fastaNames[grepl(contigName, fastaNames, fixed=TRUE)]
            ### ----------------------------------------------------------------
            ### Among them, find the forward names
            ### ----------------------------------------------------------------
            forwardSelectNames <-
                contigSubGroupNames[grepl(REGEX_SuffixForward,
                                          contigSubGroupNames)]
            ### ----------------------------------------------------------------
            ### Among them, find the reverse names
            ### ----------------------------------------------------------------
            reverseSelectNames <-
                contigSubGroupNames[grepl(REGEX_SuffixReverse,
                                          contigSubGroupNames)]
            warnings <- checkGreplForward(forwardSelectNames,
                                          warnings[[1]], warnings[[2]])
            warnings <- checkGreplReverse(reverseSelectNames,
                                          warnings[[1]], warnings[[2]])
            # lapply to create SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(forwardSelectNames, function(forwardName){
                newSangerRead <- new("SangerRead",
                                     printLevel           = printLevel,
                                     inputSource          = inputSource,
                                     readFeature          = "Forward Read",
                                     readFileName         = FASTA_File,
                                     fastaReadName        = forwardName,
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
            names(forwardReadList) <- forwardSelectNames
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseReadList <- lapply(reverseSelectNames, function(reverseName){
                newSangerRead <- new("SangerRead",
                                     printLevel           = printLevel,
                                     inputSource          = inputSource,
                                     readFeature          = "Reverse Read",
                                     readFileName         = FASTA_File,
                                     fastaReadName        = reverseName,
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
            names(reverseReadList) <- reverseSelectNames
        } else if (inputSource == "FASTA" && processMethod == "CSV") {
            if (printLevel == "SangerContig") {
                log_info("  >> You are using CSV Name Conversion Method ",
                         "to group reads in FASTA file (with Csv file)!")
                log_info(">> Your contig name is ", contigName)
            }
            # Read FASTA file
            readFasta <- read.fasta(FASTA_File, as.string = TRUE)
            fastaNames <- names(readFasta)
            csvFile <- read.csv(CSV_NamesConversion, header = TRUE)
            selectedCsvFile <- csvFile[csvFile$contig == contigName, ]
            # forward CSV matching
            forwardCsv <- selectedCsvFile[selectedCsvFile$direction == "F", ]
            forwardCsvReads <- as.character(forwardCsv$reads)
            forwardReads <- intersect(fastaNames, forwardCsvReads)
            warnings <- checkCSVConvForward(forwardReads, 
                                            warnings[[1]], warnings[[2]])
            # lapply to create SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- lapply(forwardReads, function(forwardName){
                newSangerRead <- new("SangerRead",
                                     printLevel         = printLevel,
                                     inputSource        = inputSource,
                                     readFeature        = "Forward Read",
                                     readFileName       = FASTA_File,
                                     fastaReadName      = forwardName,
                                     geneticCode        = geneticCode,
                                     TrimmingMethod       = TrimmingMethod,
                                     M1TrimmingCutoff     = M1TrimmingCutoff,
                                     M2CutoffQualityScore = M2CutoffQualityScore,
                                     M2SlidingWindowSize  = M2SlidingWindowSize,
                                     baseNumPerRow        = baseNumPerRow,
                                     heightPerRow         = heightPerRow,
                                     signalRatioCutoff    = signalRatioCutoff,
                                     showTrimmed          = showTrimmed)
            })
            names(forwardReadList) <- forwardReads
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            # reverse CSV matching
            reverseCsv <- selectedCsvFile[selectedCsvFile$direction == "R", ]
            reverseCsvReads <- as.character(reverseCsv$reads)
            reverseReads <- intersect(fastaNames, reverseCsvReads)
            warnings <- checkCSVConvReverse(reverseReads, 
                                            warnings[[1]], warnings[[2]])
            reverseReadList <- lapply(reverseReads, function(reverseName){
                newSangerRead <- new("SangerRead",
                                     printLevel           = printLevel,
                                     inputSource          = inputSource,
                                     readFeature          = "Reverse Read",
                                     readFileName         = FASTA_File,
                                     fastaReadName        = reverseName,
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
            names(reverseReadList) <- reverseReads
        }
        
        
        ########################################################################
        ### ** 1. Post-creation check: minReadLength check 
        ### seqLen >= minReadLength
        ########################################################################
        readNum <- length(forwardReadList) + length(reverseReadList)
        log_info("   >> The number of reads detected: ", readNum)
        forwardReadListFilter <- lapply(forwardReadList, function(read) {
            if (read@objectResults@creationResult) {
                if (inputSource == "ABIF") {
                    seqLen <- read@QualityReport@trimmedFinishPos -
                        read@QualityReport@trimmedStartPos
                } else if (inputSource == "FASTA") {
                    seqLen <- length(read@primarySeq)
                }
                if (seqLen >= minReadLength) {
                    ### --------------------------------------------------------
                    ### Success: readResultTable (SangerContig Level)
                    ### --------------------------------------------------------
                    read
                } else {
                    msg <- paste0("  >> ", read@readFileName, 
                                  " is shorter than 'minReadLength' ",
                                  minReadLength,". This read is created but skipped!\n")
                    warnings[[1]] <<- c(warnings[[1]], msg)
                    warnings[[2]] <<- c(warnings[[2]], "MIN_READ_LENGTH_SKIPPED_WARN")
                    ### --------------------------------------------------------
                    ### Failed: readResultTable (SangerContig Level)
                    ### --------------------------------------------------------
                    if (inputSource == "ABIF") {
                        row <- data.frame(basename(read@readFileName), 
                                          FALSE, "MIN_READ_LENGTH_ERROR", msg, 
                                          read@inputSource,  read@readFeature)
                    } else if (inputSource == "FASTA") {
                        row <- data.frame(basename(read@fastaReadName), 
                                          FALSE, "MIN_READ_LENGTH_ERROR", msg, 
                                          read@inputSource,  read@readFeature)
                    }
                    names(row) <- readResultTableName
                    readResultTable <<- rbind(readResultTable, row)
                    NULL
                }
            } else {
                ### ------------------------------------------------------------
                ### Failed: readResultTable (SangerRead Level)
                ### ------------------------------------------------------------
                if (inputSource == "ABIF") {
                    row <- data.frame(basename(read@readFileName), 
                                      read@objectResults@creationResult, 
                                      read@objectResults@errorTypes, 
                                      read@objectResults@errorMessages, 
                                      read@inputSource, read@readFeature)
                } else if (inputSource == "FASTA") {
                    row <- data.frame(basename(read@fastaReadName), 
                                      read@objectResults@creationResult, 
                                      read@objectResults@errorTypes, 
                                      read@objectResults@errorMessages, 
                                      read@inputSource, read@readFeature)
                }
                names(row) <- readResultTableName
                readResultTable <<- rbind(readResultTable, row)
                NULL
            }
        })
        reverseReadListFilter <- lapply(reverseReadList, function(read) {
            if (read@objectResults@creationResult) {
                if (inputSource == "ABIF") {
                    seqLen <- read@QualityReport@trimmedFinishPos -
                        read@QualityReport@trimmedStartPos
                } else if (inputSource == "FASTA") {
                    seqLen <- length(read@primarySeq)
                }
                if (seqLen >= minReadLength) {
                    ### --------------------------------------------------------
                    ### Success: readResultTable (SangerContig Level)
                    ### --------------------------------------------------------
                    read
                } else {
                    msg <- paste0("  >> ", read@readFileName, 
                                  " is shorter than 'minReadLength' ",
                                  minReadLength,". This read is created but skipped!\n")
                    warnings[[1]] <<- c(warnings[[1]], msg)
                    warnings[[2]] <<- c(warnings[[2]], "MIN_READ_LENGTH_SKIPPED_WARN")
                    ### --------------------------------------------------------
                    ### Failed: readResultTable (SangerContig Level)
                    ### --------------------------------------------------------
                    if (inputSource == "ABIF") {
                        row <- data.frame(basename(read@readFileName), 
                                          FALSE, "MIN_READ_LENGTH_ERROR", msg, 
                                          read@inputSource,  read@readFeature)
                    } else if (inputSource == "FASTA") {
                        row <- data.frame(basename(read@fastaReadName), 
                                          FALSE, "MIN_READ_LENGTH_ERROR", msg, 
                                          read@inputSource,  read@readFeature)
                    }
                    names(row) <- readResultTableName
                    readResultTable <<- rbind(readResultTable, row)
                    NULL
                }
            } else {
                ### ------------------------------------------------------------
                ### Failed: readResultTable (SangerRead Level)
                ### ------------------------------------------------------------
                if (inputSource == "ABIF") {
                    row <- data.frame(basename(read@readFileName), 
                                      read@objectResults@creationResult, 
                                      read@objectResults@errorTypes, 
                                      read@objectResults@errorMessages, 
                                      read@inputSource, read@readFeature)
                } else if (inputSource == "FASTA") {
                    row <- data.frame(basename(read@fastaReadName), 
                                      read@objectResults@creationResult, 
                                      read@objectResults@errorTypes, 
                                      read@objectResults@errorMessages, 
                                      read@inputSource, read@readFeature)
                }
                names(row) <- readResultTableName
                readResultTable <<- rbind(readResultTable, row)
                NULL
            }                 
        })            
        
        ########################################################################
        ### ** 2. Post-creation check: minReadsNum check 
        ### 'forwardNumber' + 'reverseNumber' number >= minReadsNum && 2
        ########################################################################
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
            ### ----------------------------------------------------------------
            ### Success: readResultTable (SangerContig Level)
            ### ----------------------------------------------------------------
            log_success("==========================================================")
            log_success("======== 'SangerContig' S4 instance is created !! ========")
            log_success("==========================================================")
            log_info("   >> ", readNumber, " read(s) created from ", inputSource, " file.")
            if (processMethod == "REGEX") {
                log_info("     >> ", forwardNumber, " reads assigned to 'forward reads' according to 'regular expression'.")
                log_info("     >> ", reverseNumber, " reads assigned to 'reverse reads' according to 'regular expression'.")
            } else if (processMethod == "CSV") {
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
            ### ----------------------------------------------------------------
            ### Failed: readResultTable (SangerContig Level)
            ### ----------------------------------------------------------------
            msg <- paste0("The number of your total reads is ", readNumber, ".",
                          "\nNumber of total reads has to be equal or more than ",
                          minReadsNum, " ('minReadsNum' that you set)")
            errors[[1]] <- c(errors[[1]], msg)
            errors[[2]] <- c(errors[[2]], "READ_NUMBER_ERROR")
            
            lapply(forwardReadListFilter, function(read) {
                if (inputSource == "ABIF") {
                    row <- data.frame(basename(read@readFileName), FALSE,
                                      "READ_NUMBER_ERROR", msg,
                                      read@inputSource, read@readFeature)
                } else if (inputSource == "FASTA") {
                    row <- data.frame(basename(read@fastaReadName), FALSE,
                                      "READ_NUMBER_ERROR", msg,
                                      read@inputSource, read@readFeature)
                }
                names(row) <- readResultTableName
                readResultTable <<- rbind(readResultTable, row)
            })
            lapply(reverseReadListFilter, function(read) {
                if (inputSource == "ABIF") {
                    row <- data.frame(basename(read@readFileName), FALSE,
                                      "READ_NUMBER_ERROR", msg,
                                      read@inputSource, read@readFeature)
                } else if (inputSource == "FASTA") {
                    row <- data.frame(basename(read@fastaReadName), FALSE,
                                      "READ_NUMBER_ERROR", msg,
                                      read@inputSource, read@readFeature)
                }
                names(row) <- readResultTableName
                readResultTable <<- rbind(readResultTable, row)
            })
        }
        ########################################################################
        ### ** 3. Passing post-creation check
        ### Add all successfully created reads into 'readResultTable'
        ########################################################################
        lapply(forwardReadListFilter, function(read) {
            if (!basename(read@readFileName) %in% readResultTable$readName) {
                if (inputSource == "ABIF") {
                    row <- data.frame(basename(read@readFileName),
                                      read@objectResults@creationResult,
                                      "None", "None",
                                      read@inputSource, read@readFeature)
                } else if (inputSource == "FASTA") {
                    row <- data.frame(basename(read@fastaReadName),
                                      read@objectResults@creationResult,
                                      "None", "None",
                                      read@inputSource, read@readFeature)
                }
                names(row) <- readResultTableName
                readResultTable <<- rbind(readResultTable, row)
            }
        })
        lapply(reverseReadListFilter, function(read) {
            if (!basename(read@readFileName) %in% readResultTable$readName) {
                if (inputSource == "ABIF") {
                    row <- data.frame(basename(read@readFileName),
                                      read@objectResults@creationResult,
                                      "None", "None",
                                      read@inputSource, read@readFeature)
                } else if (inputSource == "FASTA") {
                    row <- data.frame(basename(read@fastaReadName),
                                      read@objectResults@creationResult,
                                      "None", "None",
                                      read@inputSource, read@readFeature)
                }
                names(row) <- readResultTableName
                readResultTable <<- rbind(readResultTable, row)
            }
        })
    }
    if (length(warnings[[1]]) != 0) {
        sapply(paste0(warnings[[2]], '\n', warnings[[1]], '\n') , 
               log_warn, simplify = FALSE)
    }
    if (length(errors[[1]]) != 0) {
        creationResult <- FALSE
        sapply(paste0(errors[[2]], '\n', errors[[1]], '\n') , 
               log_error, simplify = FALSE)
        inputSource            = ""
        processMethod          = ""
        FASTA_File          = NULL
        CSV_NamesConversion     = NULL
        ABIF_Directory        = NULL
        contigName             = NULL
        REGEX_SuffixForward    = NULL
        REGEX_SuffixReverse    = NULL
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
        }
        log_debug("   >> For more information, please run 'object'")
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
                   FASTA_File             = FASTA_File,
                   CSV_NamesConversion    = CSV_NamesConversion,
                   ABIF_Directory         = ABIF_Directory,
                   contigName             = contigName,
                   REGEX_SuffixForward    = REGEX_SuffixForward,
                   REGEX_SuffixReverse    = REGEX_SuffixReverse,
                   forwardReadList        = forwardReadListFilter,
                   reverseReadList        = reverseReadListFilter,
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
