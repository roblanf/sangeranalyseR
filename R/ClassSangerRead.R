#' @title SangerRead
#'
#' @description  An S4 class extending sangerseq S4 class which corresponds to a single ABIF file in Sanger sequencing.
#'
#' @slot objectResults
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot readFeature The direction of the Sanger read. The value must be \code{"Forward Read"} or \code{"Reverse Read"}.
#' @slot readFileName The filename of the target input file.
#' @slot fastaReadName If \code{inputSource} is \code{"FASTA"}, then this value has to be the name of the read inside the FASTA file; if \code{inputSource} is \code{"ABIF"}, then this value is \code{NULL} by default.
#' @slot abifRawData An S4 class containing all fields in the ABIF file. It is the abif class defined in sangerseqR package.
#' @slot QualityReport A S4 class containing quality trimming related inputs and trimming results.
#' @slot ChromatogramParam A S4 class containing chromatogram inputs.
#' @slot primaryAASeqS1 A polypeptide translated from primary DNA sequence starting from the first nucleic acid.
#' @slot primaryAASeqS2 A polypeptide translated from primary DNA sequence starting from the second nucleic acid.
#' @slot primaryAASeqS3 A polypeptide translated from primary DNA sequence starting from the third nucleic acid.
#' @slot geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @slot primarySeqRaw The raw primary sequence from sangerseq class in sangerseqR package before base calling.
#' @slot secondarySeqRaw The raw secondary sequence from sangerseq class in sangerseqR package before base calling.
#' @slot peakPosMatrixRaw The raw peak position matrix from sangerseq class in sangerseqR package before base calling.
#' @slot peakAmpMatrixRaw The raw peak amplitude matrix from sangerseq class in sangerseqR package before base calling.
#'
#' @name SangerRead-class
#' @exportClass SangerRead
#'
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R
#' @import sangerseqR
#' @examples
#' ## Input From ABIF file format
#' # Forward Read
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFFN <- file.path(inputFilesPath,
#'                              "Allolobophora_chlorotica",
#'                              "ACHLO",
#'                              "Achl_ACHLO006-09_1_F.ab1")
#' sangerReadF <- new("SangerRead",
#'                     inputSource           = "ABIF",
#'                     readFeature           = "Forward Read",
#'                     readFileName          = A_chloroticaFFN,
#'                     geneticCode           = GENETIC_CODE,
#'                     TrimmingMethod        = "M1",
#'                     M1TrimmingCutoff      = 0.0001,
#'                     M2CutoffQualityScore  = NULL,
#'                     M2SlidingWindowSize   = NULL,
#'                     baseNumPerRow         = 100,
#'                     heightPerRow          = 200,
#'                     signalRatioCutoff     = 0.33,
#'                     showTrimmed           = TRUE)
#'
#' # Reverse Read
#' A_chloroticaRFN <- file.path(inputFilesPath,
#'                              "Allolobophora_chlorotica",
#'                              "ACHLO",
#'                              "Achl_ACHLO006-09_2_R.ab1")
#' sangerReadR <- new("SangerRead",
#'                     inputSource           = "ABIF",
#'                     readFeature           = "Reverse Read",
#'                     readFileName          = A_chloroticaRFN,
#'                     geneticCode           = GENETIC_CODE,
#'                     TrimmingMethod        = "M1",
#'                     M1TrimmingCutoff      = 0.0001,
#'                     M2CutoffQualityScore  = NULL,
#'                     M2SlidingWindowSize   = NULL,
#'                     baseNumPerRow         = 100,
#'                     heightPerRow          = 200,
#'                     signalRatioCutoff     = 0.33,
#'                     showTrimmed           = TRUE)
#'
#'
#' ## Input From FASTA file format
#' # Forward Read
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFFNfa <- file.path(inputFilesPath,
#'                                "fasta",
#'                                "SangerRead",
#'                                "Achl_ACHLO006-09_1_F.fa")
#' readNameFfa <- "Achl_ACHLO006-09_1_F"
#' sangerReadFfa <- new("SangerRead",
#'                      inputSource        = "FASTA",
#'                      readFeature        = "Forward Read",
#'                      readFileName       = A_chloroticaFFNfa,
#'                      fastaReadName      = readNameFfa,
#'                      geneticCode        = GENETIC_CODE)
#' # Reverse Read
#' A_chloroticaRFNfa <- file.path(inputFilesPath,
#'                                "fasta",
#'                                "SangerRead",
#'                                "Achl_ACHLO006-09_2_R.fa")
#' readNameRfa <- "Achl_ACHLO006-09_2_R"
#' sangerReadRfa <- new("SangerRead",
#'                      inputSource   = "FASTA",
#'                      readFeature   = "Reverse Read",
#'                      readFileName  = A_chloroticaRFNfa,
#'                      fastaReadName = readNameRfa,
#'                      geneticCode   = GENETIC_CODE)
setClass(
    "SangerRead",
    ### ------------------------------------------------------------------------
    ### Input type of each variable of 'SangerContig'.
    ###     * Inherit from 'sangerseq' from sangerseqR.
    ### ------------------------------------------------------------------------
    contains="sangerseq",
    slots=c(objectResults       = "ObjectResults",
            inputSource         = "character",
            readFeature         = "character",
            readFileName        = "character",
            fastaReadName       = "characterORNULL",
            geneticCode         = "character",
            abifRawData         = "abifORNULL",
            QualityReport       = "QualityReportORNULL",
            ChromatogramParam   = "ChromatogramParamORNULL",
            primaryAASeqS1      = "AAString",
            primaryAASeqS2      = "AAString",
            primaryAASeqS3      = "AAString",
            primarySeqRaw       = "DNAString",
            secondarySeqRaw     = "DNAString",
            peakPosMatrixRaw    = "matrix",
            peakAmpMatrixRaw    = "matrix")
)

### ============================================================================
### Overwrite initialize for SangerRead (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerRead",
          function(.Object,
                   printLevel           = "SangerRead",
                   inputSource          = "ABIF",
                   readFeature          = "",
                   readFileName         = "",
                   fastaReadName        = NULL,
                   geneticCode          = GENETIC_CODE,
                   TrimmingMethod       = "M1",
                   M1TrimmingCutoff     = 0.0001,
                   M2CutoffQualityScore = NULL,
                   M2SlidingWindowSize  = NULL,
                   baseNumPerRow        = 100,
                   heightPerRow         = 200,
                   signalRatioCutoff    = 0.33,
                   showTrimmed          = TRUE) {
    creationResult <- TRUE
    errors <- list(character(0), character(0))
    readResultTableName <- c("readName","creationResult", "errorType", 
                             "errorMessage", "inputSource", "direction")
    ############################################################################
    ### First layer of pre-checking: filename exists (Must)
    ############################################################################
    errors <- checkReadFileNameExist (readFileName, errors[[1]], errors[[2]])
    if (length(errors[[1]]) == 0) {
        ########################################################################
        ### Second layer of checking: Check parameters for both ABIF and FASTA
        ########################################################################
        errors <- checkReadFileName (readFileName, inputSource, errors[[1]], errors[[2]])
        if (printLevel == "SangerRead") {
            errors <- checkInputSource (inputSource, errors[[1]], errors[[2]])
            errors <- checkReadFeature (readFeature, errors[[1]], errors[[2]])
            errors <- checkGeneticCode (geneticCode, errors[[1]], errors[[2]])
        }
        if (length(errors[[1]]) == 0) {
            if (inputSource == "ABIF") {
                ################################################################
                ### Third layer of checking: Check parameters for ABIF only
                ################################################################
                ##### ----------------------------------------------------------
                ##### Input parameter prechecking for TrimmingMethod.
                ##### ----------------------------------------------------------
                if (printLevel == "SangerRead") {
                    errors <- checkTrimParam(TrimmingMethod,
                                             M1TrimmingCutoff,
                                             M2CutoffQualityScore,
                                             M2SlidingWindowSize,
                                             errors[[1]], errors[[2]])
                    ##### ----------------------------------------------------------
                    ##### Input parameter prechecking for ChromatogramParam.
                    ##### ----------------------------------------------------------
                    errors <- checkBaseNumPerRow(baseNumPerRow, errors[[1]], errors[[2]])
                    errors <- checkHeightPerRow(heightPerRow, errors[[1]], errors[[2]])
                    errors <- checkSignalRatioCutoff(signalRatioCutoff, errors[[1]], errors[[2]])
                    errors <- checkShowTrimmed(showTrimmed, errors[[1]], errors[[2]])
                }
                if(length(errors[[1]]) == 0) {
                    if (printLevel == "SangerRead") {
                        log_info('------------------------------------------------')
                        log_info("-------- Creating 'SangerRead' instance --------")
                        log_info('------------------------------------------------')
                        log_info('>> ', readFeature, ": Creating abif & sangerseq ...")
                        log_info("    >> Creating ", readFeature , " raw abif ...")
                    }
                    abifRawData = read.abif(readFileName)
                    if (printLevel == "SangerRead") {
                        log_info("    >> Creating ", readFeature , " raw sangerseq ...")
                    }
                    readSangerseq <- sangerseq(abifRawData)
                    primarySeqID <- readSangerseq@primarySeqID
                    secondarySeqID <- readSangerseq@secondarySeqID
                    
                    ### --------------------------------------------------------
                    ### With non-raw & raw primarySeq / secondarySeq
                    ### --------------------------------------------------------
                    if (readFeature == "Forward Read") {
                        primarySeqRaw    <- readSangerseq@primarySeq
                        primarySeq       <- readSangerseq@primarySeq
                        secondarySeqRaw  <- readSangerseq@secondarySeq
                        secondarySeq     <- readSangerseq@secondarySeq
                        traceMatrix      <- readSangerseq@traceMatrix
                        peakPosMatrixRaw <- readSangerseq@peakPosMatrix
                        peakPosMatrix    <- readSangerseq@peakPosMatrix
                        peakAmpMatrixRaw <- readSangerseq@peakAmpMatrix
                        peakAmpMatrix    <- readSangerseq@peakAmpMatrix
                    } else if (readFeature == "Reverse Read") {
                        primarySeqRaw    <- readSangerseq@primarySeq
                        primarySeq       <- readSangerseq@primarySeq
                        secondarySeqRaw  <- readSangerseq@secondarySeq
                        secondarySeq     <- readSangerseq@secondarySeq
                        traceMatrix      <- readSangerseq@traceMatrix
                        peakPosMatrixRaw <- readSangerseq@peakPosMatrix
                        peakPosMatrix    <- readSangerseq@peakPosMatrix
                        peakAmpMatrixRaw <- readSangerseq@peakAmpMatrix
                        peakAmpMatrix    <- readSangerseq@peakAmpMatrix
                    }
                    
                    ### --------------------------------------------------------
                    ### Definition of 'PCON.1' & 'PCON.2'
                    ##### PCON.1: char => Per-base quality values (edited)
                    ##### PCON.2: char => Per-base quality values
                    ### --------------------------------------------------------
                    ### --------------------------------------------------------
                    ### 1. Running 'MakeBaseCall'!
                    ### --------------------------------------------------------
                    ## Reverse the 'traceMatrix' and 'peakPosMatrixRaw' before 
                    ##   running MakeBaseCallsInside function.
                    MBCResult <-
                        MakeBaseCallsInside (traceMatrix, peakPosMatrixRaw,
                                             abifRawData@data$PCON.2,
                                             signalRatioCutoff, readFeature, 
                                             printLevel)
                    ### ========================================================
                    ### 2. Update Once (Only during creation)
                    ###    Basecall primary seq length will be same !
                    ###    Quality Score is reversed in MakeBaseCallsInside 
                    ###    function !!
                    ### ========================================================
                    qualityPhredScores <- MBCResult[["qualityPhredScores"]]
                    ### --------------------------------------------------------
                    ##### 'QualityReport' creation
                    ### --------------------------------------------------------
                    QualityReport <-
                        new("QualityReport",
                            qualityPhredScores    = qualityPhredScores,
                            TrimmingMethod        = TrimmingMethod,
                            M1TrimmingCutoff      = M1TrimmingCutoff,
                            M2CutoffQualityScore  = M2CutoffQualityScore,
                            M2SlidingWindowSize   = M2SlidingWindowSize)
                    
                    ### ========================================================
                    ### 3. Update everytime
                    ### ========================================================
                    primarySeq <- MBCResult[["primarySeq"]]
                    secondarySeq <- MBCResult[["secondarySeq"]]
                    peakPosMatrix <- MBCResult[["peakPosMatrix"]]
                    peakAmpMatrix <- MBCResult[["peakAmpMatrix"]]
                    ### --------------------------------------------------------
                    ##### 'QualityReport' & 'ChromatogramParam' creation
                    ### --------------------------------------------------------
                    ChromatogramParam <-
                        new("ChromatogramParam",
                            baseNumPerRow     = baseNumPerRow,
                            heightPerRow      = heightPerRow,
                            signalRatioCutoff = signalRatioCutoff,
                            showTrimmed       = showTrimmed)
                    trimmedStartPos <- QualityReport@trimmedStartPos
                    trimmedFinishPos <- QualityReport@trimmedFinishPos
                }
            } else if (inputSource == "FASTA") {
                abifRawData <- NULL
                primarySeqID <- "From fasta file"
                secondarySeqID <- ""
                primarySeqRaw <- DNAString()
                if (printLevel == "SangerRead") {
                    log_info('------------------------------------------------')
                    log_info("-------- Creating 'SangerRead' instance --------")
                    log_info('------------------------------------------------')
                    log_info(readFeature, ": Creating SangerRead from FASTA ...")
                }
                readFasta <- read.fasta(readFileName, as.string = TRUE)
                ### ------------------------------------------------------------
                ### Get the Target Filename !!
                ### ------------------------------------------------------------
                fastaNames <- names(readFasta)
                targetFastaName <- fastaNames[fastaNames == fastaReadName]
                if(isEmpty(targetFastaName)) {
                    log_error(paste0("The name '", fastaReadName, 
                                     "' is not in the '", 
                                     basename(readFileName),
                                     "' FASTA file"))
                }
                primarySeq <- 
                    DNAString(as.character(readFasta[[targetFastaName]]))
                secondarySeqRaw   <- DNAString()
                secondarySeq      <- DNAString()
                traceMatrix       <- matrix()
                peakPosMatrixRaw  <- matrix()
                peakPosMatrix     <- matrix()
                peakAmpMatrixRaw  <- matrix()
                peakAmpMatrix     <- matrix()
                QualityReport     <- NULL
                ChromatogramParam <- NULL
                trimmedStartPos <- 0
                trimmedFinishPos <- length(primarySeq)
            }
            
            ## Double check again that there are no errors
            if(length(errors[[1]]) != 0) {
                primaryAASeqS1 <- AAString("")
                primaryAASeqS2 <- AAString("")
                primaryAASeqS3 <- AAString("")
            } else {
                AASeqResult    <- calculateAASeq (primarySeq, trimmedStartPos,
                                                  trimmedFinishPos, geneticCode)
                primaryAASeqS1 <- AASeqResult[["primaryAASeqS1"]]
                primaryAASeqS2 <- AASeqResult[["primaryAASeqS2"]]
                primaryAASeqS3 <- AASeqResult[["primaryAASeqS3"]]
                log_success("--------------------------------------------------------")
                log_success("-------- 'SangerRead' S4 instance is created !! --------")
                log_success("--------------------------------------------------------")
                log_success("   >> '", basename(readFileName), "' is created (", readFeature, "; ", inputSource, ").")
                if (printLevel == "SangerRead") {
                    if (TrimmingMethod == "M1" && printLevel == "SangerRead") {
                        log_info("   >> Read is trimmed by 'M1 - Mott’s trimming algorithm'.")
                    } else if (TrimmingMethod == "M2" && printLevel == "SangerRead") {
                        log_info("   >> Read is trimmed by 'M2 - sliding window method'.")
                    }   
                }
            }
        }
    }
    if (length(errors[[1]]) != 0) {
        creationResult <- FALSE
        readResultTable <- data.frame(basename(readFileName), 
                                      creationResult, errors[[2]], errors[[1]], 
                                      inputSource, readFeature)
        sapply(paste0(errors[[2]], errors[[1]], '\n') , 
               log_error, simplify = FALSE)
        # Create df to store reads that failed to be created
        inputSource         = ""
        fastaReadName       = ""
        readFeature         = ""
        readFileName        = ""
        geneticCode         = ""
        primarySeqID        = ""
        primarySeqRaw       = DNAString("")
        primarySeq          = DNAString("")
        secondarySeqID      = ""
        secondarySeqRaw     = DNAString("")
        secondarySeq        = DNAString("")
        primaryAASeqS1      = AAString("")
        primaryAASeqS2      = AAString("")
        primaryAASeqS3      = AAString("")
        traceMatrix         = matrix()
        peakPosMatrix       = matrix()
        peakPosMatrixRaw    = matrix()
        peakAmpMatrix       = matrix()
        peakAmpMatrixRaw    = matrix()
        abifRawData         = NULL
        QualityReport       = NULL
        ChromatogramParam   = NULL
    } else {
        readResultTable <- data.frame(basename(readFileName), 
                                      creationResult, "None", "None", 
                                      inputSource, readFeature)
    }
    
    if (printLevel == "SangerRead") {
        log_debug("   >> For more information, please run 'object' or 'readTable(object)'.")
        # Print the readResultTable clue
        if (nrow(readResultTable) != 0 && ncol(readResultTable) != 0) {
            names(readResultTable) <- readResultTableName
            log_debug("   >> Run 'object@objectResults@readResultTable' ",
                      "to check the result of the Sanger read")
        }
    }
    names(readResultTable) <- readResultTableName
    objectResults <- new("ObjectResults", creationResult = creationResult,
                         errorMessages = errors[[1]], errorTypes = errors[[2]],
                         warningMessages = character(0), 
                         warningTypes = character(0), printLevel = printLevel, 
                         readResultTable = readResultTable)
    callNextMethod(.Object,
                   objectResults       = objectResults,
                   inputSource         = inputSource,
                   fastaReadName       = fastaReadName,
                   readFeature         = readFeature,
                   readFileName        = readFileName,
                   geneticCode         = geneticCode,
                   primarySeqID        = primarySeqID,
                   primarySeqRaw       = primarySeqRaw,
                   primarySeq          = primarySeq,
                   secondarySeqID      = secondarySeqID,
                   secondarySeqRaw     = secondarySeqRaw,
                   secondarySeq        = secondarySeq,
                   primaryAASeqS1      = primaryAASeqS1,
                   primaryAASeqS2      = primaryAASeqS2,
                   primaryAASeqS3      = primaryAASeqS3,
                   traceMatrix         = traceMatrix,
                   peakPosMatrix       = peakPosMatrix,
                   peakPosMatrixRaw    = peakPosMatrixRaw,
                   peakAmpMatrix       = peakAmpMatrix,
                   peakAmpMatrixRaw    = peakAmpMatrixRaw,
                   abifRawData         = abifRawData,
                   QualityReport       = QualityReport,
                   ChromatogramParam   = ChromatogramParam)
})
