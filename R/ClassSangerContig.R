#' @title SangerContig
#'
#' @description  An S4 class containing forward and reverse SangerRead lists and alignment, consensus read results which corresponds to a contig in Sanger sequencing.
#'
#' @slot inputSource The input source of the raw file. It must be \code{"ABIF"} or \code{"FASTA"}. The default value is \code{"ABIF"}.
#' @slot fastaReadName If \code{inputSource} is \code{"FASTA"}, then this value has to be the FASTA file; if \code{inputSource} is \code{"ABIF"}, then this value is \code{""} by default.
#' @slot namesConversionCSV The file path to the CSV file that provides read names that follow the naming regulation. If \code{inputSource} is \code{"FASTA"}, then users need to prepare the csv file or make sure the original names inside FASTA file are valid; if \code{inputSource} is \code{"ABIF"}, then this value is \code{NULL} by default.
#' @slot parentDirectory The parent directory of all of the reads contained in ABIF format you wish to analyse. In SangerContig, all reads must be in the first layer in this directory.
#' @slot contigName The contig name of all the reads in \code{parentDirectory}.
#' @slot suffixForwardRegExp The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be \code{"_F.ab1"}.
#' @slot suffixReverseRegExp The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be \code{"_R.ab1"}.
#' @slot forwardReadList The list of SangerRead S4 instances which are all forward reads.
#' @slot reverseReadList The list of SangerRead S4 instances which are all reverse reads.
#' @slot trimmingMethodSC The read trimming method for all SangerRead S4 instances in SangerContig. The value must be \code{"M1"} (the default) or \code{'M2'}. All SangerRead must have the same trimming method.
#' @slot minReadsNum The minimum number of reads required to make a consensus sequence, must be 2 or more. The default value is \code{2}.
#' @slot minReadLength Reads shorter than this will not be included in the readset. The default \code{20} means that all reads with length of 20 or more will be included. Note that this is the length of a read after it has been trimmed.
#' @slot refAminoAcidSeq An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is \code{""}.
#' @slot minFractionCall Minimum fraction of the sequences required to call a consensus sequence for SangerContig at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
#' @slot maxFractionLost Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerContig (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
#' @slot geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
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
#' @rdname SangerContig-class
#' @exportClass SangerContig
#'
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerRead.R
#' @examples
#' ## Input From ABIF file format
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
#'                      showTrimmed           = TRUE)
#'
#' ## Input From FASTA file format
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' fastaFN <- file.path(rawDataDir, "fasta",
#'                      "SangerContig", "Achl_ACHLO006-09.fa")
#' namesConversionCSV <- file.path(rawDataDir, "fasta", "SangerContig", "names_conversion_1.csv")
#' contigName <- "Achl_ACHLO006-09"
#' suffixForwardRegExpFa <- "_[0-9]*_F$"
#' suffixReverseRegExpFa <- "_[0-9]*_R$"
#' sangerContigFa <- new("SangerContig",
#'                       inputSource           = "FASTA",
#'                       fastaFileName         = fastaFN,
#'                       namesConversionCSV    = namesConversionCSV,
#'                       contigName            = contigName,
#'                       suffixForwardRegExp   = suffixForwardRegExpFa,
#'                       suffixReverseRegExp   = suffixReverseRegExpFa,
#'                       refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN"
#'                       )
setClass("SangerContig",
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerContig'
         ### -------------------------------------------------------------------
         representation(inputSource               = "character",
                        fastaFileName             = "character",
                        namesConversionCSV        = "characterORNULL",
                        parentDirectory           = "character",
                        contigName                = "character",
                        suffixForwardRegExp       = "character",
                        suffixReverseRegExp       = "character",
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
#' Constructor method of SangerContig Class.
#'
#' @name SangerContig
#' @rdname SangerContig-class
setMethod("initialize",
          "SangerContig",
          function(.Object, ...,
                   inputSource            = "ABIF",
                   fastaFileName          = "",
                   namesConversionCSV     = NULL,
                   parentDirectory        = "",
                   contigName             = "",
                   suffixForwardRegExp    = "_F.ab1",
                   suffixReverseRegExp    = "_R.ab1",
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
                   processorsNum          = NULL) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()
    errors <- checkInputSource (inputSource, errors)
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking for contigSeq parameter
    ### ------------------------------------------------------------------------    errors <- checkMinReadsNum(minReadsNum, errors)
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

    if (length(errors) == 0) {
        processorsNum <- getProcessors (processorsNum)
        if (inputSource == "ABIF") {
            errors <- checkParentDirectory (parentDirectory, errors)
            errors <- checkNamesConversionCSV(namesConversionCSV,
                                              inputSource, errors)
            ### ----------------------------------------------------------------
            ### 'forwardAllReads' & 'reverseAllReads' files prechecking
            ### ----------------------------------------------------------------
            parentDirFiles <- list.files(parentDirectory)
            contigSubGroupFiles <-
                parentDirFiles[grepl(contigName,
                                     parentDirFiles, fixed=TRUE)]
            forwardSelectInputFiles <-
                contigSubGroupFiles[grepl(suffixForwardRegExp,
                                          contigSubGroupFiles)]
            reverseSelectInputFiles <-
                contigSubGroupFiles[grepl(suffixReverseRegExp,
                                          contigSubGroupFiles)]
            forwardAllReads <- lapply(parentDirectory, file.path,
                                      forwardSelectInputFiles)
            reverseAllReads <- lapply(parentDirectory, file.path,
                                      reverseSelectInputFiles)

            ### ----------------------------------------------------------------
            ### 'forwardNumber' + 'reverseNumber' number > 2
            ### ----------------------------------------------------------------
            forwardNumber <- length(forwardAllReads[[1]])
            reverseNumber <- length(reverseAllReads[[1]])
            if ((forwardNumber + reverseNumber) < 2) {
                msg <- "\n'Number of total reads has to be more than two."
                errors <- c(errors, msg)
            }

            ### ----------------------------------------------------------------
            ### 'forwardAllReads'  files prechecking (must exist)
            ### ----------------------------------------------------------------
            forwardAllErrorMsg <- sapply(c(forwardAllReads[[1]]),
                                         function(filePath) {
                if (!file.exists(filePath)) {
                    msg <- paste("\n'", filePath, "' forward read file does ",
                                 "not exist.\n", sep = "")
                    return(msg)
                }
                return()
            })
            errors <- c(errors, unlist(forwardAllErrorMsg), use.names = FALSE)

            ### ----------------------------------------------------------------
            ### 'reverseAllReads'  files prechecking (must exist)
            ### ----------------------------------------------------------------
            reverseAllErrorMsg <- sapply(c(reverseAllReads[[1]]),
                                         function(filePath) {
                if (!file.exists(filePath)) {
                    msg <- paste("\n'", filePath, "'",
                                 " reverse read file does not exist.\n",
                                 sep = "")
                    return(msg)
                }
                return()
            })
            errors <- c(errors, unlist(reverseAllErrorMsg), use.names = FALSE)

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

            ### ----------------------------------------------------------------
            ### Prechecking success. Start to create multiple reads.
            ### ----------------------------------------------------------------
            if (length(errors) != 0) {
                stop(errors)
            }
            trimmingMethodSC = TrimmingMethod
            # sapply to create SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- sapply(forwardAllReads[[1]], function(forwardN){
                new("SangerRead",
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
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseReadList <- sapply(reverseAllReads[[1]], function(reverseN){
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
        } else if (inputSource == "FASTA") {
            errors <- checkFastaFileName(fastaFileName, errors)
            errors <- checkNamesConversionCSV(namesConversionCSV,
                                              inputSource, errors)
            if(length(errors) != 0) {
                stop(errors)
            }
            readFasta <- read.fasta(fastaFileName, as.string = TRUE)
            if (!is.null(namesConversionCSV)) {
                message("    * Reading CSV file and matching names !!")
                fastaNames <- names(readFasta)
                csvFile <- read.csv(namesConversionCSV, header = TRUE)
                tmpFastaNames <- sapply(fastaNames, function(fastaName) {
                    as.character(csvFile[csvFile$original_read_name %in%
                                             fastaName, ]$analysis_read_name)
                })
                names(tmpFastaNames) <- NULL
                names(readFasta) <- tmpFastaNames
            }
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
                contigSubGroupNames[grepl(suffixForwardRegExp,
                                          contigSubGroupNames)]
            ### ----------------------------------------------------------------
            ### Among them, find the reverse names
            ### ----------------------------------------------------------------
            reverseSelectNames <-
                contigSubGroupNames[grepl(suffixReverseRegExp,
                                          contigSubGroupNames)]
            trimmingMethodSC <- ""
            # sapply to create SangerRead list.
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (forward list)
            ### ----------------------------------------------------------------
            forwardReadList <- sapply(forwardSelectNames, function(forwardName){
                new("SangerRead",
                    inputSource        = inputSource,
                    readFeature        = "Forward Read",
                    readFileName       = fastaFileName,
                    fastaReadName      = forwardName,
                    namesConversionCSV = namesConversionCSV,
                    geneticCode        = geneticCode)
            })
            ### ----------------------------------------------------------------
            ### "SangerRead" S4 class creation (reverse list)
            ### ----------------------------------------------------------------
            reverseReadList <- sapply(reverseSelectNames, function(reverseName){
                new("SangerRead",
                    inputSource        = inputSource,
                    readFeature        = "Reverse Read",
                    readFileName       = fastaFileName,
                    fastaReadName      = reverseName,
                    namesConversionCSV = namesConversionCSV,
                    geneticCode        = geneticCode)
            })
        }
        CSResult <- calculateContigSeq (inputSource      = inputSource,
                                        forwardReadList  = forwardReadList,
                                        reverseReadList  = reverseReadList,
                                        refAminoAcidSeq  = refAminoAcidSeq,
                                        minFractionCall  = minFractionCall,
                                        maxFractionLost  = maxFractionLost,
                                        geneticCode      = geneticCode,
                                        acceptStopCodons = acceptStopCodons,
                                        readingFrame     = readingFrame,
                                        processorsNum    = processorsNum)
        contigGapfree <- CSResult$consensusGapfree
        diffsDf <- CSResult$diffsDf
        aln2 <- CSResult$aln2
        dist <- CSResult$dist
        dend <- CSResult$dend
        indels <- CSResult$indels
        stopsDf <- CSResult$stopsDf
        spDf <- CSResult$spDf
        message("  >> 'SangerContig' S4 instance is created !!")
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   inputSource            = inputSource,
                   fastaFileName          = fastaFileName,
                   namesConversionCSV     = namesConversionCSV,
                   parentDirectory        = parentDirectory,
                   contigName             = contigName,
                   suffixForwardRegExp    = suffixForwardRegExp,
                   suffixReverseRegExp    = suffixReverseRegExp,
                   forwardReadList        = forwardReadList,
                   reverseReadList        = reverseReadList,
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
