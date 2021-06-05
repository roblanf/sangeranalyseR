checkInputSource <- function(inputSource, errors, errorTypes) {
    if (typeof(inputSource) != "character") {
        msg<- "'inputSource' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if (inputSource != "ABIF" && inputSource != "FASTA") {
            msg <- "'inputSource' must be 'ABIF' or 'FASTA'"
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkProcessMethod <- function(inputSource, processMethod, errors, errorTypes) {
    if (typeof(processMethod) != "character") {
        msg<- "'processMethod' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if (processMethod != "REGEX" && processMethod != "CSV") {
            msg <- "'processMethod' must be 'REGEX' or 'CSV'"
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkContigName <- function(contigName, errors, errorTypes) {
    if (is.null(contigName)) {
        msg<- "'contigName' must not be NULL. 'contigName' is missing."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
    }
    return(list(errors, errorTypes))
}

checkGeneticCode <- function(geneticCode, errors, errorTypes) {
    if (typeof(geneticCode) != "character") {
        msg<- "'geneticCode' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if(!("*" %in% geneticCode)) {
            msg <- "'geneticCode' does not specify any stop codons."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkRefAAS <- function(refAminoAcidSeq, errors, errorTypes) {
    if (typeof(refAminoAcidSeq) != "character") {
        msg<- "'refAminoAcidSeq' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
    }
    return(list(errors, errorTypes))
}

### ============================================================================
### ConsensusRead related: 'minReadsNum', 'minReadLength', 'minFractionCall'
###                        'maxFractionLost' prechecking
### ============================================================================
checkMinReadsNum <- function(minReadsNum, errors, errorTypes) {
    if (!is.numeric(minReadsNum)) {
        msg <- "'minReadsNum' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (minReadsNum%%1!=0) {
            msg <- "'minReadsNum' must be integer."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        }
        if (minReadsNum == 0) {
            msg <- "'minReadsNum' cannot be zero."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkMinReadLength <- function(minReadLength, errors, errorTypes) {
    if (!is.numeric(minReadLength)) {
        msg <- "'minReadLength' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (minReadLength%%1!=0) {
            msg <- "'minReadLength' must be integer."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkMinFractionCall <- function(minFractionCall, errors, errorTypes) {
    if (!is.numeric(minFractionCall)) {
        msg <- "'minFractionCall' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (minFractionCall > 1 || minFractionCall < 0) {
            msg <- "'minFractionCall' must be between 0 and 1."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
        } 
    }
    return(list(errors, errorTypes))
}

checkMaxFractionLost <- function(maxFractionLost, errors, errorTypes) {
    if (!is.numeric(maxFractionLost)) {
        msg <- "'maxFractionLost' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (maxFractionLost > 1 || maxFractionLost < 0) {
            msg <- "'maxFractionLost' must be between 0 and 1."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkAcceptStopCodons <- function(acceptStopCodons, errors, errorTypes) {
    if (!is.logical(acceptStopCodons)) {
        msg <- "'acceptStopCodons' must be 'TRUE' or 'FALSE'"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    }
    return(list(errors, errorTypes))
}

checkReadingFrame <- function(readingFrame, errors, errorTypes) {
    if (!is.numeric(readingFrame)) {
        msg <- "'readingFrame' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if(!readingFrame %in% c(1,2,3)) {
            msg <- "'readingFrame' must be 1, 2, or 3."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkProcessorsNum <- function(processorsNum, errors, errorTypes) {
    if (!is.numeric(processorsNum)) {
        msg <- "'processorsNum' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (!(processorsNum %% 1 == 0) && !is.null(processorsNum)) {
            msg <- "'processorsNum' must be integer."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

### ============================================================================
### 'ABIF_Directory' prechecking
### ============================================================================
checkABIF_Directory <- function(ABIF_Directory, errors, errorTypes) {
    if (typeof(ABIF_Directory) != "character") {
        msg<- "'ABIF_Directory' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else if (is.null(ABIF_Directory)) {
        msg<- "'ABIF_Directory' cannot be NULL."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if (!dir.exists(ABIF_Directory)) {
            msg <- paste("'", ABIF_Directory, "'",
                         " parent directory does not exist.", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "DIRECTORY_NOT_EXIST_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

### ============================================================================
### Quality trimming related: 'TrimmingMethod', 'M1TrimmingCutoff',
###                           'M2CutoffQualityScore', 'M2SlidingWindowSize'
### ============================================================================
checkTrimParam <- function(TrimmingMethod, M1TrimmingCutoff,
                           M2CutoffQualityScore, M2SlidingWindowSize, 
                           errors, errorTypes) {
    if (TrimmingMethod == "M1") {
        if (!is.numeric(M1TrimmingCutoff)) {
            msg<- "'M1TrimmingCutoff' must be numeric (You choose M1)."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        } else {
            # Ristriction about M1TrimmingCutoff !
            if (M1TrimmingCutoff > 1 || M1TrimmingCutoff < 0) {
                msg <- paste("Your input M1TrimmingCutoff is: '",
                             M1TrimmingCutoff, "' is invalid.",
                             "'M1TrimmingCutoff' should",
                             "be between 0 and 1.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
            }
        }
    } else if (TrimmingMethod == "M2") {
        if (!is.numeric(M2CutoffQualityScore)) {
            msg<- "'M2CutoffQualityScore' must be numeric (You choose M2)."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        } else {
            if (M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
                M2CutoffQualityScore%%1!=0) {
                msg <- paste("Your input M2CutoffQualityScore is: '",
                             M2CutoffQualityScore, "' is invalid.",
                             "'M2CutoffQualityScore' should",
                             "be between 0 and 60.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
            }
        }
        if (!is.numeric(M2SlidingWindowSize)) {
            msg<- "'M2SlidingWindowSize' must be numeric (You choose M2)."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        } else {
            if (M2SlidingWindowSize > 40 || M2SlidingWindowSize < 0 ||
                M2SlidingWindowSize%%1!=0) {
                msg <- paste("Your input M2SlidingWindowSize is: '",
                             M2SlidingWindowSize, "' is invalid.",
                             "'M2SlidingWindowSize' should",
                             "be between 0 and 40.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
            }
        }
    } else {
        msg <- "'TrimmingMethod' must be 'M1' or 'M2'."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    }
    return(list(errors, errorTypes))
}

### ============================================================================
### 'baseNumPerRow', 'heightPerRow', signalRatioCutoff', 'showTrimmed' prechecking
### ============================================================================
checkBaseNumPerRow <- function(baseNumPerRow, errors, errorTypes) {
    if (!is.numeric(baseNumPerRow)) {
        msg <- "'baseNumPerRow' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (baseNumPerRow%%1!=0) {
            msg <- "'baseNumPerRow' must be integer."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        }
        if (baseNumPerRow < 0 || baseNumPerRow > 200) {
            msg <- "'baseNumPerRow' must be between 0 and 200."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkHeightPerRow <- function(heightPerRow, errors, errorTypes) {
    if (!is.numeric(heightPerRow)) {
        msg <- "'heightPerRow' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (heightPerRow%%1!=0) {
            msg <- "'heightPerRow' must be integer."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
        }
        if (heightPerRow < 50 || heightPerRow > 600) {
            msg <- "'heightPerRow' must be between 50 and 600."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

### ============================================================================
### MakeBaseCalls Utilities function
### ============================================================================
checkSignalRatioCutoff <- function(signalRatioCutoff, errors, errorTypes) {
    if (!is.numeric(signalRatioCutoff)) {
        msg <- "'signalRatioCutoff' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (signalRatioCutoff < 0 || signalRatioCutoff > 1) {
            msg <- "'signalRatioCutoff' must be between 0 and 1."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_RANGE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkShowTrimmed <- function(showTrimmed, errors, errorTypes) {
    if (!is.logical(showTrimmed)) {
        msg <- "'showTrimmed' must be between TRUE and FALSE."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    }
    return(list(errors, errorTypes))
}

checkFASTA_File <- function(inputSource, FASTA_File, errors, errorTypes) {
    if (typeof(FASTA_File) != "character") {
        msg<- "'FASTA_File' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
    } else if (is.null(FASTA_File)) {
        msg<- "'FASTA_File' cannot be NULL."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if (!file.exists(FASTA_File)) {
            cat ("FASTA_File", FASTA_File)
            msg <- paste("'", FASTA_File, "'",
                         " file does not exist.", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "FILE_NOT_EXIST_ERROR")
        }
        if (is.na(str_extract(basename(FASTA_File), ".fa$")) &&
            is.na(str_extract(basename(FASTA_File), ".fasta$"))) {
            msg <- paste("'", FASTA_File, "'",
                         " file extension must be '.fa' or '.fasta'.",
                         sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "FILE_TYPE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkREGEX_SuffixForward <- function(REGEX_SuffixForward, errors, errorTypes){
    if (typeof(REGEX_SuffixForward) != "character") {
        msg<- "'REGEX_SuffixForward' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
    } else if (is.null(REGEX_SuffixForward)) {
        msg<- "'REGEX_SuffixForward' cannot be NULL."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    }
    return(list(errors, errorTypes))
}

checkREGEX_SuffixReverse <- function(REGEX_SuffixReverse, errors, errorTypes) {
    if (typeof(REGEX_SuffixReverse) != "character") {
        msg<- "'REGEX_SuffixReverse' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
    } else if (is.null(REGEX_SuffixReverse)) {
        msg<- "'REGEX_SuffixReverse' cannot be NULL."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    }
    return(list(errors, errorTypes))
}

checkCSV_NamesConversion <- function(CSV_NamesConversion, errors, errorTypes) {
    if (typeof(CSV_NamesConversion) != "character") {
        msg<- "'CSV_NamesConversion' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
    } else if (is.null(CSV_NamesConversion)) {
        msg<- "'CSV_NamesConversion' cannot be NULL."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if (!file.exists(CSV_NamesConversion)) {
            msg <- paste("'", CSV_NamesConversion, "'",
                         " file does not exist.", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "CSV_FILE_NOT_EXIST_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkReadFileNameExist <- function(readFileName, errors, errorTypes) {
    if (!file.exists(readFileName)) {
        msg <- paste("'", readFileName, "'",
                     " file does not exist.", sep = "")
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "FILE_NOT_EXIST_ERROR")
    }
    return(list(errors, errorTypes))
}

checkReadFileName <- function(readFileName, inputSource, errors, errorTypes) {
    if (inputSource == "ABIF") {
        if (is.na(str_extract(basename(readFileName), ".ab1$"))) {
            msg <- paste("'", readFileName, "'",
                         " file extension must be '.ab1'.", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "FILE_TYPE_ERROR")
        }
    } else if (inputSource == "FASTA") {
        if (is.na(str_extract(basename(readFileName), ".fa$")) &&
            is.na(str_extract(basename(readFileName), ".fasta$"))) {
            msg <- paste("'", readFileName, "'",
                         " file extension must be '.fa' or '.fasta'.",
                         sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "FILE_TYPE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkTargetFastaName <- function(targetFastaName, fastaReadName, 
                                 readFileName, errors, errorTypes) {
    if(isEmpty(targetFastaName)) {
        msg <- paste0("The name '", fastaReadName, 
                      "' is not in the '", 
                      basename(readFileName),
                      "' FASTA file")
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "FASTA_NAME_NOT_EXIST")
    }
    return(list(errors, errorTypes))
}

checkReadFeature <- function(readFeature, errors, errorTypes) {
    if (typeof(readFeature) != "character") {
        msg<- "'readFeature' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if (readFeature != "Forward Read" && readFeature != "Reverse Read") {
            msg <- "'readFeature' must be 'Forward Read' or 'Reverse Read'"
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

checkQualityPhredScores <- function(qualityPhredScores, errors, errorTypes) {
    if (length(qualityPhredScores) == 0) {
        msg <- paste("'qualityPhredScores' length cannot be zero.")
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    }
    if (!all(qualityPhredScores%%1 == 0)) {
        msg <- "All elements in 'qualityPhredScores' vector must be integer."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    }
    return(list(errors, errorTypes))
}

checkAb1FastaCsv <- function(ABIF_Directory, FASTA_File,
                             CSV_NamesConversion, inputSource, 
                             errors, errorTypes) {
    if (!file.exists(CSV_NamesConversion)) {
        msg <- paste("CSV_NamesConversion: '",  
                     CSV_NamesConversion, "'",
                     " file does not exist.", sep = "")
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "FILE_NOT_EXIST_ERROR")
    } else {
        warnings <- character()
        if (inputSource == "ABIF") {
            csvFile <- read.csv(CSV_NamesConversion, header = TRUE)
            csvReads <- as.character(csvFile$reads)
            parentDirFiles <- list.files(ABIF_Directory)
            sourceReads <- parentDirFiles[grepl("\\.ab1$", parentDirFiles)]
        } else if (inputSource == "FASTA") {
            csvFile <- read.csv(CSV_NamesConversion, header = TRUE)
            csvReads <- as.character(csvFile$reads)
            readFasta <- read.fasta(FASTA_File, as.string = TRUE)
            sourceReads <- names(readFasta)
        }
        ########################################################################
        ### Check 1: all reads in the read folder are listed in the csv
        ########################################################################
        readInCsvWarningMsg <-
            lapply(sourceReads,
                   function(sourceRead) {
                       if (!(sourceRead %in% csvReads)) {
                           msg <- paste("'", sourceRead, "' is not in the ",
                                        "csv file (", CSV_NamesConversion, ")",
                                        sep = "")
                           return(msg)}
                       return()})
        warnings <- c(warnings, unlist(readInCsvWarningMsg), use.names = FALSE)
        
        ########################################################################
        ### Check 2: all reads listed in the csv file are in the reads folder
        ########################################################################
        readInSourceWarningMsg <-
            lapply(csvReads,
                   function(csvRead) {
                       if (!(csvRead %in% sourceReads)) {
                           msg <- paste("'", csvRead, "' is not in the ",
                                        "parent directory.", sep = "")
                           return(msg)}
                       return()})
        warnings <- c(warnings, 
                      unlist(readInSourceWarningMsg), use.names = FALSE)
        
        # csv file has all of the columns, no extra columns
        ########################################################################
        ### Check 3: csv file has all of the columns, no extra columns
        ########################################################################
        if (!("contig" %in% colnames(csvFile))) {
            msg <- paste("'contig' is not in the csv file (",
                         CSV_NamesConversion, ")", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "CSV_MISMATCH_ERROR")
        } else if (!("direction" %in% colnames(csvFile))) {
            msg <- paste("'direction' is not in the csv file (",
                         CSV_NamesConversion, ")", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "CSV_MISMATCH_ERROR")
        } else if (!("reads" %in% colnames(csvFile))) {
            msg <- paste("'reads' is not in the csv file (",
                         CSV_NamesConversion, ")", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "CSV_MISMATCH_ERROR")
        }
        
        ########################################################################
        ### Check 4: F/R column has only F's and R's.
        ########################################################################
        csvFileDirLen <- length(unique(as.character(csvFile$direction)))
        if (!((csvFileDirLen == 1 || csvFileDirLen == 2) && 
              ('F' %in% unique(as.character(csvFile$direction)) || 
               'R' %in% unique(as.character(csvFile$direction))))) {
            msg <- paste("In the 'direction' column of your CSV file, 
                     you can only have 'F' and 'R' ", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "CSV_VALUE_ERROR")
        }
        if (length(warnings) != 0) {
            invisible(lapply(warnings, log_warn))
        }
    }
    return(list(errors, errorTypes))
}

checkGreplForward <- function(forwardSelectInputFiles, warnings, warningsType) {
    if (length(forwardSelectInputFiles) == 0) {
        msg <- paste0("Your 'contigName' and 'REGEX_SuffixForward' regular ", 
                      "expression parameters can not match any forward reads.")
        warnings <- c(warnings, msg)
        warningsType <- c(warningsType, "REGEX_MATCH_WARN")   
    }
    return(list(warnings, warningsType))
}

checkGreplReverse <- function(reverseSelectInputFiles, warnings, warningsType) {
    if (length(reverseSelectInputFiles) == 0) {
        msg <- paste0("Your 'contigName' and 'REGEX_SuffixReverse' regular ", 
                      "expression parameters can not match any reverse reads.")
        warnings <- c(warnings, msg)
        warningsType <- c(warningsType, "REGEX_MATCH_WARN")   
    }
    return(list(warnings, warningsType))
}

checkCSVConvForward <- function(forwardReads, warnings, warningsType) {
    if (length(forwardReads) == 0) {
        msg <- paste0("The names of forward reads in your 'CSV_NamesConversion'",
                      " CSV file do not match any forward reads.")
        warnings <- c(warnings, msg)
        warningsType <- c(warningsType, "CSV_MATCH_WARN")   
    }
    return(list(warnings, warningsType))
}

checkCSVConvReverse <- function(reverseReads, warnings, warningsType) {
    if (length(reverseReads) == 0) {
        msg <- paste0("The names of reverse reads in your 'CSV_NamesConversion'", 
                      " CSV file do not match any reverse reads.")
        warnings <- c(warnings, msg)
        warningsType <- c(warningsType, "CSV_MATCH_WARN")   
    }
    return(list(warnings, warningsType))
}
