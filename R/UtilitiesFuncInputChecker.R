checkReadFileNameExist <- function(readFileName, errors, errorTypes) {
    if (!file.exists(readFileName)) {
        msg <- paste("'", readFileName, "'",
                     " file does not exist.", sep = "")
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "FILE_NOT_EXISTS_ERROR")
    }
    return(list(errors, errorTypes))
}

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

checkFastaFileName <- function(inputSource, fastaFileName, errors, errorTypes) {
    if (typeof(fastaFileName) != "character" && !is.null(fastaFileName)) {
        msg<- "'fastaFileName' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
    } else {
        if (inputSource == "FASTA") {
            if (!file.exists(fastaFileName)) {
                cat ("fastaFileName", fastaFileName)
                msg <- paste("'", fastaFileName, "'",
                             " file does not exist.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "FILE_NOT_EXISTS_ERROR")
            }
            if (is.na(str_extract(basename(fastaFileName), ".fa$")) &&
                is.na(str_extract(basename(fastaFileName), ".fasta$"))) {
                msg <- paste("'", fastaFileName, "'",
                             " file extension must be '.fa' or '.fasta'.",
                             sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "FILE_TYPE_ERROR")
            }
        } else if (inputSource == "ABIF") {
            if (!is.null(fastaFileName)) {
                msg<- "'fastaFileName' must be null (inputSource is 'ABIF')."
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
            }
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
        if (!is.null(M2CutoffQualityScore)) {
            msg<- "'M2CutoffQualityScore' must be null (You choose M1)."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
        if (!is.null(M2SlidingWindowSize)) {
            msg<- "'M2SlidingWindowSize' must be null (You choose M1)."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
    } else if (TrimmingMethod == "M2") {
        if (!is.null(M1TrimmingCutoff)) {
            msg<- "'M1TrimmingCutoff' must be null (You choose M2)."
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
        }
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

checkAcceptStopCodons <- function(acceptStopCodons, errors, errorTypes) {
    if (!is.logical(acceptStopCodons)) {
        msg <- "'acceptStopCodons' must be 'TRUE' or 'FALSE'"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
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
### 'parentDirectory' prechecking
### ============================================================================
checkParentDirectory <- function(parentDirectory, errors, errorTypes) {
    if (typeof(parentDirectory) != "character" && !is.null(parentDirectory)) {
        msg<- "'parentDirectory' must be character type."
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
    } else {
        if (!dir.exists(parentDirectory)) {
            msg <- paste("'", parentDirectory, "'",
                         " parent directory does not exist.", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "DIRECTORY_NOT_EXISTS_ERROR")
        }
    }
    return(list(errors, errorTypes))
}

### ============================================================================
### 'baseNumPerRow', 'signalRatioCutoff', 'showTrimmed' prechecking
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


checkNamesConversionCSV <- function(nameConvMethod, namesConversionCSV, 
                                    contigName, suffixForwardRegExp,
                                    suffixReverseRegExp, inputSource, 
                                    errors, errorTypes) {
    if (typeof(nameConvMethod) != "character" && !is.null(nameConvMethod)) {
        msg <- "'nameConvMethod' must be numeric"
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_TYPE_ERROR")
    } else {
        if (nameConvMethod == "REGEX") {
            # CSV file needs to be NULL
            if (!is.null(namesConversionCSV)) {
                msg <- paste("namesConversionCSV: '", namesConversionCSV, "'",
                             " needs to be null.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
            }
            # contigName, suffixForwardRegExp, suffixReverseRegExp cannot be NULL
            if (is.null(contigName)) {
                msg <- paste("contigName cannot be null.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
            }
            if (is.null(suffixForwardRegExp)) {
                msg <- paste("suffixForwardRegExp cannot be null.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
            }
            if (is.null(suffixReverseRegExp)) {
                msg <- paste("suffixReverseRegExp cannot be null.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")   
            }
        } else if (nameConvMethod == "CSV") {
            if(is.null(namesConversionCSV)) {
                msg <- paste("namesConversionCSV cannot be null.", sep = "")
                errors <- c(errors, msg)
                errorTypes <- c(errorTypes, "PARAMETER_VALUE_ERROR")
            }
        }
    }
    return(list(errors, errorTypes))
}

checkAb1FastaCsv <- function(parentDirectory, fastaFileName,
                             namesConversionCSV, inputSource, 
                             errors, errorTypes) {
    if (!file.exists(namesConversionCSV)) {
        msg <- paste("namesConversionCSV: '",  
                     namesConversionCSV, "'",
                     " file does not exist.", sep = "")
        errors <- c(errors, msg)
        errorTypes <- c(errorTypes, "FILE_NOT_EXISTS_ERROR")
    } else {
        warnings <- character()
        if (inputSource == "ABIF") {
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            csvReads <- as.character(csvFile$reads)
            parentDirFiles <- list.files(parentDirectory)
            sourceReads <- parentDirFiles[grepl("\\.ab1$", parentDirFiles)]
        } else if (inputSource == "FASTA") {
            csvFile <- read.csv(namesConversionCSV, header = TRUE)
            csvReads <- as.character(csvFile$reads)
            readFasta <- read.fasta(fastaFileName, as.string = TRUE)
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
                                        "csv file (", namesConversionCSV, ")",
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
                         namesConversionCSV, ")", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "CSV_MISMATCH_ERROR")
        } else if (!("direction" %in% colnames(csvFile))) {
            msg <- paste("'direction' is not in the csv file (",
                         namesConversionCSV, ")", sep = "")
            errors <- c(errors, msg)
            errorTypes <- c(errorTypes, "CSV_MISMATCH_ERROR")
        } else if (!("reads" %in% colnames(csvFile))) {
            msg <- paste("'reads' is not in the csv file (",
                         namesConversionCSV, ")", sep = "")
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