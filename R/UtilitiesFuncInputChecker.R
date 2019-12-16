checkFastaFileName <- function(fastaFileName, errors) {
    if (!file.exists(fastaFileName)) {
        cat ("fastaFileName", fastaFileName)
        msg <- paste("\n'", fastaFileName, "'",
                     " file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (is.na(str_extract(basename(fastaFileName), ".fa$")) &&
        is.na(str_extract(basename(fastaFileName), ".fasta$"))) {
        msg <- paste("\n'", fastaFileName, "'",
                     " file extension must be '.fa' or '.fasta'.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkReadFileName <- function(readFileName, inputSource, errors) {
    if (!file.exists(readFileName)) {
        cat ("readFileName", readFileName)
        msg <- paste("\n'", readFileName, "'",
                     " file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (inputSource == "ABIF") {
        if (is.na(str_extract(basename(readFileName), ".ab1$"))) {
            msg <- paste("\n'", readFileName, "'",
                         " file extension must be '.ab1'.\n", sep = "")
            errors <- c(errors, msg)
        }
    } else if (inputSource == "FASTA") {
        if (is.na(str_extract(basename(readFileName), ".fa$")) &&
            is.na(str_extract(basename(readFileName), ".fasta$"))) {
            msg <- paste("\n'", readFileName, "'",
                         " file extension must be '.fa' or '.fasta'.\n",
                         sep = "")
            errors <- c(errors, msg)
        }
    }
    return(errors)
}

### ============================================================================
### QualityReport related: 'readFeature', 'qualityPhredScores'
### ============================================================================
checkInputSource <- function(inputSource, errors) {
    if (inputSource != "ABIF" && inputSource != "FASTA") {
        msg <- "\n'inputSource' must be 'ABIF' or 'FASTA'\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkReadFeature <- function(readFeature, errors) {
    if (readFeature != "Forward Read" && readFeature != "Reverse Read") {
        msg <- "\n'readFeature' must be 'Forward Read' or 'Reverse Read'\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkQualityPhredScores <- function(qualityPhredScores, errors) {
    if (length(qualityPhredScores) == 0) {
        msg <- paste("\n'qualityPhredScores' length cannot be zero.\n")
        errors <- c(errors, msg)
    }
    if (!all(qualityPhredScores%%1 == 0)) {
        msg <- "\nAll elements in 'qualityPhredScores' vector must be integer.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### Quality trimming related: 'TrimmingMethod', 'M1TrimmingCutoff',
###                           'M2CutoffQualityScore', 'M2SlidingWindowSize'
### ============================================================================
checkTrimParam <- function(TrimmingMethod, M1TrimmingCutoff,
                           M2CutoffQualityScore, M2SlidingWindowSize, errors) {
    if (TrimmingMethod == "M1") {
        if (!is.numeric(M1TrimmingCutoff)) {
            msg<- "\n'M1TrimmingCutoff' must be numeric (You choose M1).\n"
            errors <- c(errors, msg)
        } else {
            # Ristriction about M1TrimmingCutoff !
            if (M1TrimmingCutoff > 1 || M1TrimmingCutoff < 0) {
                msg <- paste("\nYour input M1TrimmingCutoff is: '",
                             M1TrimmingCutoff, "' is invalid.",
                             "'M1TrimmingCutoff' should",
                             "be between 0 and 1.\n", sep = "")
                errors <- c(errors, msg)
            }
        }
        if (!is.null(M2CutoffQualityScore)) {
            msg<- "\n'M2CutoffQualityScore' must be null (You choose M1).\n"
            errors <- c(errors, msg)
        }
        if (!is.null(M2SlidingWindowSize)) {
            msg<- "\n'M2SlidingWindowSize' must be null (You choose M1).\n"
            errors <- c(errors, msg)
        }
    } else if (TrimmingMethod == "M2") {
        if (!is.null(M1TrimmingCutoff)) {
            msg<- "\n'M1TrimmingCutoff' must be null (You choose M2).\n"
            errors <- c(errors, msg)
        }
        if (!is.numeric(M2CutoffQualityScore)) {
            msg<- "\n'M2CutoffQualityScore' must be numeric (You choose M2).\n"
            errors <- c(errors, msg)
        } else {
            if (M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
                M2CutoffQualityScore%%1!=0) {
                msg <- paste("\nYour input M2CutoffQualityScore is: '",
                             M2CutoffQualityScore, "' is invalid.",
                             "'M2CutoffQualityScore' should",
                             "be between 0 and 60.\n", sep = "")
                errors <- c(errors, msg)
            }
        }
        if (!is.numeric(M2SlidingWindowSize)) {
            msg<- "\n'M2SlidingWindowSize' must be numeric (You choose M2).\n"
            errors <- c(errors, msg)
        } else {
            if (M2SlidingWindowSize > 40 || M2SlidingWindowSize < 0 ||
                M2SlidingWindowSize%%1!=0) {
                msg <- paste("\nYour input M2SlidingWindowSize is: '",
                             M2SlidingWindowSize, "' is invalid.",
                             "'M2SlidingWindowSize' should",
                             "be between 0 and 40.\n", sep = "")
                errors <- c(errors, msg)
            }
        }
    } else {
        msg <- "\n'TrimmingMethod' must be 'M1' or 'M2'.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### ConsensusRead related: 'minReadsNum', 'minReadLength', 'minFractionCall'
###                        'maxFractionLost' prechecking
### ============================================================================
checkMinReadsNum <- function(minReadsNum, errors) {
    if (minReadsNum%%1!=0) {
        msg <- "\n'minReadsNum' must be integer.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMinReadLength <- function(minReadLength, errors) {
    if (minReadLength%%1!=0) {
        msg <- "\n'minReadLength' must be integer.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}



checkMinFractionCall <- function(minFractionCall, errors) {
    if (minFractionCall > 1 || minFractionCall < 0) {
        msg <- "\n'minFractionCall' must be between 0 and 1.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMaxFractionLost <- function(maxFractionLost, errors) {
    if (maxFractionLost > 1 || maxFractionLost < 0) {
        msg <- "\n'maxFractionLost' must be between 0 and 1.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}


checkGeneticCode <- function(geneticCode, errors) {
    if(!("*" %in% geneticCode)) {
        msg <- "\n'geneticCode' does not specify any stop codons.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}


checkReadingFrame <- function(readingFrame, errors) {
    if(!readingFrame %in% c(1,2,3)) {
        msg <- "\n'readingFrame' must be 1, 2, or 3.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkAcceptStopCodons <- function(acceptStopCodons, errors) {
    if (!is.logical(acceptStopCodons)) {
        msg <- "\n'acceptStopCodons' must be 'TRUE' or 'FALSE'\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkProcessorsNum <- function(processorsNum, errors) {
    if (!(processorsNum %% 1 == 0) && !is.null(processorsNum)) {
        msg <- "\n'processorsNum' must be integer.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### 'parentDirectory' prechecking
### ============================================================================
checkParentDirectory <- function(parentDirectory, errors) {
    if (!dir.exists(parentDirectory)) {
        msg <- paste("\n'", parentDirectory, "'",
                     " parent directory does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
}

### ============================================================================
### 'baseNumPerRow', 'signalRatioCutoff', 'showTrimmed' prechecking
### ============================================================================
checkBaseNumPerRow <- function(baseNumPerRow, errors) {
    if (baseNumPerRow%%1!=0) {
        msg <- "\n'baseNumPerRow' must be integer.\n"
        errors <- c(errors, msg)
    }
    if (baseNumPerRow < 0 || baseNumPerRow > 200) {
        msg <- "\n'baseNumPerRow' must be between 0 and 200.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkHeightPerRow <- function(heightPerRow, errors) {
    if (heightPerRow%%1!=0) {
        msg <- "\n'heightPerRow' must be integer.\n"
        errors <- c(errors, msg)
    }
    if (heightPerRow < 50 || heightPerRow > 600) {
        msg <- "\n'heightPerRow' must be between 0 and 200.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### MakeBaseCalls Utilities function
### ============================================================================
checkSignalRatioCutoff <- function(signalRatioCutoff, errors) {
    if (signalRatioCutoff < 0 || signalRatioCutoff > 1) {
        msg <- "\n'signalRatioCutoff' must be between 0 and 1.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

checkShowTrimmed <- function(showTrimmed, errors) {
    if (!is.logical(showTrimmed)) {
        msg <- "\n'showTrimmed' must be between TRUE and FALSE.\n"
        errors <- c(errors, msg)
    }
    return(errors)
}

