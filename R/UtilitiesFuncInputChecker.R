checkReadFileName <- function(readFileName, errors) {
    if (!file.exists(readFileName)) {
        cat ("readFileName", readFileName)
        msg <- paste("\n'", readFileName, "'",
                     " foward read file does not exist.\n", sep = "")
        errors <- c(errors, msg)
    }
}

### ============================================================================
### QualityReport related: 'readFeature', 'qualityPhredScores'
### ============================================================================
checkReadFeature <- function(readFeature, errors) {
    if (readFeature != "Forward Read" && readFeature != "Reverse Read") {
        msg <- paste("\n'readFeature' must be
                     'Forward Read' or 'Reverse Read'\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkQualityPhredScores <- function(qualityPhredScores, errors) {
    if (length(qualityPhredScores) == 0) {
        msg <- paste("\n'qualityPhredScores'
                               length cannot be zero.\n")
        errors <- c(errors, msg)
    }
    if (!all(qualityPhredScores%%1 == 0)) {
        msg <- paste("\nAll elements in 'qualityPhredScores' vector
                               must be integer.\n")
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
            msg<- paste("\n'M1TrimmingCutoff' must be numeric",
                        "(You choose M1).\n")
            errors <- c(errors, msg)
        } else {
            # Ristriction about M1TrimmingCutoff !
            # if (M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
            #     M2CutoffQualityScore%%1!=0) {
            #     msg <- paste("\n'Your input M2CutoffQualityScore is: ",
            #                  M2CutoffQualityScore, "' is invalid.",
            #                  "'M2CutoffQualityScore' should",
            #                  "be between 0 and 60.\n", sep = "")
            #     errors <- c(errors, msg)
            # }
        }
        if (!is.null(M2CutoffQualityScore)) {
            msg<- paste("\n'M2CutoffQualityScore' must be null",
                        "(You choose M1).\n")
            errors <- c(errors, msg)
        }
        if (!is.null(M2SlidingWindowSize)) {
            msg<- paste("\n'M2SlidingWindowSize' must be null",
                        "(You choose M1).\n")
            errors <- c(errors, msg)
        }
    } else if (TrimmingMethod == "M2") {
        if (!is.null(M1TrimmingCutoff)) {
            msg<- paste("\n'M1TrimmingCutoff' must be null",
                        "(You choose M2).\n")
            errors <- c(errors, msg)
        }
        if (!is.numeric(M2CutoffQualityScore)) {
            msg<- paste("\n'M2CutoffQualityScore' must be numeric",
                        "(You choose M2).\n")
            errors <- c(errors, msg)
        } else {
            if (M2CutoffQualityScore > 60 || M2CutoffQualityScore < 0 ||
                M2CutoffQualityScore%%1!=0) {
                msg <- paste("\n'Your input M2CutoffQualityScore is: ",
                             M2CutoffQualityScore, "' is invalid.",
                             "'M2CutoffQualityScore' should",
                             "be between 0 and 60.\n", sep = "")
                errors <- c(errors, msg)
            }
        }
        if (!is.numeric(M2SlidingWindowSize)) {
            msg<- paste("\n'M2SlidingWindowSize' must be numeric",
                        "(You choose M2).\n")
            errors <- c(errors, msg)
        } else {
            if (M2SlidingWindowSize > 20 || M2SlidingWindowSize < 0 ||
                M2SlidingWindowSize%%1!=0) {
                msg <- paste("\n'Your input M2SlidingWindowSize is: ",
                             M2SlidingWindowSize, "' is invalid.",
                             "'M2SlidingWindowSize' should",
                             "be between 0 and 20.\n", sep = "")
                errors <- c(errors, msg)
            }
        }
    } else {
        msg <- paste("\n'TrimmingMethod' must be 'M1' or 'M2'.\n")
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
        msg <- paste("\n'minReadsNum' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMinReadLength <- function(minReadLength, errors) {
    if (minReadLength%%1!=0) {
        msg <- paste("\n'minReadLength' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}



checkMinFractionCall <- function(minFractionCall, errors) {
    if (minFractionCall > 1 || minFractionCall < 0) {
        msg <- paste("\n'minFractionCall' must be between 0 and 1.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMaxFractionLost <- function(maxFractionLost, errors) {
    if (maxFractionLost > 1 || maxFractionLost < 0) {
        msg <- paste("\n'maxFractionLost' must be between 0 and 1.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}


checkGeneticCode <- function(geneticCode, errors) {
    if(!("*" %in% geneticCode)) {
        msg <- paste("\n'geneticCode' does not specify any stop codons.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}


checkReadingFrame <- function(readingFrame, errors) {
    if(!readingFrame %in% c(1,2,3)) {
        msg <- paste("\n'readingFrame' must be 1, 2, or 3.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkAcceptStopCodons <- function(acceptStopCodons, errors) {
    if (!is.logical(acceptStopCodons)) {
        msg <- paste("\n'acceptStopCodons' must be 'TRUE' or 'FALSE'\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkProcessorsNum <- function(processorsNum, errors) {
    if (!(processorsNum %% 1 == 0)) {
        msg <- paste("\n'processorsNum' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### 'parentDirectory' prechecking
### ============================================================================
checkParentDirectory <- function(parentDirectory, errors) {
    if (!file.exists(parentDirectory)) {
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
        msg <- paste("\n'baseNumPerRow' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (baseNumPerRow < 0 || baseNumPerRow > 200) {
        msg <- paste("\n'baseNumPerRow' must be between 0 and 200.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkHeightPerRow <- function(heightPerRow, errors) {
    if (heightPerRow%%1!=0) {
        msg <- paste("\n'heightPerRow' must be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    if (heightPerRow < 50 || heightPerRow > 600) {
        msg <- paste("\n'heightPerRow' must be between 0 and 200.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

### ============================================================================
### MakeBaseCalls Utilities function
### ============================================================================
checkSignalRatioCutoff <- function(signalRatioCutoff, errors) {
    if (signalRatioCutoff < 0 || signalRatioCutoff > 1) {
        msg <- paste("\n'signalRatioCutoff' must be between 0 and 1.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkShowTrimmed <- function(showTrimmed, errors) {
    if (!is.logical(showTrimmed)) {
        msg <- paste("\n'showTrimmed' must be between TRUE and FALSE.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}



