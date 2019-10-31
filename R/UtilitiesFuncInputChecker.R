### ============================================================================
### Quality related: 'cutoffQualityScore' & 'slidingWindowSize' prechecking
### ============================================================================
checkCutoffQualityScore <- function(cutoffQualityScore, errors) {
    if (cutoffQualityScore > 60 || cutoffQualityScore < 0 ||
        cutoffQualityScore%%1!=0) {
        msg <- paste("\nYour input 'cutoffQualityScore' is: ",
                     cutoffQualityScore, "is invalid.",
                     "'cutoffQualityScore' should",
                     "be between 0 and 60.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkSlidingWindowSize <- function(slidingWindowSize, errors) {
    if (slidingWindowSize > 20 || slidingWindowSize < 0 ||
        slidingWindowSize%%1!=0) {
        msg <- paste("\nYour input 'slidingWindowSize' is: ",
                     slidingWindowSize, "is invalid.",
                     "'slidingWindowSize' should",
                     "be between 0 and 20.\n", sep = "")
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
        msg <- paste("\n''minReadsNum' have to be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMinReadLength <- function(minReadLength, errors) {
    if (minReadLength%%1!=0) {
        msg <- paste("\n'minReadLength' have to be integer.\n", sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}



checkMinFractionCall <- function(minFractionCall, errors) {
    if (minFractionCall > 1 || minFractionCall < 0) {
        msg <- paste("\n'minFractionCall' have to be between 0 and 1.\n",
                     sep = "")
        errors <- c(errors, msg)
    }
    return(errors)
}

checkMaxFractionLost <- function(maxFractionLost, errors) {
    if (maxFractionLost > 1 || maxFractionLost < 0) {
        msg <- paste("\n'maxFractionLost' have to be between 0 and 1.\n",
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

# Waiting list ~~
# refAminoAcidSeq
# geneticCode
# acceptStopCodons
# readingFrame
# processorsNum

# Origin value
# cutoffQualityScore     = 20,
# slidingWindowSize      = 5,
# refAminoAcidSeq        = "",
#
# minReadsNum            = 2,
# minReadLength          = 20,
#
# minFractionCall        = 0.5,
# maxFractionLost        = 0.5,
#
# geneticCode            = GENETIC_CODE,
# acceptStopCodons       = TRUE,
#
# readingFrame           = 1,
# processorsNum          = 1
