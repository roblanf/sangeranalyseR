# =============================================================================
# Input validators
#
# Each public `check*` function accepts the parallel (errors, errorTypes)
# pair and returns `list(errors, errorTypes)` after appending any new
# diagnostics. This contract is preserved exactly — the public surface is
# unchanged.
#
# Internal dot-prefixed helpers consolidate the common patterns (enum,
# range, file/dir existence, extension regex). Each public validator is now
# a thin wrapper over these helpers.
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

# Append (msg, type) to the parallel diagnostic lists. Single concat point —
# replaces the per-validator `c()` calls.
.errAppend <- function(errors, errorTypes, msg, type) {
    list(c(errors, msg), c(errorTypes, type))
}

# Type predicate check.
.requireType <- function(value, name, predicate, expected,
                          errors, errorTypes,
                          type = "PARAMETER_TYPE_ERROR") {
    if (!predicate(value)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", name, "' must be ", expected, "."),
                          type))
    }
    list(errors, errorTypes)
}

# Enum check. value must be a single character ∈ allowed.
.requireEnum <- function(value, name, allowed, errors, errorTypes,
                          type = "PARAMETER_VALUE_ERROR") {
    if (typeof(value) != "character") {
        return(.errAppend(errors, errorTypes,
                          paste0("'", name, "' must be character type."),
                          type))
    }
    if (length(value) != 1L || !(value %in% allowed)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", name, "' must be one of: ",
                                 paste(sQuote(allowed), collapse = ", "),
                                 "."),
                          type))
    }
    list(errors, errorTypes)
}

# Numeric range check. lo / hi may be NULL to skip a side. integer = TRUE
# additionally requires the value to be a whole number.
.requireRange <- function(value, name, lo, hi, errors, errorTypes,
                           integer = FALSE) {
    if (!is.numeric(value)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", name, "' must be numeric"),
                          "PARAMETER_TYPE_ERROR"))
    }
    if (integer && (value %% 1 != 0)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", name, "' must be integer."),
                          "PARAMETER_TYPE_ERROR"))
    }
    if (!is.null(lo) && value < lo) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", name, "' must be between ", lo,
                                 " and ", hi, "."),
                          "PARAMETER_RANGE_ERROR"))
    }
    if (!is.null(hi) && value > hi) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", name, "' must be between ", lo,
                                 " and ", hi, "."),
                          "PARAMETER_RANGE_ERROR"))
    }
    list(errors, errorTypes)
}

# File extension regex check. `pattern` MUST be a valid regex with escaped
# dots — e.g. "\\.fa(sta)?$" or "\\.ab1$".
.requireExt <- function(path, pattern, errors, errorTypes,
                         message,
                         type = "FILE_TYPE_ERROR") {
    if (!grepl(pattern, basename(path), ignore.case = FALSE)) {
        return(.errAppend(errors, errorTypes, message, type))
    }
    list(errors, errorTypes)
}

# -----------------------------------------------------------------------------
# Public validators (preserved signatures)
# -----------------------------------------------------------------------------

checkInputSource <- function(inputSource, errors, errorTypes) {
    .requireEnum(inputSource, "inputSource", c("ABIF", "FASTA"),
                  errors, errorTypes)
}

checkProcessMethod <- function(inputSource, processMethod, errors, errorTypes) {
    .requireEnum(processMethod, "processMethod", c("REGEX", "CSV"),
                  errors, errorTypes)
}

checkContigName <- function(contigName, errors, errorTypes) {
    if (is.null(contigName)) {
        return(.errAppend(errors, errorTypes,
                          "'contigName' must not be NULL. 'contigName' is missing.",
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

checkGeneticCode <- function(geneticCode, errors, errorTypes) {
    if (typeof(geneticCode) != "character") {
        return(.errAppend(errors, errorTypes,
                          "'geneticCode' must be character type.",
                          "PARAMETER_VALUE_ERROR"))
    }
    if (!("*" %in% geneticCode)) {
        return(.errAppend(errors, errorTypes,
                          "'geneticCode' does not specify any stop codons.",
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

checkRefAAS <- function(refAminoAcidSeq, errors, errorTypes) {
    .requireType(refAminoAcidSeq, "refAminoAcidSeq",
                  function(x) typeof(x) == "character", "character type",
                  errors, errorTypes,
                  type = "PARAMETER_VALUE_ERROR")
}

# ----- ConsensusRead-related numerics ---------------------------------------

checkMinReadsNum <- function(minReadsNum, errors, errorTypes) {
    out <- .requireRange(minReadsNum, "minReadsNum", lo = NULL, hi = NULL,
                          errors, errorTypes, integer = TRUE)
    if (length(out[[1]]) == length(errors) && is.numeric(minReadsNum) &&
        minReadsNum == 0) {
        out <- .errAppend(out[[1]], out[[2]],
                           "'minReadsNum' cannot be zero.",
                           "PARAMETER_VALUE_ERROR")
    }
    out
}

checkMinReadLength <- function(minReadLength, errors, errorTypes) {
    .requireRange(minReadLength, "minReadLength", lo = NULL, hi = NULL,
                   errors, errorTypes, integer = TRUE)
}

checkMinFractionCall <- function(minFractionCall, errors, errorTypes) {
    .requireRange(minFractionCall, "minFractionCall", lo = 0, hi = 1,
                   errors, errorTypes)
}

checkMaxFractionLost <- function(maxFractionLost, errors, errorTypes) {
    .requireRange(maxFractionLost, "maxFractionLost", lo = 0, hi = 1,
                   errors, errorTypes)
}

checkAcceptStopCodons <- function(acceptStopCodons, errors, errorTypes) {
    if (!is.logical(acceptStopCodons)) {
        return(.errAppend(errors, errorTypes,
                          "'acceptStopCodons' must be 'TRUE' or 'FALSE'",
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

checkReadingFrame <- function(readingFrame, errors, errorTypes) {
    if (!is.numeric(readingFrame)) {
        return(.errAppend(errors, errorTypes,
                          "'readingFrame' must be numeric",
                          "PARAMETER_TYPE_ERROR"))
    }
    if (!(readingFrame %in% c(1, 2, 3))) {
        return(.errAppend(errors, errorTypes,
                          "'readingFrame' must be 1, 2, or 3.",
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

checkProcessorsNum <- function(processorsNum, errors, errorTypes) {
    if (is.null(processorsNum)) return(list(errors, errorTypes))
    if (!is.numeric(processorsNum)) {
        return(.errAppend(errors, errorTypes,
                          "'processorsNum' must be numeric",
                          "PARAMETER_TYPE_ERROR"))
    }
    if (processorsNum %% 1 != 0) {
        return(.errAppend(errors, errorTypes,
                          "'processorsNum' must be integer.",
                          "PARAMETER_TYPE_ERROR"))
    }
    list(errors, errorTypes)
}

# ----- Directory / file existence -------------------------------------------

checkABIF_Directory <- function(ABIF_Directory, errors, errorTypes) {
    if (is.null(ABIF_Directory)) {
        return(.errAppend(errors, errorTypes,
                          "'ABIF_Directory' cannot be NULL.",
                          "PARAMETER_VALUE_ERROR"))
    }
    if (typeof(ABIF_Directory) != "character") {
        return(.errAppend(errors, errorTypes,
                          "'ABIF_Directory' must be character type.",
                          "PARAMETER_VALUE_ERROR"))
    }
    if (!dir.exists(ABIF_Directory)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", ABIF_Directory, "'",
                                 " parent directory does not exist."),
                          "DIRECTORY_NOT_EXIST_ERROR"))
    }
    list(errors, errorTypes)
}

# ----- Trimming parameters --------------------------------------------------

checkTrimParam <- function(TrimmingMethod, M1TrimmingCutoff,
                           M2CutoffQualityScore, M2SlidingWindowSize,
                           errors, errorTypes) {
    if (TrimmingMethod == "M1") {
        out <- .requireRange(M1TrimmingCutoff, "M1TrimmingCutoff",
                              lo = 0, hi = 1, errors, errorTypes)
        return(out)
    }
    if (TrimmingMethod == "M2") {
        out <- .requireRange(M2CutoffQualityScore, "M2CutoffQualityScore",
                              lo = 0, hi = 60, errors, errorTypes,
                              integer = TRUE)
        out <- .requireRange(M2SlidingWindowSize, "M2SlidingWindowSize",
                              lo = 0, hi = 40, out[[1]], out[[2]],
                              integer = TRUE)
        return(out)
    }
    .errAppend(errors, errorTypes,
                "'TrimmingMethod' must be 'M1' or 'M2'.",
                "PARAMETER_VALUE_ERROR")
}

# ----- Chromatogram parameters ----------------------------------------------

checkBaseNumPerRow <- function(baseNumPerRow, errors, errorTypes) {
    .requireRange(baseNumPerRow, "baseNumPerRow", lo = 0, hi = 200,
                   errors, errorTypes, integer = TRUE)
}

checkHeightPerRow <- function(heightPerRow, errors, errorTypes) {
    .requireRange(heightPerRow, "heightPerRow", lo = 50, hi = 600,
                   errors, errorTypes, integer = TRUE)
}

checkSignalRatioCutoff <- function(signalRatioCutoff, errors, errorTypes) {
    .requireRange(signalRatioCutoff, "signalRatioCutoff", lo = 0, hi = 1,
                   errors, errorTypes)
}

checkShowTrimmed <- function(showTrimmed, errors, errorTypes) {
    if (!is.logical(showTrimmed)) {
        return(.errAppend(errors, errorTypes,
                          "'showTrimmed' must be between TRUE and FALSE.",
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

# ----- FASTA / REGEX / CSV path checks --------------------------------------

# NOTE: Regex patterns below escape the leading dot. The pre-Phase-4 versions
# used unescaped dots (`".fa$"`, `".ab1$"`) which incorrectly matched any
# 3-character extension ending in "fa"/"ab1" (e.g. ".Xfa", ".Xab1").
.FASTA_EXT_REGEX <- "\\.(fa|fasta)$"
.AB1_EXT_REGEX   <- "\\.ab1$"

checkFASTA_File <- function(inputSource, FASTA_File, errors, errorTypes) {
    if (is.null(FASTA_File)) {
        return(.errAppend(errors, errorTypes,
                          "'FASTA_File' cannot be NULL.",
                          "PARAMETER_VALUE_ERROR"))
    }
    if (typeof(FASTA_File) != "character") {
        return(.errAppend(errors, errorTypes,
                          "'FASTA_File' must be character type.",
                          "PARAMETER_VALUE_ERROR"))
    }
    if (!file.exists(FASTA_File)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", FASTA_File, "' file does not exist."),
                          "FILE_NOT_EXIST_ERROR"))
    }
    .requireExt(FASTA_File, .FASTA_EXT_REGEX, errors, errorTypes,
                 message = paste0("'", FASTA_File,
                                  "' file extension must be '.fa' or '.fasta'."))
}

## Issue #92 fix:
## REGEX_SuffixForward / REGEX_SuffixReverse can be NULL (or NA) to indicate
## a single-direction dataset (e.g. forward-only 16S barcoding). In that
## case the construction code substitutes an internal sentinel that never
## matches any filename, so the existing "no reads detected" warning path
## handles the missing direction gracefully without crashing the build.
.NEVER_MATCH_REGEX <- "NEVER_MATCH_FORWARD_REVERSE_ONLY_SENTINEL"

checkREGEX_SuffixForward <- function(REGEX_SuffixForward, errors, errorTypes) {
    if (is.null(REGEX_SuffixForward) || (length(REGEX_SuffixForward) == 1L &&
                                          is.na(REGEX_SuffixForward))) {
        return(list(errors, errorTypes))   # accept NULL/NA: reverse-only run
    }
    if (typeof(REGEX_SuffixForward) != "character") {
        return(.errAppend(errors, errorTypes,
                          "'REGEX_SuffixForward' must be character type.",
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

checkREGEX_SuffixReverse <- function(REGEX_SuffixReverse, errors, errorTypes) {
    if (is.null(REGEX_SuffixReverse) || (length(REGEX_SuffixReverse) == 1L &&
                                          is.na(REGEX_SuffixReverse))) {
        return(list(errors, errorTypes))   # accept NULL/NA: forward-only run
    }
    if (typeof(REGEX_SuffixReverse) != "character") {
        return(.errAppend(errors, errorTypes,
                          "'REGEX_SuffixReverse' must be character type.",
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

checkCSV_NamesConversion <- function(CSV_NamesConversion, errors, errorTypes) {
    if (is.null(CSV_NamesConversion)) {
        return(.errAppend(errors, errorTypes,
                          "'CSV_NamesConversion' cannot be NULL.",
                          "PARAMETER_VALUE_ERROR"))
    }
    if (typeof(CSV_NamesConversion) != "character") {
        return(.errAppend(errors, errorTypes,
                          "'CSV_NamesConversion' must be character type.",
                          "PARAMETER_VALUE_ERROR"))
    }
    if (!file.exists(CSV_NamesConversion)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", CSV_NamesConversion,
                                 "' file does not exist."),
                          "CSV_FILE_NOT_EXIST_ERROR"))
    }
    list(errors, errorTypes)
}

checkReadFileNameExist <- function(readFileName, errors, errorTypes) {
    if (!file.exists(readFileName)) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", readFileName, "' file does not exist."),
                          "FILE_NOT_EXIST_ERROR"))
    }
    list(errors, errorTypes)
}

checkReadFileName <- function(readFileName, inputSource, errors, errorTypes) {
    if (inputSource == "ABIF") {
        return(.requireExt(readFileName, .AB1_EXT_REGEX, errors, errorTypes,
                            message = paste0("'", readFileName,
                                             "' file extension must be '.ab1'.")))
    }
    if (inputSource == "FASTA") {
        return(.requireExt(readFileName, .FASTA_EXT_REGEX, errors, errorTypes,
                            message = paste0("'", readFileName,
                                             "' file extension must be ",
                                             "'.fa' or '.fasta'.")))
    }
    list(errors, errorTypes)
}

checkTargetFastaName <- function(targetFastaName, fastaReadName,
                                 readFileName, errors, errorTypes) {
    if (isEmpty(targetFastaName)) {
        return(.errAppend(errors, errorTypes,
                          paste0("The name '", fastaReadName,
                                 "' is not in the '", basename(readFileName),
                                 "' FASTA file"),
                          "FASTA_NAME_NOT_EXIST"))
    }
    list(errors, errorTypes)
}

checkReadFeature <- function(readFeature, errors, errorTypes) {
    .requireEnum(readFeature, "readFeature",
                  c("Forward Read", "Reverse Read"),
                  errors, errorTypes)
}

checkQualityPhredScores <- function(qualityPhredScores, errors, errorTypes) {
    if (length(qualityPhredScores) == 0) {
        return(.errAppend(errors, errorTypes,
                          "'qualityPhredScores' length cannot be zero.",
                          "PARAMETER_VALUE_ERROR"))
    }
    # Use is.integer fast-path — falls back to whole-number check if double.
    is_whole <- if (is.integer(qualityPhredScores)) TRUE else
                isTRUE(all.equal(qualityPhredScores,
                                  as.integer(qualityPhredScores)))
    if (!is_whole) {
        return(.errAppend(errors, errorTypes,
                          paste0("All elements in 'qualityPhredScores' ",
                                 "vector must be integer."),
                          "PARAMETER_VALUE_ERROR"))
    }
    list(errors, errorTypes)
}

# -----------------------------------------------------------------------------
# Cross-cutting structural / heavy I/O checks (kept procedural)
# -----------------------------------------------------------------------------

checkAb1FastaCsv <- function(ABIF_Directory, FASTA_File,
                             CSV_NamesConversion, inputSource,
                             errors, errorTypes) {
    if (!file.exists(CSV_NamesConversion)) {
        return(.errAppend(errors, errorTypes,
                          paste0("CSV_NamesConversion: '",
                                 CSV_NamesConversion,
                                 "' file does not exist."),
                          "FILE_NOT_EXIST_ERROR"))
    }
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

    # Vectorized membership checks — replaces the earlier per-element lapply.
    missing_in_csv <- setdiff(sourceReads, csvReads)
    if (length(missing_in_csv) > 0L) {
        warnings <- c(warnings,
                       paste0("'", missing_in_csv,
                              "' is not in the csv file (",
                              CSV_NamesConversion, ")"))
    }
    missing_in_src <- setdiff(csvReads, sourceReads)
    if (length(missing_in_src) > 0L) {
        warnings <- c(warnings,
                       paste0("'", missing_in_src,
                              "' is not in the parent directory."))
    }

    required_cols <- c("contig", "direction", "reads")
    missing_cols <- setdiff(required_cols, colnames(csvFile))
    if (length(missing_cols) > 0L) {
        return(.errAppend(errors, errorTypes,
                          paste0("'", missing_cols[[1]],
                                 "' is not in the csv file (",
                                 CSV_NamesConversion, ")"),
                          "CSV_MISMATCH_ERROR"))
    }

    dirVals <- unique(as.character(csvFile$direction))
    if (!(length(dirVals) %in% c(1L, 2L) &&
          all(dirVals %in% c("F", "R")))) {
        return(.errAppend(errors, errorTypes,
                          "In the 'direction' column of your CSV file, you can only have 'F' and 'R'",
                          "CSV_VALUE_ERROR"))
    }

    if (length(warnings) != 0L) invisible(lapply(warnings, log_warn))
    list(errors, errorTypes)
}

# ----- Warning-only diagnostics (never produce errors) -----------------------

checkGreplForward <- function(forwardSelectInputFiles, warnings, warningsType) {
    if (length(forwardSelectInputFiles) == 0L) {
        return(.errAppend(warnings, warningsType,
                          paste0("Your 'contigName' and 'REGEX_SuffixForward' ",
                                 "regular expression parameters can not match ",
                                 "any forward reads."),
                          "REGEX_MATCH_WARN"))
    }
    list(warnings, warningsType)
}

checkGreplReverse <- function(reverseSelectInputFiles, warnings, warningsType) {
    if (length(reverseSelectInputFiles) == 0L) {
        return(.errAppend(warnings, warningsType,
                          paste0("Your 'contigName' and 'REGEX_SuffixReverse' ",
                                 "regular expression parameters can not match ",
                                 "any reverse reads."),
                          "REGEX_MATCH_WARN"))
    }
    list(warnings, warningsType)
}

checkCSVConvForward <- function(forwardReads, warnings, warningsType) {
    if (length(forwardReads) == 0L) {
        return(.errAppend(warnings, warningsType,
                          paste0("The names of forward reads in your ",
                                 "'CSV_NamesConversion' CSV file do not ",
                                 "match any forward reads."),
                          "CSV_MATCH_WARN"))
    }
    list(warnings, warningsType)
}

checkCSVConvReverse <- function(reverseReads, warnings, warningsType) {
    if (length(reverseReads) == 0L) {
        return(.errAppend(warnings, warningsType,
                          paste0("The names of reverse reads in your ",
                                 "'CSV_NamesConversion' CSV file do not ",
                                 "match any reverse reads."),
                          "CSV_MATCH_WARN"))
    }
    list(warnings, warningsType)
}
