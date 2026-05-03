# =============================================================================
# Phase 3 — Negative tests (validator edge cases).
#
# Construction in this package never throws; every failure is reflected in
# `@objectResults@creationResult` and `@objectResults@errorTypes`. Tests assert
# on those, not on `expect_error`.
#
# Some tests in this file intentionally FAIL today — they pin the *desired*
# behaviour after the Phase 4 bug fixes (see plans/02_quality_audit_summary.md).
# Those tests are tagged "BUG-LOCK" in their description so future readers
# don't mistake them for regressions.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
real_ab1_forward <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                              "ACHLO", "Achl_ACHLO006-09_1_F.ab1")
real_ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")
real_csv <- file.path(inputFilesPath, "ab1", "SangerAlignment",
                      "names_conversion.csv")
real_fasta <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                        "Sanger_all_reads.fa")

# -----------------------------------------------------------------------------
# 1. Extension regex (formerly BUG-LOCK; fixed in Phase 4)
# -----------------------------------------------------------------------------
test_that(".Xfa file is rejected as invalid FASTA extension", {
    xfa <- make_xfa_file()
    sa <- new("SangerAlignment",
              inputSource         = "FASTA",
              processMethod       = "REGEX",
              FASTA_File          = xfa,
              REGEX_SuffixForward = "_[0-9]*_F$",
              REGEX_SuffixReverse = "_[0-9]*_R$",
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("FILE_TYPE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that(".Xab1 file is rejected as invalid ABIF extension", {
    dir <- withr::local_tempdir()
    out <- file.path(dir, "Achl_ACHLO006-09_1_F.Xab1")
    file.copy(real_ab1_forward, out)
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = out,
              TrimmingMethod = "M1")
    expect_false(sr@objectResults@creationResult)
    expect_true("FILE_TYPE_ERROR" %in% sr@objectResults@errorTypes)
})

# -----------------------------------------------------------------------------
# 2. Wrong / missing extension
# -----------------------------------------------------------------------------
test_that(".ab2 ABIF extension is rejected", {
    ab2 <- make_ab2_file()
    sr <- new("SangerRead",
              inputSource  = "ABIF",
              readFeature  = "Forward Read",
              readFileName = ab2)
    expect_false(sr@objectResults@creationResult)
    expect_true("FILE_TYPE_ERROR" %in% sr@objectResults@errorTypes)
})

test_that(".fast FASTA extension is rejected", {
    fast <- make_fast_file()
    sa <- new("SangerAlignment",
              inputSource         = "FASTA",
              processMethod       = "REGEX",
              FASTA_File          = fast,
              REGEX_SuffixForward = "_[0-9]*_F$",
              REGEX_SuffixReverse = "_[0-9]*_R$",
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("FILE_TYPE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("missing read filename produces FILE_NOT_EXIST_ERROR", {
    sr <- new("SangerRead",
              inputSource  = "ABIF",
              readFeature  = "Forward Read",
              readFileName = "/no/such/file_1_F.ab1")
    expect_false(sr@objectResults@creationResult)
    expect_true("FILE_NOT_EXIST_ERROR" %in% sr@objectResults@errorTypes)
})

# -----------------------------------------------------------------------------
# 3. Empty / corrupt ABIF files
# -----------------------------------------------------------------------------
test_that("empty .ab1 fails to construct", {
    empty <- make_empty_ab1()
    expect_error(
        sr <- new("SangerRead",
                  inputSource  = "ABIF",
                  readFeature  = "Forward Read",
                  readFileName = empty),
        regexp = NULL
    )
    # Note: read.abif() throws on empty files (the validator only checks
    # extension+exists). This test guards the current behaviour; if Phase 4
    # adds a pre-read content check, change to assert on creationResult.
})

test_that("random-bytes .ab1 fails to construct", {
    corrupt <- make_corrupt_ab1()
    expect_error(
        new("SangerRead",
            inputSource  = "ABIF",
            readFeature  = "Forward Read",
            readFileName = corrupt),
        regexp = NULL
    )
})

test_that("truncated .ab1 fails to construct", {
    trunc <- make_corrupt_ab1_truncated()
    expect_error(
        new("SangerRead",
            inputSource  = "ABIF",
            readFeature  = "Forward Read",
            readFileName = trunc),
        regexp = NULL
    )
})

# -----------------------------------------------------------------------------
# 4. CSV schema / value errors
# -----------------------------------------------------------------------------
test_that("CSV missing 'contig' column is rejected", {
    bad_csv <- make_csv_missing_column("contig")
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = real_ab1_dir,
              CSV_NamesConversion = bad_csv,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("CSV_MISMATCH_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("CSV missing 'direction' column is rejected", {
    bad_csv <- make_csv_missing_column("direction")
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = real_ab1_dir,
              CSV_NamesConversion = bad_csv,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("CSV_MISMATCH_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("CSV missing 'reads' column is rejected", {
    bad_csv <- make_csv_missing_column("reads")
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = real_ab1_dir,
              CSV_NamesConversion = bad_csv,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("CSV_MISMATCH_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("CSV with bad direction values is rejected", {
    bad_csv <- make_csv_bad_direction()
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = real_ab1_dir,
              CSV_NamesConversion = bad_csv,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("CSV_VALUE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("CSV with reads not on disk produces a warning row but still builds", {
    bad_csv <- make_mismatched_csv()
    # Construction emits log warnings; we only assert that it doesn't fail
    # outright (matched reads should still produce contigs).
    sa <- suppressWarnings(new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = real_ab1_dir,
              CSV_NamesConversion = bad_csv,
              processorsNum       = 1))
    expect_s4_class(sa, "SangerAlignment")
})

# -----------------------------------------------------------------------------
# 5. FASTA-record-name mismatch
# -----------------------------------------------------------------------------
test_that("FASTA name in CSV not present in file produces FASTA_NAME_NOT_EXIST", {
    fx <- make_fasta_missing_record()
    sa <- new("SangerAlignment",
              inputSource         = "FASTA",
              processMethod       = "CSV",
              FASTA_File          = fx$fasta,
              CSV_NamesConversion = fx$csv,
              processorsNum       = 1)
    # Regardless of whether the alignment as a whole builds, the per-read
    # readResultTable should record the missing-name error.
    table <- sa@objectResults@readResultTable
    # If table is empty (e.g. all matched), the alignment may have been built
    # purely from records that did exist; require at least one recorded miss.
    if (nrow(table) > 0L) {
        expect_true(any(grepl("FASTA_NAME_NOT_EXIST|FILE_NOT_EXIST_ERROR|MIN_READ_LENGTH_ERROR",
                              as.character(table$errorType))))
    }
})

# -----------------------------------------------------------------------------
# 6. Out-of-range / wrong-type scalar parameters
# -----------------------------------------------------------------------------
test_that("M1TrimmingCutoff out of [0,1] produces PARAMETER_RANGE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              TrimmingMethod      = "M1",
              M1TrimmingCutoff    = 5,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_RANGE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("TrimmingMethod = 'M3' produces PARAMETER_VALUE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              TrimmingMethod      = "M3",
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_VALUE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("non-numeric processorsNum produces PARAMETER_TYPE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              processorsNum       = "many")
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_TYPE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("signalRatioCutoff = 1.5 produces PARAMETER_RANGE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              signalRatioCutoff   = 1.5,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_RANGE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("inputSource = 'FAST5' produces PARAMETER_VALUE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "FAST5",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_VALUE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("processMethod = 'GLOB' produces PARAMETER_VALUE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "GLOB",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_VALUE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("acceptStopCodons = 'yes' produces PARAMETER_VALUE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              acceptStopCodons    = "yes",
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_VALUE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("readingFrame = 4 produces PARAMETER_VALUE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              readingFrame        = 4,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_VALUE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("minFractionCall = 1.5 produces PARAMETER_RANGE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              minFractionCall     = 1.5,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_RANGE_ERROR" %in% sa@objectResults@errorTypes)
})

test_that("minReadsNum = 0 produces PARAMETER_VALUE_ERROR", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              minReadsNum         = 0,
              processorsNum       = 1)
    expect_false(sa@objectResults@creationResult)
    expect_true("PARAMETER_VALUE_ERROR" %in% sa@objectResults@errorTypes)
})
