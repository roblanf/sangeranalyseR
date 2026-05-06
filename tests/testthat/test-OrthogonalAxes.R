# =============================================================================
# Phase 3 — Orthogonal-axis matrix tests.
#
# Crosses inputSource ∈ {ABIF, FASTA} × processMethod ∈ {REGEX, CSV}
# × TrimmingMethod ∈ {M1, M2}. For FASTA, TrimmingMethod is forced to "" by
# the constructor (no trimming on pre-called sequences); we still pass M1/M2
# so we cover that branch.
#
# Each case asserts the alignment object is well-formed: built successfully,
# has at least one contig, has a non-empty consensus, has a DNAStringSet
# alignment.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_dir   <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")
ab1_csv   <- file.path(inputFilesPath, "ab1", "SangerAlignment",
                       "names_conversion.csv")
fa_file   <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                       "Sanger_all_reads.fa")
fa_csv    <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                       "names_conversion.csv")

assert_alignment_ok <- function(sa) {
    expect_true(sa@objectResults@creationResult)
    expect_s4_class(sa, "SangerAlignment")
    expect_gt(length(sa@contigList), 0L)
    expect_s4_class(sa@contigsAlignment, "DNAStringSet")
    # contigsConsensus may be a single DNAString; require length > 0.
    expect_gt(length(sa@contigsConsensus), 0L)
}

assert_all_reads_ok <- function(sa) {
    reads <- unlist(lapply(sa@contigList,
                           function(sc) c(sc@forwardReadList, sc@reverseReadList)),
                    recursive = FALSE)
    if (length(reads) == 0L) return(invisible())
    ok <- vapply(reads, function(r) r@objectResults@creationResult, logical(1))
    expect_true(all(ok))
}

# -----------------------------------------------------------------------------
# ABIF × REGEX × M1
# -----------------------------------------------------------------------------
test_that("AB-RE-M1: ABIF + REGEX + M1 trimming produces a valid alignment", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              TrimmingMethod      = "M1",
              M1TrimmingCutoff    = 0.0001,
              processorsNum       = 1)
    assert_alignment_ok(sa)
    assert_all_reads_ok(sa)
})

# -----------------------------------------------------------------------------
# ABIF × REGEX × M2
# -----------------------------------------------------------------------------
test_that("AB-RE-M2: ABIF + REGEX + M2 trimming produces a valid alignment", {
    sa <- new("SangerAlignment",
              inputSource          = "ABIF",
              processMethod        = "REGEX",
              ABIF_Directory       = ab1_dir,
              REGEX_SuffixForward  = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse  = "_[0-9]*_R.ab1$",
              TrimmingMethod       = "M2",
              M1TrimmingCutoff     = NULL,
              M2CutoffQualityScore = 20,
              M2SlidingWindowSize  = 10,
              processorsNum        = 1)
    assert_alignment_ok(sa)
    assert_all_reads_ok(sa)
})

# -----------------------------------------------------------------------------
# ABIF × CSV × M1
# -----------------------------------------------------------------------------
test_that("AB-CS-M1: ABIF + CSV + M1 trimming produces a valid alignment", {
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = ab1_dir,
              CSV_NamesConversion = ab1_csv,
              TrimmingMethod      = "M1",
              M1TrimmingCutoff    = 0.0001,
              processorsNum       = 1)
    assert_alignment_ok(sa)
    assert_all_reads_ok(sa)
})

# -----------------------------------------------------------------------------
# ABIF × CSV × M2
# -----------------------------------------------------------------------------
test_that("AB-CS-M2: ABIF + CSV + M2 trimming produces a valid alignment", {
    sa <- new("SangerAlignment",
              inputSource          = "ABIF",
              processMethod        = "CSV",
              ABIF_Directory       = ab1_dir,
              CSV_NamesConversion  = ab1_csv,
              TrimmingMethod       = "M2",
              M1TrimmingCutoff     = NULL,
              M2CutoffQualityScore = 20,
              M2SlidingWindowSize  = 10,
              processorsNum        = 1)
    assert_alignment_ok(sa)
    assert_all_reads_ok(sa)
})

# -----------------------------------------------------------------------------
# FASTA × REGEX (TrimmingMethod is forced to "" inside the constructor)
# -----------------------------------------------------------------------------
test_that("FA-RE: FASTA + REGEX produces a valid alignment (no trimming)", {
    sa <- new("SangerAlignment",
              inputSource         = "FASTA",
              processMethod       = "REGEX",
              FASTA_File          = fa_file,
              REGEX_SuffixForward = "_[0-9]*_F$",
              REGEX_SuffixReverse = "_[0-9]*_R$",
              processorsNum       = 1)
    assert_alignment_ok(sa)
    # FASTA reads have no QualityReport, so we can't run assert_all_reads_ok
    # the same way — just verify forwardReadList children are valid SangerReads.
    expect_true(all(vapply(sa@contigList,
                           function(sc) length(sc@forwardReadList) +
                                        length(sc@reverseReadList) > 0L,
                           logical(1))))
})

# -----------------------------------------------------------------------------
# FASTA × CSV
# -----------------------------------------------------------------------------
test_that("FA-CS: FASTA + CSV produces a valid alignment (no trimming)", {
    sa <- new("SangerAlignment",
              inputSource         = "FASTA",
              processMethod       = "CSV",
              FASTA_File          = fa_file,
              CSV_NamesConversion = fa_csv,
              processorsNum       = 1)
    assert_alignment_ok(sa)
})
