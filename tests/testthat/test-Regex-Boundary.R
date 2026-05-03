# =============================================================================
# Phase 4 — Regex-boundary regression tests.
#
# The Phase 2 audit found that `checkFASTA_File`'s regex was `".fa$"` and
# `checkReadFileName`'s ABIF branch used `".ab1$"`. Both have an unescaped `.`
# that matches any character, so:
#
#   "Sanger.Xfa"   matched the .fa regex   (BUG)
#   "Sanger.Xab1"  matched the .ab1 regex  (BUG)
#
# Phase 4 escaped the dot. These tests pin every adjacent boundary so the
# pattern can't be loosened in a future refactor without flagging.
# =============================================================================

ns <- asNamespace("sangeranalyseR")

# -----------------------------------------------------------------------------
# FASTA pattern: \\.(fa|fasta)$
# -----------------------------------------------------------------------------
test_that("FASTA pattern accepts canonical extensions", {
    expect_true(grepl(ns$.FASTA_EXT_REGEX, "x.fa"))
    expect_true(grepl(ns$.FASTA_EXT_REGEX, "x.fasta"))
    expect_true(grepl(ns$.FASTA_EXT_REGEX, "Sanger_all_reads.fa"))
})

test_that("FASTA pattern rejects regex-fooling extensions", {
    expect_false(grepl(ns$.FASTA_EXT_REGEX, "x.Xfa"))
    expect_false(grepl(ns$.FASTA_EXT_REGEX, "x.Xfasta"))
    expect_false(grepl(ns$.FASTA_EXT_REGEX, "x.fast"))
    expect_false(grepl(ns$.FASTA_EXT_REGEX, "x.faa"))
    expect_false(grepl(ns$.FASTA_EXT_REGEX, "x.fa.bak"))
    expect_false(grepl(ns$.FASTA_EXT_REGEX, "fasta_only_filename"))
    expect_false(grepl(ns$.FASTA_EXT_REGEX, ""))
})

# -----------------------------------------------------------------------------
# ABIF pattern: \\.ab1$
# -----------------------------------------------------------------------------
test_that("ABIF pattern accepts canonical .ab1 extensions", {
    expect_true(grepl(ns$.AB1_EXT_REGEX, "Achl_ACHLO006-09_1_F.ab1"))
    expect_true(grepl(ns$.AB1_EXT_REGEX, "x.ab1"))
})

test_that("ABIF pattern rejects regex-fooling extensions", {
    expect_false(grepl(ns$.AB1_EXT_REGEX, "x.Xab1"))
    expect_false(grepl(ns$.AB1_EXT_REGEX, "x.ab2"))
    expect_false(grepl(ns$.AB1_EXT_REGEX, "x.ab1.bak"))
    expect_false(grepl(ns$.AB1_EXT_REGEX, "ab1_only_filename"))
    expect_false(grepl(ns$.AB1_EXT_REGEX, ""))
})

# -----------------------------------------------------------------------------
# End-to-end: the public `check*` wrappers reject all the boundary cases.
# -----------------------------------------------------------------------------
test_that("checkReadFileName ABIF rejects .Xab1 / .ab2 / .ab1.bak", {
    for (bad in c("/tmp/x.Xab1", "/tmp/x.ab2", "/tmp/x.ab1.bak")) {
        out <- checkReadFileName(bad, "ABIF", character(0), character(0))
        expect_equal(out[[2]], "FILE_TYPE_ERROR",
                     info = paste("input:", bad))
    }
})

test_that("checkReadFileName FASTA rejects .Xfa / .fast / .faa", {
    for (bad in c("/tmp/x.Xfa", "/tmp/x.fast", "/tmp/x.faa")) {
        out <- checkReadFileName(bad, "FASTA", character(0), character(0))
        expect_equal(out[[2]], "FILE_TYPE_ERROR",
                     info = paste("input:", bad))
    }
})

test_that("checkReadFileName FASTA accepts .fa and .fasta", {
    for (good in c("/tmp/x.fa", "/tmp/x.fasta")) {
        out <- checkReadFileName(good, "FASTA", character(0), character(0))
        expect_length(out[[1]], 0L)
    }
})

test_that("checkReadFileName ABIF accepts .ab1", {
    out <- checkReadFileName("/tmp/x.ab1", "ABIF",
                              character(0), character(0))
    expect_length(out[[1]], 0L)
})
