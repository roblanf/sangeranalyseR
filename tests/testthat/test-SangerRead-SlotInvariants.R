# =============================================================================
# Phase 3 — Slot-duplication invariants.
#
# `SangerRead` extends `sangerseq` (sangerseqR) and additionally stores five
# "Raw" copies of slots that already exist on the parent class:
#   primarySeqRaw, secondarySeqRaw, peakPosMatrixRaw, peakAmpMatrixRaw,
#   plus traceMatrix (inherited as-is).
#
# Phase 4 §4 plans to drop the redundant "*Raw" slots and rely on the
# inherited parent slots. These tests guard byte-identical equivalence today
# so the slot collapse can't silently regress data.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
forward_file <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                          "ACHLO", "Achl_ACHLO006-09_1_F.ab1")
reverse_file <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                          "ACHLO", "Achl_ACHLO006-09_2_R.ab1")

assert_raw_slots_match_parent <- function(file, feature) {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = feature,
              readFileName   = file,
              TrimmingMethod = "M1")
    expect_true(sr@objectResults@creationResult)

    # Reload the source via sangerseqR directly. This is exactly what the
    # SangerRead initializer does internally (R/ClassSangerRead.R:200,204).
    abif_obj <- sangerseqR::read.abif(file)
    ref      <- sangerseqR::sangerseq(abif_obj)

    # Sequence equality — compared as character to be format-agnostic.
    expect_identical(as.character(sr@primarySeqRaw),
                     as.character(ref@primarySeq))
    expect_identical(as.character(sr@secondarySeqRaw),
                     as.character(ref@secondarySeq))

    # Numeric matrix equality.
    expect_identical(sr@peakPosMatrixRaw, ref@peakPosMatrix)
    expect_identical(sr@peakAmpMatrixRaw, ref@peakAmpMatrix)

    # Inherited slot — verifies the inheritance contract is intact.
    expect_identical(sr@traceMatrix, ref@traceMatrix)

    # Sanity: object reports as both SangerRead AND sangerseq.
    expect_true(is(sr, "SangerRead"))
    expect_true(is(sr, "sangerseq"))
}

test_that("Forward SangerRead: *Raw slots match sangerseq parent slots", {
    assert_raw_slots_match_parent(forward_file, "Forward Read")
})

test_that("Reverse SangerRead: *Raw slots match sangerseq parent slots", {
    assert_raw_slots_match_parent(reverse_file, "Reverse Read")
})

test_that("primarySeqID and secondarySeqID are propagated from sangerseq", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = forward_file,
              TrimmingMethod = "M1")
    ref <- sangerseqR::sangerseq(sangerseqR::read.abif(forward_file))
    # These are inherited slots — they live on the sangerseq parent.
    expect_identical(sr@primarySeqID,   ref@primarySeqID)
    expect_identical(sr@secondarySeqID, ref@secondarySeqID)
})

test_that("primarySeq (post-MakeBaseCalls) length matches qualityPhredScores", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = forward_file,
              TrimmingMethod = "M1")
    # MakeBaseCallsInside re-runs base calling and may produce a different
    # number of called bases than the raw primarySeq. The length of the
    # post-basecall primarySeq must equal the length of qualityPhredScores
    # (this is the invariant `MakeBaseCallsInside` is supposed to maintain).
    expect_equal(length(sr@primarySeq),
                 length(sr@QualityReport@qualityPhredScores))
})
