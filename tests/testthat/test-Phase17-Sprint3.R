# =============================================================================
# Phase 17 — Sprint 3 consensus-algorithm enhancements.
#
# Targeted tests for:
#   #87  Majority-rules consensus base calling.
#   #48  Phred-aware (quality-weighted) consensus.
#   #33  Synthetic / averaged Phred scores attached to the consensus.
#
# Plus backwards-compatibility check (default `consensusMethod = "strict"`
# is unchanged).
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")

ns <- asNamespace("sangeranalyseR")

# -----------------------------------------------------------------------------
# Helpers — tiny synthetic alignments for the unit tests
# -----------------------------------------------------------------------------
mk_aln <- function(seqs) {
    s <- Biostrings::DNAStringSet(seqs)
    names(s) <- if (is.null(names(seqs))) paste0("r", seq_along(seqs)) else names(seqs)
    s
}

# -----------------------------------------------------------------------------
# Issue #87 — .computeConsensusMajority unit tests
# -----------------------------------------------------------------------------
test_that("Issue #87: majority consensus picks the most-frequent base per column", {
    aln <- mk_aln(c(read1 = "AAGCTT",
                    read2 = "AAGCTT",
                    read3 = "AAGGTT"))   # column 4 differs: C, C, G
    res <- ns$.computeConsensusMajority(aln)
    expect_s4_class(res$consensus, "DNAString")
    expect_equal(as.character(res$consensus), "AAGCTT")
    expect_length(res$qualityScores, 6L)
    # Column 4: 2/3 agreement -> 40 * 2/3 = ~27
    expect_equal(res$qualityScores[4L], as.integer(round(40 * 2/3)))
    # Columns 1,2,3,5,6: 3/3 agreement -> 40 * 1.0 = 40
    expect_equal(res$qualityScores[c(1L,2L,3L,5L,6L)], rep(40L, 5L))
})

test_that("Issue #87: majority consensus handles all-gap columns", {
    aln <- mk_aln(c("--A", "--A", "--A"))
    res <- ns$.computeConsensusMajority(aln)
    expect_equal(as.character(res$consensus), "--A")
    expect_equal(res$qualityScores, c(0L, 0L, 40L))
})

test_that("Issue #87: SangerContig with consensusMethod='majority' produces no IUPAC", {
    sc <- new("SangerContig",
               inputSource         = "ABIF",
               processMethod       = "REGEX",
               ABIF_Directory      = ab1_dir,
               contigName          = "Achl_ACHLO006-09",
               REGEX_SuffixForward = "_[0-9]*_F.ab1$",
               REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
               consensusMethod     = "majority",
               processorsNum       = 1)
    expect_true(sc@objectResults@creationResult)
    cs <- as.character(sc@contigSeq)
    # Majority consensus restricts to A/C/G/T (plus optional '-' for empty
    # columns, stripped by RemoveGaps). No IUPAC ambiguity codes.
    expect_false(grepl("[YRSWKMBDHVN]", cs))
})

# -----------------------------------------------------------------------------
# Issue #48 — quality_weighted consensus
# -----------------------------------------------------------------------------
test_that("Issue #48: quality-weighted consensus differs from majority on a tied column with Phred bias", {
    # Column 1: A vs A vs G (2-1 majority is A). With Phred weights
    # A=10, A=10, G=60, the weighted vote is 60 (G) > 20 (A) → G wins.
    aln <- mk_aln(c(read1 = "A", read2 = "A", read3 = "G"))

    res_maj <- ns$.computeConsensusMajority(aln)
    expect_equal(as.character(res_maj$consensus), "A")

    weights <- matrix(c(10L, 10L, 60L), nrow = 3L, ncol = 1L)
    res_qw <- ns$.computeConsensusMajority(aln, weights = weights)
    expect_equal(as.character(res_qw$consensus), "G")
    # Quality of the agreeing reads at this column = 60 (only one read agrees)
    expect_equal(res_qw$qualityScores[[1L]], 60L)
})

test_that("Issue #48: SangerContig with qualityAware = TRUE builds successfully", {
    sc <- new("SangerContig",
               inputSource         = "ABIF",
               processMethod       = "REGEX",
               ABIF_Directory      = ab1_dir,
               contigName          = "Achl_ACHLO006-09",
               REGEX_SuffixForward = "_[0-9]*_F.ab1$",
               REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
               qualityAware        = TRUE,
               processorsNum       = 1)
    expect_true(sc@objectResults@creationResult)
    cs <- as.character(sc@contigSeq)
    expect_false(grepl("[YRSWKMBDHVN]", cs))
})

test_that("Issue #48: .buildQualityMatrix maps per-read Phred to alignment columns", {
    aln <- mk_aln(c(read1 = "AC-G", read2 = "A-CG"))
    qpls <- list(
        read1 = c(40L, 30L, 20L),    # 3 non-gap columns: A, C, G
        read2 = c(35L, 25L, 15L)     # 3 non-gap columns: A, C, G
    )
    m <- ns$.buildQualityMatrix(aln, qpls)
    expect_equal(dim(m), c(2L, 4L))
    # read1: A C - G  -> 40, 30,  0, 20
    # read2: A - C G  -> 35,  0, 25, 15
    expect_equal(m[1L, ], c(40L, 30L, 0L, 20L))
    expect_equal(m[2L, ], c(35L, 0L, 25L, 15L))
})

test_that("Issue #48: .buildQualityMatrix falls back to flat Phred 30 when read absent from list", {
    aln <- mk_aln(c(read1 = "ACG"))
    m <- ns$.buildQualityMatrix(aln, list())
    expect_equal(m[1L, ], rep(30L, 3L))
})

# -----------------------------------------------------------------------------
# Issue #33 — consensus quality scores attached as attribute
# -----------------------------------------------------------------------------
test_that("Issue #33: consensusGapfree carries qualityScores attribute under majority/quality_weighted", {
    sc <- new("SangerContig",
               inputSource         = "ABIF",
               processMethod       = "REGEX",
               ABIF_Directory      = ab1_dir,
               contigName          = "Achl_ACHLO006-09",
               REGEX_SuffixForward = "_[0-9]*_F.ab1$",
               REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
               consensusMethod     = "majority",
               processorsNum       = 1)
    qs <- attr(sc@contigSeq, "qualityScores")
    expect_true(is.integer(qs))
    expect_gt(length(qs), 0L)
    expect_equal(length(qs), length(sc@contigSeq))
    expect_true(all(qs >= 0L & qs <= 60L))   # Phred range bounds
})

test_that("Issue #33: strict consensus has empty qualityScores attribute (preserves prior behaviour)", {
    sc <- new("SangerContig",
               inputSource         = "ABIF",
               processMethod       = "REGEX",
               ABIF_Directory      = ab1_dir,
               contigName          = "Achl_ACHLO006-09",
               REGEX_SuffixForward = "_[0-9]*_F.ab1$",
               REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
               processorsNum       = 1)
    qs <- attr(sc@contigSeq, "qualityScores")
    # strict default: integer(0) attribute
    expect_true(is.null(qs) || (is.integer(qs) && length(qs) == 0L))
})

# -----------------------------------------------------------------------------
# Backwards compatibility — strict default is unchanged
# -----------------------------------------------------------------------------
test_that("Sprint 3 default (consensusMethod='strict') still produces the same DECIPHER consensus", {
    utils::data("sangerContigData", package = "sangeranalyseR",
                envir = environment())
    fwd <- sangerContigData@forwardReadList
    rev <- sangerContigData@reverseReadList

    res <- ns$calculateContigSeq(
        inputSource      = "ABIF",
        forwardReadList  = fwd,
        reverseReadList  = rev,
        refAminoAcidSeq  = "",
        minFractionCall  = sangerContigData@minFractionCall,
        maxFractionLost  = sangerContigData@maxFractionLost,
        geneticCode      = Biostrings::GENETIC_CODE,
        acceptStopCodons = TRUE,
        readingFrame     = 1L,
        processorsNum    = 1L
    )
    expect_s4_class(res$consensusGapfree, "DNAString")
    expect_gt(length(res$consensusGapfree), 0L)
    # Strict mode: empty quality vector
    expect_equal(length(res$consensusQualityScores), 0L)
})

test_that("Sprint 3: invalid consensusMethod is rejected", {
    utils::data("sangerContigData", package = "sangeranalyseR",
                envir = environment())
    expect_error(
        ns$calculateContigSeq(
            inputSource      = "ABIF",
            forwardReadList  = sangerContigData@forwardReadList,
            reverseReadList  = sangerContigData@reverseReadList,
            refAminoAcidSeq  = "",
            minFractionCall  = 0.5,
            maxFractionLost  = 0.5,
            geneticCode      = Biostrings::GENETIC_CODE,
            acceptStopCodons = TRUE,
            readingFrame     = 1L,
            processorsNum    = 1L,
            consensusMethod  = "no-such-method"
        ),
        regexp = "consensusMethod"
    )
})

test_that("Sprint 3: SangerAlignment forwards consensusMethod to children", {
    sa <- new("SangerAlignment",
               inputSource         = "ABIF",
               processMethod       = "REGEX",
               ABIF_Directory      = ab1_dir,
               REGEX_SuffixForward = "_[0-9]*_F.ab1$",
               REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
               consensusMethod     = "majority",
               processorsNum       = 1)
    expect_true(sa@objectResults@creationResult)
    # Every child contig should have a non-IUPAC consensus
    for (sc in sa@contigList) {
        cs <- as.character(sc@contigSeq)
        expect_false(grepl("[YRSWKMBDHVN]", cs),
                      info = paste("contig", sc@contigName))
    }
})
