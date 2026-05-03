# =============================================================================
# Phase 6 — lazy amino-acid translation + BiocParallel plumbing.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")
ab1_fwd <- file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1")

# -----------------------------------------------------------------------------
# Lazy AA: default behaviour
# -----------------------------------------------------------------------------
test_that("lazyAA = TRUE leaves primaryAASeqS1/S2/S3 slots empty", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1",
              lazyAA         = TRUE)
    expect_true(sr@objectResults@creationResult)
    expect_equal(length(sr@primaryAASeqS1), 0L)
    expect_equal(length(sr@primaryAASeqS2), 0L)
    expect_equal(length(sr@primaryAASeqS3), 0L)
})

test_that("lazyAA = FALSE eagerly populates AA slots", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1",
              lazyAA         = FALSE)
    expect_true(sr@objectResults@creationResult)
    expect_gt(length(sr@primaryAASeqS1), 0L)
    expect_gt(length(sr@primaryAASeqS2), 0L)
    expect_gt(length(sr@primaryAASeqS3), 0L)
})

test_that("primaryAASeqS{1,2,3}() accessors return the same value as eager slots", {
    sr_lazy  <- new("SangerRead", inputSource="ABIF",
                    readFeature="Forward Read", readFileName=ab1_fwd,
                    TrimmingMethod="M1", lazyAA=TRUE)
    sr_eager <- new("SangerRead", inputSource="ABIF",
                    readFeature="Forward Read", readFileName=ab1_fwd,
                    TrimmingMethod="M1", lazyAA=FALSE)
    expect_identical(as.character(primaryAASeqS1(sr_lazy)),
                     as.character(sr_eager@primaryAASeqS1))
    expect_identical(as.character(primaryAASeqS2(sr_lazy)),
                     as.character(sr_eager@primaryAASeqS2))
    expect_identical(as.character(primaryAASeqS3(sr_lazy)),
                     as.character(sr_eager@primaryAASeqS3))
})

test_that("primaryAASeqS1() returns the cached slot when populated eagerly", {
    sr <- new("SangerRead", inputSource="ABIF",
              readFeature="Forward Read", readFileName=ab1_fwd,
              TrimmingMethod="M1", lazyAA=FALSE)
    expect_identical(primaryAASeqS1(sr), sr@primaryAASeqS1)
})

# -----------------------------------------------------------------------------
# .resolveBPPARAM helper
# -----------------------------------------------------------------------------
ns <- asNamespace("sangeranalyseR")

test_that(".resolveBPPARAM honours an explicit BPPARAM", {
    bp <- BiocParallel::SerialParam()
    expect_identical(ns$.resolveBPPARAM(processorsNum = 4, BPPARAM = bp), bp)
})

test_that(".resolveBPPARAM maps processorsNum = 1 to SerialParam", {
    bp <- ns$.resolveBPPARAM(processorsNum = 1)
    expect_s4_class(bp, "SerialParam")
})

test_that(".resolveBPPARAM maps processorsNum >= 2 to a parallel param", {
    bp <- ns$.resolveBPPARAM(processorsNum = 2)
    if (.Platform$OS.type == "windows") {
        expect_s4_class(bp, "SnowParam")
    } else {
        expect_s4_class(bp, "MulticoreParam")
    }
    expect_equal(BiocParallel::bpnworkers(bp), 2L)
})

test_that(".resolveBPPARAM falls back to SerialParam on non-numeric input", {
    expect_s4_class(ns$.resolveBPPARAM(processorsNum = "many"), "SerialParam")
})

# -----------------------------------------------------------------------------
# BPPARAM end-to-end: SerialParam vs MulticoreParam(2) produce identical SA
# -----------------------------------------------------------------------------
test_that("SangerAlignment under SerialParam vs MulticoreParam yields equal consensus", {
    sa_s <- SangerAlignment(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = ab1_dir,
        REGEX_SuffixForward = "_[0-9]*_F.ab1$",
        REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
        BPPARAM             = BiocParallel::SerialParam()
    )
    bp_par <- if (.Platform$OS.type == "windows")
                  BiocParallel::SnowParam(workers = 2)
              else
                  BiocParallel::MulticoreParam(workers = 2)
    sa_p <- SangerAlignment(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = ab1_dir,
        REGEX_SuffixForward = "_[0-9]*_F.ab1$",
        REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
        BPPARAM             = bp_par
    )
    expect_true(sa_s@objectResults@creationResult)
    expect_true(sa_p@objectResults@creationResult)
    expect_identical(as.character(sa_s@contigsConsensus),
                     as.character(sa_p@contigsConsensus))
    expect_identical(sort(names(sa_s@contigList)),
                     sort(names(sa_p@contigList)))
})

# -----------------------------------------------------------------------------
# Backwards compat: existing SR/SC/SA construction paths still produce
# valid objects with the new defaults (lazy=TRUE, BPPARAM=NULL).
# -----------------------------------------------------------------------------
test_that("SangerAlignment with default new-Phase-6 params still builds", {
    sa <- SangerAlignment(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = ab1_dir,
        REGEX_SuffixForward = "_[0-9]*_F.ab1$",
        REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
        processorsNum       = 1
    )
    expect_true(sa@objectResults@creationResult)
    expect_gt(length(sa@contigList), 0L)
})
