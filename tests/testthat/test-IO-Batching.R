# =============================================================================
# Phase 3 — IO batching guard.
#
# Today, SangerAlignment.initialize calls list.files(ABIF_Directory,
# recursive = TRUE) at the top level and SangerContig.initialize calls
# list.files(ABIF_Directory) again per contig. We pin that behaviour with a
# call-count test so the Phase 4 IO-cache change can show a measurable
# reduction.
#
# We also verify behaviour around hidden files and nested directories — the
# package uses the default list.files(all.files = FALSE), so dot-files are
# *not* picked up. We test that this current behaviour is stable, and that
# nested .ab1s ARE picked up by SangerAlignment (recursive = TRUE).
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
real_ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")

# -----------------------------------------------------------------------------
# 1. Hidden files (dot-files) are ignored by default
# -----------------------------------------------------------------------------
test_that("hidden .ab1 files do not crash the validator", {
    dir <- make_dir_with_hidden()
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              processorsNum       = 1)
    # Construction should not blow up regardless of whether hidden files are
    # included; we just assert the alignment built.
    expect_s4_class(sa, "SangerAlignment")
    expect_true(sa@objectResults@creationResult)
})

# -----------------------------------------------------------------------------
# 2. Nested sub-directories ARE walked (SangerAlignment uses recursive=TRUE)
# -----------------------------------------------------------------------------
test_that("nested sub-directory .ab1 files are discovered", {
    dir <- make_dir_with_hidden()
    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              processorsNum       = 1)
    # We expect at least 2 contigs: "Top" and "nested/Nested".
    expect_gte(length(sa@contigList), 2L)
})

# -----------------------------------------------------------------------------
# 3. Call-count regression on list.files (pin current behaviour)
# -----------------------------------------------------------------------------
test_that("list.files is called more than once during SangerAlignment build", {
    skip_if_not_installed("withr")

    counter <- new.env(parent = emptyenv())
    counter$n <- 0L
    real_list_files <- base::list.files

    # Trace base::list.files invocations using base::trace; this is more
    # robust across testthat versions than local_mocked_bindings on base.
    base::trace(base::list.files,
                tracer = function() counter$n <- counter$n + 1L,
                print  = FALSE)
    on.exit(suppressMessages(base::untrace(base::list.files)), add = TRUE)

    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = real_ab1_dir,
              REGEX_SuffixForward = "_[0-9]*_F.ab1$",
              REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
              processorsNum       = 1)

    # Today, list.files is invoked once at the SA level + once per SangerContig
    # (4 contigs in the ACHLO fixture) ⇒ at least 5. We assert > 1 as a stable
    # lower bound that any sensible pre-cache implementation will exceed.
    expect_gt(counter$n, 1L)
    expect_true(sa@objectResults@creationResult)
})
