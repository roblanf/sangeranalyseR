# =============================================================================
# Phase 8 — globalTrimApp gadget: server logic must be wired correctly and
# the entry point must validate its inputs. We can't actually open the
# Shiny gadget under testthat, but we can call the function with synthetic
# inputs and assert the expected error paths fire.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")
fa_file <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                     "Sanger_all_reads.fa")

test_that("globalTrimApp rejects non-SangerAlignment input", {
    expect_error(globalTrimApp(list()),
                  regexp = "SangerAlignment")
    expect_error(globalTrimApp("not an SA"),
                  regexp = "SangerAlignment")
})

test_that("globalTrimApp rejects FASTA-derived SangerAlignments", {
    sa_fa <- new("SangerAlignment",
                  inputSource         = "FASTA",
                  processMethod       = "REGEX",
                  FASTA_File          = fa_file,
                  REGEX_SuffixForward = "_[0-9]*_F$",
                  REGEX_SuffixReverse = "_[0-9]*_R$",
                  processorsNum       = 1)
    expect_error(globalTrimApp(sa_fa),
                  regexp = "ABIF")
})

test_that("globalTrimApp accepts an ABIF SangerAlignment (entry-point validation)", {
    # We can't launch the gadget under testthat, but we can confirm that the
    # input-validation phase passes and the function reaches runGadget. We
    # mock runGadget to bypass the actual UI launch.
    sa_ab <- new("SangerAlignment",
                  inputSource         = "ABIF",
                  processMethod       = "REGEX",
                  ABIF_Directory      = ab1_dir,
                  REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                  REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
                  processorsNum       = 1)
    expect_true(sa_ab@objectResults@creationResult)

    # Mock shiny::runGadget so the function returns immediately with the SA.
    out <- testthat::with_mocked_bindings(
        runGadget = function(ui, server, ...) "MOCK_GADGET_LAUNCHED",
        .package  = "shiny",
        code      = globalTrimApp(sa_ab)
    )
    expect_equal(out, "MOCK_GADGET_LAUNCHED")
})
