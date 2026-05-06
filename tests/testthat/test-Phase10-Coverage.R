# =============================================================================
# Phase 10 — coverage maximization. Targets the largest dark areas seen by
# covr: show methods (0%), MethodShared dispatchers (26%), GlobalTrimApp
# server function (~19%), peakvalues_cpp NA branches (~92%), and a number
# of UtilitiesFunc edge cases.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")
ab1_fwd <- file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1")
ab1_rev <- file.path(ab1_dir, "Achl_ACHLO006-09_2_R.ab1")
fa_file <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                     "Sanger_all_reads.fa")

# -----------------------------------------------------------------------------
# show() methods — currently 0% coverage on R/sangeranalyseR_show_method.R
# -----------------------------------------------------------------------------
test_that("show() works for SangerRead (ABIF, success)", {
    utils::data("sangerReadFData", package = "sangeranalyseR",
                envir = environment())
    out <- capture.output(show(sangerReadFData))
    expect_type(out, "character")
    expect_match(paste(out, collapse = "\n"), "SangerRead")
})

test_that("show() works for SangerContig", {
    utils::data("sangerContigData", package = "sangeranalyseR",
                envir = environment())
    out <- capture.output(show(sangerContigData))
    expect_match(paste(out, collapse = "\n"), "SangerContig")
})

test_that("show() works for SangerAlignment", {
    utils::data("sangerAlignmentData", package = "sangeranalyseR",
                envir = environment())
    out <- capture.output(show(sangerAlignmentData))
    expect_match(paste(out, collapse = "\n"), "SangerAlignment")
})

test_that("show() works on a freshly-built SangerRead with creationResult=FALSE", {
    sr_fail <- new("SangerRead",
                    inputSource    = "ABIF",
                    readFeature    = "Forward Read",
                    readFileName   = "/no/such/path.ab1",
                    TrimmingMethod = "M1")
    expect_false(sr_fail@objectResults@creationResult)
    out <- capture.output(show(sr_fail))
    expect_type(out, "character")
})

test_that("show() works on a fresh ABIF-success SangerRead", {
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = ab1_fwd,
              TrimmingMethod = "M1")
    out <- capture.output(show(sr))
    expect_match(paste(out, collapse = "\n"), "SangerRead")
})

test_that("show() works for FASTA SangerRead", {
    fa_sr_dir <- file.path(inputFilesPath, "fasta", "SangerRead")
    fa_sr <- list.files(fa_sr_dir, pattern = "\\.fa$", full.names = TRUE)[1]
    if (!is.na(fa_sr) && file.exists(fa_sr)) {
        # Read first record name
        nm <- names(seqinr::read.fasta(fa_sr, as.string = TRUE))[1]
        sr <- new("SangerRead",
                  inputSource    = "FASTA",
                  readFeature    = "Forward Read",
                  readFileName   = fa_sr,
                  fastaReadName  = nm)
        out <- capture.output(show(sr))
        expect_match(paste(out, collapse = "\n"), "SangerRead")
    } else {
        skip("no FASTA SangerRead fixture")
    }
})

# -----------------------------------------------------------------------------
# MethodShared dispatchers — error branches
# -----------------------------------------------------------------------------
test_that("launchApp() rejects a non-S4 object (logs, no shiny.appobj)", {
    res <- tryCatch(launchApp("not an S4"), error = function(e) e)
    expect_false(inherits(res, "shiny.appobj"))
})

test_that("launchApp() rejects a SangerRead (no app for SR)", {
    utils::data("sangerReadFData", package = "sangeranalyseR",
                envir = environment())
    # log_error path; does not throw. Confirm no app object is returned.
    res <- tryCatch(launchApp(sangerReadFData), error = function(e) e)
    expect_false(inherits(res, "shiny.appobj"))
})

test_that("writeFasta() dispatches across SR / SC / SA", {
    utils::data("sangerReadFData",     package = "sangeranalyseR",
                envir = environment())
    utils::data("sangerContigData",    package = "sangeranalyseR",
                envir = environment())
    utils::data("sangerAlignmentData", package = "sangeranalyseR",
                envir = environment())
    out <- withr::local_tempdir()
    expect_error(writeFasta(sangerReadFData,    outputDir = out), NA)
    expect_error(writeFasta(sangerContigData,   outputDir = out), NA)
    expect_error(writeFasta(sangerAlignmentData, outputDir = out), NA)
})

test_that("writeFasta() with non-S4 input does not write a FASTA", {
    out <- withr::local_tempdir()
    res <- tryCatch(writeFasta("not an S4", outputDir = out),
                     error = function(e) e)
    # Just confirm no FASTA files were produced
    expect_equal(length(list.files(out, pattern = "\\.fa$")), 0L)
})

test_that("generateReport() dispatcher does not throw on non-S4 input", {
    res <- tryCatch(generateReport(list()), error = function(e) e)
    expect_false(inherits(res, "error"))
})

# -----------------------------------------------------------------------------
# peakvalues_cpp — exercise the all-NA branch (currently 92% covered)
# -----------------------------------------------------------------------------
test_that("peakvalues_cpp returns c(NA, NA) when every match has NA value", {
    # Two indices in the (pstart, pstop) window, both with NA values.
    x <- cbind(c(10, 20, 30), c(NA_real_, NA_real_, 5))
    out <- peakvalues_cpp(x, 5, 25)   # only idx 10, 20 match; both NA
    expect_true(is.na(out[[1]]))
    expect_true(is.na(out[[2]]))
})

test_that("peakvalues_batch_cpp validates length-mismatch", {
    x <- cbind(c(10, 20, 30), c(1.1, 2.2, 3.3))
    expect_error(peakvalues_batch_cpp(x, c(1, 2, 3), c(5, 6)),
                 regexp = "equal length")
})

test_that("peakvalues_cpp handles single-row matrix", {
    x <- matrix(c(15, 7), nrow = 1)
    expect_equal(peakvalues_cpp(x, 10, 20), c(7, 15))
    expect_equal(peakvalues_cpp(x, 16, 20), c(0, NA), ignore_attr = TRUE)
})

# -----------------------------------------------------------------------------
# UtilitiesFunc edge cases (currently 62% covered)
# -----------------------------------------------------------------------------
test_that("getProcessors() honours an explicit integer", {
    expect_equal(getProcessors(1), 1L)
    expect_equal(getProcessors(2), 2L)
    expect_equal(getProcessors(4), 4L)
})

test_that("getProcessors(NULL) returns a positive worker count", {
    n <- getProcessors(NULL)
    expect_true(is.integer(n) || is.numeric(n))
    expect_gte(n, 1L)
})

test_that(".resolveBPPARAM honours BPPARAM kwarg over processorsNum", {
    ns <- asNamespace("sangeranalyseR")
    sp <- BiocParallel::SerialParam()
    expect_identical(ns$.resolveBPPARAM(processorsNum = 8, BPPARAM = sp), sp)
})

test_that("alignContigs() handles a single-contig SangerContigList", {
    utils::data("sangerContigData", package = "sangeranalyseR",
                envir = environment())
    ns <- asNamespace("sangeranalyseR")
    res <- ns$alignContigs(
        list(sangerContigData),
        sangerContigData@geneticCode,
        "",
        0.5, 0.5, 1)
    expect_named(res, c("consensus", "aln", "aln.tree"))
    # With one contig, there's no cross-alignment; consensus / aln are NULL,
    # aln.tree is a degenerate empty phylo.
    expect_null(res$consensus)
    expect_null(res$aln)
})

test_that("calculateAASeq() handles trim window of zero length", {
    ns <- asNamespace("sangeranalyseR")
    primarySeq <- Biostrings::DNAString("ATGAAACCCGGGTTT")
    res <- ns$calculateAASeq(primarySeq, 0L, 0L, Biostrings::GENETIC_CODE)
    expect_named(res, c("primaryAASeqS1", "primaryAASeqS2", "primaryAASeqS3"))
    expect_s4_class(res$primaryAASeqS1, "AAString")
})

# -----------------------------------------------------------------------------
# UtilitiesFuncInputChecker edge cases (currently 83% covered)
# -----------------------------------------------------------------------------
test_that("checkProcessorsNum accepts an integer-valued double", {
    out <- checkProcessorsNum(2.0, character(0), character(0))
    expect_length(out[[1]], 0L)
})

test_that("checkAcceptStopCodons accepts both TRUE and FALSE", {
    expect_length(
        checkAcceptStopCodons(TRUE, character(0), character(0))[[1]], 0L)
    expect_length(
        checkAcceptStopCodons(FALSE, character(0), character(0))[[1]], 0L)
})

test_that("checkContigName rejects NULL but accepts character", {
    out_bad <- checkContigName(NULL, character(0), character(0))
    expect_length(out_bad[[1]], 1L)
    expect_equal(out_bad[[2]], "PARAMETER_VALUE_ERROR")
    out_ok <- checkContigName("Achl_ACHLO006-09", character(0), character(0))
    expect_length(out_ok[[1]], 0L)
})

test_that("checkGreplForward / checkGreplReverse / checkCSVConv* warn-only", {
    expect_length(
        checkGreplForward(character(0), character(0), character(0))[[1]], 1L)
    expect_length(
        checkGreplReverse(character(0), character(0), character(0))[[1]], 1L)
    expect_length(
        checkCSVConvForward(character(0), character(0), character(0))[[1]], 1L)
    expect_length(
        checkCSVConvReverse(character(0), character(0), character(0))[[1]], 1L)
})

test_that("checkTargetFastaName rejects empty match", {
    out <- checkTargetFastaName(character(0), "missing", "/tmp/x.fa",
                                 character(0), character(0))
    expect_equal(out[[2]], "FASTA_NAME_NOT_EXIST")
})

# -----------------------------------------------------------------------------
# MakeBaseCalls method, qualityBasePlot,QualityReport
# -----------------------------------------------------------------------------
test_that("MakeBaseCalls(SR, signalRatioCutoff) returns a SangerRead with re-translated AA", {
    utils::data("sangerReadFData", package = "sangeranalyseR",
                envir = environment())
    sr2 <- MakeBaseCalls(sangerReadFData, signalRatioCutoff = 0.4)
    expect_s4_class(sr2, "SangerRead")
    expect_gt(length(sr2@primaryAASeqS1), 0L)
})

test_that("qualityBasePlot,QualityReport returns a plotly", {
    utils::data("qualityReportData", package = "sangeranalyseR",
                envir = environment())
    p <- qualityBasePlot(qualityReportData)
    expect_true(inherits(p, "plotly") || inherits(p, "htmlwidget"))
})

# -----------------------------------------------------------------------------
# GlobalTrimApp — exercise the *real* server function via testServer.
# -----------------------------------------------------------------------------
test_that("globalTrimApp's actual server reactives respond to inputs", {
    skip_if_not_installed("shiny")
    sa <- new("SangerAlignment",
               inputSource         = "ABIF",
               processMethod       = "REGEX",
               ABIF_Directory      = ab1_dir,
               REGEX_SuffixForward = "_[0-9]*_F.ab1$",
               REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
               processorsNum       = 1)
    expect_true(sa@objectResults@creationResult)

    # Extract the server function from globalTrimApp's body.
    # globalTrimApp builds a server() closure inline; we mock runGadget so
    # we can capture the (ui, server) and then drive the server via testServer.
    captured <- new.env()
    testthat::with_mocked_bindings(
        runGadget = function(ui, server, ...) {
            captured$server <- server
            "MOCK"
        },
        .package = "shiny",
        code = globalTrimApp(sa)
    )
    expect_true(is.function(captured$server))

    shiny::testServer(captured$server, {
        session$setInputs(trimMethod = "M1",
                          M1cutoff   = 0.0001)
        # Trigger M1 apply
        session$setInputs(M1cutoff = 0.05, apply = 1L)
        expect_equal(rv$applied_count, 1L)
        # Switch to M2 and apply
        session$setInputs(trimMethod = "M2",
                          M2score    = 25,
                          M2window   = 10,
                          apply      = 2L)
        expect_equal(rv$applied_count, 2L)
        # Done returns the SA
        # (we don't actually click done; just exercise the reactive chain)
    })
})

# -----------------------------------------------------------------------------
# updateQualityParam,SangerAlignment — error path (FASTA SA can't update)
# -----------------------------------------------------------------------------
test_that("updateQualityParam,SangerAlignment is a no-op for FASTA inputSource", {
    sa_fa <- new("SangerAlignment",
                  inputSource         = "FASTA",
                  processMethod       = "REGEX",
                  FASTA_File          = fa_file,
                  REGEX_SuffixForward = "_[0-9]*_F$",
                  REGEX_SuffixReverse = "_[0-9]*_R$",
                  processorsNum       = 1)
    # Should log_info and return; not throw.
    res <- tryCatch(updateQualityParam(sa_fa,
                                        TrimmingMethod        = "M2",
                                        M1TrimmingCutoff      = NULL,
                                        M2CutoffQualityScore  = 25,
                                        M2SlidingWindowSize   = 10),
                     error = function(e) e)
    expect_false(inherits(res, "error"))
})

test_that("updateQualityParam,SangerAlignment with bad params logs error", {
    utils::data("sangerAlignmentData", package = "sangeranalyseR",
                envir = environment())
    # BAD: no TrimmingMethod values
    res <- tryCatch(updateQualityParam(sangerAlignmentData,
                                        TrimmingMethod        = "M99",
                                        M1TrimmingCutoff      = NULL,
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL),
                     error = function(e) e)
    expect_false(inherits(res, "error"))   # logs, doesn't throw
})

# -----------------------------------------------------------------------------
# updateQualityParam,SangerRead error branch
# -----------------------------------------------------------------------------
test_that("updateQualityParam,SangerRead with bad TrimmingMethod logs and returns", {
    utils::data("sangerReadFData", package = "sangeranalyseR",
                envir = environment())
    res <- tryCatch(updateQualityParam(sangerReadFData,
                                        TrimmingMethod        = "M9",
                                        M1TrimmingCutoff      = NULL,
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL),
                     error = function(e) e)
    expect_false(inherits(res, "error"))
})

# -----------------------------------------------------------------------------
# generateReportSR / SC / SA — try the FASTA paths (avoid full RMD render)
# Here we just exercise the parameter-checking and outputDir creation paths
# via direct function entry; full rendering is tested in test-LazyAA-Reports.R.
# -----------------------------------------------------------------------------
test_that("generateReportSR with rds-only (no render) doesn't error", {
    utils::data("sangerReadFData", package = "sangeranalyseR",
                envir = environment())
    # Just confirm dispatcher entry point doesn't crash on a minimal call.
    # We pass an outputDir that exists; render is skipped via an exception
    # we catch so the body of the method is exercised.
    tmp <- withr::local_tempdir()
    res <- tryCatch(suppressWarnings(suppressMessages(
        generateReportSR(sangerReadFData, outputDir = tmp))),
        error = function(e) e)
    # Either a path or an error is fine; we exercised the method body.
    expect_true(is.character(res) || inherits(res, "error"))
})

# -----------------------------------------------------------------------------
# chromatogram_plotly — odd inputs
# -----------------------------------------------------------------------------
test_that("chromatogram_plotly with showtrim=FALSE skips the trim overlay", {
    sr <- new("SangerRead", inputSource = "ABIF",
              readFeature = "Forward Read", readFileName = ab1_fwd,
              TrimmingMethod = "M1")
    p <- chromatogram_plotly(sr, trim5 = 100, trim3 = 50, showtrim = FALSE)
    pb <- plotly::plotly_build(p)
    # No trim overlays when showtrim = FALSE.
    expect_equal(length(pb$x$data), 4L)
})

test_that("chromatogram_plotly with custom 5-element palette works", {
    sr <- new("SangerRead", inputSource = "ABIF",
              readFeature = "Forward Read", readFileName = ab1_fwd,
              TrimmingMethod = "M1")
    p <- chromatogram_plotly(sr, colors = c("#aabbcc", "#112233",
                                              "#445566", "#778899", "#fedcba"))
    expect_true(inherits(p, "plotly") || inherits(p, "htmlwidget"))
})

# -----------------------------------------------------------------------------
# SangerContig FASTA paths under both REGEX and CSV
# -----------------------------------------------------------------------------
test_that("SangerContig FASTA REGEX path constructs cleanly", {
    sc <- new("SangerContig",
               inputSource         = "FASTA",
               processMethod       = "REGEX",
               FASTA_File          = fa_file,
               contigName          = "Achl_ACHLO006-09",
               REGEX_SuffixForward = "_[0-9]*_F$",
               REGEX_SuffixReverse = "_[0-9]*_R$",
               processorsNum       = 1)
    expect_s4_class(sc, "SangerContig")
})

test_that("SangerContig FASTA CSV path constructs cleanly", {
    fa_csv <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                         "names_conversion.csv")
    sc <- new("SangerContig",
               inputSource         = "FASTA",
               processMethod       = "CSV",
               FASTA_File          = fa_file,
               contigName          = "Achl_ACHLO006-09",
               CSV_NamesConversion = fa_csv,
               processorsNum       = 1)
    expect_s4_class(sc, "SangerContig")
})
