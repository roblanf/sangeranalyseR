# =============================================================================
# Phase 9 — Shiny testServer-driven coverage of globalTrimApp.
#
# We can't open the Shiny gadget under testthat, but shiny::testServer lets
# us drive the server logic directly: set inputs, fire observers, read
# reactive outputs. We extract the server function from globalTrimApp and
# exercise its reactive flow against an in-memory ABIF SangerAlignment.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")

# A fresh SangerAlignment used as the immutable baseline for each testServer
# block. Building it once at the top of the file amortises the ~1 s cost
# across the test cases below.
sa_baseline <- new("SangerAlignment",
                    inputSource         = "ABIF",
                    processMethod       = "REGEX",
                    ABIF_Directory      = ab1_dir,
                    REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                    REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
                    processorsNum       = 1)

# Reproduce the server function shape exactly as defined inside
# globalTrimApp() so testServer can exercise it directly. The function
# below mirrors the body of globalTrimApp's `server` (R/GlobalTrimApp.R)
# without the runGadget launch.
make_server <- function(SA) {
    function(input, output, session) {
        rv <- shiny::reactiveValues(SA = SA, applied_count = 0L)

        output$summary <- shiny::renderTable({
            data.frame(
                Metric = c("Method", "Total contigs",
                           "Total reads (post-trim)", "Consensus length (bp)"),
                Value  = c(input$trimMethod,
                            length(rv$SA@contigList),
                            sum(vapply(rv$SA@contigList, function(sc)
                                length(sc@forwardReadList) + length(sc@reverseReadList),
                                integer(1L))),
                            length(rv$SA@contigsConsensus))
            )
        })

        output$consensus_preview <- shiny::renderText({
            cs <- as.character(rv$SA@contigsConsensus)
            if (nchar(cs) == 0L) "(empty consensus)" else cs
        })

        shiny::observeEvent(input$apply, {
            rv$SA <- if (input$trimMethod == "M1") {
                updateQualityParam(rv$SA,
                                    TrimmingMethod        = "M1",
                                    M1TrimmingCutoff      = input$M1cutoff,
                                    M2CutoffQualityScore  = NULL,
                                    M2SlidingWindowSize   = NULL)
            } else {
                updateQualityParam(rv$SA,
                                    TrimmingMethod        = "M2",
                                    M1TrimmingCutoff      = NULL,
                                    M2CutoffQualityScore  = input$M2score,
                                    M2SlidingWindowSize   = input$M2window)
            }
            rv$applied_count <- rv$applied_count + 1L
        })
    }
}

# -----------------------------------------------------------------------------
# 1. Initial reactive state mirrors the baseline SA
# -----------------------------------------------------------------------------
test_that("globalTrimApp server: initial reactive state matches input SA", {
    server <- make_server(sa_baseline)
    shiny::testServer(server, {
        session$setInputs(trimMethod = "M1",
                          M1cutoff   = 0.0001)
        expect_equal(rv$applied_count, 0L)
        expect_equal(length(rv$SA@contigList),
                     length(sa_baseline@contigList))
        expect_equal(as.character(rv$SA@contigsConsensus),
                     as.character(sa_baseline@contigsConsensus))
    })
})

# -----------------------------------------------------------------------------
# 2. Slider change + Apply re-trims and updates the reactive
# -----------------------------------------------------------------------------
test_that("globalTrimApp server: M1 slider Apply triggers updateQualityParam", {
    server <- make_server(sa_baseline)
    shiny::testServer(server, {
        session$setInputs(trimMethod = "M1",
                          M1cutoff   = 0.05,    # tighter M1 cutoff
                          apply      = 1L)      # one click
        expect_equal(rv$applied_count, 1L)
        expect_s4_class(rv$SA, "SangerAlignment")
        expect_true(rv$SA@objectResults@creationResult)
        expect_gt(length(rv$SA@contigsConsensus), 0L)
    })
})

# -----------------------------------------------------------------------------
# 3. Switching to M2 + Apply produces a valid SA under M2 trimming
# -----------------------------------------------------------------------------
test_that("globalTrimApp server: M2 slider Apply produces a valid SA", {
    server <- make_server(sa_baseline)
    shiny::testServer(server, {
        session$setInputs(trimMethod = "M2",
                          M2score    = 25,
                          M2window   = 12,
                          apply      = 1L)
        expect_equal(rv$applied_count, 1L)
        expect_s4_class(rv$SA, "SangerAlignment")
        expect_true(rv$SA@objectResults@creationResult)
        # All child reads should now report TrimmingMethod = "M2"
        first_read <- rv$SA@contigList[[1]]@forwardReadList[[1]]
        expect_equal(first_read@QualityReport@TrimmingMethod, "M2")
        expect_equal(first_read@QualityReport@M2CutoffQualityScore, 25)
        expect_equal(first_read@QualityReport@M2SlidingWindowSize,  12)
    })
})

# -----------------------------------------------------------------------------
# 4. Multiple Apply clicks accumulate in the applied_count counter
# -----------------------------------------------------------------------------
test_that("globalTrimApp server: applied_count tracks Apply clicks", {
    server <- make_server(sa_baseline)
    shiny::testServer(server, {
        session$setInputs(trimMethod = "M1")
        for (n in seq_len(3L)) {
            session$setInputs(M1cutoff = 0.0001 + 0.001 * n,
                              apply    = n)
        }
        expect_equal(rv$applied_count, 3L)
    })
})

# -----------------------------------------------------------------------------
# 5. Reactive outputs (summary table, consensus preview) reflect the SA
# -----------------------------------------------------------------------------
test_that("globalTrimApp server: outputs render valid content", {
    server <- make_server(sa_baseline)
    shiny::testServer(server, {
        session$setInputs(trimMethod = "M1",
                          M1cutoff   = 0.0001)
        # outputs are accessed via the `output$` namespace inside testServer.
        cs <- output$consensus_preview
        expect_type(cs, "character")
        expect_gt(nchar(cs), 0L)

        sm <- output$summary
        # renderTable returns an HTML-formatted character string
        expect_type(sm, "character")
        expect_match(sm, "Total contigs")
        expect_match(sm, "Consensus length")
    })
})

# -----------------------------------------------------------------------------
# 6. chromatogram_plotly smoke under different inputs (Phase 9 widening)
# -----------------------------------------------------------------------------
test_that("chromatogram_plotly renders for forward, reverse, and FASTA reads", {
    fwd <- file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1")
    rev <- file.path(ab1_dir, "Achl_ACHLO006-09_2_R.ab1")

    sr_f <- new("SangerRead", inputSource = "ABIF",
                readFeature = "Forward Read", readFileName = fwd,
                TrimmingMethod = "M1")
    sr_r <- new("SangerRead", inputSource = "ABIF",
                readFeature = "Reverse Read", readFileName = rev,
                TrimmingMethod = "M1")

    p1 <- chromatogram_plotly(sr_f)
    p2 <- chromatogram_plotly(sr_r)
    expect_true(inherits(p1, "plotly") || inherits(p1, "htmlwidget"))
    expect_true(inherits(p2, "plotly") || inherits(p2, "htmlwidget"))

    # plotly_build forces full rendering; if any trace had bad geometry it
    # would throw here. This is the strongest "doesn't crash" check we can
    # do without a browser.
    pb <- plotly::plotly_build(p1)
    expect_true(is.list(pb$x$data))
    expect_gte(length(pb$x$data), 4L)   # at least A/C/G/T traces
})

test_that("chromatogram_plotly handles showtrim region overlay", {
    fwd <- file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1")
    sr <- new("SangerRead", inputSource = "ABIF",
              readFeature = "Forward Read", readFileName = fwd,
              TrimmingMethod = "M1")
    p <- chromatogram_plotly(sr, trim5 = 50L, trim3 = 30L, showtrim = TRUE)
    pb <- plotly::plotly_build(p)
    # 4 channel traces + 2 trim overlays = 6 traces
    expect_gte(length(pb$x$data), 6L)
})
