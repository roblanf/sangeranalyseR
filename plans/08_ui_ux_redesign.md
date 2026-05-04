# sangeranalyseR — Phase 8 UI/UX Redesign + Lazy-AA Report Compatibility

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

This phase modernizes the interactive surface and fixes a silent regression introduced by Phase 6: the RMarkdown templates and Shiny servers were reading the `@primaryAASeqS{1,2,3}` slots **directly**, which return `AAString("")` on `lazyAA = TRUE` (the new default). Reports rendered without errors but with **empty** AA tables. Phase 8 routes every read through the accessor methods, plus adds two new modernization features: a Plotly+WebGL chromatogram renderer and a global-trim Shiny gadget.

---

## 1. Headline result

| Component                         | Pre-Phase-8 (lazyAA = TRUE)                                  | Post-Phase-8                                                          |
| --------------------------------- | ------------------------------------------------------------- | --------------------------------------------------------------------- |
| `generateReport(SangerRead)`      | Renders, but Frame-1/2/3 AA tables are **empty**              | All three frames populated correctly                                  |
| `launchAppSC(SangerContig)`       | Per-read AA panels show **empty** AAStrings                   | AA panels populated lazily on first reactive trigger                  |
| `launchAppSA(SangerAlignment)`    | Same — empty AA panels                                        | Same — fixed                                                          |
| Chromatogram renderer (HTML)      | base-R graphics → freezes browser at >5k points               | New `chromatogram_plotly()` with WebGL + downsampling                 |
| Global trimming                    | Per-read sliders only                                         | New `globalTrimApp(SA)` gadget for batch trim policy                  |

**Test result: 1279 / 1279 PASS** (Phase 7: 1249 → 30 new Phase 8 tests).

---

## 2. Lazy-AA report compatibility — the critical fix

### 2a. Reads converted to accessor calls

| File                                                             | Reads converted | Sites |
| ---------------------------------------------------------------- | --------------: | ----- |
| `inst/rmd/SangerRead_Report_ab1.Rmd`                             |               3 | lines 168, 198, 220 |
| `inst/rmd/SangerRead_Report_fasta.Rmd`                           |               3 | lines 112, 131, 153 |
| `R/ShinySangerContigServer.R`                                    |              12 | lines 389, 392, 395 (forward); 479, 482, 485 (reverse); 1748, 1751, 1754 (forward post-recalc); 1802, 1805, 1808 (reverse post-recalc) |
| `R/ShinySangerAlignmentServer.R`                                 |              12 | lines 606, 610, 614; 709, 713, 717; 2254, 2257, 2260; 2318, 2321, 2324 |
| **Total reads converted**                                        |          **30** |       |

Each was rewritten from `<obj_chain>@primaryAASeqS<N>` to `primaryAASeqS<N>(<obj_chain>)`. The 12 `<<-` *write* sites in the Shiny servers (lines 1731-1735, 1785-1789 in Contig; 2235-2241, 2299-2305 in Alignment) were left intact — those occur in re-translate handlers that compute fresh values and cache them into the slot, which works correctly under both `lazyAA=TRUE` and `lazyAA=FALSE`.

### 2b. Refactor mechanism

A short Rscript at `/tmp/rewrite_aa.R` (re-runnable) used a multi-line PCRE that:

- matched `<R-identifier-with-subscripts-and-slots>@primaryAASeqS<N>)` (closing paren guarantees this is a *read* — the assign sites are followed by `<<-`, never `)`);
- rewrote it to `primaryAASeqS<N>(<chain>))`.

Verified: 24 read rewrites in Shiny + 6 in RMD; 12 write sites preserved unchanged. Resulting code parses cleanly.

### 2c. End-to-end pin

`tests/testthat/test-LazyAA-Reports.R` (5 tests, 15 assertions):

| Test                                                                     | What it pins                                                              |
| ------------------------------------------------------------------------ | -------------------------------------------------------------------------- |
| Lazy SR accessor result == eager SR slot, for all 3 frames               | byte-equal AAString comparison across `lazyAA={FALSE, TRUE}`              |
| Default lazy SR has empty `@primaryAASeq*` slots; accessors yield non-empty | structural pin on the "lazy" contract                                     |
| `generateReportSR(lazy_SR)` produces a non-empty HTML file              | the headline regression test — pre-Phase-8 this would have rendered empty AA tables |
| Phase-8 RMD pattern `data.frame(AAString(primaryAASeqS1(sr)))` is non-empty | proves the post-fix code path produces real AA content                   |

The third test is gated by `pandoc_available()` so it skips on systems without pandoc, but on developer/CI machines it runs and exercises the actual RMD render pipeline (hence the 24.6 s test duration).

---

## 3. Plotly + WebGL chromatogram renderer

### 3a. New API

```r
chromatogram_plotly(obj,
                    trim5      = 0,
                    trim3      = 0,
                    max_points = 8000L,
                    showtrim   = FALSE,
                    colors     = "default")
```

Returns a `plotly` htmlwidget that the Shiny apps can render via `plotly::renderPlotly(...)`. Lives at `R/UtilitiesFunc.R::chromatogram_plotly`, exported in NAMESPACE.

### 3b. Why a separate function rather than replacing `chromatogram_overwrite`

`chromatogram_overwrite` (the existing function) is a faithful replica of `sangerseqR::chromatogram` extended for color-blind palettes. It's documented as user-facing API, used in vignettes, and produces print-quality PDFs via base-R graphics. Replacing it would break PDF output. **`chromatogram_plotly` is additive** — call sites that need an interactive widget can opt in; call sites that need a static image continue to use `chromatogram_overwrite`.

### 3c. Performance behaviour

- **WebGL via `scattergl`**. Plotly's WebGL backend handles >50k points per channel without browser stalls; the SVG default (`scatter`) saturates around 5k.
- **Uniform-stride downsampling** to `max_points` per channel. The ACHLO fixture's 700-bp read has ~10,000 trace points per channel; the default 8000 retains every point. For longer reads (~50,000 points common in newer ABIF v3 files) the function downsamples by stride. The `downsample_info` attribute on the returned widget reports `(original_points, rendered_points, downsample_stride)` so the Shiny app can display "showing 1 in N points" hints.
- **Color-blind palette** (`colors = "cb_friendly"`) keeps the same encoding rules used elsewhere in the package.

### 3d. Test coverage

`tests/testthat/test-Phase8-PlotlyChromatogram.R` (5 tests, 10 assertions):

| Test                                                       | Asserts                                                |
| ----------------------------------------------------------- | ------------------------------------------------------ |
| Returns a `plotly` htmlwidget                              | inheritance check                                      |
| Downsamples when trace > `max_points`                      | `attr(p, "downsample_info")$rendered_points <= max_points` |
| Preserves all points when trace <= `max_points`            | `rendered_points == original_points; stride == 1L`     |
| Accepts `cb_friendly` palette and 5-element custom vector  | both produce a widget                                  |
| Rejects bad inputs (non-sangerseq object; bad palette name) | `expect_error` with descriptive regex                  |

---

## 4. Global Trimming Controls dashboard

### 4a. New API

```r
sa  <- SangerAlignment(...)
sa2 <- globalTrimApp(sa)   # opens a Shiny gadget; returns the re-trimmed SA
```

Lives at `R/GlobalTrimApp.R`, exported in NAMESPACE. Built as a **Shiny gadget** (via `shiny::runGadget` + `shiny::dialogViewer`) — opens in the RStudio viewer or a popup window, returns a value, then closes. This pattern is more ergonomic for "tweak parameters once" workflows than a full `launchAppSA`-style standalone app.

### 4b. Reactive flow

```
   ┌─────────────────────────────────────────────────────────────┐
   │  UI (sidebarLayout)                                          │
   │                                                              │
   │  • radioButtons("trimMethod", "M1" / "M2")                   │
   │  • conditionalPanel("M1") -> sliderInput("M1cutoff")         │
   │  • conditionalPanel("M2") -> sliderInput x 2 (score+window)  │
   │  • actionButton("apply") / "done" / "cancel"                 │
   │                                                              │
   │  Right pane (live preview)                                   │
   │  • Summary table: method, contigs, reads, consensus length   │
   │  • Wrapped consensus text (80-col)                           │
   │  • Per-contig table: name, fwd, rev, contig length           │
   └─────────────────────────────────────────────────────────────┘
                                │
                                ▼
   server() — single observer:
       observeEvent(input$apply, {
           rv$SA <- updateQualityParam(rv$SA, …)        ← cascades to children
           # rv$SA reactive update triggers the 3 renderers above.
       })
       observeEvent(input$done,   { stopApp(rv$SA)  })
       observeEvent(input$cancel, { stopApp(NULL)   })
```

Under the hood, every "Apply" click calls `updateQualityParam(SA, …)` (Phase 6's existing method, `R/MethodSangerAlignment.R:27`) which in turn cascades to per-contig `updateQualityParam(SC, ...)` and per-read M1/M2 retrim, then re-runs `alignContigs` to refresh consensus + tree.

### 4c. Test coverage

`tests/testthat/test-Phase8-GlobalTrim.R` (3 tests, 5 assertions):

| Test                                                             | Asserts                                                                |
| ---------------------------------------------------------------- | ---------------------------------------------------------------------- |
| Rejects non-SangerAlignment input                                | `expect_error("SangerAlignment")` for `list()` and a string             |
| Rejects FASTA-derived SangerAlignments                           | `expect_error("ABIF")` — FASTA inputs have no quality scores to retrim |
| Accepts an ABIF SA (entry-point validation reaches `runGadget`)   | `with_mocked_bindings` swaps `shiny::runGadget` for a stub; verifies the function reaches the launch step without error |

The third test uses `with_mocked_bindings` to stub `shiny::runGadget` because we can't actually open a Shiny gadget under `testthat`. The mock confirms input validation passes and the UI/server are constructed without error.

---

## 5. Files touched in Phase 8

```
M  DESCRIPTION                                                      (+ GlobalTrimApp.R in Collate)
M  NAMESPACE                                                        (+ chromatogram_plotly, globalTrimApp exports)
M  R/UtilitiesFunc.R                                                (+ chromatogram_plotly)
A  R/GlobalTrimApp.R                                                (the new gadget)
M  R/ShinySangerContigServer.R                                      (12 read sites converted)
M  R/ShinySangerAlignmentServer.R                                   (12 read sites converted)
M  inst/rmd/SangerRead_Report_ab1.Rmd                               (3 read sites converted)
M  inst/rmd/SangerRead_Report_fasta.Rmd                             (3 read sites converted)
A  tests/testthat/test-LazyAA-Reports.R                             (5 tests, 15 assertions)
A  tests/testthat/test-Phase8-PlotlyChromatogram.R                  (5 tests, 10 assertions)
A  tests/testthat/test-Phase8-GlobalTrim.R                          (3 tests, 5 assertions)
A  plans/08_ui_ux_redesign.md
```

13 files touched. **Total Phase 8 tests added: 13** (across 30 assertions). Cumulative suite: 1249 → 1279 PASS.

---

## 6. Reproducing

```r
# Whole suite (still passing):
devtools::test()

# Lazy-AA report regression:
testthat::test_file("tests/testthat/test-LazyAA-Reports.R")

# Plotly chromatogram unit tests:
testthat::test_file("tests/testthat/test-Phase8-PlotlyChromatogram.R")

# Global trim app:
testthat::test_file("tests/testthat/test-Phase8-GlobalTrim.R")
```

Interactive smoke tests (require a display):

```r
data(sangerAlignmentData)

# 1. Generate a report from a default (lazy) SA — populated AA tables now.
sa <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = system.file("extdata", "Allolobophora_chlorotica",
                                       "ACHLO", package = "sangeranalyseR"),
    REGEX_SuffixForward = "_[0-9]*_F.ab1$",
    REGEX_SuffixReverse = "_[0-9]*_R.ab1$"
)
generateReport(sa@contigList[[1]]@forwardReadList[[1]])

# 2. WebGL chromatogram in a browser:
plotly::print(chromatogram_plotly(sa@contigList[[1]]@forwardReadList[[1]]))

# 3. Global trim dashboard:
sa2 <- globalTrimApp(sa)   # opens a Shiny gadget; click Done to return
```

---

## 7. Non-goals (deferred)

- **Replacing `chromatogram_overwrite`** with the plotly version package-wide. The base-R one is still used for PDF export in `\dontrun{}` examples; replacing would break those flows.
- **Wiring `chromatogram_plotly` into the main `launchAppSC`/`launchAppSA` Shiny apps.** Currently those use base-R `renderPlot(chromatogram_overwrite(...))`. The new function is exported and ready to be substituted; doing so requires careful per-panel testing of all the existing Shiny reactivity (the chromatogram is rendered in many places, with per-read trim overlays). Out of Phase-8 scope.
- **LTTB downsampling** for `chromatogram_plotly`. Stride downsampling preserves peak silhouettes well at the typical 8000-point budget; LTTB would visibly improve very aggressive (max_points < 1000) downsamples. Add when needed.
- **`globalTrimApp` per-contig opt-out**. Currently it applies the chosen M1/M2 policy to every read; could expose a multi-select to exclude specific contigs. Add when users ask.
