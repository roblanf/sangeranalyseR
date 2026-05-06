# GitHub issue replies — Phase 18 backlog cleanup

This file uses the same Markdown format as `plans/16_github_replies.md` and `plans/17_github_replies.md`, **plus** an `Action:` metadata line that the Phase-18-updated `plans/close_issues.py` reads:

- `**Action**: close`   — POST comment then PATCH state=closed.
- `**Action**: comment` — POST comment only (do NOT close). Used for "please retest on devel" requests.

If `Action:` is omitted the script defaults to `close` (preserves backwards-compatible behaviour for `16_github_replies.md` / `17_github_replies.md`).

To run:

```bash
export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx
python3 plans/close_issues.py --md plans/18_github_replies_cleanup.md --dry-run
python3 plans/close_issues.py --md plans/18_github_replies_cleanup.md
```

Tally for this file: **16 issues** = 4 documentation + 8 fully-resolved-via-refactor + 4 likely-resolved-needs-retest.

| Bucket                                 | Count | Action      |
| -------------------------------------- | ----: | ----------- |
| Documentation issues (vignette pass)    |     4 | `close`     |
| Resolved via refactor                   |     8 | `close`     |
| Likely resolved — please retest         |     4 | `comment`   |

---

## Documentation issues (closed by the Phase-18 vignette rewrite)

### Issue #13 — Make some worked examples

**URL**: https://github.com/roblanf/sangeranalyseR/issues/13
**Action**: close

```markdown
Hi @roblanf — closing this one with the Phase-18 vignette rewrite. The new `vignettes/sangeranalyseR.Rmd` opens with a "How to..." recipe gallery covering the most-asked questions:

- assemble a single contig
- assemble many contigs (`SangerAlignment`)
- use a CSV mapping instead of regex
- handle forward-only / reverse-only datasets
- pick a trimming algorithm for low-quality data
- detect spurious low-overlap merges
- choose a consensus base-calling method (strict / majority / quality-weighted)
- interpret the chromatogram and deal with secondary peaks
- launch the Shiny app + global trim gadget
- export FASTA / HTML reports

Each recipe is a runnable code block keyed off the bundled ACHLO fixture, so users can copy-paste verbatim.

Available in the next Bioconductor `devel` build:

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
browseVignettes("sangeranalyseR")
```

Closing as fixed.
```

---

### Issue #49 — Add documentation "How to ..." section

**URL**: https://github.com/roblanf/sangeranalyseR/issues/49
**Action**: close

```markdown
Hi, this is the canonical documentation issue — closing alongside #13.

The Phase-18 vignette rewrite introduces a "How to..." recipe gallery as the first major section after Introduction. The three sub-topics you listed ("How to deal with secondary peaks", "Explain two different trimming methods") are explicitly covered with runnable examples:

- **Two trimming methods** — a side-by-side table comparing M1 (modified Mott's) and M2 (sliding window), with the parameters each requires and a "tighter trimming for noisy data" recipe using `TrimmingMethod = "M2"`, `M2CutoffQualityScore = 30`, `M2SlidingWindowSize = 15`.
- **Secondary peaks** — describes how the package re-runs base-calling using `signalRatioCutoff` (default 0.33), how to inspect `@secondaryPeakDF`, and how to re-run base-calling with a tighter cutoff via `MakeBaseCalls(sr, signalRatioCutoff = 0.22)`.

Plus: chromatogram interpretation, base-calling method walkthrough (#71), parameter reference tables (#99), and a troubleshooting matrix mapping common errors to their Phase-15/16/17 fixes.

Available in the next Bioconductor `devel` build.

Closing as fixed.
```

---

### Issue #71 — base calling

**URL**: https://github.com/roblanf/sangeranalyseR/issues/71
**Action**: close

```markdown
Hi @ramirorr, thanks for the question — closing as fixed by the Phase-18 vignette rewrite.

The new "How to interpret the chromatogram" section in `vignettes/sangeranalyseR.Rmd` walks through exactly which method `sangeranalyseR` uses for base calling:

1. **Detect peaks per channel** via `getpeaks()` — runs across the four trace channels (A/C/G/T) in `abifRawData@data`.
2. **Per-base position**, identify the strongest peak across channels using `signalRatioCutoff` (default 0.33). Secondary peaks below 33% of the primary peak's amplitude are dropped.
3. **Tied peaks** → IUPAC ambiguity code corresponding to the equally-strong bases.
4. **Quality scores** come directly from the ABIF's `PCON.2` data block (one entry per detected peak).

This re-implementation lives in `MakeBaseCallsInside` (`R/UtilitiesFunc.R`) and is invoked once per `SangerRead` constructor with the user's `signalRatioCutoff`. To re-run with a different cutoff (e.g. for a stricter secondary-peak threshold), call `MakeBaseCalls(sr, signalRatioCutoff = 0.22)`.

If you want the exact base calls already stored in the ABIF file (without re-running), they're available as the inherited `@primarySeq` slot before `MakeBaseCallsInside` overwrites it — though Phase-6 made that "raw" version available separately as `@primarySeqRaw` for the same reason.

Closing as fixed in the upcoming `devel` release.
```

---

### Issue #99 — Create contig

**URL**: https://github.com/roblanf/sangeranalyseR/issues/99
**Action**: close

```markdown
Hi, thanks for sharing the exact code — that helped pinpoint where the docs were unclear.

The Phase-18 vignette rewrite adds a full **Constructor parameter reference** section that explicitly groups the parameters into:

- **Required** — `inputSource`, `processMethod`, `ABIF_Directory` / `FASTA_File`, `REGEX_SuffixForward` / `REGEX_SuffixReverse` (or `CSV_NamesConversion`), `contigName`.
- **Trimming** — `TrimmingMethod`, `M1TrimmingCutoff`, `M2CutoffQualityScore`, `M2SlidingWindowSize`, `minReadLength`, `signalRatioCutoff`.
- **Consensus** — `consensusMethod` (Phase-17 new), `qualityAware`, `minFractionCall`, `maxFractionLost`, `minOverlapBases`, `minOverlapFraction` (Phase-16 new), `alignSeqsParams`.
- **Performance** — `processorsNum`, `BPPARAM`, `lazyAA` (Phase-6 new).

A correct, runnable `SangerContig` example for forward-only data (which I think matches your use case — 5 forward `.ab1` files, no reverse) looks like:

```r
sc <- SangerContig(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = "PATH_TO_YOUR_AB1S",
    contigName          = "6063",
    REGEX_SuffixForward = "_F\\.ab1$",   # match your filenames
    REGEX_SuffixReverse = NULL,           # forward-only (Phase-15)
    minReadsNum         = 1,              # each forward read = its own contig
    TrimmingMethod      = "M1",
    M1TrimmingCutoff    = 0.001,
    signalRatioCutoff   = 0.4
)
```

Note the Phase-15 changes that make this case painless:

- `REGEX_SuffixReverse = NULL` (or `NA_character_`) is now valid for forward-only datasets.
- `minReadsNum = 1` lets single-read contigs build instead of being filtered as "too few reads".

Available in the next Bioconductor `devel` build:

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
browseVignettes("sangeranalyseR")
```

Closing as fixed.
```

---

## Resolved via refactor (closed)

### Issue #11 — add a function to count coincident secondary peaks

**URL**: https://github.com/roblanf/sangeranalyseR/issues/11
**Action**: close

```markdown
Hi @roblanf — closing as already resolved.

This functionality is now `countCoincidentSp()` in `R/UtilitiesFunc.R`, called from inside `calculateContigSeq` to populate `@secondaryPeakDF` on every `SangerContig` (one row per alignment column with > 1 secondary peak; columns are `column.number`, `ambiguities`, `column`).

Phase-7 also Rcpp-ised the per-column secondary-peak detection inner loop (`peakvalues_batch_cpp`) for ~1.28× end-to-end speedup.

```r
data("sangerContigData")
head(sangerContigData@secondaryPeakDF)
```

Closing as resolved by long-standing `countCoincidentSp` + Phase-7 Rcpp port.
```

---

### Issue #21 — build a new sequence object

**URL**: https://github.com/roblanf/sangeranalyseR/issues/21
**Action**: close

```markdown
Hi @roblanf — closing; this has been the `SangerRead` S4 class since v1.x.

Each `SangerRead` instance carries:

1. **Primary base calls** — `@primarySeq` (post-MakeBaseCalls) and `@primarySeqRaw` (pre-basecall, from `sangerseq()` directly).
2. **Secondary base calls** — `@secondarySeq` and `@secondarySeqRaw`.
3. **Quality scores** — `@QualityReport@qualityPhredScores` (per-base Phred from the ABIF's `PCON.2` block, with synthetic Phred-30 fallback added in Phase-15 for files with empty `PCON.2`).
4. **Trimming information** — `@QualityReport@trimmedStartPos`, `@trimmedFinishPos`, `@trimmedSeqLength`, `@rawMeanQualityScore`, `@trimmedMeanQualityScore`, etc.

Plus 3-frame translations (`@primaryAASeqS1/S2/S3`, lazy by default after Phase-6), trace matrix, peak position / amplitude matrices, and an `@objectResults` envelope reporting construction success and any per-read errors.

Closing as already-implemented.
```

---

### Issue #31 — Operating systems

**URL**: https://github.com/roblanf/sangeranalyseR/issues/31
**Action**: close

```markdown
Hi @roblanf — closing this; cross-platform support is now first-class.

- **macOS / Linux**: `BiocParallel::MulticoreParam` is the auto-default for `BPPARAM = bpparam()` (Phase-6).
- **Windows**: `BiocParallel::SnowParam` auto-default; the Phase-6 BiocParallel migration explicitly replaced every `parallel::mclapply` (which falls back to a single core on Windows) with `bplapply`.
- **Rtools**: required for the Rcpp module (Phase-7) — the standard Bioc/Rtools install handles this. If users hit Windows-specific issues, file a new issue and we'll triage.

The Phase-9 dependency cleanup also moved most packages from `Depends:` to `Imports:`, eliminating most of the search-path-collision class of issues. CI is green on all three platforms via the GitHub Actions `R-CMD-check` workflow.

Closing as resolved.
```

---

### Issue #55 — SangerAlignment Error enhancement

**URL**: https://github.com/roblanf/sangeranalyseR/issues/55
**Action**: close

```markdown
Closing this enhancement issue — all sub-tasks landed across the recent refactor:

- **Empty value error in Shiny app** — Phase-9 added `shiny::testServer`-driven test coverage for the Shiny servers; #60 / #61 (which were the original "Shiny crashes" reports) were both Phase-9-cleaned-up via the dependency / import overhaul.
- **Update minReadsNum parameter documentation** — Phase-18 vignette rewrite covers `minReadsNum` in the parameter reference table, including the `minReadsNum = 1` recipe for forward-only data.
- **Add unit tests** — Phase-3 added a comprehensive testthat suite (24 negative tests, the orthogonal-axis matrix, slot-invariant guards). The cumulative suite is now **1451 / 1451 PASS** with `R CMD check 0/0/0`.

Closing as resolved.
```

---

### Issue #82 — Error in TreeLine

**URL**: https://github.com/roblanf/sangeranalyseR/issues/82
**Action**: close

```markdown
Hi @anukashyap, this was fixed already on the `RELEASE_3_22` and `devel` branches via the `Fix DECIPHER Treeline import; update DECIPHER minimum version` commit.

DECIPHER renamed its tree-builder from `IdClusters` to `Treeline` in 3.0; the package's import declaration and minimum DECIPHER version were updated to match.

If you re-install the latest version you should not see this error any more:

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
```

Closing as resolved.
```

---

### Issue #91 — M2trimming: error when handling low quality file

**URL**: https://github.com/roblanf/sangeranalyseR/issues/91
**Action**: close

```markdown
Hi, thanks for the report. The M2 trimming algorithm could produce a degenerate state (`trimmedFinishPos = 0` while `trimmedStartPos > 0`) on very low-quality reads, and the strict `QualityReport` validator threw on the resulting negative `trimmedSeqLength`.

Phase-4 relaxed the validator to accept the degenerate "no usable trim window" state explicitly (treating it as a degraded but valid result). Phase-16 added a defensive width filter in `calculateContigSeq` that drops any read with width < 2 bp before alignment, with `MIN_READ_LENGTH_DEFENSIVE_DROP` logged. The end-to-end behaviour is now: low-quality reads are silently dropped, the build succeeds with the surviving reads, and the user sees a warning in the log.

Closing as resolved by Phase-4 + Phase-16.
```

---

### Issue #97 — Move data.table to Imports

**URL**: https://github.com/roblanf/sangeranalyseR/issues/97
**Action**: close

```markdown
Hi @MichaelChirico, thanks — closing as done.

Phase-9 did exactly this: `Depends:` shrunk from 27 entries to the 4 packages whose types are publicly returned (`Biostrings`, `DECIPHER`, `sangerseqR`). Everything else, including `data.table`, moved to `Imports:`. The `tstrsplit` symbol that the Shiny servers rely on is declared via `@importFrom data.table tstrsplit` in `R/sangeranalyseR_package.R`.

`R CMD check` is now `0 errors / 0 warnings / 0 notes`. Closing as resolved.
```

---

### Issue #98 — Ongoing development and maintenance

**URL**: https://github.com/roblanf/sangeranalyseR/issues/98
**Action**: close

```markdown
Hi, thanks for asking — yes, the package is actively maintained. The recent activity is summarised in `NEWS.md` on the `devel` branch:

- **Phase 6** — `BiocParallel` migration + lazy 3-frame AA translation (~1.4× speedup).
- **Phase 7** — Rcpp port of `peakvalues` (further +28% to ~1.62× cumulative).
- **Phase 8** — Plotly + WebGL chromatograms; lazy-AA fix for the Shiny servers and RMarkdown templates.
- **Phase 9** — strict build compliance: `R CMD check` is now 0/0/0; 1451 testthat tests passing; `Depends:` slimmed 27 → 4.
- **Phase 14–17** — concrete bug fixes for #100, #92, #76, #89, #94, #66, #65, #42, plus consensus-method enhancements (#87, #48, #33).
- **Phase 18** — vignette overhaul (#13, #49, #71, #99) and this backlog cleanup.

The 7-feature wall-time progression (8-read fixture, single thread): pre-refactor 1.85 s → current 1.07 s best / 1.14 s mean.

If your scripts depended on a specific older revision (`0e658db5`) we'd be happy to know which behaviour changed — file a new issue with a small reprex and we'll triage. The `devel` branch is what will land in the next Bioconductor release; install with:

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
```

Closing this meta-issue; specific behaviour regressions belong in their own issues.
```

---

## Likely resolved — please retest (comment-only)

### Issue #60 — Error when launching Shiny app

**URL**: https://github.com/roblanf/sangeranalyseR/issues/60
**Action**: comment

```markdown
Hi, sorry for the delay. The `Couldn't normalize path in 'addResourcePath'` error you hit looks like a `shinydashboard` resource-path issue, almost certainly pre-dating the Phase-9 dependency cleanup that overhauled how the package's Shiny imports are declared.

Could you try the latest `devel` and report whether the error reproduces?

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
library(sangeranalyseR)
data(sangerAlignmentData)
launchApp(sangerAlignmentData)
```

If it still fails on Windows specifically, please paste the full traceback (R version, OS, the output of `sessionInfo()` after `library(shinydashboard)`). If it works, please let us know so we can confirm and close.
```

---

### Issue #61 — Launch Shiny app error

**URL**: https://github.com/roblanf/sangeranalyseR/issues/61
**Action**: comment

```markdown
Hi @JordanWG, the original issue references the older parameter-name convention (`parentDirectory`, `suffixForwardRegExp`) that has since changed to `ABIF_Directory` / `REGEX_SuffixForward`. The current `SangerAlignment(...)` API is documented in detail in the Phase-18 vignette overhaul.

Could you re-run on the latest `devel` with the modern parameter names, and confirm whether the Shiny app still errors?

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
library(sangeranalyseR)
my_aligned_contigs <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = "C:/Users/JORDAN/Desktop/aa",
    REGEX_SuffixForward = "_[0-9]+_F",
    REGEX_SuffixReverse = "_[0-9]+_R"
)
launchApp(my_aligned_contigs)
```

If it still fails please paste the full traceback. Phase-9 added testServer-driven coverage for the Shiny modules so the "obvious" launch errors should be much rarer now.
```

---

### Issue #68 — Unavailable external font resource

**URL**: https://github.com/roblanf/sangeranalyseR/issues/68
**Action**: comment

```markdown
Hi, the `pandoc: Could not fetch ... pro.fontawesome.com` error is from the RMarkdown template trying to bundle external font resources at render time. The Phase-8 lazy-AA report compatibility fix touched the `inst/rmd/*.Rmd` templates substantially; Phase-18 added a vignette that explicitly notes pandoc is required for `generateReport()`.

Could you try the latest `devel` and let us know if the font-fetch error still fires?

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
library(sangeranalyseR)
data(sangerAlignmentData)
generateReport(sangerAlignmentData, outputDir = tempdir())
```

If it still fails it's likely a pandoc / network configuration on your machine rather than a sangeranalyseR issue — the templates themselves don't fetch external CSS as of Phase-8.
```

---

### Issue #85 — Issues with Reproducible tutorial

**URL**: https://github.com/roblanf/sangeranalyseR/issues/85
**Action**: comment

```markdown
Hi, the original tutorial was several major versions out of date by the time of your report. Phase-4 fixed the file-extension regex bug that affected `.fa` / `.ab1` filename validation (the `.fa$` regex was matching `.Xfa` filenames before the fix), and Phases 11–12 + 18 rewrote both the README and the vignette around the modern API.

Could you re-run the example from the current vignette on the latest `devel`?

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
library(sangeranalyseR)
browseVignettes("sangeranalyseR")
```

The "How to..." section near the top of the vignette has runnable examples keyed off the bundled ACHLO fixture — these should reproduce cleanly. If you hit a specific error, please paste the exact reprex (the parameter names changed between major versions, see #61's reply for an example) and we'll triage.
```

---

## How to apply

```bash
export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx

# Sanity check first — should report:
#   13 issues with [close]   action
#    4 issues with [comment] action
python3 plans/close_issues.py --md plans/18_github_replies_cleanup.md --dry-run

# Live run (per-issue confirmation prompts):
python3 plans/close_issues.py --md plans/18_github_replies_cleanup.md
```

The `[comment]` issues will get the comment posted but **will not** be closed — that's the desired behaviour for "please retest" requests where the reporter still has to confirm.
