# sangeranalyseR 1.21.1 (development)

## Performance

- **`BiocParallel`-backed per-`SangerRead` construction loop.** Replaces all `parallel::mclapply` call sites with `BiocParallel::bplapply`. The constructor now accepts a `BPPARAM = bpparam()` argument and a `.resolveBPPARAM(processorsNum, BPPARAM)` helper maps the legacy integer `processorsNum` argument onto the right backend (`SerialParam` for 1 worker, `MulticoreParam` on Unix and `SnowParam` on Windows for ≥ 2 workers). Cross-platform parallel is now first-class.
- **Lazy 3-frame amino-acid translation (new `lazyAA = TRUE` default).** When no `refAminoAcidSeq` is supplied, `Biostrings::translate` is no longer called eagerly during `SangerRead` construction. The slots `primaryAASeqS{1,2,3}` start as empty `AAString`s; the new exported accessors `primaryAASeqS1()`, `primaryAASeqS2()`, `primaryAASeqS3()` compute on demand. Eliminates ~35 % of construction wall time on protein-coding ABIF reads.
- **Rcpp port of the peak-detection inner loop.** New `src/peakvalues.cpp` exports `peakvalues_cpp` and `peakvalues_batch_cpp`. The latter processes all peak windows for one channel in a single `.Call`, eliminating ~22,400 per-window R-to-C++ marshalling round-trips per `SangerAlignment` build. Saves ~320 ms per build on the bundled fixture (≈ 1.28× end-to-end).

  Cumulative end-to-end timings on the bundled `Allolobophora_chlorotica/ACHLO` fixture (8 reads, 4 contigs, mean of 5 repetitions, single thread):

  | Milestone                                        | Wall time | vs. baseline |
  | ------------------------------------------------ | --------: | -----------: |
  | Pre-refactor baseline (eager AA, R-only)         | 1.85 s    | 1.00×        |
  | Lazy AA + BiocParallel plumbing                  | 1.33 s    | 1.39×        |
  | **+ Rcpp `peakvalues_batch_cpp`**                | **1.07 s (best) / 1.14 s (mean)** | **~1.62× / 1.62×** |

  Full methodology and raw artifacts under `plans/05_e2e_validation_report.md`, `plans/06_scaling_summary.md`, and `plans/07_rcpp_optimization_log.md`.

- Removed redundant `mclapply` parallelism around `nPairwiseDiffs` and `oneAmbiguousColumn` — the per-element work was sub-millisecond and fork-setup overhead dominated. They are now serial `lapply` calls and ~5 % faster on small alignments.

## New features

- `chromatogram_plotly(obj, max_points = 8000, showtrim = FALSE, colors = "default")` — interactive Plotly htmlwidget rendering of Sanger chromatograms via `scattergl` (WebGL). Uniform-stride downsampling caps points-per-channel; the returned widget carries a `downsample_info` attribute reporting the original / rendered counts.
- `globalTrimApp(SA)` — Shiny gadget that exposes M1 / M2 trim sliders across an entire `SangerAlignment`. Each "Apply" click calls `updateQualityParam(SA, ...)` which cascades to every child read; live previews of consensus length, contig count, and per-contig stats update reactively. Returns the re-trimmed `SangerAlignment` on "Done".
- `primaryAASeqS1(sr)`, `primaryAASeqS2(sr)`, `primaryAASeqS3(sr)` — lazy AA accessors that return the cached slot when populated, otherwise compute on demand via the existing `calculateAASeq` helper.
- `BPPARAM` argument added to `SangerRead()`, `SangerContig()`, `SangerAlignment()`. Defaults to `NULL` (derived from `processorsNum`); pass any `BiocParallelParam` to override.
- `lazyAA` argument (default `TRUE`) on the same three constructors. Pass `lazyAA = FALSE` for the legacy direct-slot eager-translation behaviour.

## Robustness

- **File-extension regex bug fixes.** The pre-Phase-4 patterns `".fa$"` / `".fasta$"` / `".ab1$"` (in `checkFASTA_File` and `checkReadFileName`) used unescaped `.`, so files like `Sanger_all_reads.Xfa` and `Achl_006_F.Xab1` slipped through as valid. Now `\\.fa(sta)?$` and `\\.ab1$`. Centralised as `.FASTA_EXT_REGEX` and `.AB1_EXT_REGEX` constants.
- **Validator framework refactor.** All 31 `check*` functions in `R/UtilitiesFuncInputChecker.R` now route through 5 internal helpers (`.errAppend`, `.requireType`, `.requireEnum`, `.requireRange`, `.requireExt`). Public signatures and error-type tags (`PARAMETER_RANGE_ERROR`, `FILE_TYPE_ERROR`, etc.) are preserved exactly — no consumer change required. File shrunk 603 → 523 lines (-13 %).
- **S4 `setValidity` invariants** added on `QualityReport`, `ChromatogramParam`, and `ObjectResults`. Catch any code path (including `slot<-` mutations after construction) that would land an out-of-range value in the slots. Sanger* user-facing classes intentionally have no validity (preserves the "construction never throws" contract).
- **Lazy-AA report compatibility.** All 30 direct `@primaryAASeqS{1,2,3}` slot reads in the Shiny servers and RMarkdown report templates were converted to the accessor functions, so reports rendered against `lazyAA = TRUE` objects no longer produce empty AA tables. (The 12 `<<-` write sites are intentionally preserved — they're reactive caches.)
- Removed stray `cat()` / `print()` / `message()` debug statements from `R/UtilitiesFunc.R`, `R/UtilitiesFuncInputChecker.R`, `R/ShinyServerModule.R`, and `R/ShinySangerContigServer.R`. User-facing reporting now uses `log_info` consistently.

## Build / compliance

- `R CMD check` is **fully clean: 0 errors / 0 warnings / 0 notes** (was 0 / 5 / 7 before Phase 9).
- 1360 testthat tests, all passing. New test files added across phases:
  - `test-Validator-EdgeCases.R`, `test-Validator-Helpers.R`, `test-Regex-Boundary.R` (validator + regex regression).
  - `test-OrthogonalAxes.R` (input-source × process-method × trim-method matrix).
  - `test-SangerRead-SlotInvariants.R`, `test-S4Validity.R` (S4 invariants).
  - `test-Rcpp-peakvalues.R` (R / C++ equivalence + 200-trial fuzz).
  - `test-LazyAA-BiocParallel.R`, `test-LazyAA-Reports.R` (Phase 6 / 8).
  - `test-Phase8-PlotlyChromatogram.R`, `test-Phase8-GlobalTrim.R`, `test-Phase9-GlobalTrim-testServer.R` (UI).
  - `test-Phase10-Coverage.R` (coverage maximisation).
- Coverage measured by `covr::package_coverage()`: **35.7 %** overall, **> 87 %** on every non-Shiny R file (the three Shiny server files at 0 % require a real browser harness; deferred).
- `DESCRIPTION` modernised: `Authors@R` (replacing the deprecated `Author:` + `Maintainer:` pair), `URL`, `BugReports`, `License: GPL-2 | file LICENSE`, `LinkingTo: Rcpp`, R version dependency `>= 4.0.0` (intentionally permissive).
- `Depends:` slimmed from 27 entries to the 4 packages whose types are publicly returned (Biostrings, DECIPHER, sangerseqR); the rest moved to `Imports:` (or `Suggests:` for vignette-only deps).
- ASCII-only source: replaced curly apostrophes, em-dashes, arrows, and multiplication signs across `R/UtilitiesFunc.R`, `R/Class*.R` (was a `R CMD check` warning).
- Re-saved `data/*.RData` with xz compression (largest file 1.5 MB → 698 KB).

## Bug fixes

GitHub-issue cleanup across three resolution sprints (Phases 15–17):

- **#100 (CSV substring contig-name match):** `processCSV.R` now matches contig names exactly instead of via `grepl(name, …)`, so a contig named `good` is no longer accidentally absorbed into a contig named `good_extra`.
- **#92 (forward-only `NULL` handling):** `SangerContig()` accepts `REGEX_SuffixReverse = NULL` (or `NA_character_`) and `minReadsNum = 1` for forward-only / reverse-only datasets. Previously the constructor errored with `argument is of length zero`.
- **#76 (missing `PCON.2` block):** ABIF reads with empty `PCON.2` quality data now get a synthetic Phred-30 score per base with a `MISSING_QUALITY_SCORES_WARN`, instead of silently constructing an unusable read.
- **#89 (single-read `writeFasta`):** `writeFastaSC` no longer errors on `SangerContig`s built from a single read.
- **#94 (low-overlap detection):** new `minOverlapBases` (default 20) and `minOverlapFraction` (default 0.4) post-alignment guards in `calculateContigSeq`. Spurious low-overlap merges emit `LOW_OVERLAP_WARN` and the contig is rejected before it propagates to the alignment.
- **#66 (degenerate consensus):** IUPAC ambiguity-code handling in `ConsensusSequence(ambiguity = TRUE)` is correctly preserved; consumers that called `as.character()` on the consensus once again see ambiguity codes rather than `N` collapse.
- **#65 (multi-contig duplicate reads):** when a read is matched into more than one contig (CSV or REGEX), `SangerAlignment` now logs `READ_ASSIGNED_MULTIPLE_CONTIGS_WARN` and assigns the read to the first matching contig only.
- **#42 (length-1 reads):** reads of width < 2 bp are dropped at alignment time with `MIN_READ_LENGTH_DEFENSIVE_DROP`, allowing the contig to build from the surviving reads instead of failing the whole alignment.
- **#91 (M2 trimming on degraded reads):** the `QualityReport` validator now accepts the degenerate "no usable trim window" state (`trimmedFinishPos == 0` while `trimmedStartPos > 0`) on extremely low-quality reads.

## New features (Phase 17 — consensus algorithms)

- **`consensusMethod` argument** on `SangerContig()` and `SangerAlignment()` with three options:
  - `"strict"` (default; pre-Phase-17 behaviour) — IUPAC ambiguity codes preserved at disagreeing columns.
  - `"majority"` — most-frequent base wins per column; ties break alphabetically. No IUPAC codes in the output.
  - `"quality_weighted"` — votes weighted by per-base Phred from the source reads; falls back to flat Phred 30 (with a warning) when scores are missing or for FASTA inputs.
- **`qualityAware = TRUE`** is a shorthand for `consensusMethod = "quality_weighted"`.
- **Per-position consensus quality scores.** `attr(@contigSeq, "qualityScores")` is now an integer vector of length `length(contigSeq)` under `"majority"` and `"quality_weighted"` modes (empty `integer(0)` under `"strict"` for backwards compatibility). Closes the long-standing #87 / #48 / #33 cluster.

## Documentation (Phase 18)

- **Vignette overhaul.** `vignettes/sangeranalyseR.Rmd` rewritten end-to-end with a "How to..." recipe gallery (10 recipes covering single-contig assembly, CSV mapping, forward-only data, low-quality trimming, low-overlap detection, consensus methods, secondary peaks, Shiny launch, FASTA / HTML export), a constructor parameter reference (4 tables), a troubleshooting matrix mapping common errors to the Phases 15–17 fixes, and a `sessionInfo()` block. Closes #13, #49, #71, #99.
- **`R CMD check --run-donttest` hardening.** Multiple pre-existing latent example bugs were fixed:
  - `inst/rmd/SangerContig_Report.Rmd` now loads `library(knitr)` so `kable()` resolves during report rendering.
  - `readTable.SangerRead` and `readTable.SangerContig` examples no longer call `readTable(sangerAlignmentData)` (no method exists for `SangerAlignment`).
  - `globalTrimApp`, `launchApp`, `launchAppSC`, `launchAppSA` examples switched from `\donttest{}` to `\dontrun{}` so `runGadget()` / auto-printed `shiny.appobj` no longer hang the example runner.
- **Maintainer tooling.** `plans/close_issues.py` (Phase 16.5) now parses an `Action: close|comment` metadata flag from each issue's reply Markdown — comment-only entries (used for "please retest on devel" responses) skip the `state=closed` PATCH. Backwards-compatible with the existing Phase-16 / Phase-17 reply files.

# sangeranalyseR 1.20.0 (current Bioconductor release)

Maintenance release on the `RELEASE_3_22` branch:

- DECIPHER `Treeline` import fix (replaces the older API).
- Bumped DECIPHER minimum version.
- Standard Bioc release-cycle version bumps (even `y` on release, odd `y` on devel).

# sangeranalyseR 1.6.1

- Fix chromatogram colour issue.

# sangeranalyseR 0.99.1

- Base class `SangerReads` designed to store each forward / reverse read.

# sangeranalyseR 0.1.0

- Project starts.
