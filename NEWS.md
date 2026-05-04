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

# sangeranalyseR 1.20.0 (current Bioconductor release)

Maintenance release on the `RELEASE_3_22` branch:

- DECIPHER `Treeline` import fix (replaces the older API).
- Bumped DECIPHER minimum version.
- Standard Bioc release-cycle version bumps (even `y` on release, odd `y` on devel).

(Pre-1.20.0 history is in the legacy `NEWS` file in DCF format.)
