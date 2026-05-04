<!-- badges: start -->
[![R-CMD-check](https://github.com/roblanf/sangeranalyseR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/roblanf/sangeranalyseR/actions/workflows/R-CMD-check.yaml)
[![BioC release](https://bioconductor.org/shields/build/release/bioc/sangeranalyseR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/sangeranalyseR/)
[![BioC devel](https://bioconductor.org/shields/build/devel/bioc/sangeranalyseR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/sangeranalyseR/)
[![codecov](https://codecov.io/gh/roblanf/sangeranalyseR/branch/devel/graph/badge.svg)](https://codecov.io/gh/roblanf/sangeranalyseR)
[![License: GPL-2](https://img.shields.io/badge/License-GPL--2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
![platform](https://img.shields.io/badge/platform-macOS_%7C_Linux_%7C_Windows-green.svg)
[![R-version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

# sangeranalyseR

**Fast, flexible, and reproducible workflows for assembling Sanger sequencing data into contigs in R.** A free and open-source alternative to Geneious, CodonCode Aligner, and Phred-Phrap-Consed.

For full documentation see **📒 [sangeranalyseR Documentation](https://sangeranalyser.readthedocs.io/en/latest/)**.

---

## What's new (devel)

The current development branch ships with a major rewrite of the assembly engine. **End-to-end `SangerAlignment(...)` builds are now ~1.7× faster** on the bundled fixtures, with full backward compatibility on the public API.

| Feature                              | What it does                                                                                                                                | Impact                                                                                  |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| **`BiocParallel` support**           | The per-`SangerRead` construction loop now uses `bplapply` instead of `parallel::mclapply`. Pass `BPPARAM = bpparam()` (or any `BiocParallelParam`) to choose your backend. | Multicore on Linux / macOS, automatic `SnowParam` on Windows — first-class parallelism cross-platform. |
| **Lazy AA translation (`lazyAA`)**   | The 3-frame `Biostrings::translate` step is skipped at construction time when no `refAminoAcidSeq` is supplied. AA frames compute lazily via `primaryAASeqS1/S2/S3()` accessors. | Removes ~35% of construction wall time on protein-coding reads.                          |
| **Rcpp `peakvalues` port**           | The peak-detection inner loop in `MakeBaseCallsInside` is now C++ via `peakvalues_batch_cpp` (one `.Call` per channel instead of one per peak). | ~1.28× end-to-end speedup on the ACHLO fixture; ~320 ms saved per `SangerAlignment` build. |
| **Plotly + WebGL chromatograms**     | New exported `chromatogram_plotly()` renders Sanger traces as Plotly htmlwidgets using `scattergl` (WebGL) with automatic stride downsampling. | Embeds chromatograms in Shiny dashboards without browser-freeze on > 50 k point traces. |
| **Global Trimming gadget**           | New exported `globalTrimApp(SA)` opens a Shiny gadget with M1 / M2 sliders that re-trim every read in a `SangerAlignment` and live-preview the consensus. | Faster batch parameter tuning than per-read sliders.                                    |
| **S4 `setValidity` invariants**      | Post-construction sanity checks on `QualityReport`, `ChromatogramParam`, and `ObjectResults`.                                                 | Catches silent slot mutations.                                                          |
| **Strict build compliance**          | `R CMD check` is fully clean (0 errors / 0 warnings / 0 notes). 1360 testthat tests pass, with coverage > 87% on every non-Shiny file.       | Production-quality build.                                                               |

Cumulative wall-time progression (8-read ACHLO fixture, mean of 5 reps, single thread):

| Milestone                                | Wall time | vs. baseline |
| ---------------------------------------- | --------: | -----------: |
| Pre-refactor baseline (eager AA, R-only) |  1.85 s   | 1.0×         |
| Lazy AA + BiocParallel plumbing          |  1.33 s   | 1.39×        |
| **Rcpp `peakvalues_batch_cpp`**          | **1.07 s (best)** / 1.14 s (mean) | **~1.62× best / 1.62× mean** |

See `plans/05_e2e_validation_report.md`, `plans/06_scaling_summary.md`, and `plans/07_rcpp_optimization_log.md` for the full benchmark methodology and raw artifacts.

---

## Installation

### From Bioconductor (recommended)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Stable release
BiocManager::install("sangeranalyseR")

# Development branch (latest features)
BiocManager::install(version = "devel")
BiocManager::install("sangeranalyseR")
```

### From GitHub

```r
# install.packages("devtools")    # if needed
devtools::install_github("roblanf/sangeranalyseR", ref = "devel")
```

### System requirements

- **R ≥ 4.0.0** (kept intentionally permissive for institutional installs).
- macOS, Linux, or Windows.
- C++17 toolchain for the Rcpp module — pre-installed on macOS (Xcode CLT), Linux (`build-essential`), and Windows (Rtools).
- Optional: `pandoc` for HTML report rendering.

---

## Quick start

A four-step end-to-end example using the bundled `Allolobophora chlorotica` fixture (8 ABIF files, 4 contigs).

### 1. Load and assemble

```r
library(sangeranalyseR)

# Locate the bundled fixture
ab1_dir <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                       package = "sangeranalyseR")

# Build the alignment — uses lazy AA + BiocParallel by default
sa <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = ab1_dir,
    REGEX_SuffixForward = "_[0-9]*_F.ab1$",
    REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
    TrimmingMethod      = "M1",
    M1TrimmingCutoff    = 0.0001,
    BPPARAM             = BiocParallel::bpparam()
)

sa@objectResults@creationResult       # TRUE
length(sa@contigList)                  # 4 contigs
length(sa@contigsConsensus)            # cross-contig consensus length
```

### 2. Tweak trimming parameters interactively

```r
# Open the global trim dashboard — sliders for M1 / M2 with live preview.
# Returns the re-trimmed SangerAlignment when the user clicks "Done".
sa2 <- globalTrimApp(sa)
```

### 3. Inspect a chromatogram in WebGL (no browser freeze on long traces)

```r
sr <- sa@contigList[[1]]@forwardReadList[[1]]

chromatogram_plotly(sr, max_points = 8000, showtrim = TRUE)
```

### 4. Export and report

```r
# FASTA: dispatches across SangerRead / SangerContig / SangerAlignment
writeFasta(sa)

# HTML report (requires pandoc) — works on lazy-AA objects out of the box
generateReport(sa)
```

### Lazy AA accessors

When no AA reference is supplied, `primaryAASeqS{1,2,3}()` compute on demand:

```r
sr <- sa@contigList[[1]]@forwardReadList[[1]]
length(sr@primaryAASeqS1)    # 0  — slot is empty under lazyAA = TRUE
primaryAASeqS1(sr)            # AAString returned by accessor (computed on call)
```

To restore eager translation (e.g. for the legacy direct-slot pattern):

```r
sr_eager <- SangerRead(
    inputSource    = "ABIF",
    readFeature    = "Forward Read",
    readFileName   = file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1"),
    TrimmingMethod = "M1",
    lazyAA         = FALSE
)
sr_eager@primaryAASeqS1       # populated at construction time
```

### Cross-platform parallel

```r
# macOS / Linux: forks via MulticoreParam
sa <- SangerAlignment(..., BPPARAM = BiocParallel::MulticoreParam(workers = 4))

# Windows: cluster-of-processes via SnowParam
sa <- SangerAlignment(..., BPPARAM = BiocParallel::SnowParam(workers = 4))

# Or just register a default once and forget about it:
BiocParallel::register(BiocParallel::SerialParam())
sa <- SangerAlignment(...)    # picks up the registered backend
```

---

## Citation

If `sangeranalyseR` is useful in your published work, please cite:

> **Kuan-Hao Chao, Kirston Barton, Sarah Palmer, and Robert Lanfear (2021).** *sangeranalyseR: simple and interactive processing of Sanger sequencing data in R.* Genome Biology and Evolution. DOI: [10.1093/gbe/evab028](https://doi.org/10.1093/gbe/evab028).

Available on [Genome Biology and Evolution (GBE)](https://academic.oup.com/gbe/advance-article/doi/10.1093/gbe/evab028/6137837?guestAccessKey=a28b32d6-ffab-41f2-8132-9c2dd28b99fe) and [Bioconductor](https://bioconductor.org/packages/release/bioc/html/sangeranalyseR.html).

---

## Maintainers

- **Kuan-Hao Chao** &lt;ntueeb05howard@gmail.com&gt; (creator, maintainer)
- **Rob Lanfear** &lt;rob.lanfear@gmail.com&gt; (author)

License: **GPL-2** (see `LICENSE`).

Issues and feature requests: [github.com/roblanf/sangeranalyseR/issues](https://github.com/roblanf/sangeranalyseR/issues).
