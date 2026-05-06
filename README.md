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

## What's new

- **~1.7× faster `SangerAlignment(...)`** thanks to a C++ peak-detection inner loop, parallel per-read construction (`BiocParallel`), and lazy 3-frame amino-acid translation that only runs when you ask for it.
- **Interactive Plotly + WebGL chromatograms** via the new `chromatogram_plotly()` — smooth scrolling and zoom on Sanger traces with tens of thousands of points.
- **Global trimming dashboard** via the new `globalTrimApp(sa)` — adjust M1 / M2 trimming parameters across an entire `SangerAlignment` with a live consensus preview.

For the full per-version changelog see [`NEWS.md`](NEWS.md).

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

### 3. Explore your data

Open the per-read Shiny app for a `SangerContig` (or use `launchApp(sa)` on a full `SangerAlignment`):

```r
launchApp(sa)
```

<img src="https://i.imgur.com/gwY6AqB.png" alt="sangeranalyseR Shiny app — interactive contig browser" style="width:100%">

You can also pull up an interactive WebGL chromatogram for a single read without launching the full app:

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
- **Rob Lanfear** &lt;rob.lanfear@gmail.com&gt; (creator, author)

License: **GPL-2** (see `LICENSE`).

Issues and feature requests: [github.com/roblanf/sangeranalyseR/issues](https://github.com/roblanf/sangeranalyseR/issues).
