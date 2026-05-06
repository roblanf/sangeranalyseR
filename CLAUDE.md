# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

This is a Bioconductor R package. All commands assume the working directory is the package root.

Build / check (from a shell):
- `R CMD build .` — build the source tarball.
- `R CMD check sangeranalyseR_*.tar.gz` — standard R CMD check.
- `R CMD BiocCheck sangeranalyseR_*.tar.gz` — Bioconductor-specific lint (required before submitting to Bioc).

Interactive development (from R, with `devtools` installed):
- `devtools::load_all()` — load all R/ files for iteration.
- `devtools::document()` — regenerate `man/*.Rd` and `NAMESPACE` from roxygen2 comments. **Always run after editing roxygen comments**; `NAMESPACE` is auto-generated ("do not edit by hand").
- `devtools::test()` — run the full testthat suite.
- `testthat::test_file("tests/testthat/test-Constructors.R")` — run a single test file. Most test files share setup loaded from `tests/testthat/helper-*.R` (testthat auto-sources these).
- `devtools::check()` — full local check.

Note: `tests/testthat.R` has `test_check("sangeranalyseR")` **commented out**, so simply running that script will not execute the suite — use `devtools::test()` or `testthat::test_dir("tests/testthat")`.

Branch convention (Bioconductor):
- `master` is devel; `RELEASE_3_2X` branches are stable Bioconductor releases.
- The `Version:` field in `DESCRIPTION` follows Bioc parity rules: odd `y` on devel, even `y` on release. Don't bump versions unless mirroring a Bioc cycle event — recent commits (`bump x.y.z version to even y prior to creation of RELEASE_3_22 branch`, etc.) show the pattern.

## Architecture

### S4 class hierarchy (the data model)

There are three primary user-facing S4 classes forming a containment hierarchy:

- **`SangerRead`** (`R/ClassSangerRead.R`) — one Sanger read (one ABIF file or one FASTA entry). Extends `sangerseq` from the `sangerseqR` package. Holds the raw chromatogram, base-called sequences, and a nested `QualityReport` and `ChromatogramParam`.
- **`SangerContig`** (`R/ClassSangerContig.R`) — a contig assembled from one or more forward + reverse `SangerRead`s. Holds two lists (`forwardReadList`, `reverseReadList`), the consensus sequence, and the alignment.
- **`SangerAlignment`** (`R/ClassSangerAlignment.R`) — a set of `SangerContig`s plus their cross-contig alignment, consensus, and `phylo` tree.

Helper classes (`R/ClassQualityReport.R`, `R/ClassChromatogramParam.R`, `R/ClassObjectResults.R`) are nested slots, not standalone user objects. `ObjectResults` records construction success/errors so that creation never throws — it returns an object whose `objectResults` slot describes what failed. Many `setClassUnion` "...ORNULL" types are used so slots can be nullable.

Class definition order matters: `DESCRIPTION`'s `Collate:` field controls source-file load order so that referenced classes exist before they're used (e.g., `ClassQualityReport.R` before `ClassSangerRead.R`). When adding a new class file, update `Collate:`.

### Construction flow

User-facing constructors `SangerRead()`, `SangerContig()`, `SangerAlignment()` live in `R/Constructors.R` and are thin wrappers around `new(...)`. The real work happens in each class's `setMethod("initialize", ...)`:

1. Validate every argument by calling `check*` helpers from `R/UtilitiesFuncInputChecker.R` — these append to `errors` / `warnings` lists rather than `stop()`-ing. Construction always returns an object.
2. Branch on **two orthogonal axes**:
   - `inputSource`: `"ABIF"` (read raw `.ab1` files via `sangerseqR`) vs `"FASTA"` (read pre-called sequences via `seqinr::read.fasta`).
   - `processMethod` (Contig/Alignment only): `"REGEX"` (group reads by filename suffix patterns like `_F.ab1$` / `_R.ab1$`) vs `"CSV"` (a names-conversion CSV maps reads → contigs).
3. For ABIF input, two trimming algorithms are selectable via `TrimmingMethod`: `"M1"` (modified Mott's, uses `M1TrimmingCutoff`) or `"M2"` (sliding-window like Trimmomatic, uses `M2CutoffQualityScore` + `M2SlidingWindowSize`). Setting one nulls the other's parameters.
4. `SangerAlignment.initialize` recursively builds child `SangerContig`s, which recursively build child `SangerRead`s — pre-checks at the parent level are skipped on nested calls via the `printLevel` argument (so each parameter is validated exactly once).

### Generic dispatch — three exported "façade" functions

`R/MethodShared.R` defines `launchApp()`, `writeFasta()`, and `generateReport()`. Each inspects `class(object)[1]` and dispatches to the appropriate class-specific exported method:

- `launchApp` → `launchAppSC` / `launchAppSA` (no `*SR` — there is no Shiny app for a single read).
- `writeFasta` → `writeFastaSR` / `writeFastaSC` / `writeFastaSA`.
- `generateReport` → `generateReportSR` / `generateReportSC` / `generateReportSA`.

Generics are declared in `R/AllGenerics.R`; per-class implementations live in `R/MethodSangerRead.R`, `R/MethodSangerContig.R`, `R/MethodSangerAlignment.R`. When adding a new operation that varies by class, follow the same pattern: generic in `AllGenerics.R`, three `setMethod` implementations, optional façade in `MethodShared.R`.

### Shiny apps and reports

- Shiny: each Shiny-capable class has a paired UI/Server file — `R/ShinySangerContigUI.R` + `R/ShinySangerContigServer.R`, and the same for `SangerAlignment`. Shared modules live in `R/ShinyServerModule.R`. The apps are quite large (~1.8k–2.3k lines each).
- Reports: `generateReport*` renders RMarkdown templates from `inst/rmd/` (`SangerRead_Report_ab1.Rmd`, `SangerRead_Report_fasta.Rmd`, `SangerContig_Report.Rmd`, `SangerAlignment_Report.Rmd`) via `rmarkdown::render`.

### External dependencies that drive structure

- `sangerseqR` — base ABIF parsing; `SangerRead` extends `sangerseq`.
- `Biostrings` / `pwalign` — DNA/AA strings, alignment primitives.
- `DECIPHER` — multiple alignment, consensus, frameshift correction, distance matrices, `Treeline`. The minimum DECIPHER version was recently bumped because `Treeline` replaced an older API (see most recent commit `Fix DECIPHER Treeline import; update DECIPHER minimum version`).
- `ape` — `phylo` tree objects (`bionjs`, `dist.dna`).
- `logger` — used throughout for `log_info` / `log_error` / `log_debug`. Many code paths log an error and return rather than `stop()`-ing, consistent with the never-throw construction model.

### Test data

Reusable fixtures live in `data/` (`*.RData`, loadable with `data(sangerReadFData)` etc.), and raw input files live under `inst/extdata/Allolobophora_chlorotica/`, `inst/extdata/Drosophila_melanogaster/`, `inst/extdata/ab1/`, `inst/extdata/fasta/`. Tests reach them via `system.file("extdata/", package = "sangeranalyseR")`.

### Documentation site

`docs/` is a built Sphinx/ReadTheDocs site (sources in `docs/source/`, built HTML in `docs/build/`). The committed `docs/build/` artifacts change frequently — only edit `docs/source/` and rebuild via `docs/Makefile`. For R-level help, edit roxygen comments in `R/` and run `devtools::document()` to regenerate `man/*.Rd`.
