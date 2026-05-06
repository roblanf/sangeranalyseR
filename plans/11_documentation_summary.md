# sangeranalyseR — Phase 11 Documentation Overhaul

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

---

## 1. Headline result

| Item                               | Before Phase 11                                                                                        | After Phase 11                                                                                                                                                          |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Project README**                 | `readme.md` (lowercase), 113 lines, 2021-era Travis CI badge, MIT license badge (which contradicted the `License: GPL-2` field), no mention of any post-2021 work. | **`README.md`** (canonical case), full rewrite covering Phase 6 → 10 work: BiocParallel, lazy AA, Rcpp port, Plotly chromatograms, globalTrimApp, build compliance, coverage. |
| **Badges**                         | License (wrong), Travis (defunct), generic ReadTheDocs, OS                                              | R-CMD-check, BioC release, BioC devel, Codecov, License (correct GPL-2), platform, R version                                                                             |
| **`devtools::document()`**         | Already up to date as of Phase 10                                                                      | Re-run; NAMESPACE / `man/` confirmed in sync; 0 added/removed.                                                                                                            |
| **`rcmdcheck`**                    | 0 / 0 / 0 (Phase 9)                                                                                    | **0 / 0 / 0** (no regression).                                                                                                                                            |
| **Tests**                          | 1360 PASS                                                                                              | **1360 PASS** (Phase 11 is documentation-only; no test changes).                                                                                                          |

---

## 2. New `README.md` structure

The file has 8 top-level sections, sized for a 2-minute read while still surfacing every refactor milestone:

```
README.md
├── Badges
│     ├── R-CMD-check (GitHub Actions)
│     ├── BioC release / BioC devel build status
│     ├── Codecov coverage
│     ├── License: GPL-2  (corrected from the old "MIT" badge)
│     ├── platform (macOS / Linux / Windows)
│     └── R-version (>= 4.0.0)
├── Tagline + ReadTheDocs link
├── What's new (devel)
│     ├── Feature table — 7 rows covering BiocParallel, lazy AA, Rcpp, Plotly,
│     │   globalTrimApp, setValidity, strict build compliance.
│     │   Each row: feature name | what it does | impact.
│     └── Cumulative wall-time progression table (1.85 s → 1.07 s, 1.62× speedup)
├── Installation
│     ├── From Bioconductor (release + devel paths)
│     ├── From GitHub (devtools::install_github)
│     └── System requirements (R >= 4.0.0; C++17 toolchain; pandoc optional)
├── Quick start
│     ├── 1. Load and assemble — full SangerAlignment(...) call with BPPARAM
│     ├── 2. Tweak trimming interactively — globalTrimApp(SA)
│     ├── 3. WebGL chromatogram — chromatogram_plotly(sr, max_points = 8000)
│     ├── 4. Export and report — writeFasta + generateReport
│     ├── Lazy AA accessors — primaryAASeqS1/S2/S3 + lazyAA = FALSE escape hatch
│     └── Cross-platform parallel — Multicore / Snow / Serial examples
├── Citation (GBE 2021 paper, BibTeX-friendly format)
└── Maintainers + License (GPL-2) + Issues link
```

### Key emphasised numbers (sourced from prior phase artifacts)

| Number                          | Source                                                                  |
| ------------------------------- | ------------------------------------------------------------------------ |
| **1.7× end-to-end speedup**     | Cumulative milestone: 1.85 s → 1.07 s (best). `plans/05`, `06`, `07`.  |
| **~35% wall-time savings (lazy AA)** | Phase-6 profile diff: `Biostrings::translate` was 35.4% self-time.   |
| **~1.28× from Rcpp port**        | Phase-7 microbench `bench_e2e.csv` (1.461 s → 1.142 s mean over 5 reps). |
| **~320 ms saved per build**      | Same source.                                                            |
| **1360 testthat tests pass**     | Phase-10 final run.                                                     |
| **>87% coverage on every non-Shiny file** | Phase-10 covr report.                                            |

### Errors corrected from the old `readme.md`

- **License badge fixed**: was `License: MIT`, but `DESCRIPTION` has `GPL-2 | file LICENSE`. Now correctly shows `GPL-2`.
- **Travis CI badge removed**: Travis-CI for OSS was discontinued in 2021. Replaced with the GitHub Actions `R-CMD-check` workflow badge.
- **Generic "documentation build" badge removed**: pointed at a generic pip badge; replaced with codecov + Bioc build badges.

---

## 3. File rename: `readme.md` → `README.md`

The repo previously tracked the file as lowercase `readme.md`. Phase 11 promotes it to the canonical `README.md` which:

1. Matches GitHub's automatic README detection.
2. Matches `R CMD build` / Bioconductor expectations.
3. Matches the `tools::Rd2HTML` cross-reference convention.

Performed via two-step `git mv` (`readme.md → tmp__readme → README.md`) because the underlying APFS filesystem is case-insensitive — a direct same-name rename would be a no-op for git. The git index now shows `README.md`; the lowercase form is gone.

---

## 4. `devtools::document()` re-sync

Re-ran `roxygen2::roxygenise()` after the README rewrite. Roxygen had no work to do (all `@importFrom`, `@export`, `@rdname`, `@aliases` from Phase 9 are still consistent with the source) — confirming Phase 9's documentation hygiene is preserved.

`NAMESPACE` and the 31 generated `man/*.Rd` files were already in lock-step with the R source from Phase 10. Phase 11 changes nothing here.

`DESCRIPTION` confirmed to carry:
- `Authors@R` (modern format, replacing the deprecated `Author:` + `Maintainer:` pair) — `Rob Lanfear (aut)`, `Kuan-Hao Chao (aut, cre)`.
- `URL: https://github.com/roblanf/sangeranalyseR`
- `BugReports: https://github.com/roblanf/sangeranalyseR/issues`
- `License: GPL-2 | file LICENSE`
- `Depends: R (>= 4.0.0)` (rolled back in Phase 10 for user compatibility).

---

## 5. Files touched in Phase 11

```
A  README.md                                       (rewrite, replacing readme.md)
D  readme.md                                       (renamed via git mv)
A  plans/11_documentation_summary.md
```

No code changes; no tests added; no DESCRIPTION / NAMESPACE / man/ changes. Phase 11 is intentionally documentation-only.

---

## 6. Reproducing

```r
# Verify documentation is in sync (no diff expected on man/ / NAMESPACE):
roxygen2::roxygenise()
# -> nothing rewritten

# Confirm clean build:
rcmdcheck::rcmdcheck(".", args = "--no-manual")
# -> STATUS 0 / 0 errors / 0 warnings / 0 notes

# Confirm full suite still passes:
devtools::test()
# -> 1360 / 1360 PASS
```

Render the README locally (e.g. via VSCode markdown preview, `pandoc README.md -o README.html`, or just view on github.com after push) to confirm badges resolve.

---

## 7. Non-goals (deferred)

- **Vignette refresh**. The vignette in `vignettes/sangeranalyseR.Rmd` still points at the old (Bioc 3.13-era) workflow. Bringing it into line with the new BPPARAM / lazyAA / Rcpp story is a separate, much larger documentation task.
- **`pkgdown` site**. Not currently part of the build; would complement the existing ReadTheDocs site.
- **NEWS.md**. The existing `NEWS` file uses an old DCF format and stops at version 1.6.1. A modern `NEWS.md` covering Phases 6–10 would help downstream packagers follow upgrades.
