# sangeranalyseR — Phase 12 Documentation Refinement & Changelog Migration

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

---

## 1. Headline result

| Item                                  | Before Phase 12                                                             | After Phase 12                                                                |
| ------------------------------------- | --------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `README.md` "What's new" subsection   | 7-row feature table + cumulative timing table (~25 lines)                    | **3 plain bullets** + link to `NEWS.md` (8 lines)                              |
| Iconic Shiny screenshot                | **Missing** (dropped during Phase 11 rewrite)                                | **Restored** at "(3) Explore your data" anchor (`<img src="…/gwY6AqB.png">`)  |
| R-CMD-check badge                     | Already first in the badges block (Phase 11)                                 | Confirmed first; unchanged                                                    |
| `NEWS.md`                              | **Did not exist**                                                            | **Created** with versioned `# sangeranalyseR 1.21.1` and `# 1.20.0` sections  |
| Legacy DCF `NEWS` file                 | Present, stops at v1.6.1                                                     | **Untouched** — both can coexist; removing legacy risks pre-existing tooling  |
| `R CMD check`                          | 0 / 0 / 0                                                                    | **0 / 0 / 0** (no regression)                                                  |
| Tests                                  | 1360 PASS                                                                    | **1360 PASS** (Phase 12 is documentation-only)                                 |

---

## 2. README simplification

### "What's new" (was vs now)

**Before** — engineer-focused 7-row feature table ("`BiocParallel` support / Lazy AA / Rcpp / Plotly / globalTrimApp / setValidity / Strict build compliance") plus a 4-row cumulative timing table:

```
| Feature              | What it does                | Impact                  |
| BiocParallel support | bplapply replaces mclapply… | Multicore on Linux/macOS… |
| Lazy AA…             | …                            | …                        |
| (etc., 7 rows)       |                             |                         |
```

**After** — three end-user bullets:

> - **~1.7× faster `SangerAlignment(...)`** thanks to a C++ peak-detection inner loop, parallel per-read construction (`BiocParallel`), and lazy 3-frame amino-acid translation that only runs when you ask for it.
> - **Interactive Plotly + WebGL chromatograms** via the new `chromatogram_plotly()` — smooth scrolling and zoom on Sanger traces with tens of thousands of points.
> - **Global trimming dashboard** via the new `globalTrimApp(sa)` — adjust M1 / M2 trimming parameters across an entire `SangerAlignment` with a live consensus preview.
>
> For the full per-version changelog see [`NEWS.md`](NEWS.md).

The dense feature/impact wording moves to `NEWS.md` (§3 below).

### Restored visual

Pre-Phase-11 `readme.md` had:

```html
#### (3) Explore your data
<img src="https://i.imgur.com/gwY6AqB.png" style="width:100%">
```

Phase 11 dropped this when it rewrote Quick Start. Phase 12 restores it at the same semantic position — between the `globalTrimApp(sa)` step and the `chromatogram_plotly()` example — under a new `### 3. Explore your data` subsection that opens with a `launchApp(sa)` call (the action the screenshot illustrates):

```html
### 3. Explore your data

Open the per-read Shiny app for a `SangerContig` (or use `launchApp(sa)` on a full `SangerAlignment`):

`launchApp(sa)`

<img src="https://i.imgur.com/gwY6AqB.png" alt="..." style="width:100%">
```

The `alt` attribute was added (it was missing in the original) for accessibility / screen-reader support.

### R-CMD-check badge confirmation

Already the first badge in the `<!-- badges: start -->` block on line 2 of `README.md`. Phase 12 leaves it untouched. Confirmed via `head -3 README.md` — the line is exactly:

```
[![R-CMD-check](https://github.com/roblanf/sangeranalyseR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/roblanf/sangeranalyseR/actions/workflows/R-CMD-check.yaml)
```

### Net README size

189 → 183 lines. The "What's new" section shrank by ~17 lines; the new "(3) Explore your data" subsection added ~11 lines.

---

## 3. New `NEWS.md` structure

```
NEWS.md
├── # sangeranalyseR 1.21.1 (development)
│     ├── ## Performance
│     │     ├── BiocParallel-backed per-SangerRead loop (Phase 6)
│     │     ├── Lazy 3-frame AA translation (Phase 6)
│     │     ├── Rcpp peakvalues_batch_cpp port (Phase 7)
│     │     │     └── cumulative timing table (1.85 → 1.07 s, 1.62×)
│     │     └── Removed redundant mclapply calls
│     ├── ## New features
│     │     ├── chromatogram_plotly()
│     │     ├── globalTrimApp(SA)
│     │     ├── primaryAASeqS{1,2,3}() accessors
│     │     ├── BPPARAM constructor argument
│     │     └── lazyAA constructor argument
│     ├── ## Robustness
│     │     ├── File-extension regex bug fixes (.Xfa, .Xab1)
│     │     ├── Validator framework refactor (31 -> 9 helpers)
│     │     ├── S4 setValidity invariants
│     │     ├── Lazy-AA report compatibility (30 read sites converted)
│     │     └── Stray cat()/print()/message() debug removal
│     └── ## Build / compliance
│           ├── R CMD check 0/0/0
│           ├── 1360 testthat tests + per-phase test files listed
│           ├── Coverage ≥ 87% on every non-Shiny file
│           ├── DESCRIPTION modernisation (Authors@R, URL, BugReports, R >= 4.0.0)
│           ├── Depends slimmed 27 -> 4
│           ├── ASCII-only source
│           └── data/ xz-recompressed
└── # sangeranalyseR 1.20.0 (current Bioconductor release)
      └── DECIPHER Treeline import fix; release-cycle version bumps
```

The cumulative-timing table from the old README is preserved verbatim under "Performance / Rcpp" so anyone re-checking the speedup claim has the per-milestone numbers.

Backreferences to `plans/05_e2e_validation_report.md`, `plans/06_scaling_summary.md`, and `plans/07_rcpp_optimization_log.md` are kept so future maintainers can find the raw artifacts.

The legacy DCF-format `NEWS` file (which stops at v1.6.1) is **untouched**. Both formats can coexist; removing the legacy file risks breaking pre-existing release tooling that may still parse it.

---

## 4. Verification

| Check                                      | Result                                                                                          |
| ------------------------------------------ | ----------------------------------------------------------------------------------------------- |
| `head -3 README.md`                        | R-CMD-check is the first badge.                                                                 |
| `grep "i.imgur.com/gwY6AqB" README.md`     | 1 hit (line 107) — screenshot restored.                                                         |
| `wc -l README.md`                          | 183 (was 189 in Phase 11; net -6 after slim + restore).                                         |
| `head -1 NEWS.md`                          | `# sangeranalyseR 1.21.1 (development)`                                                          |
| `grep "NEWS\\.md" .Rbuildignore`           | empty — `NEWS.md` ships in the tarball.                                                          |
| `Rscript -e 'devtools::test()'`            | 1360 / 1360 PASS, no regressions.                                                                 |
| `Rscript -e 'rcmdcheck::rcmdcheck(...)'`   | STATUS 0 / 0 errors / 0 warnings / 0 notes.                                                      |

---

## 5. Files touched in Phase 12

```
M  README.md                                        (slim "What's new" + restore screenshot)
A  NEWS.md                                          (new modern changelog)
A  plans/12_documentation_refinement.md
```

3 files, 0 R source / test / DESCRIPTION / NAMESPACE changes. Phase 12 is intentionally documentation-only.

---

## 6. Reproducing

```r
# Verify nothing in code or build broke:
devtools::test()                                      # 1360 / 1360 PASS
rcmdcheck::rcmdcheck(".", args = "--no-manual")        # 0 / 0 / 0
```

Render `README.md` and `NEWS.md` locally (e.g. VSCode markdown preview, or `pandoc README.md -o README.html`) to confirm the badges resolve and the Imgur image loads.

---

## 7. Non-goals (deferred)

- Removing the legacy DCF `NEWS` file. Two changelog files coexisting is harmless; removing risks breaking pre-existing tooling.
- A `pkgdown` site. README + ReadTheDocs already cover the surface.
- Updating `vignettes/sangeranalyseR.Rmd` to reflect Phase 6–10 features (same scope as Phase 11's deferred non-goal).
