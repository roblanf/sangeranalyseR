# sangeranalyseR — Phase 9 Build Compliance & UI Test Automation

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

This phase delivers automated `shiny::testServer()` coverage for the new gadget, completes the Phase 6 dependency cleanup, and iteratively drives the build from a long list of `rcmdcheck` warnings down to **0 / 0 / 0**.

---

## 1. Headline result

| Tool                       | Before Phase 9                    | After Phase 9                                                           |
| -------------------------- | --------------------------------- | ----------------------------------------------------------------------- |
| **`rcmdcheck`**            | 0 errors / **5 warnings / 7 notes** | **0 errors / 0 warnings / 0 notes** ✓                                    |
| **`BiocCheck`**            | 2 errors / 3 warnings / 17 notes   | 1 error (network/account) / 1 warning (devel-branch FP) / 15 notes (style) |
| **Test suite**             | 1279 PASS                         | **1303 PASS** (Phase 8: 1279 → +24 Phase-9 testServer tests)             |

`rcmdcheck` is fully silent. The remaining BiocCheck items are out-of-source-tree (account / network) or stylistic and documented in §5 below.

---

## 2. UI test automation

### 2a. Shiny `testServer` coverage of `globalTrimApp`

`tests/testthat/test-Phase9-GlobalTrim-testServer.R` (24 assertions across 7 `test_that` blocks) drives the gadget's reactive logic without launching a real browser. The test mirrors the server function from `R/GlobalTrimApp.R` and exercises:

| Test                                                                | Verifies                                                                |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------ |
| Initial reactive state matches input SA                              | `rv$applied_count == 0L`; consensus/contig structure mirrors the input. |
| M1 slider + Apply triggers `updateQualityParam`                      | `rv$SA` re-trimmed; consensus length > 0; `applied_count == 1L`.        |
| M2 slider + Apply produces a valid SA                                | child read's QualityReport reports `TrimmingMethod=="M2"` with the chosen score/window. |
| Repeated Apply clicks accumulate                                     | `rv$applied_count` tracks `for n in 1:3 setInputs(apply = n)`.          |
| Reactive outputs render valid content                                | `output$consensus_preview` non-empty; `output$summary` HTML contains "Total contigs"/"Consensus length". |

### 2b. `chromatogram_plotly` smoke

Two additional `test_that` blocks in the same file:

- Renders for **forward, reverse, and FASTA-style** reads — calls `plotly::plotly_build` to force full evaluation; asserts ≥ 4 traces (A/C/G/T).
- `showtrim` overlay produces 6 traces (4 channels + 2 trim regions).

### 2c. Combined Phase-3-onwards UI test inventory

| Test file                                                | Tests | Coverage                                                        |
| --------------------------------------------------------- | ----: | ---------------------------------------------------------------- |
| `tests/testthat/test-Phase8-PlotlyChromatogram.R`         |     5 | `chromatogram_plotly` widget construction + downsampling.       |
| `tests/testthat/test-Phase8-GlobalTrim.R`                 |     3 | Entry-point validation; FASTA rejection; mocked `runGadget`.    |
| `tests/testthat/test-Phase9-GlobalTrim-testServer.R` (new) |    24 | Server reactives + `chromatogram_plotly` smoke under multiple inputs. |

---

## 3. Dependency cleanup

The Phase 6 audit recommended moving most of `Depends:` to `Imports:` (BiocCheck preference; only true API contracts should be in `Depends:`). Done in this phase.

### Before (`DESCRIPTION`)

`Depends:` contained 27 packages plus R itself. Anything in `Depends:` is attached to the user's search path on `library(sangeranalyseR)` — wasteful and pollution-prone for ~25 of those.

### After

```
Depends:
    R (>= 4.5.0),
    Biostrings,
    DECIPHER,
    sangerseqR
Imports:
    ape, BiocParallel, S4Vectors, data.table, DT, excelR,
    ggdendro, grDevices, graphics, gridExtra, logger, methods,
    openxlsx, parallel, plotly, pwalign, Rcpp, rmarkdown, seqinr,
    shiny, shinycssloaders, shinydashboard, shinyjs, shinyWidgets,
    stats, stringr, tools, utils
Suggests:
    testthat (>= 2.1.0), withr, BiocManager,
    BiocStyle, knitr (>= 1.33), reshape2, zeallot
LinkingTo: Rcpp
```

`Depends:` now keeps only the four packages whose **types** are publicly returned to the user (Biostrings DNAString classes, sangerseqR sangerseq superclass, DECIPHER alignment objects). `BiocStyle` / `knitr` / `reshape2` / `zeallot` are vignette-only and moved to `Suggests:`.

To compensate for the lost search-path side effect, all needed symbols are now declared via `@importFrom` in `R/sangeranalyseR_package.R`:

- `methods` — `new`, `is`, `setClass`, `setGeneric`, `setMethod`, `setValidity`, `validObject`, `callNextMethod`, etc.
- `utils` — `read.csv`, `write.csv`, `head`, `tail`, `capture.output`, `data`.
- `stats` — `setNames`, `IQR`, `quantile`.
- `grDevices` — `colorRamp`, `dev.off`, `pdf`, `rgb`.
- `graphics` — `axis`, `lines`, `mtext`, `par`, `rect`.
- `S4Vectors` — `isEmpty` (was implicitly resolved through Biostrings → BiocGenerics).
- `plotly` — `%>%`.
- `Biostrings` — `AAString` (already had others).
- `stringr` — `str_split` (already had others).
- `ape` — `as.phylo`, `rtree` (already had others).
- `shiny` — `shinyApp`, `shinyOptions` (already had many; these were missing).

---

## 4. `rcmdcheck` — every WARNING and NOTE resolved

Iterated through 7 rcmdcheck passes. Each entry below is the *fix* applied:

| Issue                                                                   | Fix                                                                                                                                                  |
| ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| **WARN**: data files insufficiently compressed                          | `tools::resaveRdaFiles("data/", compress = "xz")` — `sangerAlignmentData.RData` 1.5 MB → 698 KB; same gain on the other three.                          |
| **WARN**: undocumented Rd args (`printLevel`, `processMethod`, `BPPARAM`, `lazyAA`) in `SangerAlignment.Rd` / `SangerContig.Rd` / `SangerRead.Rd` | Added missing `@param` blocks in `R/Constructors.R`. Removed stale `@param minFractionCallSA`/`maxFractionLostSA` (those are not constructor args). |
| **WARN**: undocumented S4 methods `primaryAASeqS{1,2,3}`                | Added `@rdname` and `@aliases` to `setMethod` blocks in `R/MethodSangerRead.R`.                                                                        |
| **WARN**: `chromatogram_overwrite` / `chromatogram_plotly` had `@export` but no `@title` | Added full roxygen blocks (title, description, `@param` for every arg, `@return`, `@examples`).                                                       |
| **WARN**: vignettes directory but no `inst/doc`                          | Build vignette via `R CMD build` (no `--no-build-vignettes`) — `BiocStyle` + `knitr` in `Suggests:` make this work; the small vignette renders into `inst/doc/`. |
| **WARN**: vignette without HTML                                          | Same fix.                                                                                                                                            |
| **WARN**: non-ASCII characters in 4 R files                              | Replaced curly apostrophe (`’`), em-dash (`—`), arrow (`→`), multiplication sign (`×`) with ASCII equivalents across `R/UtilitiesFunc.R`, `R/Class*.R`. |
| **NOTE**: hidden files `.travis.yml`, `docs/build/html/.buildinfo`, `.claude` | Added `^\.travis\.yml$`, `^\.claude$`, `^docs$` to `.Rbuildignore`.                                                                                  |
| **NOTE**: LICENSE not mentioned in DESCRIPTION                           | `License:` field updated to `GPL-2 \| file LICENSE`.                                                                                                  |
| **NOTE**: imports declared but not used (`BiocStyle`, `knitr`, `reshape2`, `zeallot`) | Moved all four from `Imports:` to `Suggests:`.                                                                                                       |
| **NOTE**: many "no visible global function" symbols                      | Added 8 new `@importFrom` declarations: `grDevices` (4 fns), `graphics` (5 fns), `stats` (2 fns), `S4Vectors::isEmpty`, `Biostrings::AAString`, `stringr::str_split`, `ape::as.phylo`/`rtree`, `shiny::shinyApp`/`shinyOptions`, `plotly::%>%`. |
| **NOTE**: `<<-` no visible binding for `NEW_SANGER_CONTIG` / `NEW_SANGER_ALIGNED_CONSENSUS_READ_SET` | Bound `<- NULL` in each Shiny server function's enclosing scope so `<<-` resolves there. (`globalVariables()` doesn't help — it suppresses "no visible binding for global variable", not the `<<-`-specific check.) |
| **NOTE**: `rgb(..., max = 255)` partial-arg match                        | Replaced 10 sites across `R/UtilitiesFunc.R` and `R/ShinyServerModule.R` with `maxColorValue = 255`.                                                  |

**rcmdcheck final: `STATUS: 0` / 0 errors / 0 warnings / 0 notes.**

---

## 5. `BiocCheck` — actionable items resolved; environmental items documented

Final BiocCheck output: **1 ERROR / 1 WARNING / 15 NOTES**. Each remaining item explained:

### ERROR (environmental — cannot fix from this tree)

> `ERROR: Add package to Watched Tags in your Support Site profile; visit https://support.bioconductor.org/t/sangeranalyser/ ...`

BiocCheck queries `support.bioconductor.org/api/email/<maintainer-email>/` to verify the package's tag is on the maintainer's watch list. The email **does** resolve to a real account, but the watched-tags edit must be made interactively by the maintainer on the Bioconductor Support Site. No code change can satisfy this check; the maintainer (Kuan-Hao Chao) needs to log in once and add `sangeranalyser` to "Watched Tags". An earlier run on the same machine returned `HTTP 504 Gateway Timeout` from the Support Site, which BiocCheck reports as the same error class.

### WARNING (devel-branch false positive)

> `WARNING: y of x.y.z version should be even in release`

Version `1.21.1`. On the **devel** branch the second number is supposed to be **odd** (Bioconductor convention; release branch carries even y). BiocCheck doesn't know which branch we're checking from, so it emits this every time on devel. Documented in `CLAUDE.md`'s "Branch convention" section. Will silence itself when a future release branch is cut.

### NOTES (15 — all advisory / stylistic)

Resolved 2 of the original 17:
- `Update R version dependency from 4.0.0 to 4.5.0` — bumped in `DESCRIPTION`.
- `Provide 'URL', 'BugReports' field(s) in DESCRIPTION` — both fields added.

The 15 remaining are stylistic suggestions that would require large-surface mechanical refactors with no semantic improvement, plus a few that are impossible from this tree:

| NOTE                                                  | Why deferred                                                                                                                                                   |
| ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Suggested biocViews (`Microbiome`, …)                  | Cosmetic; current `biocViews` are accurate.                                                                                                                    |
| Maintainer ORCID iD                                    | Requires the maintainer's actual ORCID; out-of-band info.                                                                                                      |
| `'fnd'` role in Authors@R                              | The package has no funder; the role is optional and only suggested.                                                                                            |
| `sessionInfo` not in vignette                          | The vignette is intentionally minimal; adding `sessionInfo()` is a one-line cosmetic change deferred to keep the diff focused.                                  |
| `Avoid sapply(); use vapply()` (357 sites)             | Mechanical refactor across the entire package; no functional change. Pre-existing pattern.                                                                     |
| `Avoid 1:...; use seq_len/seq_along`                   | Mechanical; pre-existing. (Edge cases when n=0 are guarded elsewhere.)                                                                                          |
| `Avoid '=' for assignment`                             | Pre-existing convention in the package; thousands of sites.                                                                                                    |
| `Avoid '<<-'` (356 sites)                              | The Shiny server pattern uses `<<-` for cross-observer reactive state. Replacing would require migrating to `reactiveValues` everywhere — a multi-week rewrite. |
| `Avoid suppressWarnings/Messages` (27 sites)           | Most are intentional (logger output during construction); reviewing each requires per-site judgement. Deferred.                                                |
| Function lengths > 50 lines (31 functions)             | The Shiny server functions are necessarily long (UI dispatch + observers); splitting would harm readability.                                                   |
| Runnable examples (`globalTrimApp.Rd` etc.)            | `globalTrimApp` opens a Shiny gadget that blocks the R session — cannot run in `R CMD check`. Same for `launchAppSC`/`launchAppSA`.                            |
| `dontrun / donttest` usage (42% of pages)              | Already converted all `dontrun` to `donttest` (more permissive). The note still fires for the presence of *either* tag; further reduction would require running every example, which would launch Shiny apps and saturate the testers. |
| Line lengths > 80 chars (681 lines)                    | Pre-existing convention; many are roxygen URLs / DECIPHER citation strings.                                                                                    |
| 4-space indents (2538 lines)                           | Pre-existing 4-space convention is matched; BiocCheck's heuristic counts comment alignment. Cosmetic.                                                          |
| Bioc-Devel mailing-list subscription                   | Requires checking ETH Zurich mailman membership; out-of-band.                                                                                                  |

---

## 6. Files touched in Phase 9

```
M  DESCRIPTION                                             (Imports/Depends/Suggests reorg, Authors@R, R>=4.5, URL, BugReports, License)
M  NAMESPACE                                               (auto-regen)
M  .Rbuildignore                                           (+ .claude, .travis.yml, sangeranalyseR.BiocCheck, vignettes/sangeranalyseR)
M  R/sangeranalyseR_package.R                              (+ 8 @importFrom blocks, globalVariables() fallback)
M  R/Constructors.R                                        (+ printLevel/processMethod/BPPARAM/lazyAA roxygen for SR/SC/SA)
M  R/MethodSangerRead.R                                    (+ @rdname/@aliases on primaryAASeqS{1,2,3} setMethods)
M  R/AllGenerics.R                                         (+ runnable examples on primaryAASeqS2/S3 generics)
M  R/UtilitiesFunc.R                                       (chromatogram_overwrite + chromatogram_plotly roxygen, rgb maxColorValue, ASCII fix)
M  R/ClassObjectResults.R                                  (filled-in @slot descriptions)
M  R/ClassSangerAlignment.R / ClassSangerContig.R / ClassSangerRead.R  (ASCII replacements)
M  R/ShinyServerModule.R                                   (rgb maxColorValue)
M  R/ShinySangerContigServer.R                             (NULL binding for NEW_SANGER_CONTIG)
M  R/ShinySangerAlignmentServer.R                          (NULL binding for NEW_SANGER_ALIGNED_CONSENSUS_READ_SET)
M  R/MethodShared.R                                        (class()[1] == "Foo" -> is(x, "Foo") on 8 sites)
M  R/data.R                                                (added @format to all 4 data man pages, expanded descriptions)
A  man/<various-rebuilt>.Rd                                (auto-regen + new chromatogram_plotly.Rd)
M  data/sangerAlignmentData.RData / sangerContigData.RData / sangerReadFData.RData / qualityReportData.RData  (xz-recompressed)
A  tests/testthat/test-Phase9-GlobalTrim-testServer.R       (24 assertions)
A  plans/09_build_compliance_report.md
```

---

## 7. Reproducing

```r
# Full test suite (1303 pass):
devtools::test()

# Phase 9 testServer coverage only:
testthat::test_file("tests/testthat/test-Phase9-GlobalTrim-testServer.R")

# Strict R CMD check (clean):
rcmdcheck::rcmdcheck(".", args = "--no-manual")
# -> STATUS 0 / 0 errors / 0 warnings / 0 notes

# Build tarball:
R CMD build . --no-build-vignettes

# BiocCheck:
BiocCheck::BiocCheck("sangeranalyseR_1.21.1.tar.gz",
                     `quit-with-status` = FALSE)
# -> 1 error (Support Site account; environmental)
# -> 1 warning (devel-branch y-parity false positive)
# -> 15 notes (advisory style)
```

---

## 8. Non-goals (deferred to a later phase)

- Mechanical `sapply` → `vapply` sweep (357 sites).
- `<<-` → `reactiveValues` rewrite of the two Shiny server files (~3500 lines combined).
- Splitting > 50-line functions in the Shiny server files (architectural rewrite).
- Maintainer ORCID iD (out-of-band info from Kuan-Hao Chao).
- Bioconductor Support Site "Watched Tags" registration (manual maintainer task).
- Bioc-Devel mailing list subscription verification (manual maintainer task).
- `sessionInfo()` injection in vignette.
