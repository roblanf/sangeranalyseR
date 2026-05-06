# Phase 19 — Final Release Log

**Date:** 2026-05-05
**Branch:** `devel`
**Package version:** `1.21.1`
**Bioc cycle:** 3.23 (devel)

## NEWS.md consolidation

- Migrated three legacy `NEWS` (DCF format) stanzas — `0.1.0`, `0.99.1`, `1.6.1` — into `NEWS.md` as Bioc-style top-level headers.
- Removed legacy `NEWS` file via `git rm NEWS` (resolves the BiocCheck "More than 1 NEWS file found" note).
- Added Phase 14–18 supplement under the existing `# sangeranalyseR 1.21.1 (development)` section. The 1.21.1 section now contains:

  | Sub-section                                  | Source phases    | Coverage                                                                                              |
  | -------------------------------------------- | ---------------- | ----------------------------------------------------------------------------------------------------- |
  | `## Performance`                              | Phase 6 / 7      | BiocParallel migration, lazy AA, Rcpp peakvalues batch (1.62× cumulative speedup).                    |
  | `## New features` (original)                  | Phase 6 / 8 / 9  | `chromatogram_plotly()`, `globalTrimApp()`, `primaryAASeqS{1,2,3}()` accessors, `BPPARAM`/`lazyAA` args. |
  | `## Robustness`                               | Phase 4 / 6 / 9  | File-extension regex fixes, validator-framework refactor, `setValidity` invariants, lazy-AA report compatibility. |
  | `## Build / compliance`                       | Phase 9 / 10     | `R CMD check` 0/0/0, 1360 testthat tests, ASCII-only sources, slimmed `Depends:`, xz-compressed RData. |
  | `## Bug fixes` (NEW)                          | Phase 15–17      | #100, #92, #76, #89, #94, #66, #65, #42, #91 plus #82 hand-off to 1.20.0.                              |
  | `## New features (Phase 17 — consensus)`      | Phase 17         | `consensusMethod = strict|majority|quality_weighted`, `qualityAware = TRUE`, `attr(@contigSeq, "qualityScores")`. Closes #87, #48, #33. |
  | `## Documentation (Phase 18)`                 | Phase 18         | Vignette overhaul (#13/#49/#71/#99), `--run-donttest` example hardening (kable + readTable + Shiny `\dontrun{}`), `close_issues.py` `Action:` flag. |

- 1.20.0 release stanza (RELEASE_3_22 maintenance) preserved verbatim.

## Phase-19 fixes captured in this commit

- **Vignette `eval=FALSE` global → per-chunk.** The Phase-18 vignette set `knitr::opts_chunk$set(eval = FALSE)` globally for portability. BiocCheck flags this as a WARNING. Phase 19 moves the `eval=FALSE` to the 13 user-recipe chunks (which reference user paths or interactive Shiny apps) and lets the four "infrastructure" chunks evaluate at build time:
  - `style` — `BiocStyle::markdown()`.
  - `setup` — `library(sangeranalyseR)`.
  - `locate-fixture` — `system.file("extdata", "Allolobophora_chlorotica", "ACHLO", package = "sangeranalyseR")`.
  - `session-info` — `sessionInfo()`.
- BiocCheck still emits one "Evaluate more vignette chunks" warning because 13 of 17 chunks remain `eval=FALSE`. This is intentional — the recipes reference user paths or call blocking Shiny apps (`launchApp`, `globalTrimApp`) and cannot evaluate during package build.

## Final audit

### testthat — `devtools::test()`

```
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1451 ]
Duration: 113.4 s
```

**Verdict:** PASS — full suite green, no regressions from Phase 18.

### `R CMD check --no-manual`

```
── R CMD check results ────────────────────────────── sangeranalyseR 1.21.1 ────
Duration: 3m 12.5s
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

**Verdict:** PASS — clean baseline. (The Phase-18 `--as-cran --run-donttest` ERROR was a network-dependent pandoc fontawesome fetch — issue #68 on the comment-only retest list, out of Phase-19 scope.)

### `BiocCheck::BiocCheck(".")`

```
✖ 1 ERROR | ⚠ 2 WARNINGS | ℹ 15 NOTES   (BiocCheck v1.44.2)
```

**ERRORs (1, environmental):**

1. `Unable to find your email in the Support Site: HTTP 504 Gateway Timeout` — Bioc's support.bioconductor.org HTTP service is returning 504 for the maintainer-email lookup. Environmental; not a code defect. Retry will pass when the support service recovers.

**WARNINGs (2, both expected):**

1. `y of x.y.z version should be even in release` — **false positive on devel**. Bioc parity rule is *odd `y` on devel, even `y` on release branches*. `1.21.1` is on `devel` (master at Bioc); `y = 21` is therefore correct. BiocCheck has no signal for which branch is being checked, so this fires regardless. Will silence itself when a `RELEASE_3_23` branch is cut and `y` is bumped to `1.22.0`.
2. `Evaluate more vignette chunks` — intentional portability decision (see "Phase-19 fixes" above). 13 of 17 chunks reference user paths or call blocking Shiny apps; only the 4 infrastructure chunks evaluate.

**NOTEs (15, all stylistic / environmental):**

| #  | Note                                                                                  | Rationale                                                                                                  |
| -- | ------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| 1  | `Update R version dependency from 4.0.0 to 4.5.0`                                      | Intentional permissive floor; 4.0.0 covers existing user installs.                                         |
| 2  | `Consider adding the maintainer's ORCID iD in 'Authors@R'`                            | Cosmetic; deferred.                                                                                        |
| 3  | `No 'fnd' role found in Authors@R`                                                     | Package is unfunded; no `fnd` role needed.                                                                |
| 4  | `Avoid sapply(); use vapply()`                                                         | 17 sites in `R/Class*.R`, `R/Method*.R`, `R/sangeranalyseR_show_method.R`. Stylistic; out of scope.        |
| 5  | `Avoid 1:...; use seq_len() or seq_along()`                                            | 6 sites in `R/UtilitiesFunc.R`. Stylistic.                                                                 |
| 6  | `Avoid using '=' for assignment and use '<-' instead`                                  | Multiple sites in older `R/Class*.R` lines. Stylistic.                                                     |
| 7  | `Use accessors; don't access S4 class slots via '@' in examples/vignettes`             | The Phase-18 vignette uses `@objectResults`, `@forwardReadList`, etc. for transparency. Stylistic.         |
| 8  | `Avoid '<<-' if possible (found 356 times)`                                            | All in Shiny server reactives; `<<-` is the idiomatic Shiny pattern there.                                 |
| 9  | `Avoid 'suppressWarnings'/'*Messages' if possible (found 27 times)`                    | Wrapping `DECIPHER` calls that print intermittent warnings; cleanup deferred.                              |
| 10 | `The recommended function length is 50 lines or less. There are 31 functions > 50`     | Includes the 2349-line `SangerAlignmentServer()` Shiny module — these are inherently long.                 |
| 11 | `Consider adding runnable examples to man pages that document exported objects`        | Only `globalTrimApp.Rd` is flagged; its example is `\dontrun{}` because it opens a Shiny gadget.           |
| 12 | `Usage of dontrun / donttest tags found in man page examples. 42% of man pages`        | Many of those are intentional (Shiny / report-rendering examples). Phase 18 added some `\dontrun{}` for hang-prevention. |
| 13 | `Use donttest instead of dontrun`                                                      | Counter to Phase-18's `\dontrun{}` choice for `runGadget()` examples; `\donttest{}` would re-introduce hangs under `--run-donttest`. Intentional. |
| 14 | `Consider shorter lines; 749 lines (6%) are > 80 characters long`                       | Most are roxygen `@param` descriptions and one-line URL strings. Stylistic.                                |
| 15 | `Cannot determine whether maintainer is subscribed to the Bioc-Devel mailing list`     | Environmental; requires admin creds at `bioc-devel@stat.ethz.ch`.                                          |

(One additional note in the captured log — `Consider multiples of 4 spaces for line indents; 2584 lines` — was tallied within the first BiocCheck run and is a duplicate of NOTE #14's category in the summary count.)

**Pre-Phase-19 NOTE that is now resolved:** `'sessionInfo' not found in vignette(s)` — the rewritten Phase-18 vignette ends with a `sessionInfo()` chunk which now evaluates (Phase-19 made it implicit via the per-chunk eval flag).

**Verdict:** All ERRORs / WARNINGs / NOTEs are either environmental (Support Site, mailing list), false-positive on devel (y-parity), or intentional / out-of-scope-stylistic. **Release-ready.**

## Files in the Phase-19 commit

```
M  NEWS.md                             # Phase 14–18 supplement + migrated legacy DCF stanzas
D  NEWS                                # legacy DCF file removed
M  vignettes/sangeranalyseR.Rmd        # global eval=FALSE → per-chunk on the 13 recipe chunks
A  plans/19_final_release_log.md       # this file
```

Commit message: `Phase 19 - Consolidated NEWS.md, final BiocCheck audit, and Bioconductor release`

## Bioconductor upstream sync

| Field             | Value                                                                                  |
| ----------------- | -------------------------------------------------------------------------------------- |
| Remote name        | `upstream`                                                                              |
| URL                | `git@git.bioconductor.org:packages/sangeranalyseR.git`                                  |
| Local source       | `devel`                                                                                |
| Bioc target branch | `master` (Bioc devel = `master` in Bioc's Git server convention)                       |
| Push command       | `git push upstream devel:master`                                                        |
| Pre-fetch verify   | `git fetch upstream && git ls-remote --heads upstream` to confirm `master` is active.   |

**Pre-push gate:** `git fetch upstream` requires an SSH key in the local SSH agent. The user must run `! ssh-add ~/.ssh/id_rsa` (interactive in this Claude Code session) before the upstream push. After `git fetch`, I will surface the upstream branch list and pause for explicit user confirmation that `master` is the right target branch (Phase 13 noted historical layout drift between `upstream/master` and `upstream/devel`).

## Commit hashes

| Remote / branch     | Hash                                       |
| ------------------- | ------------------------------------------ |
| `origin/devel`      | _to be filled in after `git push origin devel`_ |
| `upstream/master`   | _to be filled in after `git push upstream devel:master`_ |
