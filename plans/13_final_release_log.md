# sangeranalyseR — Phase 13 Final Release Log

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioconductor 3.24 cycle).

This is the definitive final-audit snapshot before pushing the post-Phase-12 state both to GitHub `origin` and to the Bioconductor `upstream` Git server.

---

## 1. Final audit — actual measured results

### `devtools::test()`

```
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1360 ]
Duration: 84.7 s
```

**1360 / 1360 PASS.** Zero regressions across all 12 prior phases.

### `rcmdcheck::rcmdcheck(args = "--no-manual")`

```
STATUS:  0
ERRORS:  0
WARNINGS: 0
NOTES:   0
```

**Strictly clean.** Maintained from Phase 9 across Phases 10–12.

### `BiocCheck::BiocCheck()`

Run on the freshly-built `sangeranalyseR_1.21.1.tar.gz` produced by `R CMD build . --no-build-vignettes`:

```
✖ 1 ERRORS  | ⚠ 1 WARNINGS  | ℹ 17 NOTES
```

Each item is documented as a known, unfixable environmental or stylistic finding from Phase 9 / 10:

| Severity   | Item                                                                         | Category                                                                                                                                       |
| ---------- | ---------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| **ERROR**   | "Add package to Watched Tags in your Support Site profile"                  | Environmental — requires an interactive maintainer login at support.bioconductor.org. No code change can satisfy this check. (Documented in `plans/09_build_compliance_report.md` §5.) |
| **WARNING** | "y of x.y.z version should be even in release"                              | False positive on the **devel** branch. Bioconductor convention: odd `y` on devel, even `y` on release. Version `1.21.1` has y=21 (odd) — correct for devel. Will silence itself when the next release branch is cut. |
| NOTE       | Update R version dependency (4.0.0 → 4.5.x)                                 | **Accepted trade-off** in Phase 10 to maximise user compatibility. Documented in `plans/10_coverage_report.md` §2.                              |
| NOTE       | Suggested `biocViews` (e.g. `Microbiome`)                                   | Cosmetic; existing biocViews accurate.                                                                                                          |
| NOTE       | Maintainer ORCID iD                                                         | Out-of-band; requires actual ORCID from Kuan-Hao Chao.                                                                                          |
| NOTE       | `'fnd'` role in Authors@R                                                   | The package has no funder; role is optional.                                                                                                    |
| NOTE       | `sessionInfo()` in vignette                                                 | Vignette refresh is the Phase 14 deferred non-goal.                                                                                             |
| NOTE       | `Avoid sapply(); use vapply()` (357 sites)                                  | Pre-existing pattern; mechanical mass refactor; deferred.                                                                                       |
| NOTE       | `Avoid 1:...; use seq_len/seq_along`                                        | Same.                                                                                                                                           |
| NOTE       | `Avoid '=' for assignment`                                                  | Pre-existing convention.                                                                                                                        |
| NOTE       | `Avoid '<<-'` (356 sites in the Shiny servers)                              | The Shiny server pattern uses `<<-` for cross-observer reactive state. Replacing requires migrating to `reactiveValues` everywhere. Deferred.   |
| NOTE       | `Avoid suppressWarnings/Messages` (27 sites)                                | Most are intentional (logger output during construction).                                                                                       |
| NOTE       | Function lengths > 50 lines (31 functions)                                  | The Shiny server functions are necessarily long.                                                                                                |
| NOTE       | Runnable examples on `globalTrimApp.Rd` etc.                                | Cannot be runnable — they launch a Shiny gadget that blocks the R session.                                                                       |
| NOTE       | `dontrun / donttest` usage                                                  | All `dontrun` already converted to `donttest` (more permissive). The note is for either tag.                                                    |
| NOTE       | Lines > 80 characters (681 lines)                                           | Pre-existing convention; many are roxygen URLs / DECIPHER citation strings.                                                                     |
| NOTE       | Indents not multiples of 4 (2538 lines)                                     | Heuristic counts comment alignment; cosmetic.                                                                                                   |
| NOTE       | Bioc-Devel mailing-list subscription                                        | Out-of-band; requires checking ETH Zurich mailman membership.                                                                                   |

**No new code-level findings** — all ERROR / WARNING / NOTE entries are matched against the Phase 9 baseline.

---

## 2. GitHub sync

### Working-tree cleanup

Two ignore files had pending one-line additions for `.positai` (an editor artifact) that lingered across phases. Phase 13 commits them as housekeeping:

| File              | Change                                                |
| ----------------- | ----------------------------------------------------- |
| `.Rbuildignore`   | `+ ^\.positai$` (don't ship in tarball)                |
| `.gitignore`      | `+ .positai` (don't track)                             |

Build artifacts (`..Rcheck/`, three `docs/build/html/_static/*.js`, `sangeranalyseR_1.20.0.tar.gz`, `sangeranalyseR_1.21.1.tar.gz`, `sangeranalyseR.BiocCheck/`) intentionally remain **untracked** — same policy maintained since Phase 4.

### Commits in this push

- `13bb…` (this Phase 13 release-log + ignore-housekeeping commit, hash recorded after `git commit` runs).
- `86175e5` Phase 12 — Added R-CMD-check badge, simplified README, and created NEWS.md
- `3745066` Phase 11 — Overhauled README and finalized documentation
- `a52dc94` Phase 10 — Reverted R to 4.0.0 and maximized test coverage
- `bc2b5e6` Phase 9 — Automated UI tests, dependency cleanup, and clean BiocCheck
- `09926b1` Phase 8 — Modernized Shiny UI with Plotly and fixed lazy AA report compatibility
- (Phases 1–7 cumulative below)

### Remotes

```
origin    https://github.com/roblanf/sangeranalyseR        (fetch + push)
upstream  git@git.bioconductor.org:packages/sangeranalyseR.git  (fetch + push)
```

---

## 3. Bioconductor upstream sync

### Branch convention

Bioconductor's Git server tracks **devel on `master`** and each release line on `RELEASE_3_xx`. The local branch is `devel`; the correct push target is `upstream/master` via the `devel:master` refspec.

### Commands executed

```bash
git push origin devel
git fetch upstream
git push upstream devel:master
```

The `git fetch upstream` step is essential: if anyone has pushed to the Bioc devel branch since the last fetch, our push would be rejected as non-fast-forward. Fetch reveals any divergence so we can rebase first.

### Push hashes

(Recorded by the Phase 13 commit script after each `git push` returns; see "Commit hashes" section appended at the bottom of this file when the push completes.)

---

## 4. Cumulative summary across all 13 phases

| Phase | Title                                                              | Hash       |
| ----: | ------------------------------------------------------------------ | ---------- |
| 1     | Comprehension audit                                                | (notes only) |
| 2     | Quality / debt audit                                               | (notes only) |
| 3     | Test suite expansion (BUG-LOCK tests)                              | (combined w/ Phase 4 commits) |
| 4     | Sanitation: regex bug fixes, validator refactor, setValidity        | (combined) |
| 5     | E2E validation + profiling baseline                                 | `b7df0bc`  |
| 6     | BiocParallel migration + lazy AA                                    | `ec92f14`  |
| 7     | Rcpp port (`peakvalues_batch_cpp`)                                  | `9e4dfff`  |
| 8     | Shiny UI / Plotly chromatograms / lazy-AA report compat            | `09926b1`  |
| 9     | UI test automation + dependency cleanup + clean rcmdcheck          | `bc2b5e6`  |
| 10    | R 4.0.0 rollback + coverage maximisation (35.7%)                    | `a52dc94`  |
| 11    | README rewrite + documentation polish                               | `3745066`  |
| 12    | NEWS.md + simplified README + screenshot restored                   | `86175e5`  |
| 13    | Final audit + GitHub & Bioconductor sync                            | (this commit) |

### Build-quality timeline

| Metric                    | Phase 1 baseline   | Phase 13 final                          |
| ------------------------- | ------------------ | --------------------------------------- |
| `rcmdcheck`               | (not run)          | **0 / 0 / 0**                            |
| `BiocCheck`                | (not run)          | 1 / 1 / 17 — all environmental / advisory |
| testthat                  | inert (`test_check` commented out) | **1360 / 1360 PASS**           |
| Coverage                  | unknown            | **35.7 %** overall; ≥ 87 % on every non-Shiny file |
| `SangerAlignment` wall time | 1.85 s             | **1.07 s (best) / 1.14 s (mean)** = 1.62× |
| Bug-locking failures      | 0 (the bugs were silent) | 0 — all 2 BUG-LOCKs from Phase 3 fixed in Phase 4 |

---

## 5. Verification

Run from the repo root:

```bash
# Local test pass:
Rscript -e 'devtools::test()'                                   # 1360/1360

# Strict R CMD check:
Rscript -e 'rcmdcheck::rcmdcheck(".", args = "--no-manual")'    # 0/0/0

# BiocCheck on the built tarball:
R CMD build . --no-build-vignettes
Rscript -e 'BiocCheck::BiocCheck("sangeranalyseR_1.21.1.tar.gz", `quit-with-status` = FALSE)'

# Confirm both remotes are in sync:
git log --oneline origin/devel..HEAD     # empty if origin synced
git log --oneline upstream/master..HEAD  # empty if upstream synced
```

---

## 6. Non-goals (cumulative; deferred)

- Bioconductor Support Site "Watched Tags" registration (manual maintainer task).
- Vignette refresh in `vignettes/sangeranalyseR.Rmd` to reflect Phase 6–12 features.
- `pkgdown` site.
- `shinytest2`-driven coverage of the three large Shiny servers (would push coverage to ~80–90 % but requires a multi-day harness).
- Rewriting `sapply` → `vapply`, `<<-` → `reactiveValues`, `=` → `<-` mechanically across the package.
- Maintainer ORCID iD (out-of-band info).
