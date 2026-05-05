# sangeranalyseR — Phase 14 Issue Backlog Triage

Source: GitHub Issues at https://github.com/roblanf/sangeranalyseR/issues, snapshot taken 2026-05-04.

**Backlog totals:** 95 issues lifetime → 52 closed, **43 open**. The 43 open issues are triaged below. Phases 1–13 already implicitly closed several of them; those are flagged `Resolved via Refactor`.

---

## 1. Triage table (43 open issues)

Categories: **Critical Bug** (breaks the pipeline), **Enhancement/Feature** (new capability), **Documentation** (text-only), **Stale/Irrelevant** (no longer applicable / fixed by refactor / abandoned).

Status: `Needs Fix`, `Needs Triage` (more info from reporter required), `Resolved via Refactor` (Phases 1–13 already addressed it).

| #   | Date       | Title                                                                             | Category             | Status                  | Priority |
| --: | ---------- | --------------------------------------------------------------------------------- | -------------------- | ----------------------- | :------: |
| 100 | 2026-03-09 | SangerAlignment CSV/ABIF returns CONTIG_NUMBER_ZERO_ERROR though SangerContigs succeed | Critical Bug         | Needs Fix               | **#1**   |
| 99  | 2026-02-03 | Create contig — user can't get parameters right                                   | Documentation        | Needs Triage            |    —     |
| 98  | 2025-10-01 | Ongoing development and maintenance                                               | Stale/Irrelevant     | **Resolved via Refactor** (Phases 1–13) |    —     |
| 97  | 2025-07-11 | Move `data.table` from `Depends:` to `Imports:`                                   | Enhancement          | **Resolved via Refactor** (Phase 9 moved `data.table` to `Imports:`) |    —     |
| 95  | 2023-12-18 | SCF compatibility                                                                  | Enhancement          | Needs Fix               |    —     |
| 94  | 2023-12-06 | Minimal-overlap F + R 16S reads (~50–100 bp overlap rejected)                     | Critical Bug         | Needs Fix               | **#3**   |
| 93  | 2023-11-24 | Batch SangerRead instance generation from AB1                                     | Enhancement          | Needs Fix               |    —     |
| 92  | 2023-10-02 | Forward-only reads rejected with `REGEX_SuffixReverse must be character type`     | Critical Bug         | Needs Fix               | **#2**   |
| 91  | 2023-03-23 | M2 trimming errors on low-quality files                                           | Critical Bug         | **Resolved via Refactor** (Phase 4 relaxed `setValidity`; degraded-state guards) |    —     |
| 90  | 2023-02-27 | Support DNA reference sequences for non-coding genes (16S)                         | Enhancement          | Needs Fix               |    —     |
| 89  | 2023-02-08 | `writeFasta()` errors on contigs with only 1 read (`writeXStringSet 'x' must be an XStringSet`) | Critical Bug | Needs Fix |    —     |
| 87  | 2023-01-31 | Majority-rules consensus base calling option                                      | Enhancement          | Needs Fix               |    —     |
| 86  | 2023-01-16 | `writeFasta(outputDir = NULL)` writes to `/tmp` silently                          | Enhancement          | Needs Fix               |    —     |
| 85  | 2022-11-02 | Reproducible-tutorial errors (vignette break)                                     | Documentation        | **Likely Resolved via Refactor** (Phase 4 regex bug fixes; Phase 11/12 README rewrites) |    —     |
| 84  | 2022-08-11 | Variant calling from reference gene (codon-level annotation)                       | Enhancement          | Needs Fix               |    —     |
| 82  | 2022-05-31 | `Error in TreeLine: could not find function "TreeLine"`                           | Critical Bug         | **Resolved via Refactor** (the DECIPHER `Treeline` import fix that was already on `RELEASE_3_22` / `devel`) |    —     |
| 76  | 2021-11-03 | `qualityPhredScores length cannot be zero` on certain ABIFs                       | Critical Bug         | Needs Fix               | **#4**   |
| 75  | 2021-10-15 | Use of trimming functions outside interactive Shiny — request stable API          | Enhancement          | Needs Fix               |    —     |
| 74  | 2021-10-07 | `launchApp()` and `writeFasta` errors                                              | Critical Bug         | Needs Triage            |    —     |
| 72  | 2021-09-22 | `generateReportSA()` "subscript out of bounds" on contigList[[1]]                 | Critical Bug         | **Likely Resolved via Refactor** (Phase 8 lazy-AA RMD fixes); needs verification | —     |
| 71  | 2021-06-25 | Question: which base-calling method is used                                       | Documentation        | Needs Fix (doc clarification) | —    |
| 68  | 2021-04-19 | `pandoc` font-fetch failure during `generateReport()`                             | Critical Bug         | **Resolved via Refactor** (Phase 8 lazy-AA RMD changes; needs reporter retest) |    —     |
| 66  | 2021-01-26 | Improper merge / contig with M1 trimming default                                  | Critical Bug         | Needs Fix               | **#5**   |
| 65  | 2021-01-22 | Handle the case where a read is assigned to >1 place                              | Enhancement          | Needs Fix               |    —     |
| 61  | 2020-10-21 | Launch Shiny app error (Windows R 4.0)                                            | Critical Bug         | **Likely Resolved via Refactor** (Phase 6 BiocParallel + Phase 9 dependency cleanup); needs reporter retest | —     |
| 60  | 2020-10-13 | Shiny `addResourcePath` cannot normalize path                                     | Critical Bug         | **Likely Resolved via Refactor** (Phase 9 cleaned up shinydashboard imports); needs reporter retest |    —     |
| 55  | 2020-07-09 | SangerAlignment Error enhancement (#52 follow-up)                                 | Enhancement          | **Mostly Resolved via Refactor** (Phase 3 added unit tests; Phase 9 added Shiny tests) |    —     |
| 50  | 2020-06-11 | Output regex matching table to user as part of logging                            | Enhancement          | Needs Fix               |    —     |
| 49  | 2020-05-26 | Add documentation "How to ..." section                                            | Documentation        | Needs Fix               |    —     |
| 48  | 2020-05-25 | Consider Phred scores when consensus building                                     | Enhancement          | Needs Fix               |    —     |
| 43  | 2020-04-16 | `readTable()` function spec                                                        | Enhancement          | **Partly Resolved** (Phase 3+ `readTable,SR/SC` exists; SA-level still missing) |    —     |
| 42  | 2020-04-15 | `minReadLength` doesn't filter very-low-quality (length-1) reads                  | Critical Bug         | Needs Fix               |    —     |
| 41  | 2020-04-14 | Low-quality reads trimming default (keep 0 bp instead of 1)                       | Enhancement          | Needs Fix               |    —     |
| 38  | 2020-04-14 | Trees should be unrooted                                                           | Enhancement          | Needs Fix               |    —     |
| 37  | 2020-04-13 | Spurious indels — does package handle them                                         | Documentation/Bug    | Needs Triage            |    —     |
| 34  | 2020-03-18 | Detect 16S indels without reference                                                | Enhancement          | Needs Fix               |    —     |
| 33  | 2020-03-18 | Add Phred quality to consensus output                                             | Enhancement          | Needs Fix               |    —     |
| 31  | 2020-02-10 | OS testing (Mac / Linux / Windows)                                                | Documentation        | **Resolved via Refactor** (Phase 6 BiocParallel + Phase 9 cross-platform CI) |    —     |
| 21  | 2017-05-12 | Build a new sequence object (primary, secondary, quality, trim)                   | Stale/Irrelevant     | **Resolved via Refactor** (the modern `SangerRead` S4 class already does this) |    —     |
| 13  | 2016-05-06 | Make some worked examples                                                          | Documentation        | Needs Fix               |    —     |
| 12  | 2016-05-06 | Exclude reads based on quality scores                                              | Enhancement          | Needs Fix               |    —     |
| 11  | 2016-05-06 | Count coincident secondary peaks                                                   | Enhancement          | **Resolved via Refactor** (`countCoincidentSp` exists in `R/UtilitiesFunc.R`) |    —     |
| 10  | 2016-05-06 | Add a function to trim to reference                                                | Enhancement          | Needs Fix               |    —     |

### Quick category counts

| Category               | Open  | Of which "Resolved via Refactor" | Net "Needs Fix"  |
| ---------------------- | ----: | -------------------------------: | ---------------: |
| Critical Bug            |  14  | 4 fully + 3 likely               | **7 confirmed**  |
| Enhancement / Feature   |  20  | 5 fully + 1 partly               | 14               |
| Documentation           |   6  | 1 fully                          | 5                |
| Stale / Irrelevant      |   3  | 3                                | 0                |
| **Total**               |  43  | 13 fully + 4 likely              | 26 actionable    |

---

## 2. Top 5 priority issues (immediate-action shortlist)

Ranked by **user-impact severity** (does it break the pipeline?) **× breadth** (how many users hit it?) **× actionability under the new architecture**.

### #1 — Issue **#100**: `SangerAlignment` CSV/ABIF returns `CONTIG_NUMBER_ZERO_ERROR` even though every individual `SangerContig` succeeds *(2026-03-09)*

**Root cause (suspected):** The CSV-driven path in `SangerAlignment.initialize` (`R/ClassSangerAlignment.R:358–413`) collapses `SangerContig` results into a list, but a per-contig failure during `alignContigs` (Phase 7's batch paths or `BPPARAM` propagation) may be filtering them all out. The Phase-3 `test-Validator-EdgeCases.R::"CSV with reads not on disk produces a warning"` proves the warning path; this looks like the missing **error-recovery** path on the success-but-aggregation-fails branch.

**Why #1:** Reported on the `devel` branch we just modernized — this is the most recent regression-shaped report, against the post-Phase-7 architecture we own. High community visibility (most recent issue).

**Estimated effort:** 1–2 days. Investigate the `SangerContigList` filter at `R/ClassSangerAlignment.R:518` (`Filter(Negate(is.null), SangerContigList)`) plus the `contigNum == 0` branch trigger.

### #2 — Issue **#92**: Forward-reads-only rejected with `REGEX_SuffixReverse must be character type` *(2023-10-02; reproduced 2024 on Bioc forum)*

**Root cause:** `checkREGEX_SuffixReverse` (`R/UtilitiesFuncInputChecker.R:285`) is unconditionally called in `SangerAlignment.initialize`. There's no opt-out for "forward only" datasets — common for 16S short-read sequencing. The fix is to allow `REGEX_SuffixReverse = NA` (or a sentinel like `""`) and skip the reverse-grouping branch when it's set.

**Why #2:** **Many** users hit this — 16S, COI-with-only-F-primer, and other half-coverage datasets are common in microbial / barcode pipelines. Closed issue **#78** is the same complaint, and #92 references that. Multiple repeat reports.

**Estimated effort:** 0.5–1 day. Add a `forwardOnly = FALSE` argument or accept `NA`/empty `REGEX_SuffixReverse`, gate the reverse construction loops on it. Plus a Phase-3-style negative test.

### #3 — Issue **#94**: Minimal-overlap F + R 16S reads forced to full-length overlap *(2023-12-06)*

**Root cause:** `calculateContigSeq` in `R/UtilitiesFunc.R:182` calls `DECIPHER::AlignSeqs` directly. That function aligns globally; without a `gapOpening` / `gapExtension` knob exposed at the user level, short-overlap pairs are squeezed into a forced alignment. Issue **#90** reports the same DECIPHER-pass-through limitation.

**Why #3:** The dominant Sanger-amplicon use case for 16S rRNA. Users are silently getting *wrong* contigs (degenerate consensus); this is worse than a hard error.

**Estimated effort:** 1–2 days. Either expose `AlignSeqs` knobs (e.g. `iterations`, `gapOpening`) as `SangerContig()` arguments, or detect minimal-overlap pairs ahead of alignment and use a primer-aware merge instead.

### #4 — Issue **#76**: `qualityPhredScores length cannot be zero` on certain ABIFs *(2021-11-03; corroborated 2023+)*

**Root cause:** Some ABIF files contain unusual data blocks (`unimplemented legacy type found in file`) where `sangerseqR::read.abif` succeeds but the `PCON.2` (per-base quality) data block is empty. The Phase-4 `checkQualityPhredScores` validator now catches this, but it's a hard error — there's no fallback path that uses the raw trace's primary base calls when quality scores are missing.

**Why #4:** Affects **batches** of files from particular sequencing instruments (Beckman / older 3500 firmware). Same problem family as the SCF compatibility request (#95).

**Estimated effort:** 1–2 days. Add an opt-in `acceptMissingQuality = FALSE` flag that, when TRUE, synthesises a flat Phred score (e.g. 30) for every base so the rest of the pipeline can run. Document that contigs from such reads have no quality-aware trimming.

### #5 — Issue **#66**: Improper merge / contig with default M1 cutoff *(2021-01-26)*

**Root cause:** With the default `M1TrimmingCutoff = 0.0001`, low-quality reads can be trimmed so aggressively that their remaining bases don't overlap. The current pipeline silently merges them anyway, producing a degenerate consensus (lots of IUPAC ambiguity codes). Reported as: "I verified that there is no similarity between the forward and reverse trimmed reads … why is a merge happening? none should happen."

**Why #5:** Silent quality-of-output regression. Users blame the package even when the trim parameters are the real issue. Documentation alone won't fix it; the pipeline should detect and warn.

**Estimated effort:** 1–2 days. Add an "alignment quality" check after `calculateContigSeq` — if the forward / reverse trimmed reads share < N bp of pairwise alignment match, log a `LOW_OVERLAP_WARN` and refuse to merge (or merge with a special flag set on the `SangerContig` `objectResults`). Issue **#65** ("Handle a read assigned to >1 place") is in the same family and could be folded into the same robustness pass.

---

## 3. "Resolved via Refactor" — issues to close immediately

These 13 + 4 likely-resolved issues have already been addressed by Phases 1–13. **Recommended action: post a closing comment referencing the relevant phase and close.**

| #   | Title                                                  | Resolved by                                                                  |
| --: | ------------------------------------------------------ | ----------------------------------------------------------------------------- |
| 11  | Count coincident secondary peaks                        | `countCoincidentSp` already exists in `R/UtilitiesFunc.R`.                    |
| 21  | Build a new sequence object (primary/secondary/quality/trim) | `SangerRead` S4 class delivers exactly this.                              |
| 31  | OS testing                                             | Phase 6 BiocParallel makes Windows first-class; Phase 9 cleaned cross-platform deps. |
| 55  | SangerAlignment Error enhancement                      | Phase 3 unit tests + Phase 9 Shiny tests cover the original sub-tasks.        |
| 82  | `Error in TreeLine`                                    | Already fixed by the Phase-pre-1 DECIPHER `Treeline` import bump on `devel`.  |
| 91  | M2 trimming errors on low-quality files                | Phase 4 relaxed `QualityReport setValidity` to accept the degenerate-window case. |
| 97  | Move `data.table` to `Imports:`                        | Phase 9 dependency cleanup did this for `data.table` and 23 others.            |
| 98  | Ongoing development and maintenance                    | Phases 1–13 themselves answer this issue; close with a pointer to `NEWS.md`.  |

Likely resolved (request reporter to retest on `devel`):

| 60, 61 | Shiny `addResourcePath` / launch errors           | Phase 9 cleaned shinydashboard / shiny imports; new Phase-9 testServer coverage. |
| 68     | pandoc font-fetch failure                          | Phase 8 RMD fixes + lazy-AA compat.                                              |
| 72     | `subscript out of bounds` in generateReportSA      | Phase 8 RMD fixes — `params$SangerAlignment@contigList[[1]]` is no longer broken with lazy-AA. |
| 85     | Reproducible-tutorial errors                        | Phase 4 regex bug fixes + Phase 11/12 README rewrites.                          |

---

## 4. Recommended sprint plan

**Sprint 1 (this iteration, 1–2 weeks):** Issues **#100, #92, #76, #89** — the 4 confirmed pipeline-breaking bugs, all individually small fixes, all in scope under the modern architecture.

**Sprint 2 (next iteration, 2–3 weeks):** Issues **#94, #66, #66's siblings #65 and #42** — the contig-merge / minimum-overlap robustness pass. These need the same "alignment-quality post-check" infrastructure, so they're naturally bundled.

**Sprint 3 (longer, 3–4 weeks):** Enhancements #87 (majority-rules consensus), #48 (Phred-aware consensus), #33 (quality in consensus output) — all touch the same `ConsensusSequence` integration point in `calculateContigSeq` and would benefit from a single rewrite.

**Documentation pass (parallel):** Close out #13, #49, #71, #99 with a single PR to the vignette + a "How to" section in the existing ReadTheDocs site.

---

## 5. Methodology + sources

- Issue list pulled from the GitHub REST API on 2026-05-04: `GET /repos/roblanf/sangeranalyseR/issues?state=open&per_page=100`. 43 open items returned (excluding pull requests, of which there are currently 0).
- Each open issue's body inspected to identify root cause and assess against post-Phase-13 codebase state.
- "Resolved via Refactor" verdicts cross-referenced against the Phase 1–13 logs in `plans/01_…` through `plans/13_…`.
- This file is the deliverable; full per-issue body extracts archived in `/tmp/issues_open.txt` during triage.
