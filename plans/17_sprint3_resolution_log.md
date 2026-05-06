# sangeranalyseR — Phase 17 Sprint 3 Resolution Log

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

Resolves GitHub Issues **#87**, **#48**, **#33** as scheduled in the Phase 14 triage. All three touch `calculateContigSeq`'s consensus phase and bundle naturally into one architectural change.

---

## 1. Headline result

| Metric                   | Before Phase 17    | After Phase 17                     |
| ------------------------ | ------------------ | ---------------------------------- |
| Open issues addressed     | n/a                | **3 (#87, #48, #33)**               |
| Test count               | 1419 PASS          | **1451 PASS** (+32 new; 0 regressions) |
| `rcmdcheck`               | 0 / 0 / 0          | **0 / 0 / 0** (no regression)       |
| Backwards compatibility  | n/a                | **Default `consensusMethod = "strict"` is unchanged** — pre-Phase-17 callers see no behaviour difference |

---

## 2. Architecture overview

Two new helpers + one new branch in `calculateContigSeq`, plus pluggable `consensusMethod` / `qualityAware` constructor arguments threaded all the way out to the `SangerAlignment()` and `SangerContig()` user wrappers.

### New helpers in `R/UtilitiesFunc.R`

| Helper                          | Purpose                                                                                                    |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `.computeConsensusMajority(aln, weights = NULL)` | Per-column plurality vote; ties break alphabetically. Synthesises per-position Phred = `40 * winner_count / total_count`, or mean Phred of agreeing reads when `weights` is supplied. Returns `list(consensus = DNAString, qualityScores = integer)`. |
| `.buildQualityMatrix(aln, qualityPhredScoresList)` | Maps each read's per-base Phred onto the alignment columns. Returns an `nreads × ncols` integer matrix; gap columns get 0; unknown reads fall back to flat Phred 30 with a warning. |

### `calculateContigSeq` branching

```r
if (consensusMethod == "strict") {
    consensus <- ConsensusSequence(aln, ...)[[1L]]   # pre-Phase-17 path
    consensusQualityScores <- integer(0)
} else {
    weights_mat <- if (consensusMethod == "quality_weighted")
        .buildQualityMatrix(aln, qualityPhredScoresList) else NULL
    res  <- .computeConsensusMajority(aln, weights = weights_mat)
    consensus              <- res$consensus
    consensusQualityScores <- res$qualityScores
}
```

After gap-stripping, the quality vector is subset to non-gap positions so `length(qualityScores) == length(consensusGapfree)` and stored both:

- as **`attr(consensusGapfree, "qualityScores")`** — survives serialisation, accessible from `attributes(sc@contigSeq)$qualityScores`.
- as **`CSResult$consensusQualityScores`** — top-level list element of the `calculateContigSeq` return.

### Constructor plumbing

`SangerContig.initialize` builds a `qualityPhredScoresList` from the filtered SangerReads (pulling `@QualityReport@qualityPhredScores`, sliced to `[trimmedStartPos+1, trimmedFinishPos]`, reversed for reverse reads) only when the user requests `qualityAware = TRUE` or `consensusMethod = "quality_weighted"`. ABIF only — FASTA inputs have no Phred scores; the `.buildQualityMatrix` fallback emits flat Phred 30 with a warning.

`SangerAlignment.initialize` forwards both new args to all 4 nested `new("SangerContig", ...)` call sites.

`Constructors.R` wrappers expose both args with documented roxygen blocks.

---

## 3. Issue-by-issue technical resolution

### Issue #87 — Majority-rules consensus base calling

**Root cause**: pre-Phase-17 the only path was `DECIPHER::ConsensusSequence(ambiguity = TRUE)`, which encodes per-column disagreements as IUPAC codes. Users wanting plain ACGT had no escape hatch.

**Fix**: new `consensusMethod = "majority"` mode. The `.computeConsensusMajority` helper does:

```r
for (j in 1:ncols) {
    bases <- aln_matrix[, j]
    bases <- bases[bases != "-"]
    if (length(bases) == 0L) cons[j] <- "-"
    else cons[j] <- names(table(bases))[which.max(table(bases))]
}
```

Tie-breaking is whatever `which.max` does on a `table()` — alphabetical because `table()` returns sorted names. Documented in the user reply.

**Tests**:

- Synthetic 3-row × 6-col alignment with one disagreeing column (`AAGCTT`/`AAGCTT`/`AAGGTT`) → consensus is `AAGCTT`, with 4/6 columns at Phred 40 and the disagreeing column at Phred 27 (= 40 × 2/3).
- All-gap column → consensus character is `-`, quality is 0.
- End-to-end `SangerContig(consensusMethod = "majority")` on real ACHLO data → consensus contains no IUPAC ambiguity codes (`expect_false(grepl("[YRSWKMBDHVN]", cs))`).

### Issue #48 — Phred-aware consensus

**Root cause**: even with majority voting, low-quality calls and high-quality calls counted equally, which is wrong when one read is much more reliable than another.

**Fix**: new `consensusMethod = "quality_weighted"` mode (alias `qualityAware = TRUE`). `.buildQualityMatrix` constructs an `nreads × ncols` integer matrix mapping each read's `@QualityReport@qualityPhredScores` (post-trim) to alignment columns. `.computeConsensusMajority` with non-NULL `weights` then sums Phred-weighted votes per base; the base with the highest total wins.

**Tests**:

- 3-read column with bases `A, A, G` and weights `10, 10, 60`: majority says `A` (2/3); quality-weighted says `G` (60 > 20). The two paths produce different answers — the test asserts both behaviours separately.
- `.buildQualityMatrix` unit test — for a 2-read alignment with gaps in different columns, the matrix correctly encodes per-column Phred values, with 0 in gap columns.
- `.buildQualityMatrix` fallback to flat Phred 30 when the read isn't in the supplied list.
- End-to-end `SangerContig(qualityAware = TRUE)` on real ACHLO data builds successfully.

### Issue #33 — Add quality to consensus

**Root cause**: the package returned only the consensus DNA, not any per-position confidence signal — making downstream filtering (e.g., "mask positions where < 80% of reads agree") impossible without re-deriving from the alignment.

**Fix**: under `consensusMethod = "majority"` or `"quality_weighted"`, attach per-position quality scores as `attr(consensusGapfree, "qualityScores")` (gap-aligned to the gap-stripped consensus). Strict mode preserves the pre-Phase-17 behaviour: `attr(...)` is absent or `integer(0)`.

This is **not** a deep per-base re-calibration like the OverlapPER paper referenced in #32 — it's a synthetic confidence score from the consensus voting itself. Documented as such in the user reply.

**Tests**:

- Majority mode: `length(attr(sc@contigSeq, "qualityScores")) == length(sc@contigSeq)`, all values in `[0, 60]`.
- Strict mode: attribute is either NULL or `integer(0)`.

---

## 4. Backwards compatibility

The default `consensusMethod = "strict"` is **byte-for-byte unchanged** from Phase 16. Verified by the regression test:

> `Sprint 3 default (consensusMethod='strict') still produces the same DECIPHER consensus`

All 1419 pre-Phase-17 tests continue to pass. Cumulative suite: 1419 → **1451 PASS**.

---

## 5. Files touched in Phase 17

```
M  R/UtilitiesFunc.R                       (+ .computeConsensusMajority, .buildQualityMatrix; calculateContigSeq branching; consensusQualityScores in return list + as attr())
M  R/ClassSangerContig.R                    (propagate consensusMethod / qualityAware; build qualityPhredScoresList)
M  R/ClassSangerAlignment.R                 (forward to 4 nested new() sites)
M  R/Constructors.R                         (add args + roxygen for SA + SC wrappers)
M  man/SangerAlignment.Rd, SangerContig.Rd  (auto-regen)
A  tests/testthat/test-Phase17-Sprint3.R    (11 test_thats / 32 assertions)
A  plans/17_github_replies.md               (3 issue replies; same Markdown format as Phase 16)
A  plans/17_sprint3_resolution_log.md
```

---

## 6. Reply format compatibility

The new replies file `plans/17_github_replies.md` uses the **same fence structure** as `plans/16_github_replies.md`. The Phase-16.5 closure script handles it without modification:

```bash
$ python3 plans/close_issues.py --md plans/17_github_replies.md --dry-run
Parsed 3 issue(s) from plans/17_github_replies.md:
  • #87   (1756 chars)  first line: 'Hi, thanks for the clear write-up...'
  • #48   (1389 chars)  first line: 'Hi, finally addressing this one...'
  • #33   (1507 chars)  first line: 'Thanks for the cross-link to #32...'
```

To post and close:

```bash
export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx
python3 plans/close_issues.py --md plans/17_github_replies.md
```

---

## 7. Verification

```r
# Targeted Sprint 3 tests:
testthat::test_file("tests/testthat/test-Phase17-Sprint3.R")
# -> 32 / 32 PASS

# Full suite:
devtools::test()
# -> 1451 / 1451 PASS

# Strict R CMD check:
roxygen2::roxygenise()
rcmdcheck::rcmdcheck(".", args = "--no-manual")
# -> 0 / 0 / 0
```

---

## 8. Sprint 4 preview

Per Phase 14 triage, the remaining open issues bundle into a documentation pass:

- **#13** Make worked examples.
- **#49** Add documentation "How to ..." section.
- **#71** Base-calling method explanation.
- **#99** Parameters tutorial / contig-creation guide.

Pure docs work — likely a single PR to `vignettes/sangeranalyseR.Rmd` plus the existing ReadTheDocs site, no R changes.
