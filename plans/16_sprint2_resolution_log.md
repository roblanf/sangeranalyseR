# sangeranalyseR — Phase 16 Sprint 2 Resolution Log

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

Resolves GitHub Issues **#94**, **#66**, **#65**, **#42** as scheduled in the Phase 14 triage report. Also pivots the Sprint 1 reply mechanism from `gh` CLI bash script to a copy-paste Markdown doc (the local `gh` is too old to support `gh issue comment`).

---

## 1. Headline result

| Metric                  | Before Phase 16            | After Phase 16                  |
| ----------------------- | -------------------------- | ------------------------------- |
| Open issues addressed    | (none in Sprint 2 fixed)   | **4 (#94, #66, #65, #42)**       |
| Test count              | 1406 PASS                  | **1419 PASS** (+13 new Phase-16 tests; 0 regressions) |
| `rcmdcheck`              | 0 / 0 / 0                  | **0 / 0 / 0** (no regression after roxygen + import additions) |
| Reply format             | `.sh` (broken on local gh 0.11) | **Markdown copy-paste** (`plans/16_github_replies.md`) |

---

## 2. Issue-by-issue technical resolution

### Issue #94 + #66 — Low-overlap merges produce degenerate IUPAC consensus

These two issues bundle: #94 reports forced full-length overlap on minimal-overlap 16S, #66 reports the resulting silently-degenerate consensus. The fix is one consolidated change to `calculateContigSeq` in `R/UtilitiesFunc.R`:

**1. Forward DECIPHER alignment params via a new `alignSeqsParams = list()` argument.**

Pre-Phase-16, `AlignSeqs(frReadSet, processors = N, verbose = FALSE)` was called with a fixed argument set. Users had no knob to tune `iterations`, `refinements`, `gapOpening`, etc. Phase 16 introduces:

```r
base_align_args <- list(myXStringSet = frReadSet,
                        processors    = processorsNum,
                        verbose       = FALSE)
extra_align_args <- alignSeqsParams[
    !names(alignSeqsParams) %in% names(base_align_args)]
aln <- do.call(AlignSeqs, c(base_align_args, extra_align_args))
```

The `setdiff` on `names()` prevents user overrides of the three core args (`myXStringSet`, `processors`, `verbose`).

**2. Post-alignment overlap-quality check** (new `minOverlapFraction` and `minOverlapBases` arguments, both default 0 to preserve old behaviour):

```r
nongap <- as.matrix(aln) != "-"
ovl <- crossprod(t(nongap))     # nreads x nreads pairwise non-gap counts
diag(ovl) <- NA_integer_
overlap_min_obs <- min(ovl, na.rm = TRUE)
threshold <- max(minOverlapFraction * shorter_pair_len,
                  as.numeric(minOverlapBases))
if (overlap_min_obs < threshold) {
    log_warn(">> LOW_OVERLAP_WARN: smallest pairwise overlap is ",
             overlap_min_obs, " bp ...")
}
```

The `crossprod(t(nongap))` is a single C-level matrix multiplication that yields the full pairwise non-gap-overlap matrix in one call.

Both new arguments thread through `SangerContig.initialize`, `SangerAlignment.initialize`, and the `Constructors.R` wrappers.

**Tests** (in `tests/testthat/test-Phase16-Sprint2.R`):

- Two synthetic 200-bp reads (one all-A, one all-T) with `minOverlapFraction = 0.5` and `minOverlapBases = 50L` triggers `LOW_OVERLAP_WARN`.
- Real ACHLO data (~600bp overlap) with the same thresholds does **not** trigger the warning (regression check).
- Passing a bogus `alignSeqsParams = list(bogus_argument_xyz = 42L)` triggers a DECIPHER-level error (proves the args reach `AlignSeqs`).
- End-to-end `SangerContig(..., minOverlapBases = 50L)` builds successfully.
- Default thresholds (0 / 0L) leave the consensus unchanged — backward compatibility regression test.

### Issue #42 — `minReadLength` doesn't filter length-1 reads

The Phase-3 `minReadLength` filter in `SangerContig.initialize` handles the typical case, but degenerate-trim states (M1 producing `trimmedFinishPos = 0` while `trimmedStartPos > 0`) could yield a length-1 entry in `frReadSet` after the SangerContig-level filter. That entry then crashed DECIPHER's `AlignSeqs` (which requires width ≥ 2).

**Fix**: defensive pre-alignment width filter in `calculateContigSeq`:

```r
too_short <- BiocGenerics::width(frReadSet) < 2L
if (any(too_short)) {
    log_warn(">> Dropping ", sum(too_short),
             " read(s) with trimmed length < 2 bp ",
             "(MIN_READ_LENGTH_DEFENSIVE_DROP).")
    frReadSet <- frReadSet[!too_short]
}
if (length(frReadSet) < 2L) {
    log_warn(">> Fewer than 2 usable reads after defensive filter; ",
             "returning empty consensus.")
    if (length(frReadSet) == 1L) {
        consensusGapfree <- frReadSet[[1L]]
    } else {
        consensusGapfree <- DNAString()
    }
    return(list("consensusGapfree" = consensusGapfree, ...))
}
```

The early return is essential — `log_error` doesn't `stop()`, so without an early return execution would still flow into `AlignSeqs` and crash on n=1.

**Test**: surgically mutate one forward `SangerRead`'s `primarySeq` to `"A"` (length 1) and `trimmedFinishPos = 1`; assert `calculateContigSeq` still produces a non-empty consensus (the bad read is dropped silently).

### Issue #65 — Reads assigned to >1 contig in CSV

**Fix**: `checkAb1FastaCsv` now detects when a single `reads` filename is mapped to multiple distinct `contig` values:

```r
read_to_contig <- aggregate(
    as.character(csvFile$contig),
    by   = list(read = as.character(csvFile$reads)),
    FUN  = function(v) length(unique(v))
)
multi_assigned <- read_to_contig$read[read_to_contig$x > 1L]
if (length(multi_assigned) > 0L) {
    msg <- paste0("Read(s) assigned to >1 distinct contig in CSV: ",
                  paste(sQuote(multi_assigned), collapse = ", "),
                  ". Each read must belong to exactly one contig.")
    log_warn(msg)
    warnings <- c(warnings, paste0(msg, " (READ_ASSIGNED_MULTIPLE_CONTIGS_WARN)"))
}
```

Logged as a warning rather than promoted to an error — some users deliberately reuse the same read across contigs (rare but legitimate); a warning is enough to surface the issue.

**Tests**:

- CSV with one read assigned to two contigs → warning logged.
- CSV with no duplicates → no warning (regression check).

---

## 3. Reply-format pivot

The Phase 15 deliverable `plans/github_replies_sprint1.sh` used `gh issue comment <#> --body "..."` and `gh issue close <#>`. The local `gh` CLI is **v0.11.0 (2020-07-16)** — predates both the `comment` subcommand and the `--body` flag. Running the script would have failed.

**Phase 16 replacement**: `plans/16_github_replies.md` contains all 8 replies (Sprint 1 reposted + Sprint 2 new) as **copy-paste-ready Markdown blocks**. For each issue:

1. Clickable GitHub URL (one-click navigation to the issue page).
2. The reply body inside a fenced ```` ```markdown ```` code block — copy the block contents (the GitHub-rendered "Copy" icon does this), paste into the issue's comment box, click **Comment**, then click **Close issue**.

This works regardless of `gh` CLI version and is also auditable post-hoc — anyone reviewing the project history can see the exact text that was posted.

---

## 4. Files touched in Phase 16

```
M  R/UtilitiesFunc.R                     (#94/#66 overlap check + alignSeqsParams; #42 defensive filter + early return)
M  R/UtilitiesFuncInputChecker.R          (#65 duplicate-read warning)
M  R/ClassSangerContig.R                  (propagate new constructor args)
M  R/ClassSangerAlignment.R               (propagate new constructor args x4 nested new() sites)
M  R/Constructors.R                       (add new args to public wrappers + roxygen)
M  R/sangeranalyseR_package.R             (+ BiocGenerics::width; + stats::aggregate)
M  DESCRIPTION                             (+ BiocGenerics in Imports)
M  NAMESPACE                               (auto-regen)
A  tests/testthat/test-Phase16-Sprint2.R   (8 test_thats / 13 assertions)
A  plans/16_github_replies.md              (consolidated Sprint 1 + 2 reply Markdown)
A  plans/16_sprint2_resolution_log.md
```

---

## 5. Verification

```r
# Targeted Sprint 2 tests:
testthat::test_file("tests/testthat/test-Phase16-Sprint2.R")
# -> 13 / 13 PASS

# Full suite (no regressions):
devtools::test()
# -> 1419 / 1419 PASS

# Strict R CMD check (rebuild after new params + imports):
roxygen2::roxygenise()
rcmdcheck::rcmdcheck(".", args = "--no-manual")
# -> 0 / 0 / 0
```

---

## 6. Sprint 3 preview (per Phase 14 triage)

Issues remaining for the next sprint — all touch the `ConsensusSequence` integration point in `calculateContigSeq` and bundle naturally:

- **#87** Majority-rules consensus base calling option.
- **#48** Phred-aware consensus building.
- **#33** Add Phred quality to consensus output.

A separate documentation pass closes:

- **#13** Worked examples.
- **#49** "How to" docs section.
- **#71** Base-calling method explanation.
- **#99** Parameters tutorial.
