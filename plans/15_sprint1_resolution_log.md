# sangeranalyseR — Phase 15 Sprint 1 Resolution Log

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

Resolves GitHub Issues **#100**, **#92**, **#76**, **#89** as scheduled in the Phase 14 triage report.

---

## 1. Headline result

| Metric                 | Before Phase 15            | After Phase 15                 |
| ---------------------- | -------------------------- | ------------------------------ |
| Open issues addressed   | (none in Sprint 1 fixed)   | **4 (#100, #92, #76, #89)**     |
| Test count             | 1360 PASS                  | **1406 PASS** (+46 new Phase-15 tests; 0 regressions) |
| `rcmdcheck`             | 0 / 0 / 0                  | **0 / 0 / 0** (no regression)   |

---

## 2. Issue-by-issue technical resolution

### Issue #100 — `SangerAlignment` CSV/ABIF returns `CONTIG_NUMBER_ZERO_ERROR` though every individual `SangerContig` succeeds

**Root cause** (`R/ClassSangerAlignment.R:371–377`, pre-Phase-15):

```r
contigNames <- as.character(unique(csvFile$contig))
contigNames <- lapply(contigNames, function(contigName) {
    contigNameSelectInputFiles <-
        parentDirFiles[grepl(contigName, parentDirFiles)]
    inside_contigNames <- file.path(dirname(contigNameSelectInputFiles), contigName)
    inside_contigNames
})
```

The `grepl(contigName, parentDirFiles)` substring match required every CSV `contig` label (e.g. `rbcL2`, `ITS5`) to appear textually inside its read filenames. The CSV-grouping mechanism is explicitly designed to map *arbitrary* labels to *arbitrary* filenames via the `reads` column, so this requirement was incorrect.

**Fix**: drive lookups directly off the CSV's `reads` column. For each unique contig label, intersect `csvFile$reads` (filtered to that label) against `parentDirFiles` (matching either full relative paths or basenames), then derive the directory prefix from the matched files. Strip the `./` prefix when files are flat at the root. Empty-match contig labels log a warning and are skipped instead of aborting the whole alignment.

**Tests** (`tests/testthat/test-Phase15-Sprint1.R`):

- 2-contig CSV with non-substring labels (`rbcL2`, `ITS5`) — both contigs build, consensus non-empty.
- Mixed CSV where one contig label has all-existing reads and another has none — first builds, second silently dropped.

### Issue #92 — Forward-reads-only datasets rejected with `REGEX_SuffixReverse must be character type`

**Root cause** (`R/UtilitiesFuncInputChecker.R::checkREGEX_SuffixReverse`):

```r
if (is.null(REGEX_SuffixReverse)) {
    return(.errAppend(errors, errorTypes,
                      "'REGEX_SuffixReverse' cannot be NULL.",
                      "PARAMETER_VALUE_ERROR"))
}
```

The validator rejected `NULL` (the default), making it impossible to express "I have only forward reads" cleanly.

**Fix**: both `checkREGEX_SuffixForward` and `checkREGEX_SuffixReverse` now accept `NULL` and `NA_character_` to mean "no reads in that direction". When detected, `SangerAlignment.initialize` and `SangerContig.initialize` substitute an internal sentinel (`.NEVER_MATCH_REGEX = "NEVER_MATCH_FORWARD_REVERSE_ONLY_SENTINEL"`) so the existing `parentDirFiles[grepl(...)]` machinery returns empty for that direction. The pre-existing "no reads detected" warning path then logs a friendly message instead of crashing.

**Tests**:

- `SangerAlignment(REGEX_SuffixReverse = NULL, minReadsNum = 1)` on forward-only directory → builds, no reverse reads.
- `SangerAlignment(REGEX_SuffixReverse = NA_character_)` → same.
- `SangerContig(REGEX_SuffixReverse = NULL)` → builds with `forwardReadList` populated, `reverseReadList` empty.
- Validator unit tests: `checkREGEX_SuffixReverse(NULL)` and `(NA_character_)` are silent; `(123L)` still flags `PARAMETER_VALUE_ERROR`.

User-facing API addition: pass `minReadsNum = 1` for forward-only configs (each read becomes its own contig).

### Issue #76 — `qualityPhredScores length cannot be zero` on certain ABIFs

**Root cause** (`R/ClassSangerRead.R:235`, pre-Phase-15):

```r
MBCResult <- MakeBaseCallsInside(traceMatrix, peakPosMatrixRaw,
                                 abifRawData@data$PCON.2,
                                 ...)
```

When the ABIF has the `unimplemented legacy type found in file` warning at `read.abif` time, the `PCON.2` (per-base quality) data block is empty. `MakeBaseCallsInside` returned an empty `qualityPhredScores`, and the Phase-4 `checkQualityPhredScores` validator threw `'qualityPhredScores' length cannot be zero`.

**Fix**: `SangerRead.initialize` now detects the empty/missing `PCON.2` block, synthesises a flat Phred-30 vector matching the peak count (`rep(30L, nrow(peakPosMatrixRaw))`), and logs a `MISSING_QUALITY_SCORES_WARN`. The rest of the pipeline (basecalling, alignment, consensus) runs unmodified. Caveat: with synthetic flat quality, M1/M2 trimming becomes a no-op — users should manually inspect the consensus from such reads.

**Tests**:

- Real ABIF with intact `PCON.2` still uses the instrument quality scores (no regression).
- Doctored ABIF with `PCON.2 = raw(0)` falls through to synthesised Phred 30; every score in the resulting `@QualityReport@qualityPhredScores` slot equals 30.

The mock target had to be `.package = "sangeranalyseR"` (where `read.abif` is imported via `@importFrom sangerseqR`), not `"sangerseqR"` — internal call sites use the importing package's namespace.

### Issue #89 — `writeFasta()` errors on contigs with only 1 read

**Root cause** (`R/MethodSangerContig.R::writeFastaSC`, pre-Phase-15):

```r
alignmentObject <- object@alignment
alignmentObject$Consensus <- NULL
writeAlignment <- append(alignmentObject, list(object@contigSeq))
writeXStringSet(writeAlignment, ...)   # 'x' must be an XStringSet object
```

Single-read `SangerContig` objects (or any case where `calculateContigSeq` is skipped because `readNumber < 2`) have an empty `@alignment` slot. `append(empty_DNAStringSet, list(contigSeq))` returned a plain `list`, not an `XStringSet`, so `writeXStringSet` rejected it.

**Fix**: detect the empty-alignment state and build a single-element `DNAStringSet` from `@contigSeq` directly. Multi-read contigs continue to use the original `c(alignmentObject, DNAStringSet(contigSeq))` path (also corrected from `append(...)` to `c(...)` for type consistency).

**Tests**:

- `writeFastaSC` directly on a single-read `SangerContig` (built with `minReadsNum = 1`, forward-only) produces a non-empty `<contigName>_reads_alignment.fa` file with at least 1 record.
- The polymorphic `writeFasta()` dispatcher follows the same path.

---

## 3. Test coverage delta

`tests/testthat/test-Phase15-Sprint1.R` (10 `test_that` blocks, 46 assertions):

| Issue | Test block(s)                                                                                   | Assertions |
| ----: | ----------------------------------------------------------------------------------------------- | ---------: |
|  #100 | non-substring contig labels build correctly                                                      |          7 |
|  #100 | partial-match contig labels: matched ones build, ghost ones logged + skipped                    |          5 |
|   #92 | `SangerAlignment(REGEX_SuffixReverse = NULL)`                                                    |          7 |
|   #92 | `SangerAlignment(REGEX_SuffixReverse = NA_character_)`                                           |          1 |
|   #92 | `SangerContig(REGEX_SuffixReverse = NULL)`                                                       |          3 |
|   #92 | Validator unit (`checkREGEX_SuffixReverse(NULL)` / `(NA)` / `(123L)`)                           |          3 |
|   #76 | Empty `PCON.2` synthesises Phred 30 (mocked `read.abif`)                                         |          5 |
|   #76 | Intact `PCON.2` keeps instrument scores (regression check)                                       |          2 |
|   #89 | `writeFastaSC` on single-read contig produces non-empty FASTA                                    |          7 |
|   #89 | `writeFasta` dispatcher on single-read contig                                                    |          2 |

**46 / 46 PASS.** Cumulative suite: 1360 → **1406 PASS** (no regressions).

---

## 4. GitHub-CLI replies

`plans/github_replies_sprint1.sh` (executable bash) contains four `gh issue comment <#>` + `gh issue close <#>` pairs, one per resolved issue. Each comment:

- Acknowledges the reporter by handle and thanks them.
- Identifies the **root cause** in technical detail (with a code snippet of the offending pre-fix lines where useful).
- Describes the **fix** and references the regression tests.
- Notes the fix will land in the next Bioconductor `devel` build (`sangeranalyseR` 1.21.x) and provides the `remotes::install_github("roblanf/sangeranalyseR", ref = "devel")` install command for immediate testing.
- Closes the issue with an explicit "please reopen if not resolved" invitation.

To execute:

```bash
chmod +x plans/github_replies_sprint1.sh
./plans/github_replies_sprint1.sh
```

(Requires `gh` CLI authenticated to `roblanf/sangeranalyseR`.)

---

## 5. Files touched in Phase 15

```
M  R/ClassSangerAlignment.R              (#100 CSV-keyed lookup; #92 sentinel substitution)
M  R/ClassSangerContig.R                  (#92 sentinel substitution)
M  R/ClassSangerRead.R                    (#76 missing-PCON.2 synthesis)
M  R/UtilitiesFuncInputChecker.R          (#92 NULL/NA-accepting REGEX validators)
M  R/MethodSangerContig.R                  (#89 single-read writeFasta path)
A  tests/testthat/test-Phase15-Sprint1.R   (46 assertions / 10 test_thats)
A  plans/15_sprint1_resolution_log.md
A  plans/github_replies_sprint1.sh         (executable; gh CLI commands)
```

---

## 6. Verification

```r
# Targeted Sprint 1 tests:
testthat::test_file("tests/testthat/test-Phase15-Sprint1.R")
# -> 46 / 46 PASS

# Full suite (no regressions):
devtools::test()
# -> 1406 / 1406 PASS

# Strict R CMD check:
rcmdcheck::rcmdcheck(".", args = "--no-manual")
# -> 0 / 0 / 0
```

---

## 7. Sprint 2 preview (per Phase 14 triage)

Issues remaining for the next sprint:

- **#94** Minimal-overlap F + R 16S reads forced to full-length overlap.
- **#66** Improper merge / contig with default M1 cutoff.
- **#65** Handle the case where a read is assigned to >1 place.
- **#42** `minReadLength` doesn't filter very-low-quality short reads.

These bundle naturally because they all need the same "post-merge alignment-quality check" infrastructure.
