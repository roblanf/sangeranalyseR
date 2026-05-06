# GitHub issue replies — Sprint 1 + Sprint 2

**Why this file replaces `plans/github_replies_sprint1.sh`**: the local `gh` CLI is v0.11.0 (2020-07-16) and predates both the `gh issue comment` subcommand and the `--body` flag. Rather than upgrade the CLI tool inside the dev environment, we deliver the replies as **copy-paste-ready Markdown blocks** that can be pasted directly into the GitHub web UI on each issue.

For each issue below:

1. The **GitHub URL** is a clickable link — click to open the issue.
2. The **reply text** is in a fenced code block. Click the "Copy" icon at the top-right of the code block, then paste into the issue's comment box on the web UI.
3. After posting the comment, click the **"Close issue"** button on the issue page.

All four Sprint 1 fixes shipped in commit `55c63e1` (Phase 15). All four Sprint 2 fixes shipped in the current Phase 16 commit. Both are on the GitHub `devel` branch ahead of the next Bioconductor `devel` build.

---

## Sprint 1

### Issue #100 — `SangerAlignment` CSV/ABIF returns CONTIG_NUMBER_ZERO_ERROR

**URL**: https://github.com/roblanf/sangeranalyseR/issues/100

```markdown
Hi @alnusmeinata, thanks for the careful report and the reproducible context.

This was a real bug in the CSV+ABIF aggregation path of `SangerAlignment.initialize`. The pre-fix code matched contig labels to filenames via a substring `grepl()`:

```r
contigNameSelectInputFiles <- parentDirFiles[grepl(contigName, parentDirFiles)]
```

That silently required every `contig` value in your CSV to appear as a substring of its associated read filenames. The CSV-grouping mechanism is precisely meant to map *arbitrary* contig labels (e.g. `rbcL2`, `ITS5`) to *arbitrary* filenames via the explicit `reads` column, so this requirement was incorrect — and exactly why each individual `SangerContig()` call worked while `SangerAlignment()` returned `CONTIG_NUMBER_ZERO_ERROR`.

**Fix (Phase 15)**: the aggregator now drives lookups directly off the CSV's `reads` column instead of fuzzy-matching the contig label against filenames. For each unique contig label, it intersects `csvFile$reads` (filtered to that label) against the actual files under `ABIF_Directory`, then derives the directory prefix from the matched files. Labels with no matching reads now log a friendly warning and skip rather than aborting the whole alignment.

Two new regression tests cover this:

- `Issue #100: SangerAlignment CSV+ABIF works when contig labels are non-substring of filenames` — uses labels `rbcL2` / `ITS5` that don't appear in any `Achl_ACHLO*` filename, asserts both contigs build correctly.
- `Issue #100: SangerAlignment CSV+ABIF logs warning when a contig label has no matching reads` — partial-match resilience.

Will be available in the next Bioconductor `devel` build (`sangeranalyseR` 1.21.x). To test against your data right now, install from the GitHub `devel` branch:

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
```

Closing as fixed; please reopen if the new code path doesn't resolve your specific case.
```

---

### Issue #92 — Forward-reads-only datasets rejected

**URL**: https://github.com/roblanf/sangeranalyseR/issues/92

```markdown
Hi @czliubio, thanks for the report — sorry for the long wait.

This affected anyone with single-direction Sanger data (forward-only barcoding, 16S half-coverage runs, etc.). The validators `checkREGEX_SuffixForward` and `checkREGEX_SuffixReverse` rejected `NULL` outright with `PARAMETER_VALUE_ERROR`, which made it impossible to express "no reverse reads" cleanly.

**Fix (Phase 15)**: both validators now accept `NULL` and `NA_character_` to mean "no reads in that direction". When either is supplied, `SangerAlignment.initialize` and `SangerContig.initialize` substitute an internal sentinel that never matches any filename, and the existing "no reads detected" warning path takes over — no crashes.

Forward-only example, which now works as expected:

```r
sa <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = "path/to/ab1s",
    REGEX_SuffixForward = "_F\\.ab1$",
    REGEX_SuffixReverse = NULL,    # explicit forward-only
    minReadsNum         = 1        # each forward read becomes its own contig
)
```

Note that `minReadsNum = 1` is needed so single-read contigs aren't filtered out by the post-construction read-count check. (The default `minReadsNum = 2` is appropriate for paired runs.)

Three regression tests cover this (`SangerAlignment` + `SangerContig` with `NULL`/`NA` reverse) plus a unit test on `checkREGEX_SuffixReverse` directly. Available in the next Bioconductor `devel` build.

Closing as fixed.
```

---

### Issue #76 — `qualityPhredScores Error`

**URL**: https://github.com/roblanf/sangeranalyseR/issues/76

```markdown
Hi @CodyG12, thanks for the report. The `unimplemented legacy type found in file` warning was the giveaway — `sangerseqR::read.abif` was succeeding on those files but returning an empty `PCON.2` quality block, which then failed our hard `'qualityPhredScores' length cannot be zero` validator. We've seen this pattern from older Beckman / 3500 firmware ABIFs.

**Fix (Phase 15)**: `SangerRead.initialize` now detects the missing/empty `PCON.2` block, synthesises a flat Phred-30 quality vector (one entry per detected peak), and logs a `MISSING_QUALITY_SCORES_WARN` so it's clear the trimmed sequence's quality assessment is not based on the original instrument's per-base scores. The rest of the pipeline (basecalling, alignment, consensus) runs normally.

Caveats: with synthetic flat quality, M1/M2 trimming becomes a no-op (no quality differential to trim against). The `SangerRead@QualityReport` slot is still populated and consistent so Shiny / report rendering continue to work, but you should manually inspect the consensus from these reads.

Two regression tests:

- Real ABIF with intact `PCON.2` continues to use the instrument quality scores (no behaviour change for the happy path).
- Doctored ABIF with empty `PCON.2` falls through to the synthesised Phred-30 path; assertion: every score in the resulting `@QualityReport@qualityPhredScores` slot equals 30.

Available in the next Bioconductor `devel` build. Install from GitHub for immediate testing:

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
```

Closing as fixed.
```

---

### Issue #89 — `writeFasta()` errors on contigs with only 1 read

**URL**: https://github.com/roblanf/sangeranalyseR/issues/89

```markdown
Thanks for the report and the precise error trace.

Confirmed bug: `SangerContig` objects built from a single read (or any case where `calculateContigSeq` is skipped because `readNumber < 2`) have an empty `@alignment` slot. `writeFastaSC` then ran:

```r
writeAlignment <- append(alignmentObject, list(object@contigSeq))
```

…which on an empty `DNAStringSet` returned a plain `list`, not an XStringSet — so `writeXStringSet` rejected it with `'x' must be an XStringSet object`.

**Fix (Phase 15)**: `writeFastaSC` now detects the empty-alignment state and writes the contig sequence on its own as a single-record `DNAStringSet` (named `<contigName>_contig`). Single-read contigs now produce a valid `<contigName>_reads_alignment.fa` file.

Two regression tests:

- `writeFastaSC` directly on a single-read `SangerContig` produces a non-empty FASTA.
- The polymorphic `writeFasta()` dispatcher follows the same path.

Available in the next Bioconductor `devel` build.

Closing as fixed.
```

---

## Sprint 2

### Issue #94 — Minimal-overlap F + R 16S reads forced to full-length overlap

**URL**: https://github.com/roblanf/sangeranalyseR/issues/94

```markdown
Hi @aliciaastr, thanks for the careful description — this captured a real limitation.

The pre-Phase-16 pipeline called `DECIPHER::AlignSeqs` with a fixed argument list, so users with low-overlap F + R reads (~50–100 bp shared region on ~800 bp reads, common in 16S barcoding) had no way to tune the alignment, and the resulting consensus was a long IUPAC-ambiguity soup over the non-overlapping flanks.

**Fix (Phase 16)**: two changes.

1. `SangerContig()` and `SangerAlignment()` now accept `alignSeqsParams = list(...)`, a named list forwarded verbatim to `DECIPHER::AlignSeqs` (or `AlignTranslation` when `refAminoAcidSeq` is supplied). You can now tune `iterations`, `refinements`, `gapOpening`, etc.:

    ```r
    sa <- SangerAlignment(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = "path/to/16S",
        REGEX_SuffixForward = "_F\\.ab1$",
        REGEX_SuffixReverse = "_R\\.ab1$",
        alignSeqsParams     = list(iterations = 1L, refinements = 1L)
    )
    ```

2. New `minOverlapFraction` and `minOverlapBases` arguments. When > 0, after read alignment the smallest pairwise non-gap overlap is computed; if it falls below the threshold a `LOW_OVERLAP_WARN` is logged so you don't silently get a degenerate consensus. Default behaviour (thresholds 0) is unchanged.

    ```r
    sa <- SangerAlignment(
        ...,
        minOverlapBases    = 50L,    # warn if any pairwise overlap < 50 bp
        minOverlapFraction = 0.05    # or < 5% of the shorter read
    )
    ```

Regression tests cover both the warning-fires path and the no-regression-on-real-data path. Available in the next Bioconductor `devel` build.

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
```

Closing as fixed.
```

---

### Issue #66 — Improper merge / contig with default M1 cutoff

**URL**: https://github.com/roblanf/sangeranalyseR/issues/66

```markdown
Hi @joshpost, thanks for the persistence; you correctly identified that the merge was happening even when forward and reverse reads shared no signal.

This is a sibling of #94. The pre-Phase-16 alignment phase trusted DECIPHER to return *something*, and the subsequent `ConsensusSequence` call dutifully filled in IUPAC ambiguity codes for the non-overlapping flanks — producing exactly the long degenerate sequences you reported.

**Fix (Phase 16)**: the same `minOverlapFraction` / `minOverlapBases` instrumentation added for #94 now flags this case. With `minOverlapBases = 50L` (or any threshold below your expected real overlap), you'd have seen:

```
WARN: LOW_OVERLAP_WARN: smallest pairwise overlap is 12 bp (between
1_Read_AB4-13_F.ab1 <-> 2_Read_AB4-13_R.ab1); required threshold is
50 bp. Consensus may contain spurious IUPAC ambiguity codes; review
carefully or tighten trimming parameters before merging.
```

…instead of silently getting a degenerate consensus. The default thresholds (0 / 0L) preserve pre-Phase-16 behaviour for backwards compatibility — opt in by passing the new parameters.

Recommended workflow for low-overlap pairs:

```r
sa <- SangerAlignment(
    ...,
    minOverlapBases    = 50L,
    minOverlapFraction = 0.05
)
```

Closing as fixed; please reopen if the warning doesn't fire on a case you expect.
```

---

### Issue #65 — Handle the case where a read is assigned to >1 place

**URL**: https://github.com/roblanf/sangeranalyseR/issues/65

```markdown
Thanks @roblanf for filing this one — closed loop with the rest of the Sprint 2 fixes.

**Fix (Phase 16)**: `checkAb1FastaCsv` (the CSV-validation helper) now detects when a single `reads` filename is mapped to >1 distinct `contig` value and logs a `READ_ASSIGNED_MULTIPLE_CONTIGS_WARN`:

```
WARN: Read(s) assigned to >1 distinct contig in CSV: 'Achl_ACHLO006-09_1_F.ab1'.
Each read must belong to exactly one contig.
(READ_ASSIGNED_MULTIPLE_CONTIGS_WARN)
```

The validator does **not** error in this case — for users who deliberately use the same read in multiple contigs (rare but legitimate), the warning is enough to surface the issue without breaking their workflow. If we get reports of users who expected a hard failure, we can promote it to an error in a future release.

Two regression tests:

- CSV with one read assigned to two contigs triggers the warning.
- CSV with no duplicates does not trigger the warning (regression check).

Closing as fixed.
```

---

### Issue #42 — `minReadLength` doesn't filter very-low-quality length-1 reads

**URL**: https://github.com/roblanf/sangeranalyseR/issues/42

```markdown
Hi @roblanf, finally addressing this one.

The root cause was a missing defensive filter in `calculateContigSeq` itself. The Phase-3 `minReadLength` filter at the SangerContig level handles the typical case, but on aggressively-trimmed degenerate inputs (where M1 produces `trimmedFinishPos = 0` for some windows) a length-1 entry could survive into the alignment input and silently break `DECIPHER::AlignSeqs`.

**Fix (Phase 16)**: `calculateContigSeq` now applies a defensive width filter before alignment — any read with trimmed primary sequence width < 2 bp is dropped with a `MIN_READ_LENGTH_DEFENSIVE_DROP` warning. If the filter brings us below 2 reads total, the function returns a degenerate-but-well-formed result (containing the single surviving read as its own consensus, or an empty `DNAString` if none survive) so the caller can fold this into a controlled `READ_NUMBER_ERROR` instead of crashing the whole alignment.

Regression test: a synthetic `SangerContig` with a length-1 forward read mixed with normal reverse reads now builds successfully — the bad read is dropped, the consensus is well-formed.

Closing as fixed.
```

---

## How to apply the closes

1. Click the URL for the issue you want to close.
2. Scroll to the comment box at the bottom.
3. Click the "Copy" icon on the relevant code block above and paste into the comment box. The triple-backtick fences delimit the comment body — they're not part of it.
4. Click **Comment** to post.
5. Then click **Close issue**.

Repeat for each of the 8 issues. After all are closed the open backlog drops from **43 → 35**.
