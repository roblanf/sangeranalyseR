# GitHub issue replies — Sprint 3 (consensus algorithms)

This file follows the **same Markdown format** as `plans/16_github_replies.md` so that `plans/close_issues.py` can parse it without modification. Run the closure script with `--md plans/17_github_replies.md`:

```bash
export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx
python3 plans/close_issues.py --md plans/17_github_replies.md --dry-run
python3 plans/close_issues.py --md plans/17_github_replies.md
```

All three Sprint 3 fixes shipped in the current Phase 17 commit on the GitHub `devel` branch.

---

## Sprint 3

### Issue #87 — Majority-rules consensus base calling?

**URL**: https://github.com/roblanf/sangeranalyseR/issues/87

```markdown
Hi, thanks for the clear write-up and the example alignment.

You're right that the pre-Phase-17 consensus path used DECIPHER's `ConsensusSequence(ambiguity = TRUE)`, which encodes per-column disagreements as IUPAC ambiguity codes. For users who prefer a plain majority call (and accept the loss of disagreement information), Phase 17 adds an explicit `consensusMethod` parameter to both `SangerContig()` and `SangerAlignment()`:

- `consensusMethod = "strict"` (default) — pre-Phase-17 behaviour, IUPAC ambiguity codes preserved.
- `consensusMethod = "majority"` — at each alignment column the most-frequent base wins; ties break alphabetically. No IUPAC ambiguity codes ever appear in the output.
- `consensusMethod = "quality_weighted"` — same as majority but votes are weighted by source-read Phred scores (see #48).

Example:

```r
sc <- SangerContig(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = "path/to/ab1s",
    contigName          = "Achl_ACHLO006-09",
    REGEX_SuffixForward = "_F\\.ab1$",
    REGEX_SuffixReverse = "_R\\.ab1$",
    consensusMethod     = "majority"
)
as.character(sc@contigSeq)        # plain ACGT, no IUPAC codes
```

Bonus: under `"majority"` and `"quality_weighted"`, Phase 17 also attaches a synthetic per-position quality vector to the consensus (see #33):

```r
attr(sc@contigSeq, "qualityScores")    # integer vector, length == length(consensus)
```

Available in the next Bioconductor `devel` build (`sangeranalyseR` 1.21.x). To test against your data right now, install from the GitHub `devel` branch:

```r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
```

Closing as fixed; please reopen if the new majority path doesn't match your expectations.
```

---

### Issue #48 — Consider Phred scores when consensus building

**URL**: https://github.com/roblanf/sangeranalyseR/issues/48

```markdown
Hi, finally addressing this one — closed-loop with #87 since both touch the same code path.

Phase 17 adds a `consensusMethod = "quality_weighted"` mode (also accessible via the shorthand `qualityAware = TRUE`) to `SangerContig()` and `SangerAlignment()`. The implementation:

1. Build a per-read × per-alignment-column Phred matrix from the source reads' `@QualityReport@qualityPhredScores` slots, with reverse reads' scores reversed to match the reverse-complemented alignment.
2. At each column, sum the Phred-weighted votes per base. The base with the highest total wins.
3. The consensus quality at each position is the mean Phred of agreeing reads at that column.

Example:

```r
sc <- SangerContig(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = "path/to/ab1s",
    contigName          = "your-contig",
    REGEX_SuffixForward = "_F\\.ab1$",
    REGEX_SuffixReverse = "_R\\.ab1$",
    qualityAware        = TRUE          # = consensusMethod = "quality_weighted"
)
attr(sc@contigSeq, "qualityScores")    # mean Phred of agreeing reads at each position
```

Caveats:

- Quality-weighted only applies for `inputSource = "ABIF"` (FASTA reads have no Phred scores). The package falls back to flat Phred 30 when scores are missing and logs a warning.
- The default `consensusMethod = "strict"` is unchanged — opt-in only.

Closing as fixed.
```

---

### Issue #33 — Add quality to consensus

**URL**: https://github.com/roblanf/sangeranalyseR/issues/33

```markdown
Thanks for the cross-link to #32 — answering both with the same change.

Phase 17 attaches per-position consensus quality scores to `@contigSeq` whenever the new `consensusMethod` is `"majority"` or `"quality_weighted"` (see #87 / #48). The scores are stored as an attribute on the `DNAString` so they survive serialisation and don't break old code:

```r
sc <- SangerContig(
    ...,
    consensusMethod = "majority"
)
attr(sc@contigSeq, "qualityScores")
#  integer vector, length == length(sc@contigSeq)
```

How the synthetic scores are computed:

- Under `"majority"`: per-column synthetic Phred = `40 * (winner_count / total_reads_at_column)`. So a 3/3 unanimous column reports Phred 40; a 2/3 majority reports Phred ~27.
- Under `"quality_weighted"`: per-column Phred = mean of the source reads' Phred scores at that column for the reads that agreed with the consensus base.

This isn't quite the OverlapPER per-base re-calibration paper @czliubio referenced in #32 (which builds a per-base quality from overlap evidence + sequencing error model), but it gives a useful per-position confidence signal for reports / downstream filtering.

Strict-mode users see no behaviour change: the `qualityScores` attribute is empty (`integer(0)`) under the default `consensusMethod = "strict"`.

Closing as fixed.
```

---

## How to apply

Use the Phase-16.5 closure script to post all 3 replies:

```bash
export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx
python3 plans/close_issues.py --md plans/17_github_replies.md
```

Or paste each reply manually via the GitHub web UI — copy the contents of each ` ```markdown ` block above into the issue's comment box, click **Comment**, then **Close issue**.
