#!/usr/bin/env bash
# =============================================================================
# Phase 15 — GitHub-CLI replies + closes for Sprint 1 issues (#100, #92, #76, #89).
#
# Runs `gh issue comment` then `gh issue close` for each issue. Re-runnable;
# `gh` is idempotent on close (won't error if already closed).
#
# Pre-requisites:
#   - gh CLI authenticated to roblanf/sangeranalyseR.
#   - Phase 15 commit pushed to origin/devel (so the linked diff is reachable).
#
# Usage:
#   chmod +x plans/github_replies_sprint1.sh
#   ./plans/github_replies_sprint1.sh
# =============================================================================

set -euo pipefail
REPO=roblanf/sangeranalyseR

# -----------------------------------------------------------------------------
# Issue #100 — CSV/ABIF returns CONTIG_NUMBER_ZERO_ERROR
# -----------------------------------------------------------------------------
gh issue comment 100 -R "$REPO" --body "$(cat <<'EOF'
Hi @alnusmeinata, thanks for the careful report and the reproducible context.

This was a real bug in the CSV+ABIF aggregation path of \`SangerAlignment.initialize\`. The pre-fix code matched contig labels to filenames via a substring \`grepl()\`:

\`\`\`r
contigNameSelectInputFiles <- parentDirFiles[grepl(contigName, parentDirFiles)]
\`\`\`

That silently required every \`contig\` value in your CSV to appear as a substring of its associated read filenames. The CSV-grouping mechanism is precisely meant to map *arbitrary* contig labels (e.g. \`rbcL2\`, \`ITS5\`) to *arbitrary* filenames via the explicit \`reads\` column, so this requirement was incorrect — and exactly why each individual \`SangerContig()\` call worked while \`SangerAlignment()\` returned \`CONTIG_NUMBER_ZERO_ERROR\`.

**Fix (Phase 15)**: the aggregator now drives lookups directly off the CSV's \`reads\` column instead of fuzzy-matching the contig label against filenames. For each unique contig label, it intersects \`csvFile$reads\` (filtered to that label) against the actual files under \`ABIF_Directory\`, then derives the directory prefix from the matched files. Labels with no matching reads now log a friendly warning and skip rather than aborting the whole alignment.

Two new regression tests cover this:

- \`Issue #100: SangerAlignment CSV+ABIF works when contig labels are non-substring of filenames\` — uses labels \`rbcL2\` / \`ITS5\` that don't appear in any \`Achl_ACHLO*\` filename, asserts both contigs build correctly.
- \`Issue #100: SangerAlignment CSV+ABIF logs warning when a contig label has no matching reads\` — partial-match resilience.

Will be available in the next Bioconductor \`devel\` build (\`sangeranalyseR\` 1.21.x). If you'd like to test against your data right now you can install from the GitHub \`devel\` branch:

\`\`\`r
remotes::install_github("roblanf/sangeranalyseR", ref = "devel")
\`\`\`

Closing as fixed; please reopen if the new code path doesn't resolve your specific case.
EOF
)"

gh issue close 100 -R "$REPO"

# -----------------------------------------------------------------------------
# Issue #92 — Forward-reads-only datasets rejected
# -----------------------------------------------------------------------------
gh issue comment 92 -R "$REPO" --body "$(cat <<'EOF'
Hi @czliubio, thanks for the report — sorry for the long wait.

This affected anyone with single-direction Sanger data (forward-only barcoding, 16S half-coverage runs, etc.). The validators \`checkREGEX_SuffixForward\` and \`checkREGEX_SuffixReverse\` rejected \`NULL\` outright with \`PARAMETER_VALUE_ERROR\`, which made it impossible to express "no reverse reads" cleanly.

**Fix (Phase 15)**: both validators now accept \`NULL\` and \`NA_character_\` to mean "no reads in that direction". When either is supplied, \`SangerAlignment.initialize\` and \`SangerContig.initialize\` substitute an internal sentinel that never matches any filename, and the existing "no reads detected" warning path takes over — no crashes.

Forward-only example, which now works as expected:

\`\`\`r
sa <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = "path/to/ab1s",
    REGEX_SuffixForward = "_F\\\\.ab1$",
    REGEX_SuffixReverse = NULL,    # explicit forward-only
    minReadsNum         = 1        # each forward read becomes its own contig
)
\`\`\`

Note that \`minReadsNum = 1\` is needed so that single-read contigs aren't filtered out by the post-construction read-count check. (The default \`minReadsNum = 2\` is appropriate for paired runs.)

Three regression tests cover this (\`SangerAlignment\` + \`SangerContig\` with \`NULL\`/\`NA\` reverse) plus a unit test on \`checkREGEX_SuffixReverse\` directly. Available in the next Bioconductor \`devel\` build.

Closing as fixed.
EOF
)"

gh issue close 92 -R "$REPO"

# -----------------------------------------------------------------------------
# Issue #76 — qualityPhredScores Error
# -----------------------------------------------------------------------------
gh issue comment 76 -R "$REPO" --body "$(cat <<'EOF'
Hi @CodyG12, thanks for the report. The \`unimplemented legacy type found in file\` warning was the giveaway — \`sangerseqR::read.abif\` was succeeding on those files but returning an empty \`PCON.2\` quality block, which then failed our hard \`'qualityPhredScores' length cannot be zero\` validator. We've seen this pattern from older Beckman / 3500 firmware ABIFs.

**Fix (Phase 15)**: \`SangerRead.initialize\` now detects the missing/empty \`PCON.2\` block, synthesises a flat Phred-30 quality vector (one entry per detected peak), and logs a \`MISSING_QUALITY_SCORES_WARN\` so it's clear the trimmed sequence's quality assessment is not based on the original instrument's per-base scores. The rest of the pipeline (basecalling, alignment, consensus) runs normally.

Caveats: with synthetic flat quality, M1/M2 trimming becomes a no-op (no quality differential to trim against). The \`SangerRead@QualityReport\` slot is still populated and consistent so Shiny / report rendering continue to work, but you should manually inspect the consensus from these reads.

Two regression tests:

- Real ABIF with intact \`PCON.2\` continues to use the instrument quality scores (no behaviour change for the happy path).
- Doctored ABIF with empty \`PCON.2\` falls through to the synthesised Phred-30 path; assertion: every score in the resulting \`@QualityReport@qualityPhredScores\` slot equals 30.

Available in the next Bioconductor \`devel\` build.

Closing as fixed.
EOF
)"

gh issue close 76 -R "$REPO"

# -----------------------------------------------------------------------------
# Issue #89 — writeFasta errors on single-read contigs
# -----------------------------------------------------------------------------
gh issue comment 89 -R "$REPO" --body "$(cat <<'EOF'
Thanks for the report and the precise error trace.

Confirmed bug: \`SangerContig\` objects built from a single read (or any case where \`calculateContigSeq\` is skipped because \`readNumber < 2\`) have an empty \`@alignment\` slot. \`writeFastaSC\` then ran:

\`\`\`r
writeAlignment <- append(alignmentObject, list(object@contigSeq))
\`\`\`

…which on an empty \`DNAStringSet\` returned a plain \`list\`, not an XStringSet — so \`writeXStringSet\` rejected it with \`'x' must be an XStringSet object\`.

**Fix (Phase 15)**: \`writeFastaSC\` now detects the empty-alignment state and writes the contig sequence on its own as a single-record \`DNAStringSet\` (named \`<contigName>_contig\`). Single-read contigs now produce a valid \`<contigName>_reads_alignment.fa\` file.

Two regression tests:

- \`writeFastaSC\` directly on a single-read \`SangerContig\` produces a non-empty FASTA.
- The polymorphic \`writeFasta()\` dispatcher follows the same path.

Available in the next Bioconductor \`devel\` build.

Closing as fixed.
EOF
)"

gh issue close 89 -R "$REPO"

echo "Sprint 1 replies posted and issues closed."
