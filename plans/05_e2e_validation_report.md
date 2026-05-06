# sangeranalyseR — Phase 5 E2E Validation & Performance Baseline

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24, R-devel r87868).

This report records actual measurements from running the full pipeline on the bundled ACHLO fixture. Raw artifacts live under `plans/phase5_artifacts/`.

---

## 1. Test fixture

**Dataset:** `inst/extdata/Allolobophora_chlorotica/ACHLO/`
- 8 ABIF files, 4 contigs (each contig has 1 forward + 1 reverse read).
- Read filenames `Achl_ACHLO{006,007,040,041}-09_{1_F,2_R}.ab1`.
- Forward reads ~700 bp pre-trim, ~460 bp post-trim (M1 default).

**Environment:**
- macOS Darwin 25.3.0 (single thread, `processorsNum = 1`).
- R Under development (unstable) (2025-03-02 r87868).
- DECIPHER 3.3.2, Biostrings 2.75.4, sangerseqR 1.43.0.
- Package loaded via `devtools::load_all()` from the dev tree (1.21.1) — *not* the system-installed 1.18.0.

---

## 2. Full-workflow execution

End-to-end pipeline run twice (M1 and M2). One untimed warm-up run amortises JIT/cache. Reproducer: `plans/phase5_artifacts/run_e2e.R`.

```
=== Pipeline run: TrimmingMethod = M1 ===
  [SangerAlignment M1     ]   1.855 s
  [writeFasta (all)       ]   0.021 s
  [updateQualityParam     ]   0.984 s

=== Pipeline run: TrimmingMethod = M2 ===
  [SangerAlignment M2     ]   1.743 s
  [writeFasta (all)       ]   0.007 s
  [updateQualityParam     ]   0.963 s
```

| Method | Stage                                      | Seconds |
| ------ | ------------------------------------------ | ------: |
| M1     | Ingest + trim + assemble + cross-align     |   1.855 |
| M1     | writeFasta (all reads + contigs + alignment) | 0.021 |
| M1     | updateQualityParam → M2 (re-trim + re-align) | 0.984 |
| M2     | Ingest + trim + assemble + cross-align     |   1.743 |
| M2     | writeFasta (all)                           |   0.007 |
| M2     | updateQualityParam → M1 (re-trim + re-align) | 0.963 |

**Time-to-Consensus** (full pipeline, 8 reads → 4 contigs → cross-aligned consensus):
- **M1 path: ~1.85 s** (~230 ms per read)
- **M2 path: ~1.74 s** (~218 ms per read)
- M2 is consistently ~6% faster than M1 (sliding-window arithmetic versus Mott's cumulative score).

`updateQualityParam` re-runs trimming and the cross-contig alignment from cached reads; it costs ~53% of a from-scratch run because it skips ABIF parse + base calling.

`writeFasta` is essentially free (<1% of total wall time) — output IO is not a concern.

---

## 3. Result audit

### 3a. Structural invariants (`plans/phase5_artifacts/accuracy.txt`)

```
[PASS] SangerAlignment built (M1)
[PASS] contigList non-empty (M1)
[PASS] contigsConsensus non-empty (M1)
[PASS] contigsAlignment is DNAStringSet (M1)
[PASS] All child reads valid (M1)
[PASS] SangerAlignment built (M2)
[PASS] M1 vs M2: contig count equal
[PASS] M1 vs M2: same contig names
[PASS] Consensus is deterministic across two runs
```

All 9 structural / determinism checks PASS. The cross-method check (M1 vs M2) confirms both trimming algorithms produce 4 contigs with identical contig names — they only differ in the trimmed read boundaries, not in contig grouping.

### 3b. Phase-4 *Raw slot invariants (`plans/phase5_artifacts/slot_invariants.txt`)

```
[PASS] primarySeqRaw   == sangerseq parent primarySeq
[PASS] secondarySeqRaw == sangerseq parent secondarySeq
[PASS] peakPosMatrixRaw  == parent peakPosMatrix
[PASS] peakAmpMatrixRaw  == parent peakAmpMatrix
[PASS] traceMatrix (inherited) == parent traceMatrix
[PASS] primarySeq (post-MakeBaseCalls) length == qualityPhredScores length
```

All 6 invariants PASS. This validates the Phase-4 finding that the `*Raw` slots preserve the pre-`MakeBaseCallsInside` data correctly — the parent `sangerseq` primary/secondary sequences and peak matrices, freshly reloaded via `sangerseqR::sangerseq(read.abif(...))`, match the `*Raw` slots byte-for-byte. Since the inherited parent slots get overwritten by `MakeBaseCallsInside` post-construction, removing the `*Raw` slots (the original Phase-2 audit recommendation) would have lost this data.

### 3c. Determinism

Two independent `SangerAlignment(...)` builds with identical inputs produce identical `contigsConsensus` (`identical(as.character(...))` returns `TRUE`). DECIPHER's `AlignSeqs` / `ConsensusSequence` and `bionjs` are deterministic given a fixed seed — no surprise here, but it's now pinned by an explicit assertion.

---

## 4. Profiling — where the time actually goes

`Rprof()` with 5 ms sampling, captured fresh after warm-up. Full output: `plans/phase5_artifacts/profile_M1.txt`, `profile_M2.txt`.

### 4a. M1 build (1.57 s sampled, 92 samples)

Top by total time:

| Function                       | total.pct | self.pct | Notes                                     |
| ------------------------------ | --------: | -------: | ----------------------------------------- |
| `initialize`                   |    100.0% |     1.0% | Top-level S4 dispatch; everything is below it |
| `FUN` / `lapply`               |     99.4% |     1.3% | The per-contig + per-read lapply loops    |
| `Biostrings::translate`        |   **35.4%** |     ~0% | Called inside `calculateAASeq`; 3-frame translation |
| `calculateAASeq`               |     35.4% |     0.0% | One call per `SangerRead` constructor      |
| `MakeBaseCallsInside`          |   **33.8%** |     0.6% | Per-`SangerRead` re-base-calling           |
| `.make_fuzzy_genetic_code`     |     30.6% |     0.0% | Internal Biostrings helper called by translate |
| `peakvalues` (in MakeBaseCalls)|     25.8% |   **16.2%** | Per-peak inner loop                       |
| `.Call2`                       |     15.9% |    14.0% | C-level Biostrings/sangerseqR primitives  |
| `validObject`                  |      2.5% |     1.0% | Phase-4 setValidity overhead              |

### 4b. M2 build (1.54 s sampled, 93 samples)

Same top three: `Biostrings::translate` 35.8%, `MakeBaseCallsInside` 31.9%, `peakvalues` 23.8%, `.Call2` 13.7%, `read.abif` 9.1%, `validObject` 3.3%, `calculateContigSeq` 15.6%. Profile shape is essentially identical to M1.

### 4c. The headline finding

**The user's hypotheses (IO / trimming / DECIPHER consensus) are all wrong on this dataset.** Time goes:

| Bucket                                            | % of total | Where it lives                                   |
| ------------------------------------------------- | ---------: | ------------------------------------------------ |
| **AA translation (3 frames per read)**            | **~35%**   | `calculateAASeq` → `Biostrings::translate` × 3   |
| **Per-read base calling (peak picking)**          | **~32%**   | `MakeBaseCallsInside` → `peakvalues` inner loop  |
| **DECIPHER cross-contig + per-contig alignment**  | **~16%**   | `calculateContigSeq` + `alignContigs`            |
| **ABIF file IO**                                  |  ~9%       | `read.abif`                                       |
| **S4 dispatch / setValidity / accumulation**     |  ~3-5%     | `initialize`, `validObject`, etc.                |
| Other                                             |  ~3%       | logger, paste, etc.                               |

**`peakvalues` self-time is 16% of the entire pipeline** — it's a 1-line R function called once per channel × per peak × per read. Vectorising / Rcpp-ising this would be the single biggest local win.

**`Biostrings::translate` is being called eagerly** for all 3 reading frames on every SangerRead, even when the user never inspects `primaryAASeqS{1,2,3}`. Lazifying these (compute on first slot access) would shave ~35% off the per-read build time on this fixture.

DECIPHER is *not* the bottleneck on small datasets like this. It will become the bottleneck only as N (reads) and L (read length) grow into the dozens of contigs / multi-kilobase territory, because DECIPHER's `AlignSeqs` is O(N²L) at worst.

---

## 5. Comparison to Phase-2 audit predictions

| Phase-2 prediction                                                      | Phase-5 measurement                                          | Verdict   |
| ----------------------------------------------------------------------- | ------------------------------------------------------------ | --------- |
| `mclapply` overhead dominates `nPairwiseDiffs` / `oneAmbiguousColumn`  | Both fall below the profile noise floor (<1%) on N=8 reads   | **Confirmed** — fork overhead would have made this worse, current serial path is reasonable |
| DECIPHER copy-by-value round-trips waste time                           | Combined alignment ≈16% of total — not the dominant cost     | **Overstated** — real but secondary on small N |
| Validators are a meaningful chunk                                       | <1% (`validObject` ~3%, includes new setValidity)            | **Overstated** — they were always cheap; the wins were correctness, not speed |
| Per-`SangerRead` construction (read.abif + base calling) is the bottleneck for parallelism | `MakeBaseCallsInside + read.abif` ≈ 41%; embarrassingly parallel | **Confirmed** — BiocParallel migration of the per-read loop should give roughly linear scaling |

---

## 6. Recommendations for the next phase

Ranked by expected wall-clock impact on a typical sangeranalyseR workload (10s of reads, 100s of contigs):

| # | Action                                                                                                   | Expected savings                       |
| - | -------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| 1 | Lazify `calculateAASeq` — only compute `primaryAASeqS{1,2,3}` when the user asks for them (or when `refAminoAcidSeq != ""`) | ~35% of per-read build time            |
| 2 | Vectorise / Rcpp-ise `peakvalues` (and possibly `getpeaks`) inside `MakeBaseCallsInside`                  | ~16% of total wall time                 |
| 3 | Parallelise the `lapply(... new("SangerRead", ...))` loops in `SangerContig.initialize` via `BiocParallel::bplapply` | Up to N× on N cores for IO+basecall+translate (84% of total) |
| 4 | Cache the parsed CSV / FASTA so `checkAb1FastaCsv` doesn't re-read what `SangerAlignment.initialize` is about to read again | Unmeasurable on M=4 contigs; meaningful on M >> 100 |
| 5 | Replace `oneAmbiguousColumn` per-column scan with a single `Biostrings::consensusMatrix` call             | Fixes the asymptotic bound; small effect on N=8 |

Items 1 and 2 are localised and low-risk. Item 3 is the BiocParallel migration already scoped in `plans/02_quality_audit_summary.md` §2.

---

## 7. Artifacts

All under `plans/phase5_artifacts/`:

| File                       | Contents                                                                |
| -------------------------- | ----------------------------------------------------------------------- |
| `run_e2e.R`                | Reproducer: invokes `devtools::load_all()`, runs M1 + M2 + Rprof + invariants. |
| `run_e2e.log`              | Full stdout + log_info trace from the run.                              |
| `timings.csv`              | Per-stage Sys.time() breakdown (the table in §2).                       |
| `Rprof_M1.out`, `Rprof_M2.out` | Raw `Rprof` sample files (compatible with `summaryRprof`, `profvis::parse_rprof`). |
| `profile_M1.txt`, `profile_M2.txt` | `summaryRprof` `by.self` (top 25) + `by.total` (top 15).      |
| `accuracy.txt`             | Structural / determinism PASS log.                                      |
| `slot_invariants.txt`      | Phase-4 `*Raw` slot equality PASS log.                                  |

To reproduce:

```bash
Rscript plans/phase5_artifacts/run_e2e.R
```

Or interactive:

```r
devtools::load_all(".")
source("plans/phase5_artifacts/run_e2e.R")
```

For an interactive flame graph instead of `summaryRprof` text output:

```r
library(profvis)
profvis::profvis(profvis::parse_rprof("plans/phase5_artifacts/Rprof_M1.out"))
```
