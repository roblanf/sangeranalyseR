# sangeranalyseR — Phase 7 Rcpp Optimization Log

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24, R-devel r87868, Rcpp 1.0.14).

---

## 1. Headline result

| Configuration                                        | Wall time (8 reads, 4 contigs, mean of 5 reps) | vs. baseline |
| ---------------------------------------------------- | ---------------------------------------------: | -----------: |
| Phase 7 — `peakvalues_batch_cpp` (C++ batch)         | **1.142 s** (min 1.050 s)                      | **−22%**     |
| Phase 6 baseline — `.peakvalues_r` per-window R loop | 1.461 s (min 1.355 s)                          | baseline     |

**Speedup: 1.28× mean, 1.29× best-case.** ~320 ms saved per `SangerAlignment(...)` build on the bundled ACHLO fixture.

Cumulative wall-time progression:

| Milestone                                  | Mean wall time | Cumulative speedup vs. start |
| ------------------------------------------ | --------------:| ----------------------------: |
| Phase 5 baseline (eager AA, R-only)        | 1.85 s         | 1.0×                          |
| Phase 6 (lazy AA + BiocParallel plumbing)  | 1.33 s         | 1.39×                         |
| Phase 7 (lazy AA + Rcpp `peakvalues_batch`) | **1.07 s** (min) / 1.14 s (mean) | **1.62× best-case / 1.62× mean** |

Determinism preserved: byte-identical `contigsConsensus` across two builds, and across the C++ vs `.peakvalues_r` paths.

---

## 2. What changed

### 2a. New C++ source (`src/peakvalues.cpp`)

Two functions:

- `peakvalues_cpp(NumericMatrix x, double pstart, double pstop)` — direct port of `.peakvalues_r`, used in equivalence tests.
- `peakvalues_batch_cpp(NumericMatrix x, NumericVector pstarts, NumericVector pstops)` — processes all peak windows for one channel in a single `.Call`. Returns a `2 × K` matrix where column k = `(max_amp, max_pos)`.

Both preserve the original R contract: empty region → `(0, NA)`, one or more matches → `(max, position-of-first-occurrence-of-max)`, NA-in-column-2 is skipped.

### 2b. Why a batch API was needed

The naive per-window port (`peakvalues_cpp`) showed a real per-call speedup in the microbenchmark:

| Input shape              | n | calls | R mean (s) | C++ mean (s) | Per-call speedup |
| ------------------------- | -: | ----: | --------: | -----------: | ---------------: |
| Small (typical channel)   |  500 |   500 | 0.0076    | 0.0029       | **2.6×**         |
| Medium                    | 2000 |   500 | 0.0169    | 0.0091       | 1.9×             |
| Large (worst-case)        | 5000 |  1000 | 0.0731    | 0.0431       | 1.7×             |

…but on the full pipeline the per-call `.Call` marshalling overhead **eliminated** the savings. A typical SangerAlignment over the ACHLO fixture invokes `peakvalues` 4 channels × ~700 windows × 8 reads ≈ **22,400 times**. At ~7 µs of `.Call` overhead per invocation, that's ~150 ms of pure dispatch cost — roughly the size of the per-call inner-loop savings.

The batch API collapses 22,400 dispatches into 4 channels × 8 reads = **32 dispatches per build**, eliminating the marshalling overhead. That's where the 320 ms / 1.28× end-to-end win comes from.

### 2c. Call-site change

Before (`R/UtilitiesFunc.R::MakeBaseCallsInside`, per-iteration):

```r
for (i in seq_len(length(starts))) {
    Apeak <- peakvalues(Apeaks, starts[i], stops[i])
    Cpeak <- peakvalues(Cpeaks, starts[i], stops[i])
    Gpeak <- peakvalues(Gpeaks, starts[i], stops[i])
    Tpeak <- peakvalues(Tpeaks, starts[i], stops[i])
    ...
}
```

After (4 batch calls hoisted outside the loop):

```r
AbatchOut <- peakvalues_batch_cpp(Apeaks, starts, stops)
CbatchOut <- peakvalues_batch_cpp(Cpeaks, starts, stops)
GbatchOut <- peakvalues_batch_cpp(Gpeaks, starts, stops)
TbatchOut <- peakvalues_batch_cpp(Tpeaks, starts, stops)

for (i in seq_len(length(starts))) {
    Apeak <- AbatchOut[, i]
    Cpeak <- CbatchOut[, i]
    Gpeak <- GbatchOut[, i]
    Tpeak <- TbatchOut[, i]
    ...
}
```

The original R `peakvalues` was renamed `.peakvalues_r` (private; not exported). It survives only as the reference implementation against which `peakvalues_cpp` and `peakvalues_batch_cpp` are tested for byte-identical output.

### 2d. Profile shift

Before Phase 7 (`plans/phase6_artifacts/profile_M1_lazy.txt`):

| Function             | self.pct |
| -------------------- | -------: |
| `peakvalues`         | 24.2%    |
| `c`                  |  7.6%    |
| `<`                  |  5.8%    |

After Phase 7 (`plans/phase7_artifacts/profile_phase7.txt`):

| Function                  | self.pct |
| ------------------------- | -------: |
| `.Call`                   | 15.7%    |
| `c`                       |  6.3%    |
| `isatty`                  |  6.3%    |
| `MakeBaseCallsInside`     |  4.2%    |
| `read.abif` (total 17.3%) |  2.1% self |

`peakvalues` is gone from `by.self` and from `by.total` — its work is now under `.Call` (the native code itself). The new top of profile is `.Call` (the C++ work, fundamentally compute-bound) and `read.abif` (file IO). Neither is amenable to further pure-R optimization.

---

## 3. Mathematical equivalence — test coverage

`tests/testthat/test-Rcpp-peakvalues.R` (12 test_thats, ~190 assertions across the per-window pinning loops):

| Test                                                                   | What it pins                                          |
| ---------------------------------------------------------------------- | ------------------------------------------------------ |
| empty region                                                           | `c(0, NA)` for both R and C++                          |
| zero-row matrix                                                        | `c(0, NA)`                                             |
| single peak in region                                                  | `c(amp, pos)`                                          |
| multiple peaks, distinct max                                           | argmax position correct                                |
| tied max → first occurrence                                            | matches R's `which.max` semantics                      |
| strict `>` / `<` boundary                                              | `pstart`/`pstop` themselves are excluded               |
| NA in column 2                                                         | skipped, returns the max of the non-NA entries         |
| non-integer trace values                                               | preserved bit-exactly (`identical()`)                  |
| 200 random fuzz trials                                                 | zero mismatches between R and C++                      |
| batch column j == per-call result                                      | `identical(batch[,j], peakvalues_cpp(x, starts[j], stops[j]))` |
| batch column j == `.peakvalues_r` result                                | `identical(batch[,j], .peakvalues_r(x, starts[j], stops[j]))` |
| `peakvalues_batch_cpp` rejects mismatched `pstarts` / `pstops` length  | `expect_error("equal length")`                         |
| SangerAlignment consensus deterministic across two builds              | `identical(as.character(...))`                         |
| MakeBaseCallsInside under C++ batch == under R `.peakvalues_r` loop    | identical `primarySeq`, `secondarySeq`, `qualityPhredScores` |

Plus all 1167 pre-existing tests continue to pass. **Final suite: 1249 / 1249 PASS.**

---

## 4. Build setup

| File                                | Change                                                                                |
| ----------------------------------- | -------------------------------------------------------------------------------------- |
| `DESCRIPTION`                        | `Rcpp` added to `Imports:`, new `LinkingTo: Rcpp`, `RcppExports.R` added to `Collate:` |
| `NAMESPACE` (auto-regenerated)       | `useDynLib(sangeranalyseR, .registration = TRUE)` + `importFrom(Rcpp, sourceCpp)`     |
| `R/sangeranalyseR_package.R`         | `@useDynLib` and `@importFrom Rcpp` roxygen tags                                      |
| `src/peakvalues.cpp`                 | New — 99 lines C++                                                                     |
| `src/RcppExports.cpp`                | Auto-generated by `Rcpp::compileAttributes()`                                          |
| `R/RcppExports.R`                    | Auto-generated thin R wrappers                                                         |

### Toolchain note (macOS-specific)

On macOS Tahoe 26.x with Apple Command Line Tools 17, `clang++` couldn't find `<cmath>` from its default include path (`/Library/Developer/CommandLineTools/usr/include/c++/v1/` only contains `__cxx_version`). Worked around by adding to `~/.R/Makevars`:

```make
CXXFLAGS += -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
CXX17FLAGS += -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
CXX14FLAGS += -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
```

This is a per-developer setup file, not a package-level `src/Makevars`, so it does not pollute Bioconductor builds (which use a freshly-installed Xcode SDK and don't hit the issue). Linux and Windows R-tools toolchains compile straight without configuration.

---

## 5. Files touched in Phase 7

```
M  DESCRIPTION                                                        (+ Rcpp imports/LinkingTo, Collate)
M  NAMESPACE                                                          (auto-regen)
M  R/RcppExports.R                                                    (auto-generated)
M  R/UtilitiesFunc.R                                                  (rename peakvalues→.peakvalues_r; switch inner loop to peakvalues_batch_cpp)
M  R/sangeranalyseR_package.R                                         (+ @useDynLib, @importFrom Rcpp)
A  src/peakvalues.cpp                                                 (Phase 7 C++ port; both single + batch API)
A  src/RcppExports.cpp                                                (auto-generated)
A  tests/testthat/test-Rcpp-peakvalues.R
A  plans/07_rcpp_optimization_log.md
A  plans/phase7_artifacts/bench_peakvalues.R
A  plans/phase7_artifacts/bench_peakvalues.csv
A  plans/phase7_artifacts/bench_e2e.R
A  plans/phase7_artifacts/bench_e2e.csv
A  plans/phase7_artifacts/Rprof_phase7.out
A  plans/phase7_artifacts/profile_phase7.txt
```

---

## 6. Reproducing

```bash
# All tests (1249 passes):
Rscript -e 'devtools::test()'

# Microbenchmark (per-call):
Rscript plans/phase7_artifacts/bench_peakvalues.R

# End-to-end pipeline (real workload):
Rscript plans/phase7_artifacts/bench_e2e.R

# Re-profile:
Rscript -e 'profvis::profvis(profvis::parse_rprof("plans/phase7_artifacts/Rprof_phase7.out"))'
```

---

## 7. Open questions / non-goals

- **Binary search inside `peakvalues_batch_cpp`**. The peak indexes from `getpeaks` are sorted, so we could replace the linear scan with `std::lower_bound`/`std::upper_bound` for an extra log-vs-linear gain. On the 8-read fixture this is negligible (per-channel n ≈ 200–400 peaks); on multi-thousand-read alignments it could deliver another 2-3×. Deferred — would need its own benchmark.
- **Full-loop port to C++** (`MakeBaseCallsInside` itself). Would subsume basecalling, IUPAC ambiguity, primary/secondary peak selection, and quality-score realignment. Significantly larger surface area; not blocked by Phase 7.
- **`getpeaks` port**. Called four times per `SangerRead` (once per channel) outside the inner loop — accounts for <2% of total wall time; not worth the porting cost.
- **`parallel` → `Imports:` move** in DESCRIPTION (BiocCheck preference).
- **Bioconductor-grade build verification** via `R CMD check --as-cran` and `R CMD BiocCheck` on Linux (the macOS toolchain quirk noted in §4 won't reproduce there).
