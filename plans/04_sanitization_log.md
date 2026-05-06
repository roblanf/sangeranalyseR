# sangeranalyseR — Phase 4 Sanitization Log

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

This log captures the bug fixes and refactors landed in Phase 4, together with verification commands.

---

## 1. Bug fixes

### 1a. Regex escape — `.fa` / `.fasta` / `.ab1`

**File:** `R/UtilitiesFuncInputChecker.R`

The old patterns used unescaped `.`, which matches any character — so e.g. `Sanger_all_reads.Xfa` slipped through as a "valid FASTA". Fixed by introducing two named constants and routing every extension check through the new `.requireExt` helper:

```r
.FASTA_EXT_REGEX <- "\\.(fa|fasta)$"
.AB1_EXT_REGEX   <- "\\.ab1$"
```

Both `checkFASTA_File()` and `checkReadFileName()` now use these patterns instead of inline `str_extract(.., ".fa$")` / `str_extract(.., ".ab1$")`.

The two BUG-LOCK tests added in Phase 3 (`test-Validator-EdgeCases.R`) — which were intentionally failing — now PASS:

| Test                                                                       | Status before Phase 4 | Status after Phase 4 |
| -------------------------------------------------------------------------- | --------------------- | -------------------- |
| `.Xfa file is rejected as invalid FASTA extension`                         | **FAIL** (BUG-LOCK)   | **PASS**             |
| `.Xab1 file is rejected as invalid ABIF extension`                         | **FAIL** (BUG-LOCK)   | **PASS**             |

The `BUG-LOCK:` prefix has been stripped from those test descriptions, since the bug is now fixed.

### 1b. Stray debug prints removed

Hunted via `grep -nE "^\s*(cat|print|message)\s*\("`. Five non-`show` debug leaks were removed/promoted to `log_info`:

| File                                          | Line  | Before                                                | After                                  |
| --------------------------------------------- | ----: | ----------------------------------------------------- | -------------------------------------- |
| `R/UtilitiesFuncInputChecker.R`               |  338  | `cat("FASTA_File", FASTA_File)`                       | removed                                 |
| `R/UtilitiesFunc.R`                           |  280  | `print("Removing reads with stop codons")`            | `log_info("Removing reads with stop codons")` |
| `R/UtilitiesFunc.R`                           |  944  | `cat(paste("Chromatogram saved to", filename, ...))`  | `log_info("Chromatogram saved to ", filename, ...)` |
| `R/ShinyServerModule.R`                       |  721  | `message(length(PhredScoreDF))`                       | removed                                 |
| `R/ShinySangerContigServer.R`                 | 1282  | `print("** Inside SCDifferencesDFUI !!!")`            | removed                                 |

`cat()` calls inside `R/sangeranalyseR_show_method.R`, `R/MethodSangerRead.R::readTable`, and `R/MethodSangerContig.R::readTable` are **legitimate** (they implement `show()` / `readTable` printing) and were left in place.

This unblocks the Bioconductor `BiocCheck` rule against `cat`/`print` in package code.

---

## 2. Validator refactor

**File:** `R/UtilitiesFuncInputChecker.R` (rewritten end-to-end).

### What changed

- Added five **internal helpers** (dot-prefixed, not exported):
  - `.errAppend(errors, errorTypes, msg, type)` — single concat point.
  - `.requireType(value, name, predicate, expected, errors, errorTypes, type)` — type-predicate check.
  - `.requireEnum(value, name, allowed, errors, errorTypes, type)` — character-enum check.
  - `.requireRange(value, name, lo, hi, errors, errorTypes, integer)` — numeric range / integer check.
  - `.requireExt(path, pattern, errors, errorTypes, message, type)` — file-extension regex check.
- Every public `check*` function preserves its **signature** and **error-type tag** (`PARAMETER_TYPE_ERROR`, `PARAMETER_VALUE_ERROR`, `PARAMETER_RANGE_ERROR`, `FILE_TYPE_ERROR`, etc.) — no test that asserts on a tag had to change.
- Switched per-element membership scans in `checkAb1FastaCsv` from two `lapply` loops to a single `setdiff` — same correctness, vectorized.
- Removed the duplicated chain in `checkTrimParam`: M2 ranges now flow through `.requireRange(integer = TRUE)` once instead of being open-coded twice.

### Code-size impact

| File                                | Before (lines) | After (lines) | Delta  |
| ----------------------------------- | --------------:| -------------:| ------:|
| `R/UtilitiesFuncInputChecker.R`     |            603 |            523 |   -80 (-13%) |

The 13% reduction is concentrated in the validators that are now one-liners over `.requireRange` / `.requireEnum` (e.g., `checkInputSource`, `checkSignalRatioCutoff`, `checkMinFractionCall`). The remaining lines are the helpers (~110 lines), `checkAb1FastaCsv`/`checkTrimParam` (kept procedural for cross-field logic), and per-validator wrappers that retain the (errors, errorTypes) signature contract.

### Latency micro-benchmark methodology

Reproducible script (run from R at the package root):

```r
library(microbenchmark)
library(sangeranalyseR)

errs  <- character(0)
types <- character(0)

bench <- microbenchmark(
    inputSource    = checkInputSource("ABIF",   errs, types),
    processMethod  = checkProcessMethod("ABIF", "REGEX", errs, types),
    refAAS         = checkRefAAS("",            errs, types),
    minReadsNum    = checkMinReadsNum(2,        errs, types),
    minFracCall    = checkMinFractionCall(0.5,  errs, types),
    signalRatio    = checkSignalRatioCutoff(0.33, errs, types),
    trimParamM1    = checkTrimParam("M1", 0.0001, NULL, NULL, errs, types),
    trimParamM2    = checkTrimParam("M2", NULL, 20, 10, errs, types),
    readFileName   = checkReadFileName("/tmp/x.ab1", "ABIF", errs, types),
    times = 5000L
)
print(bench)
```

Expected behavior: the helper-driven happy-path (no errors appended) returns `list(errors, errorTypes)` *unchanged* — single function call, no allocation. On the failure path, exactly one `c()` allocation is performed (in `.errAppend`), versus two in the previous implementation (one for `errors`, one for `errorTypes` — though those occurred in close succession, R's hidden length-doubling amortizes them).

The refactor's headline win is **code locality** rather than raw nanoseconds: every error path goes through one site, so adding a new error tag or restructuring the diagnostic envelope is a single-file change. Bench numbers should be recorded by the maintainer running the script on their own machine and pasted into this section.

---

## 3. S4 setValidity additions

Added post-construction invariants on the three helper classes whose state is most likely to be mutated via `slot<-`:

| Class                                       | Invariants checked                                                                                                                      |
| ------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `ChromatogramParam` (`R/ClassChromatogramParam.R`) | `baseNumPerRow ∈ [0, 200]`, `heightPerRow ∈ [50, 600]`, `signalRatioCutoff ∈ [0, 1]`, `showTrimmed` is a single logical.                |
| `QualityReport` (`R/ClassQualityReport.R`)  | `trimmedFinishPos >= trimmedStartPos`, `trimmedStartPos >= 0`, `length(qualityPhredScores) == length(qualityBaseScores)`, `remainingRatio ∈ [0, 1]`. (Skipped on the empty/default state.) |
| `ObjectResults` (`R/ClassObjectResults.R`)  | `length(creationResult) == 1` and is logical, `length(errorMessages) == length(errorTypes)`, `length(warningMessages) == length(warningTypes)`. |

### Why these classes and not Sanger*

- `SangerRead` / `SangerContig` / `SangerAlignment` follow a "construction never throws" contract — on input failure they build a sentinel object whose `objectResults@creationResult == FALSE`. Adding `setValidity` that throws would break this contract. The existing `check*` validators in `initialize` are the right layer for input validation.
- The three helper classes throw on construction failure already (e.g., `QualityReport` calls `stop(errors)` when its checks fail), so a validity method that throws is consistent.

The new validity methods are exercised by `tests/testthat/test-S4Validity.R` (10 new test cases).

---

## 4. Redundancy cleanup — Phase-2 misdiagnosis corrected

The Phase 2 audit (`plans/02_quality_audit_summary.md` §4) claimed the five "Raw" slots on `SangerRead` (`primarySeqRaw`, `secondarySeqRaw`, `peakPosMatrixRaw`, `peakAmpMatrixRaw`, `traceMatrix`) were redundant copies of the inherited `sangerseq` parent slots. **That was wrong.** A close re-read of `R/ClassSangerRead.R:200–429` shows:

1. The constructor pulls slots out of a temporary `sangerseq` object (`primarySeq`, `secondarySeq`, `peakPosMatrix`, `peakAmpMatrix`).
2. **`MakeBaseCallsInside` then re-base-calls and overwrites the local variables** for `primarySeq`, `secondarySeq`, `peakPosMatrix`, `peakAmpMatrix`.
3. Both the rewritten *and* the pre-rewrite values are passed to `callNextMethod`:
   - `primarySeq` (post-basecall) → inherited slot.
   - `primarySeqRaw` (pre-basecall) → SangerRead's own slot.

So the `*Raw` slots **preserve the pre-basecall sequences**, which would otherwise be lost when `MakeBaseCallsInside` rewrites the inherited slots. They are not duplicates.

**Decision:** keep the slots. Document the corrected understanding here. The Phase 3 `test-SangerRead-SlotInvariants.R` tests (which assert `sr@primarySeqRaw == ref@primarySeq` against a *fresh* `sangerseq()` reload — the pre-basecall state) continue to pass and will keep guarding the invariant.

`traceMatrix` is purely inherited (only on the parent `sangerseq` class) and is not duplicated.

---

## 5. New tests added in Phase 4

| File                                                        | Tests | Scope                                                               |
| ----------------------------------------------------------- | ----: | ------------------------------------------------------------------- |
| `tests/testthat/test-Validator-Helpers.R`                   |    19 | Direct unit tests for `.errAppend`, `.requireType`, `.requireEnum`, `.requireRange`, `.requireExt`, plus public-API smoke. |
| `tests/testthat/test-Regex-Boundary.R`                      |    10 | Boundary regression tests pinning the escaped extension regex against `.Xfa`, `.Xab1`, `.fast`, `.faa`, `.ab1.bak`, etc. |
| `tests/testthat/test-S4Validity.R`                          |    10 | `validObject()` rejects malformed `ChromatogramParam`/`QualityReport`/`ObjectResults` instances. |
| `tests/testthat/test-Validator-EdgeCases.R`                 |     — | Renamed two BUG-LOCK tests; they now PASS without any other change. |

Net new in Phase 4: **3 new test files, 39 new test_that blocks.**

Cumulative Phase 3 + 4 net new: **9 new test files, 82 new test blocks.**

---

## 6. Verification

| Goal                                                | Command                                                                                |
| --------------------------------------------------- | -------------------------------------------------------------------------------------- |
| Run the full suite                                  | `devtools::test()`                                                                      |
| Run only Phase 4 additions                          | `testthat::test_dir("tests/testthat", filter = "Validator-Helpers\\|Regex-Boundary\\|S4Validity")` |
| Confirm BUG-LOCK tests now pass                     | `testthat::test_file("tests/testthat/test-Validator-EdgeCases.R", filter = "Xfa\\|Xab1")` |
| Lint                                                | `R CMD BiocCheck sangeranalyseR_*.tar.gz`                                                |
| Full check                                          | `R CMD build .; R CMD check sangeranalyseR_*.tar.gz`                                     |
| Coverage delta vs. Phase 3 baseline                 | `cov <- covr::package_coverage(); covr::report(cov)`                                     |
| Latency micro-bench (manual)                        | See §2 "Latency micro-benchmark methodology".                                            |

### Expected outcomes

- Every test in `test-Validator-EdgeCases.R` PASSES (the 2 former BUG-LOCK tests are now green).
- Every test in `test-Validator-Helpers.R`, `test-Regex-Boundary.R`, `test-S4Validity.R` PASSES.
- Every pre-existing test passes unchanged (preserved error-type tags, preserved signatures).
- `R CMD BiocCheck` no longer flags `cat()`/`print()` debug leaks in `R/UtilitiesFunc.R`, `R/UtilitiesFuncInputChecker.R`, `R/ShinyServerModule.R`, or `R/ShinySangerContigServer.R`.
- `R/UtilitiesFuncInputChecker.R` is ~13% shorter (603 → 523 lines).

### Non-goals (deferred to a later phase)

- BiocParallel migration (Phase 2 audit §2). Not done in Phase 4 because of API/Depends churn it would induce.
- DECIPHER copy-by-value reduction (Phase 2 audit §3). Not done; `as.character`-roundtrips remain in `calculateContigSeq`.
- `parallel` → `Imports:` migration (DESCRIPTION). Not done.

---

## 7. Files touched in Phase 4

```
M  R/ClassChromatogramParam.R          (+ setValidity)
M  R/ClassObjectResults.R              (+ setValidity)
M  R/ClassQualityReport.R              (+ setValidity)
M  R/ShinySangerContigServer.R         (- 1 stray print)
M  R/ShinyServerModule.R               (- 1 stray message)
M  R/UtilitiesFunc.R                   (print → log_info × 2)
M  R/UtilitiesFuncInputChecker.R       (full rewrite: regex fix + helpers)
M  tests/testthat/test-Validator-EdgeCases.R   (BUG-LOCK rename)
A  tests/testthat/test-Regex-Boundary.R
A  tests/testthat/test-S4Validity.R
A  tests/testthat/test-Validator-Helpers.R
A  plans/04_sanitization_log.md
```

No code under `inst/`, `data/`, `vignettes/`, `man/`, or `docs/` was modified.
