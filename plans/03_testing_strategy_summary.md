# sangeranalyseR — Phase 3 Testing Strategy Summary

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

This document records the test-suite expansion landed in Phase 3 — what was added, what it covers, which tests are deliberately failing today (locking in the desired post-fix behaviour for Phase 4), and how to measure coverage going forward.

---

## 1. Test inventory

### Files

| File                                                          | Status   | Tests | Scope                                                                                                |
| ------------------------------------------------------------- | -------- | ----: | ----------------------------------------------------------------------------------------------------- |
| `tests/testthat.R`                                             | edited   |     — | `test_check("sangeranalyseR")` un-commented so `R CMD check` actually runs the suite.                 |
| `tests/testthat/helper-Fixtures.R`                             | new      |     — | 11 synthetic-fixture builders (under `tempdir()`, no `inst/extdata/` writes).                          |
| `tests/testthat/test-Constructors.R`                           | filled (was empty) |  3  | Wrapper smoke tests for `SangerRead()`, `SangerContig()`, `SangerAlignment()`.                          |
| `tests/testthat/test-Validator-EdgeCases.R`                    | new      |    21 | Negative tests: `.Xfa`/`.Xab1`/`.ab2`/`.fast` extensions, empty/corrupt/truncated `.ab1`, malformed CSV (4 schema variants), out-of-range scalars (10 variants). |
| `tests/testthat/test-OrthogonalAxes.R`                         | new      |     6 | Cartesian matrix ABIF/FASTA × REGEX/CSV × M1/M2.                                                       |
| `tests/testthat/test-SangerRead-SlotInvariants.R`              | new      |     4 | Byte-equality of `SangerRead` raw slots against the `sangerseq` parent.                                |
| `tests/testthat/test-IO-Batching.R`                            | new      |     3 | Hidden-files / nested-dirs behaviour + `list.files` call-count regression.                            |
| `tests/testthat/test-Coverage-Smoke.R`                         | new      |     6 | Breadth: `qualityBasePlot`, `MakeBaseCalls`, `readTable`, `writeFasta`, `updateQualityParam`, `data()`. |

### Pre-existing test files (unchanged)

| File                                                      | Tests | Scope (summary)                                            |
| --------------------------------------------------------- | ----: | ----------------------------------------------------------- |
| `test-SangerRead-ABIF.R`                                  |    ~25 | Forward/reverse ABIF reads, MakeBaseCalls variants.          |
| `test-SangerRead-FASTA.R`                                 |    ~10 | FASTA reads.                                                 |
| `test-SangerContig-ABIF.R`                                |     ~6 | Contig assembly from ABIF.                                   |
| `test-SangerContig-FASTA.R`                               |     ~3 | Contig from FASTA.                                           |
| `test-SangerAlignment-ABIF.R`                             |    ~10 | Alignment from ABIF.                                         |
| `test-SangerAlignment-FASTA.R`                            |     ~3 | Alignment from FASTA.                                        |
| `test-{SangerRead,SangerContig,SangerAlignment}-Compare.R` |    ~3 | Cross-mode equivalence checks.                              |
| `test-QualityReport.R`                                    |     ~5 | `QualityReport` direct construction + slot values.          |
| `test-interal-functions.R`                                |    ~10 | Trim/base-call helpers in `R/UtilitiesFunc.R`.              |

### Helpers

| File                                                    | Purpose                                      |
| ------------------------------------------------------- | -------------------------------------------- |
| `helper-Fixtures.R` (new)                               | Synthetic file builders.                      |
| `helper-SangerRead-{ABIF,FASTA}-{forward,reverse}.R`    | Pre-existing `sangerReadF` / `sangerReadR`.   |
| `helper-SangerContig-{ABIF,FASTA}.R`                    | Pre-existing contig fixtures.                 |
| `helper-SangerAlignment-{ABIF,FASTA}.R`                 | Pre-existing alignment fixtures.              |
| `helper-SangerRead-prechecking.R`                       | Pre-existing failure-object fixture.          |

Net new: **6 test files, 1 helper file, 43 new test_that blocks.**

---

## 2. Negative-test catalogue

The package's "construction never throws" contract means each negative test asserts on `obj@objectResults@creationResult == FALSE` and on `errorTypes`, *not* on `expect_error`. Three tests use `expect_error` because the failure path bypasses the validator and propagates from `read.abif`/`sangerseq` itself (empty / random-byte / truncated `.ab1`).

| #  | Scenario                                          | Fixture builder                              | Expected `errorType`           | Status today                                                                                                         |
| --: | ------------------------------------------------ | --------------------------------------------- | ------------------------------ | -------------------------------------------------------------------------------------------------------------------- |
|  1 | `.Xfa` extension on FASTA file                  | `make_xfa_file()`                             | `FILE_TYPE_ERROR`              | **FAIL today** — BUG-LOCK. `checkFASTA_File`'s `str_extract(.., ".fa$")` has an unescaped `.` that matches "Xfa". Fixed by escaping in Phase 4. |
|  2 | `.Xab1` extension on ABIF read                   | inline (file copy)                            | `FILE_TYPE_ERROR`              | **FAIL today** — BUG-LOCK. Same bug class in `checkReadFileName`.                                                    |
|  3 | `.ab2` extension on ABIF read                    | `make_ab2_file()`                             | `FILE_TYPE_ERROR`              | PASS                                                                                                                  |
|  4 | `.fast` extension on FASTA file                  | `make_fast_file()`                            | `FILE_TYPE_ERROR`              | PASS                                                                                                                  |
|  5 | Missing read filename                            | inline path                                   | `FILE_NOT_EXIST_ERROR`         | PASS                                                                                                                  |
|  6 | Empty `.ab1` (zero bytes)                        | `make_empty_ab1()`                            | (throws from `read.abif`)      | PASS — `expect_error`                                                                                                  |
|  7 | Random-byte `.ab1`                               | `make_corrupt_ab1()`                          | (throws from `read.abif`)      | PASS — `expect_error`                                                                                                  |
|  8 | Truncated `.ab1` (32-byte prefix)                | `make_corrupt_ab1_truncated()`                | (throws from `sangerseq`)      | PASS — `expect_error`                                                                                                  |
|  9 | CSV missing `contig` column                       | `make_csv_missing_column("contig")`           | `CSV_MISMATCH_ERROR`           | PASS                                                                                                                  |
| 10 | CSV missing `direction` column                    | `make_csv_missing_column("direction")`        | `CSV_MISMATCH_ERROR`           | PASS                                                                                                                  |
| 11 | CSV missing `reads` column                        | `make_csv_missing_column("reads")`            | `CSV_MISMATCH_ERROR`           | PASS                                                                                                                  |
| 12 | CSV with `direction` values other than F/R        | `make_csv_bad_direction()`                    | `CSV_VALUE_ERROR`              | PASS                                                                                                                  |
| 13 | CSV references reads not on disk                  | `make_mismatched_csv()`                       | (warning row only)             | PASS — alignment still builds                                                                                         |
| 14 | FASTA name in CSV not present in file             | `make_fasta_missing_record()`                 | `FASTA_NAME_NOT_EXIST` (or row) | PASS (best-effort assert)                                                                                             |
| 15 | `M1TrimmingCutoff = 5`                           | bundled                                       | `PARAMETER_RANGE_ERROR`        | PASS                                                                                                                  |
| 16 | `TrimmingMethod = "M3"`                           | bundled                                       | `PARAMETER_VALUE_ERROR`        | PASS                                                                                                                  |
| 17 | `processorsNum = "many"`                          | bundled                                       | `PARAMETER_TYPE_ERROR`         | PASS                                                                                                                  |
| 18 | `signalRatioCutoff = 1.5`                         | bundled                                       | `PARAMETER_RANGE_ERROR`        | PASS                                                                                                                  |
| 19 | `inputSource = "FAST5"`                           | bundled                                       | `PARAMETER_VALUE_ERROR`        | PASS                                                                                                                  |
| 20 | `processMethod = "GLOB"`                          | bundled                                       | `PARAMETER_VALUE_ERROR`        | PASS                                                                                                                  |
| 21 | `acceptStopCodons = "yes"` (string, not logical)  | bundled                                       | `PARAMETER_VALUE_ERROR`        | PASS                                                                                                                  |
| 22 | `readingFrame = 4`                                | bundled                                       | `PARAMETER_VALUE_ERROR`        | PASS                                                                                                                  |
| 23 | `minFractionCall = 1.5`                           | bundled                                       | `PARAMETER_RANGE_ERROR`        | PASS                                                                                                                  |
| 24 | `minReadsNum = 0`                                 | bundled                                       | `PARAMETER_VALUE_ERROR`        | PASS                                                                                                                  |

**Two BUG-LOCK tests fail today** (rows 1 and 2). They are not regressions; they pin desired behaviour and will turn green once Phase 4 fixes the regex (`"\\.fa(sta)?$"` instead of `".fa$"`).

---

## 3. Orthogonal-axis matrix

`tests/testthat/test-OrthogonalAxes.R` runs all 6 cells:

| Case      | inputSource | processMethod | TrimmingMethod | Fixture                                                                       | Previously tested?                                |
| --------- | ----------- | ------------- | -------------- | ----------------------------------------------------------------------------- | -------------------------------------------------- |
| AB-RE-M1  | ABIF        | REGEX         | M1             | `inst/extdata/Allolobophora_chlorotica/ACHLO/`                                | Partially — covered by `test-SangerAlignment-ABIF.R` |
| AB-RE-M2  | ABIF        | REGEX         | M2             | same                                                                          | **No**                                              |
| AB-CS-M1  | ABIF        | CSV           | M1             | `ACHLO/` + `inst/extdata/ab1/SangerAlignment/names_conversion.csv`            | Partially                                           |
| AB-CS-M2  | ABIF        | CSV           | M2             | same                                                                          | **No**                                              |
| FA-RE     | FASTA       | REGEX         | (forced "")    | `inst/extdata/fasta/SangerAlignment/Sanger_all_reads.fa`                      | Partially — covered by `test-SangerAlignment-FASTA.R` |
| FA-CS     | FASTA       | CSV           | (forced "")    | same FASTA + `inst/extdata/fasta/SangerAlignment/names_conversion.csv`        | **No**                                              |

For each case the assertions are:

```r
expect_true(sa@objectResults@creationResult)
expect_s4_class(sa, "SangerAlignment")
expect_gt(length(sa@contigList), 0L)
expect_s4_class(sa@contigsAlignment, "DNAStringSet")
expect_gt(length(sa@contigsConsensus), 0L)
```

Plus, for ABIF cases, every child `SangerRead`'s `creationResult` is `TRUE`.

The matrix exposed three previously-untested branches (AB-RE-M2, AB-CS-M2, FA-CS); confirming all six pass means the constructor's branch table is fully exercised.

---

## 4. Slot-invariant guard (Phase-4 safety net)

`tests/testthat/test-SangerRead-SlotInvariants.R` pins five slot-equality invariants between `SangerRead` and its `sangerseq` parent:

| Slot on `SangerRead`     | Compared to                                         | Why it matters                                                                                                  |
| ------------------------ | --------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| `@primarySeqRaw`         | `sangerseq(read.abif(file))@primarySeq`              | Phase 4 §4 plans to drop this slot and rely on the inherited parent slot.                                       |
| `@secondarySeqRaw`       | `sangerseq(read.abif(file))@secondarySeq`            | Same.                                                                                                            |
| `@peakPosMatrixRaw`      | `sangerseq(read.abif(file))@peakPosMatrix`           | Same.                                                                                                            |
| `@peakAmpMatrixRaw`      | `sangerseq(read.abif(file))@peakAmpMatrix`           | Same.                                                                                                            |
| `@traceMatrix` (inherited) | `sangerseq(read.abif(file))@traceMatrix`           | Inheritance contract — must remain identical even after slot collapse.                                          |

A fifth test verifies `length(@primarySeq) == length(@QualityReport@qualityPhredScores)`, the invariant that `MakeBaseCallsInside` is supposed to maintain post-basecall.

All five tests **PASS today** and serve as a regression net before any S4 slot pruning lands.

---

## 5. IO-batching guard

`tests/testthat/test-IO-Batching.R` adds three checks on `checkABIF_Directory` and the surrounding tree-walks:

| Test                                                    | What it asserts                                                                                                             |
| ------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| Hidden `.hidden_F.ab1` does not crash construction       | The validator silently ignores dot-files (default `list.files(all.files = FALSE)`); construction succeeds.                  |
| Nested sub-directory `.ab1` files are discovered         | `list.files(..., recursive = TRUE)` at the SA level picks up reads in `nested/`; `length(@contigList) >= 2`.                  |
| `list.files` is invoked > 1 time during a SA build        | Pins today's tree-walk count. After Phase 4's IO cache, this assertion will be tightened to `expect_equal(counter$n, 1L)`. |

The third test uses `base::trace(base::list.files, ...)` rather than `local_mocked_bindings` because `local_mocked_bindings` has restrictions on mocking base functions. The trace is removed via `on.exit(untrace(...))`.

**Baseline call count today** (ACHLO fixture, 4 contigs): expected ≥ 5 (1 SA-level + 1 per contig). Phase-4 target: 1.

---

## 6. Coverage baseline

To measure:

```r
install.packages("covr")          # if needed
cov <- covr::package_coverage(type = "tests", quiet = FALSE)
print(cov)
covr::report(cov, file = "coverage.html")   # do not commit
```

Per-file coverage targets (post Phase 3):

| File                                  | Pre-Phase-3 (estimated) | Phase-3 target |
| ------------------------------------- | -----------------------: | --------------: |
| `R/UtilitiesFuncInputChecker.R`       | ~25–30%                  | **> 60%**        |
| `R/Class*.R` (3 user-facing classes)  | ~40–55%                  | > 50%           |
| `R/UtilitiesFunc.R`                   | ~35–45%                  | > 45%           |
| `R/MethodShared.R` (façade dispatch)  | ~20%                     | > 70%           |
| `R/Method{Sanger*,QualityReport}.R`   | ~30–40%                  | > 45%           |
| `R/Constructors.R`                    | 0% (test file was empty) | **> 80%**        |

Areas intentionally **not** covered by automated tests:

- `R/ShinySangerContig{UI,Server}.R`, `R/ShinySangerAlignment{UI,Server}.R`, `R/ShinyServerModule.R` — exercised manually via `launchApp()`; ~4,500 lines, all reactive and requiring a browser.
- `inst/rmd/*.Rmd` report templates — exercised via `\dontrun{}` examples, not unit tests.
- `R/sangeranalyseR_show_method.R` — only `print` formatting, low value.

---

## 7. Known failing tests (locked-in failures)

Currently failing on `devel`, intentionally:

| Test name (testthat)                                                       | File                                       | Reason                                                            | Resolution                                                          |
| -------------------------------------------------------------------------- | ------------------------------------------ | ----------------------------------------------------------------- | -------------------------------------------------------------------- |
| `BUG-LOCK: .Xfa file is rejected as invalid FASTA extension`               | `test-Validator-EdgeCases.R`               | `checkFASTA_File` regex `".fa$"` has unescaped `.`                | Phase 4: change to `"\\.fa(sta)?$"` in `R/UtilitiesFuncInputChecker.R:344` |
| `BUG-LOCK: .Xab1 file is rejected as invalid ABIF extension`               | `test-Validator-EdgeCases.R`               | `checkReadFileName` regex `".ab1$"` has unescaped `.`             | Phase 4: change to `"\\.ab1$"` in `R/UtilitiesFuncInputChecker.R:414` |

CI / `R CMD check` will report these as test failures until the Phase-4 regex fix lands. **This is desired** — they document the bug and its desired remediation in executable form.

---

## 8. How to run

| Goal                            | Command                                                                                            |
| ------------------------------- | -------------------------------------------------------------------------------------------------- |
| Run the whole suite             | `devtools::test()`                                                                                  |
| Run one file                    | `testthat::test_file("tests/testthat/test-Validator-EdgeCases.R")`                                  |
| Run one test                    | `testthat::test_file(file, filter = "Xfa")`                                                          |
| Run via `R CMD check`           | `R CMD build .; R CMD check sangeranalyseR_*.tar.gz` (now actually executes the suite)               |
| Run via Bioconductor lint       | `R CMD BiocCheck sangeranalyseR_*.tar.gz`                                                            |
| Coverage report                 | `cov <- covr::package_coverage(); covr::report(cov, file = "coverage.html")`                         |
| Just the validator coverage     | `covr::file_coverage("R/UtilitiesFuncInputChecker.R", "tests/testthat/test-Validator-EdgeCases.R")`  |

To run the existing failing BUG-LOCK tests in isolation:

```r
testthat::test_file(
    "tests/testthat/test-Validator-EdgeCases.R",
    filter = "BUG-LOCK"
)
```

Both should report **fail** until Phase 4 lands the regex fix.
