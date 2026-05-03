# sangeranalyseR — Phase 2 Technical-Debt & Performance Audit

Working tree: `/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR`, branch `devel`, version 1.21.1 (Bioc 3.24).

---

## 1. Validation overhead — the 31 `check*` functions

### Where they live and how they are invoked

All 31 validators sit in `R/UtilitiesFuncInputChecker.R`. Every one of them follows the same shape:

```r
checkX <- function(value, errors, errorTypes) {
    if (...) {
        errors     <- c(errors, msg)
        errorTypes <- c(errorTypes, "PARAMETER_*_ERROR")
    }
    return(list(errors, errorTypes))
}
```

So every call performs **two `c()` vector-grow operations** plus **constructs and returns a 2-element `list`** that the caller immediately destructures into the next call. R copies the `errors` / `errorTypes` vectors on each `c()`. With 31 validators that's 62 grow-and-copy operations on every successful construction even when nothing fails — quadratic in the number of accumulated errors.

### Are they called inside hot loops?

**Yes, indirectly — but the loop nesting is mostly hidden by `printLevel`.** Tracing call sites:

| Construction site                                           | Validators run there                                                                                          |
| ----------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `SangerAlignment.initialize` (`R/ClassSangerAlignment.R:155`) | 11 SA-level checks (input source, process method, genetic code, refAAS, minReadsNum/Length/Frac, maxFracLost, acceptStop, readingFrame, processorsNum) + chromatogram/trim checks |
| `SangerContig.initialize`   (`R/ClassSangerContig.R:160`)     | Same set + `checkContigName`, `checkAb1FastaCsv` (heavy — see below)                                            |
| `SangerRead.initialize`     (`R/ClassSangerRead.R:135`)        | `checkReadFileNameExist`, `checkReadFileName`, `checkInputSource`, `checkReadFeature`, `checkGeneticCode`, plus all chromatogram/trim checks for ABIF |

The `printLevel == "SangerRead"` guard inside `SangerRead.initialize` already skips the upper-axis checks when the read is being constructed by a parent contig. So when `SangerAlignment` builds, say, **N reads across M contigs**, the per-read validators that *do* run for each child are:

```r
checkReadFileNameExist(readFileName, ...)   # file.exists() — N stat() calls
checkReadFileName(readFileName, inputSource, ...)  # str_extract regex
```

Plus `MakeBaseCallsInside` etc. So the *per-read* validation footprint is small (two `file.exists` + a regex per ABIF file). The real waste is at the **alignment level**: every parameter that's already been checked by `SangerAlignment.initialize` is *not* re-validated by children — but `SangerAlignment.initialize` itself runs 15+ `c()` shuffles for what amount to compile-time-constant range checks.

### What can be consolidated

The 31 functions collapse cleanly into 4 buckets:

| Bucket                  | Validators                                                                                                                                    | Replacement                                                                                                              |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| Enum value              | `checkInputSource`, `checkProcessMethod`, `checkReadFeature`, `checkAcceptStopCodons`, `checkShowTrimmed`, `checkReadingFrame`               | A single helper `enumCheck(name, value, allowed, errors)` driven by a static list. Kills 6 near-duplicate functions.     |
| Numeric range / integer | `checkMinReadsNum`, `checkMinReadLength`, `checkMinFractionCall`, `checkMaxFractionLost`, `checkProcessorsNum`, `checkBaseNumPerRow`, `checkHeightPerRow`, `checkSignalRatioCutoff`, `checkRefAAS` (length 0 case implicit) | A single `numCheck(name, value, lo, hi, integer = TRUE/FALSE, errors)` driven by a small validator table. Kills 9 functions. |
| File / path             | `checkReadFileNameExist`, `checkReadFileName`, `checkABIF_Directory`, `checkFASTA_File`, `checkCSV_NamesConversion`, `checkTargetFastaName`  | Single helper `pathCheck(name, path, kind, ext = c(".ab1"), errors)`. Kills 5.                                            |
| Cross-cutting / "did the matcher actually find anything" | `checkGreplForward`, `checkGreplReverse`, `checkCSVConvForward`, `checkCSVConvReverse` | Already collapsible to `checkMatchedAny(label, vec, type = "REGEX_MATCH_WARN" | "CSV_MATCH_WARN", warnings)`. Kills 4.    |
| Special / structural    | `checkAb1FastaCsv`, `checkTrimParam`, `checkGeneticCode`, `checkContigName`, `checkQualityPhredScores`                                       | Keep as-is — they encode real cross-field logic.                                                                          |

So **31 → ~9** validators (5 specials + 4 generic helpers driven by tables).

### S4 `setValidity` — should we use it?

S4 validity is the right tool for **shape** invariants ("`minFractionCall` must be in [0,1]", "`TrimmingMethod` ∈ {M1,M2}"), and *not* for cross-object I/O concerns ("does this CSV file exist on disk and have a `direction` column"). Recommendation:

- Move all **enum** and **range** checks into `setValidity("SangerRead", ...)`, `setValidity("SangerContig", ...)`, `setValidity("QualityReport", ...)`. R's S4 system runs them automatically on `new(...)`, with proper error chaining (`validObject` returns *all* invariant violations at once instead of the home-grown `errors` list).
- Keep **path/IO/cross-field** checks (`checkABIF_Directory`, `checkAb1FastaCsv`, `checkFASTA_File`, `checkCSV_NamesConversion`, `checkTargetFastaName`, `checkGreplForward/Reverse`, `checkCSVConvForward/Reverse`) as procedural validators — they have side effects (`file.exists`, `read.csv`, `read.fasta`) and produce *warnings*, not invariant failures.
- The package's "construction never throws — it returns an object whose `objectResults@creationResult = FALSE`" contract is incompatible with vanilla `setValidity` (which throws). Reconcile by wrapping `validObject(...)` in `tryCatch` inside each `initialize` and folding the condition into the `errors` accumulator. Net: one `tryCatch` replaces 15+ `check*` calls per `initialize`.

### Specific overhead pattern: `c(errors, msg)` growth

Every validator does `errors <- c(errors, msg)`. That's a copy of the whole vector. Replace with `errors[[length(errors) + 1L]] <- msg` *or* (better) collect into `list()` and `unlist()` at the end. With error-free runs this is invisible, but on a malformed batch the cost is `O(K^2)` for K errors.

### Top 5 slowest validation checks

Ranked by *latency in the realistic call path*, not LOC:

| Rank | Function                        | Why it's slow                                                                                                                                      | Mitigation                                                                                                                                      |
| ---: | ------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **1** | `checkAb1FastaCsv` (line 475)   | Re-reads the **entire CSV** with `read.csv` and, for FASTA input, re-reads the **entire FASTA** with `read.fasta` — *just to validate*. Then runs two `lapply` membership scans of N×M. The same files are read again moments later inside the constructor. | Cache parsed CSV/FASTA on the parent object (or pass it in as an arg). Replace `lapply` membership probes with `setdiff(sourceReads, csvReads)` (single C-level call). |
| **2** | `checkFASTA_File` (line 327)    | `file.exists` + 2× `str_extract` regex on every call. Also has a leftover `cat()` debug print. The regex pattern `".fa$"` is wrong (missing `\\`), so it matches "Xfa" filenames too. | Fix regex (`"\\.fa(sta)?$"`), drop `cat`, single `tools::file_ext()` call.                                                                       |
| **3** | `checkABIF_Directory` (line 175) | Calls `dir.exists` on the user-supplied path. Cheap by itself, but `SangerAlignment` then immediately calls `list.files(ABIF_Directory, recursive = TRUE)` and **does it again** inside `SangerContig` for each contig. The directory tree gets walked repeatedly. | Walk once at the top level; pass `parentDirFiles` down through the constructor stack.                                                            |
| **4** | `checkReadFileNameExist` + `checkReadFileName` (lines 402, 412) | Run **per `SangerRead`**, i.e. once per ABIF file. For an alignment of ~100 reads this is 200 `file.exists`/`str_extract` syscalls. With network filesystems each `stat()` is non-trivial. | Validate the whole batch at the parent: one `file.exists()` on a vector; pass the boolean mask down.                                             |
| **5** | `checkQualityPhredScores` (line 461) | `all(qualityPhredScores %% 1 == 0)` allocates an N-length numeric vector and an N-length logical for every read just to verify integrality. For typical 800-bp reads × hundreds of reads this allocates GBs over a session. | `is.integer(qualityPhredScores)` would suffice (the upstream `MakeBaseCallsInside` already produces integer Phreds), or use `isTRUE(all(qualityPhredScores == as.integer(qualityPhredScores)))` short-circuited. |

Honourable mention: `checkTrimParam` (line 199) duplicates the M2 range check that `M2inside_calculate_trimming` then performs again — collapse to one site of truth.

---

## 2. Parallel efficiency — `parallel` vs `BiocParallel`

### Inventory of parallel call sites

```
R/UtilitiesFunc.R:29   processors = detectCores(all.tests = FALSE, logical = FALSE)
R/UtilitiesFunc.R:159  r = mclapply(is, oneAmbiguousColumn, aln=aln, mc.cores = processorsNum)        # countCoincidentSp
R/UtilitiesFunc.R:256  stops = as.numeric(unlist(mclapply(frReadSet, countStopSodons, ...)))         # calculateContigSeq, refAAS path
R/UtilitiesFunc.R:283  stops = as.numeric(unlist(mclapply(frReadSet, countStopSodons, ...)))         # calculateContigSeq, no-refAAS path
R/UtilitiesFunc.R:324  diffs = mclapply(aln, nPairwiseDiffs, subject = consensus, ...)               # calculateContigSeq, post-alignment
```

Plus the imports declared in `R/sangeranalyseR_package.R:33` (`@importFrom parallel mclapply detectCores`) and `DESCRIPTION` listing `parallel` in `Depends`. There are **zero** uses of `BiocParallel`, `bplapply`, `parLapply`, `makeCluster`, or `SnowParam`. There is also no `future`/`furrr`.

### Problems

1. **Windows correctness.** `parallel::mclapply` falls back to a **single core** on Windows because forking is not available. The package "documents" this implicitly in `getProcessors` (`R/UtilitiesFunc.R:4`), which forces `processorsNum = 1` whenever the OS isn't macOS/Linux. So Windows users effectively run serial today even when they pass `processorsNum = 8`. This is a correctness/UX bug, not just a performance one.
2. **Nested forking.** `calculateContigSeq` is called from inside the `lapply` over contigs in `SangerAlignment.initialize` (no parallelism at that outer loop, fortunately), and `mclapply` is called from inside `calculateContigSeq` *up to three times in a row* (lines 256, 283, 324) and once again in `countCoincidentSp` (line 159). Each `mclapply` forks; on macOS the per-fork copy of an in-memory `DNAStringSet` and the entire R session is non-trivial. For typical Sanger inputs (≤ 100 reads × ~800 bp) the **fork cost dominates** — `nPairwiseDiffs` is a `compareStrings` + two `str_count` calls, work that completes in microseconds per element. Parallelizing it via `mclapply` adds ~50–200 ms of fork overhead per call, often making the parallel version slower than serial.
3. **Wrong granularity.** The expensive work in this package is:
   - `read.abif` + `MakeBaseCallsInside` per `SangerRead` (one `.ab1` parse + per-peak basecalling). **Currently serial.**
   - `AlignSeqs` / `AlignTranslation` in `calculateContigSeq` and `alignContigs`. **DECIPHER is already internally parallel** via its own `processors=` arg.
   - `CorrectFrameshifts`. **Already DECIPHER-parallel.**
   The work that's currently `mclapply`-ed (`nPairwiseDiffs`, `oneAmbiguousColumn`, `countStopSodons`) is fast per-element. We're parallelizing the wrong things.
4. **No registered backend, no abstraction.** Bioconductor convention is that any package with parallel work accepts a `BPPARAM` argument that defaults to `bpparam()`, so users can plug `SerialParam`, `MulticoreParam`, `SnowParam`, or `BatchtoolsParam` without code changes. Today the only knob is an integer `processorsNum`.
5. **`detectCores` is discouraged in libraries.** R-core and CRAN policy: don't auto-grab all cores in a package. `getProcessors()` does exactly that, defeating user/cluster scheduler limits.

### Migration plan to `BiocParallel`

**Goal:** replace every `mclapply` and the `getProcessors` helper with `BiocParallel::bplapply(..., BPPARAM = BPPARAM)` plumbed through every entry point.

1. **DESCRIPTION**
   - Add `BiocParallel` to `Depends:` (or `Imports:`, preferred — most of the package's `Depends:` entries should be `Imports:` anyway, but that's a separate cleanup).
   - Remove `parallel` from `Depends:` (only `parallel::detectCores` would remain, and `BiocParallel` provides equivalents).
2. **Public API**
   - In `SangerRead()`, `SangerContig()`, `SangerAlignment()`, `MakeBaseCalls`, `updateQualityParam`, and `launchApp*` / `generateReport*`, add `BPPARAM = bpparam()` (with `processorsNum` retained as a deprecated alias that constructs a `MulticoreParam(processorsNum)` for one or two release cycles).
   - In `calculateContigSeq`, `alignContigs`, `countCoincidentSp`, replace the `processorsNum` argument with `BPPARAM`.
3. **Replace each call site**

   | Today                                                           | After                                                            |
   | --------------------------------------------------------------- | ---------------------------------------------------------------- |
   | `mclapply(is, oneAmbiguousColumn, aln=aln, mc.cores = processorsNum)` | `bplapply(is, oneAmbiguousColumn, aln=aln, BPPARAM = BPPARAM)`    |
   | `mclapply(frReadSet, countStopSodons, …, mc.cores = processorsNum)`   | `bplapply(frReadSet, countStopSodons, …, BPPARAM = BPPARAM)`      |
   | `mclapply(aln, nPairwiseDiffs, subject = consensus, mc.cores = …)`    | `bplapply(aln, nPairwiseDiffs, subject = consensus, BPPARAM = …)` |
   | `getProcessors(processorsNum)` (then passed to `DECIPHER::*(processors = …)`) | `bpnworkers(BPPARAM)` for DECIPHER's `processors=` argument; or pass `SerialParam()` and let DECIPHER handle threading. |
4. **Right-size the granularity.** Do this *during* the migration, not after:
   - Remove `mclapply` from `nPairwiseDiffs` and `oneAmbiguousColumn` — they're too cheap. Use plain `lapply`.
   - Keep `bplapply` for the `countStopSodons` loop (which translates each read 3× into AA — non-trivial).
   - **Add** `bplapply` to the per-read loops in `SangerContig.initialize` (the four `lapply(forwardAllReads, function(forwardN) new("SangerRead", …))` blocks, lines 311 / 332 / 370 / 398 / 451 / 472 / 511 / 538 of `R/ClassSangerContig.R`). Each iteration does a `read.abif` + `MakeBaseCallsInside` + `QualityReport` build, which is the actual bottleneck.
   - Honour `bpprogressbar(BPPARAM)` so the Shiny app can wire a progress callback.
5. **Don't double-parallelise.** When DECIPHER does its own threading via `processors=`, don't wrap *it* in `bplapply` too. Either let DECIPHER thread and use `SerialParam()` outside, or use `bplapply` outside and pass `processors = 1` to DECIPHER. Default policy: leave DECIPHER threading on; serial outside.
6. **Tests.** Add a test that runs the whole pipeline under `BiocParallel::register(SerialParam())` and one that runs it under `MulticoreParam(2)` (skipped on Windows). Fixture: existing `inst/extdata/Allolobophora_chlorotica/ACHLO/`.
7. **Windows.** `BiocParallel::SnowParam()` is the portable choice. Document `MulticoreParam` as Unix-only; let `bpparam()` pick the right default.
8. **Deprecation window.** Keep `processorsNum` accepted for two devel cycles, emit `.Deprecated("BPPARAM")` when it's non-NULL.

### Estimated wins

- **Windows:** changes from "always serial" to "actually parallel" via `SnowParam`. Single biggest correctness win.
- **macOS/Linux on small jobs (≤ 20 reads, ≤ 5 contigs):** a measurable speed-up by *removing* the `mclapply` calls on `nPairwiseDiffs` / `oneAmbiguousColumn` (fork overhead currently dominates).
- **Larger jobs (≥ 100 reads):** parallelizing `SangerRead` construction is roughly linear in worker count up to ~8 cores; that's where most wall-clock time is spent today.

---

## 3. DECIPHER integration — copy-by-value bottlenecks

### Where data is converted

`R/UtilitiesFunc.R:182` (`calculateContigSeq`) and `R/UtilitiesFunc.R:46` (`alignContigs`) are the integration surface to DECIPHER. Sequence-format conversions trace as follows for every contig:

```
SangerRead@primarySeq   (DNAString, in C-level XStringSet pool)
  └─ as.character()         ← copy 1: DNAString → character
      └─ substr(...)        ← copy 2: character window for trimming
          └─ unlist(c(fwd, rev))
              └─ DNAStringSet(...)   ← copy 3: character → DNAStringSet
                  ├─ CorrectFrameshifts(...)   → DECIPHER returns new XStringSet
                  ├─ AlignSeqs / AlignTranslation → new XStringSet
                  ├─ ConsensusSequence → DNAString
                  ├─ DistanceMatrix    → matrix
                  ├─ Treeline          → list
                  ├─ RemoveGaps(DNAStringSet(consensus))[[1]]  ← copy 4: re-wrap to drop gaps
                  └─ countCoincidentSp(aln, ...) → uses subseq + as.character per column (copy 5×width)
```

So a single contig does **3 round-trip conversions** between `DNAString[Set]` and base-R character before alignment even begins. DECIPHER's API works on `XStringSet`, so the right path is:

1. Build a `DNAStringSet` directly from `SangerRead@primarySeq` slots **without going through `as.character`**. Use `XVector::subseq` (or `Biostrings::subseq`) to apply trimming windows; `subseq` returns a view, not a copy. For reverse reads, `reverseComplement(subseq(primarySeq, start, end))`.
2. Replace `unlist(c(fwd, rev))` + `DNAStringSet(…)` with `c(fwdSet, rcRevSet)` where both are already `DNAStringSet`s — that's an O(1) concatenation of internal pointers.

### `oneAmbiguousColumn` — the real allocator

`R/UtilitiesFunc.R:170` currently does, **for every alignment column**:

```r
col = as.character(subseq(aln, i, i))
str = paste(col, sep="", collapse="")
```

That's two allocations per column × `aln@ranges@width[1]` columns × N reads. For a 1500-column alignment of 30 reads it's ~90,000 string allocations. DECIPHER provides `consensusMatrix(aln, as.prob = FALSE)` and `Biostrings::alphabetFrequency(aln, baseOnly = FALSE)` which return a matrix of counts in one C-level call. Replace the column-by-column scan with one `consensusMatrix` and a vectorised "count of IUPAC ambiguity codes per column ≥ 2" expression.

### `as.character` round-trips at module boundaries

```
R/UtilitiesFunc.R:53  as.character(SangerContig@contigSeq)             # alignContigs
R/UtilitiesFunc.R:191 as.character(forwardRead@primarySeq)             # calculateContigSeq fwd
R/UtilitiesFunc.R:201 as.character(reverseComplement(reverseRead@primarySeq))  # rev
R/UtilitiesFunc.R:789 as.character(SangerReadInst@primarySeq)          # SangerReadInnerTrimming
```

In every case the value is reassembled into a `DNAStringSet` immediately afterwards. Drop `as.character`; work on `XStringSet` directly. Estimated saving: a few ms per read but it scales with sequence length and matters more in `updateQualityParam` (called every time the user moves a Shiny slider).

### `RemoveGaps(DNAStringSet(consensus))[[1]]`

Line 343. `DNAStringSet(consensus)` wraps a single `DNAString` into a set of length 1 just so `RemoveGaps` accepts it; the `[[1]]` then unwraps. `RemoveGaps` accepts an `XStringSet`, so create the consensus as a 1-element `DNAStringSet` from the start (it already comes from `ConsensusSequence(aln, ...)[[1]]`) and avoid the wrap/unwrap.

### `frReadSet[which(frReadSetLen>0)]`

Line 262. `which()` materializes an integer vector of indices; `frReadSet[...]` on `DNAStringSet` accepts a logical vector directly. Saves an allocation, more importantly clarifies intent.

### Recommendations (DECIPHER section)

| Action                                                                                         | File:line                       | Expected effect                    |
| ---------------------------------------------------------------------------------------------- | ------------------------------- | ---------------------------------- |
| Build forward/reverse `DNAStringSet`s with `subseq` + `reverseComplement` instead of `as.character` + `substr` | `UtilitiesFunc.R:190–215`       | Eliminate 2× full-sequence copies per read |
| Replace `mclapply(is, oneAmbiguousColumn, …)` with a single `consensusMatrix` / `alphabetFrequency` call | `UtilitiesFunc.R:155–177`        | One C call instead of N×width R closures |
| Drop `DNAStringSet(consensus)` wrap/unwrap in `RemoveGaps`                                     | `UtilitiesFunc.R:343`            | Trivial, but cleaner                 |
| Use logical indexing on `DNAStringSet`                                                         | `UtilitiesFunc.R:262`            | Trivial                              |
| Same conversion fixes in `alignContigs` (`vapply(... as.character) → DNAStringSet`)            | `UtilitiesFunc.R:51–55`          | One copy avoided per `SangerAlignment` build |
| Cache trimmed sequences on the `SangerRead` (or its `QualityReport`) so re-alignment in `updateQualityParam` doesn't re-trim from raw on every Shiny slider tick | `MethodSangerContig.R`/`Alignment.R` | Big win for interactive UX        |

---

## 4. Inheritance audit — `SangerRead contains="sangerseq"`

### What we use from `sangerseqR`

- `read.abif(filename)` → returns the `abif` raw container (data block, header, raw trace channels). Used once, on construction (`R/ClassSangerRead.R:200`).
- `sangerseq(abifRawData)` → wraps `abif` into a `sangerseq` S4 with `traceMatrix`, `peakPosMatrix`, `peakAmpMatrix`, `primarySeq`, `secondarySeq`, `primarySeqID`, `secondarySeqID`. Used once, on construction (`R/ClassSangerRead.R:204`). All slots are then **pulled out** and reassigned into the new `SangerRead`'s own slots.
- `chromatogram(...)` is **shadowed** by `chromatogram_overwrite` in `R/UtilitiesFunc.R:802`, exported from this package. We don't actually call `sangerseqR::chromatogram` in user-facing flows.
- `primarySeq` is `@importFrom`-ed but is just a slot accessor — it works on `sangerseq`/`SangerRead` because of inheritance.

So the practical surface from `sangerseqR` is **two functions** (`read.abif`, `sangerseq`) plus the inheritance contract.

### Is the inheritance pulling its weight?

Probably not. The `SangerRead` initializer **copies every meaningful slot out of the temporary `sangerseq` object into its own slots** (`primarySeqID`, `secondarySeqID`, `primarySeqRaw`, `secondarySeqRaw`, `traceMatrix`, `peakPosMatrixRaw`, `peakAmpMatrixRaw`). After that, the `sangerseq` parent slots are mostly redundant duplicates. `sangerseq` also carries its own `primarySeq` / `secondarySeq` / `peakPosMatrix` / `peakAmpMatrix` slots, so the child has *both* `primarySeqRaw` (own slot) and `primarySeq` (parent slot via inheritance), which is a maintenance trap.

**Options:**

- (a) **Keep inheritance** and stop duplicating slots. Drop `primarySeqRaw`, `secondarySeqRaw`, `peakPosMatrixRaw`, `peakAmpMatrixRaw` from `SangerRead` and instead read them via the inherited slots from the parent `sangerseq` (which itself stores the values that came out of `sangerseq()`). Saves ~5 slots × N reads of memory and removes the "two sources of truth" risk.
- (b) **Switch to composition.** Replace `contains = "sangerseq"` with a `slot sangerseqRaw = "sangerseqORNULL"` and explicitly delegate. Cleaner, but a public API change for any user that was relying on `is(x, "sangerseq")`.

Recommend (a) for this release; (b) only if a 2.x major bump is on the table.

### Modernizing legacy `sangerseqR` calls with `pwalign`

`pwalign` is the new home for `Biostrings::pairwiseAlignment` and `compareStrings` (the Bioc 3.19/3.20 split moved them out of Biostrings to slim it down). The package **already imports `pwalign::compareStrings`** (line 22 of `R/sangeranalyseR_package.R`, line 150 of `R/UtilitiesFunc.R`), so the migration is partially done. There is **no use of `pairwiseAlignment` anywhere** in the package — sequence alignment is delegated to DECIPHER's `AlignSeqs`/`AlignTranslation`.

What's *not* yet modernised:

- `R/UtilitiesFunc.R:155–177` (`countCoincidentSp` / `oneAmbiguousColumn`) does ad-hoc per-column ambiguity detection in pure R. As noted in §3, this should become a single `Biostrings::consensusMatrix` or `alphabetFrequency` call — doesn't need `pwalign` but is in the same neighbourhood.
- The current `pwalign::compareStrings` use is fine; could be vectorised over the entire alignment (it accepts `XStringSet` patterns) instead of being called inside `mclapply`. That's a micro-optimisation (~2× fewer R-level calls) and also makes the code work without parallel.

There are no remaining `Biostrings::pairwiseAlignment`-style legacy calls in this codebase, so the migration cost is low.

---

## Top-5 prioritised punch list

| # | Item                                                                                                                    | Files                                                                  | Effort | Impact                                                  |
| - | ----------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- | ------ | ------------------------------------------------------- |
| 1 | Migrate `parallel::mclapply` → `BiocParallel::bplapply`; add `BPPARAM` to public API; deprecate `processorsNum`         | `R/UtilitiesFunc.R`, `R/Constructors.R`, `R/Method*.R`, `DESCRIPTION`   | M      | Fixes Windows correctness; restores meaningful parallelism |
| 2 | Cache CSV/FASTA read in `checkAb1FastaCsv` and pass parsed structures through; stop re-reading the file in validators    | `R/UtilitiesFuncInputChecker.R`, `R/Class*Alignment*.R`, `R/Class*Contig*.R` | S      | Removes the slowest validator                             |
| 3 | Eliminate `DNAString[Set] ↔ character` round-trips; replace `oneAmbiguousColumn` loop with `consensusMatrix`             | `R/UtilitiesFunc.R:155–177, 190–215, 51–55, 343`                        | S      | Trims memory churn; speeds up Shiny re-trim flow         |
| 4 | Collapse 31 validators → 4 generic helpers + 5 specials; route shape checks through `setValidity`                       | `R/UtilitiesFuncInputChecker.R`, `R/Class*.R`                            | M      | Smaller code, faster construction, better error reports  |
| 5 | Drop redundant `*Raw` slots on `SangerRead`; rely on inherited `sangerseq` slots                                        | `R/ClassSangerRead.R`                                                   | S      | Less memory per read; one source of truth                  |

---

## Compliance notes

- `parallel` should not be in `Depends:` — no symbols need to be on the search path. Move to `Imports:` (along with most other `Depends:` entries — `BiocCheck` will flag this).
- `getProcessors()` calling `detectCores()` violates Bioconductor/CRAN multicore guidelines. Replace with `BiocParallel::bpparam()` defaulting to `SerialParam()`.
- The `cat("FASTA_File", FASTA_File)` debug print at `R/UtilitiesFuncInputChecker.R:338` will trip `BiocCheck`'s "no `cat`/`print` in package code" rule. Use `log_info` (already imported).
- `R/UtilitiesFunc.R:362` (`MakeBaseCallsInside`) and the rest of `UtilitiesFunc.R` would benefit from `@noRd` roxygen tags so they don't accidentally appear as exported help pages.

---

## Verification

After implementing the punch list:

1. `devtools::test()` — full testthat suite passes.
2. `R CMD check --as-cran` — no new NOTEs/WARNINGs.
3. `R CMD BiocCheck` — `parallel` warning gone; `detectCores` warning gone; `cat()` warning gone.
4. Benchmark on `inst/extdata/Allolobophora_chlorotica/ACHLO/` (8 ABIF files): wall-clock for `SangerAlignment(...)` should drop on macOS/Linux and become non-trivially parallel on Windows for the first time.
5. `BiocParallel::register(SerialParam()); SangerAlignment(...)` produces identical output to `BiocParallel::register(MulticoreParam(2)); SangerAlignment(...)` (regression-test the consensus and tree).
