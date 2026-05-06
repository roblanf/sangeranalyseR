## Phase 5 — End-to-end verification + profiling.
##
## Loads the development tree (NOT the installed package), runs the full
## SangerAlignment pipeline against the bundled ACHLO fixture under both M1
## and M2 trimming, validates structural invariants, and produces both
## per-stage Sys.time() breakdowns and an Rprof summary.
##
## Output files (written next to this script):
##   timings.csv        Per-stage wall-clock timings for M1 and M2.
##   profile_M1.txt     summaryRprof() output for the M1 run.
##   profile_M2.txt     summaryRprof() output for the M2 run.
##   accuracy.txt       Consensus / alignment / read-result invariants.
##   slot_invariants.txt  *Raw slot vs. sangerseq parent slot equality.

suppressPackageStartupMessages({
    library(devtools)
})

repo_root <- "/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR"
stopifnot(dir.exists(file.path(repo_root, "R")))
out_dir <- file.path(repo_root, "plans", "phase5_artifacts")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Repo:    ", repo_root, "\n")
cat("Out dir: ", out_dir,   "\n")
cat("Loading dev tree ...\n")
suppressMessages(load_all(repo_root, quiet = TRUE))
cat("sangeranalyseR loaded; version:", as.character(packageVersion("sangeranalyseR")), "\n\n")

ab1_dir <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                       package = "sangeranalyseR")
stopifnot(dir.exists(ab1_dir))
ab1_files <- list.files(ab1_dir, pattern = "\\.ab1$", full.names = TRUE)
cat("Fixture: ", ab1_dir, "\n")
cat("Files:   ", length(ab1_files), "\n\n")

REGEX_F <- "_[0-9]*_F.ab1$"
REGEX_R <- "_[0-9]*_R.ab1$"

## -- helpers -----------------------------------------------------------------

stage_time <- function(label, expr) {
    gc(verbose = FALSE)
    t0 <- Sys.time()
    val <- force(eval(expr, envir = parent.frame()))
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  [%-22s] %.3f s\n", label, elapsed))
    list(value = val, elapsed = elapsed)
}

run_pipeline <- function(method = c("M1","M2")) {
    method <- match.arg(method)
    cat(sprintf("\n=== Pipeline run: TrimmingMethod = %s ===\n", method))
    args <- list(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = ab1_dir,
        REGEX_SuffixForward = REGEX_F,
        REGEX_SuffixReverse = REGEX_R,
        TrimmingMethod      = method,
        processorsNum       = 1
    )
    if (method == "M1") {
        args$M1TrimmingCutoff     <- 0.0001
        args$M2CutoffQualityScore <- NULL
        args$M2SlidingWindowSize  <- NULL
    } else {
        args$M1TrimmingCutoff     <- NULL
        args$M2CutoffQualityScore <- 20
        args$M2SlidingWindowSize  <- 10
    }
    timings <- list()
    timings$ingest_and_assemble <- stage_time(
        sprintf("SangerAlignment %s", method),
        quote(do.call(SangerAlignment, args))
    )
    sa <- timings$ingest_and_assemble$value

    timings$writeFasta <- stage_time(
        "writeFasta (all)",
        quote({
            tmp <- file.path(tempdir(), paste0("e2e_", method))
            dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
            writeFastaSA(sa, outputDir = tmp)
            tmp
        })
    )
    timings$updateQP <- stage_time(
        "updateQualityParam",
        if (method == "M1")
            quote(updateQualityParam(sa, TrimmingMethod = "M2",
                                     M1TrimmingCutoff = NULL,
                                     M2CutoffQualityScore = 25,
                                     M2SlidingWindowSize = 12))
        else
            quote(updateQualityParam(sa, TrimmingMethod = "M1",
                                     M1TrimmingCutoff = 0.001,
                                     M2CutoffQualityScore = NULL,
                                     M2SlidingWindowSize = NULL))
    )
    list(sa = sa, timings = timings)
}

## -- run M1 + M2 -------------------------------------------------------------

cat("\n----- Warming up (one untimed run to amortise JIT/cache) -----\n")
invisible(run_pipeline("M1"))

m1 <- run_pipeline("M1")
m2 <- run_pipeline("M2")

## -- Rprof on a fresh M1 run -------------------------------------------------

cat("\n=== Rprof: profiling a fresh M1 SangerAlignment build ===\n")
prof_file_M1 <- file.path(out_dir, "Rprof_M1.out")
Rprof(prof_file_M1, interval = 0.005, line.profiling = FALSE)
sa_prof <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = ab1_dir,
    REGEX_SuffixForward = REGEX_F,
    REGEX_SuffixReverse = REGEX_R,
    TrimmingMethod      = "M1",
    M1TrimmingCutoff    = 0.0001,
    processorsNum       = 1
)
Rprof(NULL)
cat("M1 profile capture complete\n")

prof_file_M2 <- file.path(out_dir, "Rprof_M2.out")
Rprof(prof_file_M2, interval = 0.005, line.profiling = FALSE)
sa_prof_m2 <- SangerAlignment(
    inputSource          = "ABIF",
    processMethod        = "REGEX",
    ABIF_Directory       = ab1_dir,
    REGEX_SuffixForward  = REGEX_F,
    REGEX_SuffixReverse  = REGEX_R,
    TrimmingMethod       = "M2",
    M2CutoffQualityScore = 20,
    M2SlidingWindowSize  = 10,
    processorsNum        = 1
)
Rprof(NULL)
cat("M2 profile capture complete\n")

prof_summary_str <- function(prof_file) {
    s <- summaryRprof(prof_file)
    out <- character()
    out <- c(out, "by.self (top 25):")
    bs <- head(s$by.self, 25L)
    bs$self.time     <- sprintf("%6.3f", bs$self.time)
    bs$total.time    <- sprintf("%6.3f", bs$total.time)
    bs$self.pct      <- sprintf("%5.1f", bs$self.pct)
    bs$total.pct     <- sprintf("%5.1f", bs$total.pct)
    out <- c(out, capture.output(print(bs)))
    out <- c(out, "", "by.total (top 15):")
    bt <- head(s$by.total, 15L)
    bt$self.time     <- sprintf("%6.3f", bt$self.time)
    bt$total.time    <- sprintf("%6.3f", bt$total.time)
    bt$self.pct      <- sprintf("%5.1f", bt$self.pct)
    bt$total.pct     <- sprintf("%5.1f", bt$total.pct)
    out <- c(out, capture.output(print(bt)))
    out <- c(out, "",
             sprintf("Total sampling time: %.3f s (%d samples)",
                     s$sampling.time, nrow(s$by.self)))
    out
}

writeLines(prof_summary_str(prof_file_M1), file.path(out_dir, "profile_M1.txt"))
writeLines(prof_summary_str(prof_file_M2), file.path(out_dir, "profile_M2.txt"))

## -- timings.csv -------------------------------------------------------------

timings_df <- rbind(
    data.frame(method = "M1",
               stage  = names(m1$timings),
               seconds = sapply(m1$timings, `[[`, "elapsed")),
    data.frame(method = "M2",
               stage  = names(m2$timings),
               seconds = sapply(m2$timings, `[[`, "elapsed"))
)
timings_df$seconds <- round(timings_df$seconds, 4)
write.csv(timings_df, file.path(out_dir, "timings.csv"), row.names = FALSE)
cat("\nTimings written to timings.csv\n")
print(timings_df)

## -- Accuracy / structural validation ---------------------------------------

cat("\n=== Accuracy + structural validation ===\n")
acc <- character()
acc_check <- function(label, expr) {
    ok <- tryCatch(isTRUE(eval(expr, envir = parent.frame())),
                    error = function(e) FALSE)
    line <- sprintf("[%s] %s", if (ok) "PASS" else "FAIL", label)
    cat(line, "\n")
    acc <<- c(acc, line)
    invisible(ok)
}

sa <- m1$sa
acc_check("SangerAlignment built (M1)",
           quote(sa@objectResults@creationResult))
acc_check("contigList non-empty (M1)",
           quote(length(sa@contigList) > 0L))
acc_check("contigsConsensus non-empty (M1)",
           quote(length(sa@contigsConsensus) > 0L))
acc_check("contigsAlignment is DNAStringSet (M1)",
           quote(is(sa@contigsAlignment, "DNAStringSet")))
acc_check("All child reads valid (M1)",
           quote({
               reads <- unlist(lapply(sa@contigList, function(sc)
                   c(sc@forwardReadList, sc@reverseReadList)),
                   recursive = FALSE)
               all(vapply(reads, function(r) r@objectResults@creationResult,
                          logical(1)))
           }))

sa2 <- m2$sa
acc_check("SangerAlignment built (M2)",
           quote(sa2@objectResults@creationResult))
acc_check("M1 vs M2: contig count equal",
           quote(length(sa@contigList) == length(sa2@contigList)))
acc_check("M1 vs M2: same contig names",
           quote(identical(sort(names(sa@contigList)),
                            sort(names(sa2@contigList)))))

## Run twice; consensus must be deterministic (same inputs / params).
cat("Determinism check: rebuilding SangerAlignment a second time ...\n")
sa_again <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = ab1_dir,
    REGEX_SuffixForward = REGEX_F,
    REGEX_SuffixReverse = REGEX_R,
    TrimmingMethod      = "M1",
    M1TrimmingCutoff    = 0.0001,
    processorsNum       = 1
)
acc_check("Consensus is deterministic across two runs",
           quote(identical(as.character(sa@contigsConsensus),
                            as.character(sa_again@contigsConsensus))))

writeLines(acc, file.path(out_dir, "accuracy.txt"))

## -- Phase-4 *Raw slot invariants -------------------------------------------

cat("\n=== Phase-4 SangerRead *Raw slot invariants ===\n")
slot_lines <- character()
chk <- function(label, expr) {
    ok <- isTRUE(eval(expr, envir = parent.frame()))
    line <- sprintf("[%s] %s", if (ok) "PASS" else "FAIL", label)
    cat(line, "\n")
    slot_lines <<- c(slot_lines, line)
    invisible(ok)
}

ab1 <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                   "Achl_ACHLO006-09_1_F.ab1", package = "sangeranalyseR")
sr  <- new("SangerRead",
            inputSource    = "ABIF",
            readFeature    = "Forward Read",
            readFileName   = ab1,
            TrimmingMethod = "M1")
ref <- sangerseqR::sangerseq(sangerseqR::read.abif(ab1))

chk("primarySeqRaw   == sangerseq parent primarySeq",
     quote(identical(as.character(sr@primarySeqRaw),
                      as.character(ref@primarySeq))))
chk("secondarySeqRaw == sangerseq parent secondarySeq",
     quote(identical(as.character(sr@secondarySeqRaw),
                      as.character(ref@secondarySeq))))
chk("peakPosMatrixRaw  == parent peakPosMatrix",
     quote(identical(sr@peakPosMatrixRaw, ref@peakPosMatrix)))
chk("peakAmpMatrixRaw  == parent peakAmpMatrix",
     quote(identical(sr@peakAmpMatrixRaw, ref@peakAmpMatrix)))
chk("traceMatrix (inherited) == parent traceMatrix",
     quote(identical(sr@traceMatrix, ref@traceMatrix)))
chk("primarySeq (post-MakeBaseCalls) length == qualityPhredScores length",
     quote(length(sr@primarySeq) == length(sr@QualityReport@qualityPhredScores)))

writeLines(slot_lines, file.path(out_dir, "slot_invariants.txt"))

cat("\n=== Phase 5 run complete. Artifacts in:", out_dir, "===\n")
