## Phase 6 — benchmark lazy AA + BiocParallel scaling.
##
## Compares:
##   (a) lazy=TRUE  + SerialParam   (the new default)
##   (b) lazy=FALSE + SerialParam   (eager AA — old behavior)
##   (c) lazy=TRUE  + MulticoreParam(workers = N) for N in 1, 2, 4
##
## All against the bundled 8-read ACHLO fixture. Output: timings.csv +
## profile_M1_lazy.txt.

suppressPackageStartupMessages({
    library(devtools)
    library(BiocParallel)
})

repo_root <- "/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR"
out_dir   <- file.path(repo_root, "plans", "phase6_artifacts")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(load_all(repo_root, quiet = TRUE))

ab1_dir <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                       package = "sangeranalyseR")
REGEX_F <- "_[0-9]*_F.ab1$"
REGEX_R <- "_[0-9]*_R.ab1$"

stage_time <- function(label, expr) {
    gc(verbose = FALSE)
    t0 <- Sys.time()
    val <- force(eval(expr, envir = parent.frame()))
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  [%-32s] %.3f s\n", label, elapsed))
    list(value = val, elapsed = elapsed)
}

run_once <- function(label, lazyAA, BPPARAM) {
    gc(verbose = FALSE)
    t0 <- Sys.time()
    sa <- SangerAlignment(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = ab1_dir,
        REGEX_SuffixForward = REGEX_F,
        REGEX_SuffixReverse = REGEX_R,
        TrimmingMethod      = "M1",
        M1TrimmingCutoff    = 0.0001,
        BPPARAM             = BPPARAM,
        lazyAA              = lazyAA
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  [%-40s] %.3f s\n", label, elapsed))
    list(elapsed = elapsed, sa = sa)
}

cat("\n=== Warm-up ===\n")
invisible(run_once("warm-up", TRUE, SerialParam()))

cat("\n=== Benchmark matrix ===\n")
res <- list()
for (rep in 1:3) {
    cat(sprintf("\n-- repetition %d --\n", rep))
    res[[length(res)+1L]] <- list(
        rep    = rep,
        config = "lazy=TRUE  + SerialParam",
        secs   = run_once("lazy=TRUE  + SerialParam", TRUE, SerialParam())$elapsed
    )
    res[[length(res)+1L]] <- list(
        rep    = rep,
        config = "lazy=FALSE + SerialParam",
        secs   = run_once("lazy=FALSE + SerialParam", FALSE, SerialParam())$elapsed
    )
    for (n in c(2, 4)) {
        bp <- if (.Platform$OS.type == "windows")
                  SnowParam(workers = n)
              else
                  MulticoreParam(workers = n)
        res[[length(res)+1L]] <- list(
            rep    = rep,
            config = sprintf("lazy=TRUE  + MulticoreParam(%d)", n),
            secs   = run_once(sprintf("lazy=TRUE + Multicore(%d)", n),
                                TRUE, bp)$elapsed
        )
    }
}

## Aggregate
df <- do.call(rbind, lapply(res, as.data.frame, stringsAsFactors = FALSE))
df$secs <- round(df$secs, 4)
write.csv(df, file.path(out_dir, "timings.csv"), row.names = FALSE)

cat("\n=== Aggregate (mean of 3 reps) ===\n")
agg <- aggregate(secs ~ config, data = df, FUN = mean)
agg$secs <- round(agg$secs, 4)
print(agg)
write.csv(agg, file.path(out_dir, "timings_summary.csv"), row.names = FALSE)

## Re-profile under the new lazy default
cat("\n=== Rprof: M1 + lazyAA=TRUE + SerialParam ===\n")
prof_file <- file.path(out_dir, "Rprof_M1_lazy.out")
Rprof(prof_file, interval = 0.005)
sa_p <- SangerAlignment(
    inputSource         = "ABIF",
    processMethod       = "REGEX",
    ABIF_Directory      = ab1_dir,
    REGEX_SuffixForward = REGEX_F,
    REGEX_SuffixReverse = REGEX_R,
    TrimmingMethod      = "M1",
    M1TrimmingCutoff    = 0.0001,
    BPPARAM             = SerialParam(),
    lazyAA              = TRUE
)
Rprof(NULL)

s <- summaryRprof(prof_file)
out <- character()
out <- c(out, "by.self (top 25):")
bs <- head(s$by.self, 25L)
for (col in c("self.time","total.time","self.pct","total.pct"))
    bs[[col]] <- sprintf("%6.3f", bs[[col]])
out <- c(out, capture.output(print(bs)))
out <- c(out, "", "by.total (top 15):")
bt <- head(s$by.total, 15L)
for (col in c("self.time","total.time","self.pct","total.pct"))
    bt[[col]] <- sprintf("%6.3f", bt[[col]])
out <- c(out, capture.output(print(bt)))
out <- c(out, "", sprintf("Total sampling time: %.3f s", s$sampling.time))
writeLines(out, file.path(out_dir, "profile_M1_lazy.txt"))

## Determinism: lazy AA accessor matches eager slot
cat("\n=== Lazy AA equivalence check ===\n")
ab1 <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                   "Achl_ACHLO006-09_1_F.ab1", package = "sangeranalyseR")
sr_lazy  <- new("SangerRead", inputSource="ABIF", readFeature="Forward Read",
                readFileName=ab1, lazyAA=TRUE)
sr_eager <- new("SangerRead", inputSource="ABIF", readFeature="Forward Read",
                readFileName=ab1, lazyAA=FALSE)
ok1 <- identical(as.character(primaryAASeqS1(sr_lazy)),
                  as.character(sr_eager@primaryAASeqS1))
ok2 <- identical(as.character(primaryAASeqS2(sr_lazy)),
                  as.character(sr_eager@primaryAASeqS2))
ok3 <- identical(as.character(primaryAASeqS3(sr_lazy)),
                  as.character(sr_eager@primaryAASeqS3))
ok_empty <- length(sr_lazy@primaryAASeqS1) == 0L &&
             length(sr_lazy@primaryAASeqS2) == 0L &&
             length(sr_lazy@primaryAASeqS3) == 0L
cat(sprintf("  lazy slot empty?       %s\n", ok_empty))
cat(sprintf("  S1 accessor == eager?  %s\n", ok1))
cat(sprintf("  S2 accessor == eager?  %s\n", ok2))
cat(sprintf("  S3 accessor == eager?  %s\n", ok3))

writeLines(c(
    sprintf("[%s] lazy slot is empty AAString",        if (ok_empty) "PASS" else "FAIL"),
    sprintf("[%s] primaryAASeqS1() == eager slot",     if (ok1)      "PASS" else "FAIL"),
    sprintf("[%s] primaryAASeqS2() == eager slot",     if (ok2)      "PASS" else "FAIL"),
    sprintf("[%s] primaryAASeqS3() == eager slot",     if (ok3)      "PASS" else "FAIL")
), file.path(out_dir, "lazy_aa_equivalence.txt"))

cat("\n=== Phase 6 bench complete. Artifacts in:", out_dir, "===\n")
