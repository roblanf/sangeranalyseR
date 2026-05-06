## Phase 7 — end-to-end SangerAlignment benchmark
##
## Compares Phase 6 baseline (R-only peakvalues) vs Phase 7 (C++ peakvalues)
## by running the full pipeline twice — once with peakvalues_cpp monkey-
## patched to the R helper, once with the native C++ version.

suppressPackageStartupMessages(library(devtools))
repo_root <- "/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR"
out_dir   <- file.path(repo_root, "plans", "phase7_artifacts")
suppressMessages(load_all(repo_root, quiet = TRUE))

ns <- asNamespace("sangeranalyseR")
ab1_dir <- system.file("extdata", "Allolobophora_chlorotica", "ACHLO",
                        package = "sangeranalyseR")

run_sa <- function() {
    SangerAlignment(
        inputSource         = "ABIF",
        processMethod       = "REGEX",
        ABIF_Directory      = ab1_dir,
        REGEX_SuffixForward = "_[0-9]*_F.ab1$",
        REGEX_SuffixReverse = "_[0-9]*_R.ab1$",
        TrimmingMethod      = "M1",
        M1TrimmingCutoff    = 0.0001,
        BPPARAM             = BiocParallel::SerialParam(),
        lazyAA              = TRUE
    )
}

time_n <- function(label, fn, reps = 5L) {
    invisible(fn())                            # warm-up
    tt <- numeric(reps)
    for (r in seq_len(reps)) {
        gc(verbose = FALSE)
        t0 <- Sys.time()
        invisible(fn())
        tt[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    cat(sprintf("[%-32s]  mean=%.3fs  min=%.3fs  reps=%d\n",
                  label, mean(tt), min(tt), reps))
    tt
}

cat("\n=== End-to-end SangerAlignment: C++ vs (monkey-patched) R ===\n")
cpp_times <- time_n("C++ peakvalues_batch (Phase 7)", run_sa)

# An R fallback that emulates the C++ batch API by looping over
# .peakvalues_r. This is the Phase-6 baseline behavior — what the package
# would do if the C++ function were unavailable.
peakvalues_batch_r <- function(x, pstarts, pstops) {
    out <- matrix(NA_real_, nrow = 2L, ncol = length(pstarts))
    for (j in seq_along(pstarts)) {
        out[, j] <- ns$.peakvalues_r(x, pstarts[j], pstops[j])
    }
    out
}

r_times <- testthat::with_mocked_bindings(
    peakvalues_batch_cpp = peakvalues_batch_r,
    .package             = "sangeranalyseR",
    code                 = time_n("R peakvalues_batch (Phase 6 baseline)", run_sa)
)

speedup <- mean(r_times) / mean(cpp_times)
cat(sprintf("\nMean speedup: %.2fx (saved %.0f ms per build)\n",
              speedup,
              (mean(r_times) - mean(cpp_times)) * 1000))

df <- data.frame(
    config  = c("phase6_R_peakvalues", "phase7_cpp_peakvalues"),
    mean_s  = c(round(mean(r_times),    4), round(mean(cpp_times),    4)),
    min_s   = c(round(min(r_times),     4), round(min(cpp_times),     4)),
    reps    = c(length(r_times), length(cpp_times))
)
write.csv(df, file.path(out_dir, "bench_e2e.csv"), row.names = FALSE)
cat("\nSaved bench_e2e.csv\n"); print(df)

## Re-profile to confirm peakvalues drops from the top of the profile.
cat("\n=== Rprof on Phase-7 SangerAlignment ===\n")
prof <- file.path(out_dir, "Rprof_phase7.out")
Rprof(prof, interval = 0.005)
sa <- run_sa()
Rprof(NULL)

s <- summaryRprof(prof)
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
writeLines(out, file.path(out_dir, "profile_phase7.txt"))

## Determinism check vs raw expected.
cat("\n=== Determinism: SangerAlignment built twice produces identical consensus ===\n")
sa2 <- run_sa()
cat(sprintf("Consensus identical?  %s\n",
              identical(as.character(sa@contigsConsensus),
                         as.character(sa2@contigsConsensus))))
