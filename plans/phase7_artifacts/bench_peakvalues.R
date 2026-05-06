## Phase 7 — peakvalues microbenchmark
## .peakvalues_r (pure R) vs peakvalues_cpp (Rcpp), on input shapes
## representative of the ABIF traceMatrix.

suppressPackageStartupMessages(library(devtools))
repo_root <- "/Users/chaokuan-hao/Documents/R_packages/sangeranalyseR"
out_dir   <- file.path(repo_root, "plans", "phase7_artifacts")
suppressMessages(load_all(repo_root, quiet = TRUE))

ns <- asNamespace("sangeranalyseR")
peakvalues_r <- ns$.peakvalues_r

mk_peaks <- function(n, seed = 42L) {
    set.seed(seed)
    idxs <- sort(sample.int(n * 4L, n, replace = FALSE))
    vals <- runif(n, min = 0, max = 1000)
    cbind(as.numeric(idxs), as.numeric(vals))
}

bench_one <- function(label, n, calls_per_run, reps) {
    x <- mk_peaks(n)
    # Sample call windows that overlap the index range.
    set.seed(7L)
    starts <- sort(sample(1:(n*4 - 50), calls_per_run, replace = TRUE))
    stops  <- starts + sample(20:80, calls_per_run, replace = TRUE)

    bench_fn <- function(fn) {
        tt <- numeric(reps)
        for (r in seq_len(reps)) {
            t0 <- Sys.time()
            for (i in seq_along(starts)) fn(x, starts[i], stops[i])
            tt[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        }
        tt
    }
    rt <- bench_fn(peakvalues_r)
    ct <- bench_fn(peakvalues_cpp)

    cat(sprintf("%-35s  N=%5d  reps=%d  R: %.4f s  C++: %.4f s  speedup: %.1fx\n",
                 label, n, reps, mean(rt), mean(ct), mean(rt) / mean(ct)))
    data.frame(
        label   = label,
        n_peaks = n,
        calls   = calls_per_run,
        reps    = reps,
        R_mean  = round(mean(rt), 6),
        cpp_mean= round(mean(ct), 6),
        R_min   = round(min(rt),  6),
        cpp_min = round(min(ct),  6),
        speedup = round(mean(rt) / mean(ct), 2)
    )
}

cat("\n=== peakvalues microbenchmark ===\n")
results <- rbind(
    bench_one("small (typical .ab1 channel)",  n =  500L, calls_per_run =  500L, reps = 5L),
    bench_one("medium",                         n = 2000L, calls_per_run =  500L, reps = 5L),
    bench_one("large (worst-case)",             n = 5000L, calls_per_run = 1000L, reps = 5L)
)

write.csv(results,
           file.path(out_dir, "bench_peakvalues.csv"),
           row.names = FALSE)
cat("\nSaved to bench_peakvalues.csv\n")
print(results)
