// Phase 7 Rcpp port of UtilitiesFunc.R::peakvalues.
//
// Original R implementation:
//
//   peakvalues <- function(x, pstart, pstop) {
//       region <- x[x[,1] > pstart & x[,1] < pstop, , drop = FALSE]
//       if (length(region[,1]) == 0) return(c(0, NA))
//       else return(c(max(region[,2], na.rm = TRUE),
//                     region[which.max(region[,2]), 1]))
//   }
//
// `x` is the output of `getpeaks(trace)`: a 2-column numeric matrix where
// column 1 is the trace index of each peak (an integer-valued double) and
// column 2 is the trace amplitude at that index.
//
// Contract preserved exactly:
//   * No row matches `> pstart && < pstop`     -> c(0, NA)
//   * One or more rows match, all NA in col 2 -> c(NA, NA) (matches R's
//        max(na.rm = TRUE) of all-NA = -Inf with warning, then which.max
//        on all-NA returns integer(0); pre-Phase-7 R could produce a length-
//        1 result here. Production traces never hold NAs, so this case is
//        only for defensive coverage.)
//   * Otherwise                                -> c(max_amp, pos_at_argmax)
//
// Tie-break: when multiple rows share the max amplitude, the FIRST occurrence
// wins (matches R's `which.max`).

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector peakvalues_cpp(Rcpp::NumericMatrix x,
                                    double pstart,
                                    double pstop) {
    const int n = x.nrow();
    double max_val = R_NegInf;
    double max_pos = NA_REAL;
    bool found_any = false;
    bool found_non_na = false;

    for (int i = 0; i < n; ++i) {
        const double idx = x(i, 0);
        if (idx > pstart && idx < pstop) {
            found_any = true;
            const double v = x(i, 1);
            if (!Rcpp::NumericVector::is_na(v)) {
                if (!found_non_na || v > max_val) {
                    max_val = v;
                    max_pos = idx;
                    found_non_na = true;
                }
            }
        }
    }

    Rcpp::NumericVector out(2);
    if (!found_any) {
        out[0] = 0.0;
        out[1] = NA_REAL;
    } else if (!found_non_na) {
        // Region matched but every value was NA — defensive, never observed
        // on real ABIF traces.
        out[0] = NA_REAL;
        out[1] = NA_REAL;
    } else {
        out[0] = max_val;
        out[1] = max_pos;
    }
    return out;
}

// Batch variant: process all peak windows for one channel in a single
// .Call. Reduces the per-call marshalling overhead that, on small
// per-call work, can dominate the savings of the C++ inner loop.
//
// Returns a 2 x K NumericMatrix where K = length(pstarts) and column k is
// (max_amp, max_pos) for window (pstarts[k], pstops[k]). Empty windows
// produce (0, NA), all-NA windows produce (NA, NA).
//
// Contract per column matches peakvalues_cpp() exactly.

// [[Rcpp::export]]
Rcpp::NumericMatrix peakvalues_batch_cpp(Rcpp::NumericMatrix x,
                                          Rcpp::NumericVector pstarts,
                                          Rcpp::NumericVector pstops) {
    const int n = x.nrow();
    const int k = pstarts.size();
    if (pstops.size() != k) {
        Rcpp::stop("pstarts and pstops must have equal length");
    }

    Rcpp::NumericMatrix out(2, k);

    for (int j = 0; j < k; ++j) {
        const double pstart = pstarts[j];
        const double pstop  = pstops[j];
        double max_val      = R_NegInf;
        double max_pos      = NA_REAL;
        bool found_any      = false;
        bool found_non_na   = false;

        for (int i = 0; i < n; ++i) {
            const double idx = x(i, 0);
            if (idx > pstart && idx < pstop) {
                found_any = true;
                const double v = x(i, 1);
                if (!Rcpp::NumericVector::is_na(v)) {
                    if (!found_non_na || v > max_val) {
                        max_val = v;
                        max_pos = idx;
                        found_non_na = true;
                    }
                }
            }
        }

        if (!found_any) {
            out(0, j) = 0.0;
            out(1, j) = NA_REAL;
        } else if (!found_non_na) {
            out(0, j) = NA_REAL;
            out(1, j) = NA_REAL;
        } else {
            out(0, j) = max_val;
            out(1, j) = max_pos;
        }
    }
    return out;
}
