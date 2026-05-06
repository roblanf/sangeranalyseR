# =============================================================================
# Synthetic-fixture builders for Phase 3 negative tests.
#
# Every builder writes to a path under tempdir() and registers cleanup with
# withr::defer so test runs leave the user's filesystem untouched. Real fixture
# files in inst/extdata/ are never modified.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")

# Path to the bundled real fixtures we copy / mutate.
.real_ab1_forward <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                               "ACHLO", "Achl_ACHLO006-09_1_F.ab1")
.real_ab1_reverse <- file.path(inputFilesPath, "Allolobophora_chlorotica",
                               "ACHLO", "Achl_ACHLO006-09_2_R.ab1")
.real_fasta_alignment <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                                   "Sanger_all_reads.fa")
.real_csv_alignment <- file.path(inputFilesPath, "ab1", "SangerAlignment",
                                 "names_conversion.csv")
.real_fasta_csv <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                             "names_conversion.csv")

# ---- FASTA / extension fixtures ---------------------------------------------

# Copy the bundled FASTA to a path with a `.Xfa` extension. The current regex
# in checkFASTA_File ('".fa$"') is buggy — the unescaped '.' matches any
# character, so 'Xfa' currently passes the type check.
make_xfa_file <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "Sanger_all_reads.Xfa")
    file.copy(.real_fasta_alignment, out)
    out
}

# Rename a real FASTA file to have a `.fast` extension.
make_fast_file <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "Sanger_all_reads.fast")
    file.copy(.real_fasta_alignment, out)
    out
}

# ---- ABIF fixtures ----------------------------------------------------------

# Rename a real .ab1 to have a `.ab2` extension.
make_ab2_file <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "Achl_ACHLO006-09_1_F.ab2")
    file.copy(.real_ab1_forward, out)
    out
}

# Zero-byte .ab1 file.
make_empty_ab1 <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "empty_1_F.ab1")
    file.create(out)
    out
}

# .ab1 file containing random bytes only.
make_corrupt_ab1 <- function(nbytes = 64L, env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "corrupt_1_F.ab1")
    set.seed(42L)
    bytes <- as.raw(sample.int(256L, size = nbytes, replace = TRUE) - 1L)
    writeBin(bytes, out)
    out
}

# Truncate a real .ab1 to `nbytes` (preserves the magic header but cuts the
# body). 32 bytes is enough to keep ABIF magic ("ABIF") and version intact.
make_corrupt_ab1_truncated <- function(nbytes = 32L, env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "truncated_1_F.ab1")
    full <- readBin(.real_ab1_forward, what = "raw",
                    n = file.info(.real_ab1_forward)$size)
    writeBin(full[seq_len(min(nbytes, length(full)))], out)
    out
}

# ---- CSV fixtures -----------------------------------------------------------

.read_real_csv <- function() utils::read.csv(.real_csv_alignment, header = TRUE,
                                             stringsAsFactors = FALSE)

# CSV where some `reads` rows reference filenames that do not exist on disk.
make_mismatched_csv <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "mismatched_names_conversion.csv")
    df <- .read_real_csv()
    df$reads[1] <- "DOES_NOT_EXIST.ab1"
    utils::write.csv(df, out, row.names = FALSE)
    out
}

# CSV missing one of {reads, direction, contig}.
make_csv_missing_column <- function(col, env = parent.frame()) {
    stopifnot(col %in% c("reads", "direction", "contig"))
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, sprintf("missing_%s.csv", col))
    df <- .read_real_csv()
    df[[col]] <- NULL
    utils::write.csv(df, out, row.names = FALSE)
    out
}

# CSV whose `direction` column has invalid values (not F/R).
make_csv_bad_direction <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out <- file.path(dir, "bad_direction.csv")
    df <- .read_real_csv()
    df$direction <- rep("X", nrow(df))
    utils::write.csv(df, out, row.names = FALSE)
    out
}

# ---- Directory layouts ------------------------------------------------------

# Build a tempdir containing two real .ab1s at the top level, a `.hidden` ab1,
# and a nested sub-directory with two more real .ab1s. Used by IO-batching
# tests to verify recursive walking and hidden-file behaviour.
make_dir_with_hidden <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    file.copy(.real_ab1_forward, file.path(dir, "Top_1_F.ab1"))
    file.copy(.real_ab1_reverse, file.path(dir, "Top_2_R.ab1"))
    file.copy(.real_ab1_forward, file.path(dir, ".hidden_1_F.ab1"))
    nested <- file.path(dir, "nested")
    dir.create(nested)
    file.copy(.real_ab1_forward, file.path(nested, "Nested_1_F.ab1"))
    file.copy(.real_ab1_reverse, file.path(nested, "Nested_2_R.ab1"))
    dir
}

# ---- FASTA-record-not-present fixtures --------------------------------------

# A FASTA with one record removed compared to the bundled CSV mapping.
make_fasta_missing_record <- function(env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    out_fa  <- file.path(dir, "Sanger_missing_record.fa")
    raw <- readLines(.real_fasta_alignment, warn = FALSE)
    # Drop the first record (header + sequence lines until the next '>').
    header_idxs <- grep("^>", raw)
    drop_start  <- header_idxs[1]
    drop_end    <- if (length(header_idxs) >= 2) header_idxs[2] - 1L else length(raw)
    keep <- raw[-seq.int(drop_start, drop_end)]
    writeLines(keep, out_fa)
    list(fasta = out_fa, csv = .real_fasta_csv)
}

# ---- Helper assertion -------------------------------------------------------

# Convenience for asserting that a row with `errorType == X` exists in a
# Sanger* object's @objectResults@readResultTable.
expect_table_has_error_type <- function(table, type) {
    testthat::expect_true(
        nrow(table) > 0L && type %in% as.character(table$errorType),
        info = sprintf("expected errorType '%s' in readResultTable; got [%s]",
                       type, paste(unique(as.character(table$errorType)),
                                   collapse = ", "))
    )
}
