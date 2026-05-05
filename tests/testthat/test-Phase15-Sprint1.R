# =============================================================================
# Phase 15 — Sprint 1 issue resolutions.
#
# Targeted tests for:
#   #100  SangerAlignment CSV/ABIF returns CONTIG_NUMBER_ZERO_ERROR even
#         though every individual SangerContig succeeds.
#   #92   Forward-reads-only datasets rejected with "REGEX_SuffixReverse
#         must be character type".
#   #76   `qualityPhredScores length cannot be zero` on certain ABIFs that
#         have no PCON.2 quality block.
#   #89   writeFasta() errors on contigs with only one read because
#         alignment slot is empty.
# =============================================================================

inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
ab1_dir <- file.path(inputFilesPath, "Allolobophora_chlorotica", "ACHLO")
fa_file <- file.path(inputFilesPath, "fasta", "SangerAlignment",
                     "Sanger_all_reads.fa")

# -----------------------------------------------------------------------------
# Issue #100: CSV labels that don't appear in filenames must still match
# -----------------------------------------------------------------------------
test_that("Issue #100: SangerAlignment CSV+ABIF works when contig labels are non-substring of filenames", {
    # Build a fixture: copy real ab1 files into a tempdir and write a CSV
    # whose `contig` labels do NOT appear in the filenames as substrings.
    tmp <- withr::local_tempdir()
    real_files <- list.files(ab1_dir, pattern = "\\.ab1$", full.names = TRUE)
    file.copy(real_files, tmp)

    # Two arbitrary contig labels ("rbcL2", "ITS5") — neither is a substring
    # of any "Achl_ACHLO*_F.ab1" / "_R.ab1" filename.
    csv_path <- file.path(tmp, "names_conversion.csv")
    fwd <- list.files(tmp, pattern = "_F\\.ab1$")
    rev <- list.files(tmp, pattern = "_R\\.ab1$")
    df <- data.frame(
        reads     = c(fwd[1:2], rev[1:2], fwd[3:4], rev[3:4]),
        direction = c("F","F","R","R","F","F","R","R"),
        contig    = c("rbcL2","rbcL2","rbcL2","rbcL2",
                      "ITS5","ITS5","ITS5","ITS5"),
        stringsAsFactors = FALSE
    )
    write.csv(df, csv_path, row.names = FALSE)

    sa <- new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = tmp,
              CSV_NamesConversion = csv_path,
              processorsNum       = 1)
    expect_true(sa@objectResults@creationResult)
    expect_equal(length(sa@contigList), 2L)
    expect_setequal(names(sa@contigList), c("rbcL2", "ITS5"))
    # Each child contig must have populated forward + reverse read lists
    for (sc in sa@contigList) {
        expect_gt(length(sc@forwardReadList), 0L)
        expect_gt(length(sc@reverseReadList), 0L)
    }
    # And the cross-contig consensus is non-empty
    expect_gt(length(sa@contigsConsensus), 0L)
})

test_that("Issue #100: SangerAlignment CSV+ABIF logs warning when a contig label has no matching reads", {
    # If the user's CSV references reads that don't exist on disk for a
    # specific contig label, the new code logs a warning and skips that
    # contig rather than aborting the whole alignment.
    tmp <- withr::local_tempdir()
    real_files <- list.files(ab1_dir, pattern = "\\.ab1$", full.names = TRUE)
    file.copy(real_files, tmp)

    csv_path <- file.path(tmp, "names_conversion.csv")
    fwd <- list.files(tmp, pattern = "_F\\.ab1$")
    rev <- list.files(tmp, pattern = "_R\\.ab1$")
    df <- data.frame(
        reads     = c(fwd[1:2], rev[1:2],
                      "NONEXISTENT_F.ab1", "NONEXISTENT_R.ab1"),
        direction = c("F","F","R","R","F","R"),
        contig    = c("good","good","good","good","ghost","ghost"),
        stringsAsFactors = FALSE
    )
    write.csv(df, csv_path, row.names = FALSE)

    sa <- suppressWarnings(new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "CSV",
              ABIF_Directory      = tmp,
              CSV_NamesConversion = csv_path,
              processorsNum       = 1))
    # `good` contig builds; `ghost` is silently dropped (no matching reads).
    expect_s4_class(sa, "SangerAlignment")
    expect_true(sa@objectResults@creationResult)
    expect_equal(length(sa@contigList), 1L)
    expect_equal(names(sa@contigList), "good")
})

# -----------------------------------------------------------------------------
# Issue #92: forward-reads-only / reverse-reads-only datasets
# -----------------------------------------------------------------------------
test_that("Issue #92: SangerAlignment accepts NULL REGEX_SuffixReverse for forward-only datasets", {
    tmp <- withr::local_tempdir()
    fwd_files <- list.files(ab1_dir, pattern = "_F\\.ab1$", full.names = TRUE)
    file.copy(fwd_files, tmp)

    sa <- suppressWarnings(new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = tmp,
              REGEX_SuffixForward = "_[0-9]*_F\\.ab1$",
              REGEX_SuffixReverse = NULL,        # explicit forward-only
              minReadsNum         = 1,           # forward-only contigs are 1-read each
              processorsNum       = 1))
    expect_s4_class(sa, "SangerAlignment")
    expect_true(sa@objectResults@creationResult)
    expect_gt(length(sa@contigList), 0L)
    # No reverse reads picked up
    for (sc in sa@contigList) {
        expect_equal(length(sc@reverseReadList), 0L)
        expect_gt(length(sc@forwardReadList), 0L)
    }
})

test_that("Issue #92: SangerAlignment accepts NA REGEX_SuffixReverse", {
    tmp <- withr::local_tempdir()
    fwd_files <- list.files(ab1_dir, pattern = "_F\\.ab1$", full.names = TRUE)
    file.copy(fwd_files, tmp)

    sa <- suppressWarnings(new("SangerAlignment",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = tmp,
              REGEX_SuffixForward = "_[0-9]*_F\\.ab1$",
              REGEX_SuffixReverse = NA_character_,
              minReadsNum         = 1,
              processorsNum       = 1))
    expect_true(sa@objectResults@creationResult)
})

test_that("Issue #92: SangerContig accepts NULL REGEX_SuffixReverse", {
    tmp <- withr::local_tempdir()
    fwd_files <- list.files(ab1_dir, pattern = "_F\\.ab1$", full.names = TRUE)
    file.copy(fwd_files, tmp)

    sc <- suppressWarnings(new("SangerContig",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = tmp,
              contigName          = "Achl_ACHLO006-09",
              REGEX_SuffixForward = "_[0-9]*_F\\.ab1$",
              REGEX_SuffixReverse = NULL,
              minReadsNum         = 1,
              processorsNum       = 1))
    expect_s4_class(sc, "SangerContig")
    expect_equal(length(sc@reverseReadList), 0L)
})

test_that("Issue #92: validators accept NULL/NA without flagging PARAMETER_VALUE_ERROR", {
    out_null <- checkREGEX_SuffixReverse(NULL,
                                          character(0), character(0))
    expect_length(out_null[[1]], 0L)

    out_na <- checkREGEX_SuffixReverse(NA_character_,
                                        character(0), character(0))
    expect_length(out_na[[1]], 0L)

    # And still rejects the truly-bad case
    out_bad <- checkREGEX_SuffixReverse(123L,
                                         character(0), character(0))
    expect_equal(out_bad[[2]], "PARAMETER_VALUE_ERROR")
})

# -----------------------------------------------------------------------------
# Issue #76: missing/empty PCON.2 quality block
# -----------------------------------------------------------------------------
test_that("Issue #76: SangerRead synthesises Phred 30 when PCON.2 is missing", {
    # Build a synthetic ABIF object whose @data$PCON.2 is empty by
    # monkey-patching read.abif locally with mocked bindings.
    real_file <- file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1")
    real_abif <- sangerseqR::read.abif(real_file)
    real_abif@data$PCON.2 <- raw(0)   # simulate missing quality block

    ## Mock the sangeranalyseR-side imported binding (not sangerseqR's own).
    ## SangerRead.initialize calls `read.abif(...)` via the namespace import,
    ## so the mock has to target the importing package's namespace.
    sr <- testthat::with_mocked_bindings(
        read.abif = function(filename) real_abif,
        .package  = "sangeranalyseR",
        code      = new("SangerRead",
                         inputSource    = "ABIF",
                         readFeature    = "Forward Read",
                         readFileName   = real_file,
                         TrimmingMethod = "M1")
    )
    expect_s4_class(sr, "SangerRead")
    expect_true(sr@objectResults@creationResult)
    # All synthesised quality scores should be 30
    qphred <- sr@QualityReport@qualityPhredScores
    expect_gt(length(qphred), 0L)
    expect_true(all(qphred == 30L))
})

test_that("Issue #76: SangerRead with intact PCON.2 still uses real quality scores", {
    real_file <- file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1")
    sr <- new("SangerRead",
              inputSource    = "ABIF",
              readFeature    = "Forward Read",
              readFileName   = real_file,
              TrimmingMethod = "M1")
    qphred <- sr@QualityReport@qualityPhredScores
    expect_gt(length(qphred), 0L)
    # Real data: not all scores are 30
    expect_false(all(qphred == 30L))
})

# -----------------------------------------------------------------------------
# Issue #89: writeFastaSC must not crash on single-read contigs
# -----------------------------------------------------------------------------
test_that("Issue #89: writeFastaSC succeeds on a single-read contig", {
    # Build a SangerContig from a single forward read by using a regex
    # that matches only that one file. The contig builder accepts
    # minReadsNum = 1 (lowering the floor) — it then skips alignment but
    # still populates @contigSeq from the lone read.
    tmp <- withr::local_tempdir()
    file.copy(file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1"), tmp)

    sc <- suppressWarnings(new("SangerContig",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = tmp,
              contigName          = "Achl_ACHLO006-09",
              REGEX_SuffixForward = "_[0-9]*_F\\.ab1$",
              REGEX_SuffixReverse = NULL,            # forward-only
              minReadsNum         = 1,
              processorsNum       = 1))
    expect_true(sc@objectResults@creationResult)
    expect_equal(length(sc@forwardReadList), 1L)
    expect_equal(length(sc@reverseReadList), 0L)
    expect_equal(length(sc@alignment), 0L)        # no alignment for n=1
    expect_gt(length(sc@contigSeq), 0L)            # but contigSeq is set

    # writeFastaSC must NOT crash; it should write a contig.fa with one
    # record.
    out_dir <- withr::local_tempdir()
    expect_error(writeFastaSC(sc, outputDir = out_dir), NA)

    written <- list.files(out_dir, pattern = "\\.fa$", full.names = TRUE)
    expect_gt(length(written), 0L)
    # The reads_alignment.fa file should exist and contain at least the
    # contig sequence record.
    aln_file <- written[grep("reads_alignment", written)]
    expect_equal(length(aln_file), 1L)
    aln_seqs <- Biostrings::readDNAStringSet(aln_file)
    expect_gte(length(aln_seqs), 1L)
})

test_that("Issue #89: writeFasta dispatcher works on a single-read contig", {
    tmp <- withr::local_tempdir()
    file.copy(file.path(ab1_dir, "Achl_ACHLO006-09_1_F.ab1"), tmp)
    sc <- suppressWarnings(new("SangerContig",
              inputSource         = "ABIF",
              processMethod       = "REGEX",
              ABIF_Directory      = tmp,
              contigName          = "Achl_ACHLO006-09",
              REGEX_SuffixForward = "_[0-9]*_F\\.ab1$",
              REGEX_SuffixReverse = NULL,
              minReadsNum         = 1,
              processorsNum       = 1))

    out_dir <- withr::local_tempdir()
    expect_error(writeFasta(sc, outputDir = out_dir), NA)
    expect_gt(length(list.files(out_dir, pattern = "\\.fa$")), 0L)
})
