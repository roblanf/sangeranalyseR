#' @name qualityReportData
#' @title QualityReport instance
#' @description A pre-built \code{QualityReport} S4 object derived from the
#'   bundled ACHLO ABIF fixture, suitable for vignette and example use
#'   without re-running the trimming pipeline.
#' @format A \code{\link{QualityReport-class}} S4 object containing
#'   per-base Phred scores, the trimmed start/finish positions, and the
#'   raw / trimmed mean / minimum quality scores.
#' @docType data
#' @keywords datasets
#' @author Kuan-Hao Chao
#' @usage data(qualityReportData)
NULL

#' @name sangerReadFData
#' @title SangerRead instance
#' @description A pre-built \code{SangerRead} S4 object for one forward
#'   ABIF read from the ACHLO fixture.
#' @format A \code{\link{SangerRead-class}} S4 object with populated
#'   \code{primarySeq}, \code{secondarySeq}, \code{traceMatrix}, and
#'   nested \code{QualityReport} / \code{ChromatogramParam}.
#' @docType data
#' @keywords datasets
#' @author Kuan-Hao Chao
#' @usage data(sangerReadFData)
NULL

#' @name sangerContigData
#' @title SangerContig instance
#' @description A pre-built \code{SangerContig} S4 object containing one
#'   forward + one reverse \code{SangerRead} from the ACHLO fixture, plus
#'   the assembled contig consensus.
#' @format A \code{\link{SangerContig-class}} S4 object.
#' @docType data
#' @keywords datasets
#' @author Kuan-Hao Chao
#' @usage data(sangerContigData)
NULL

#' @name sangerAlignmentData
#' @title SangerAlignment instance
#' @description A pre-built \code{SangerAlignment} S4 object aggregating
#'   four contigs from the ACHLO fixture.
#' @format A \code{\link{SangerAlignment-class}} S4 object containing
#'   \code{contigList}, \code{contigsAlignment}, \code{contigsConsensus},
#'   and a \code{contigsTree} \code{phylo} instance.
#' @docType data
#' @keywords datasets
#' @author Kuan-Hao Chao
#' @usage data(sangerAlignmentData)
NULL
