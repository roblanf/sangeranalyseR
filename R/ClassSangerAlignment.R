#' @export
setOldClass("phylo")

#' @title SangerAlignment
#'
#' @description  An S4 class for storing multiple single consensus reads to
#'
#' @slot parentDirectory .
#' @slot suffixForwardRegExp.
#' @slot suffixReverseRegExp .
#' @slot minFractionCallSCSet .
#' @slot maxFractionLostSCSet .
#' @slot geneticCode .
#' @slot refAminoAcidSeq .
#' @slot contigList .
#' @slot contigsConsensus .
#' @slot contigsAlignment .
#' @slot contigsTree .
#'
#' @name SangerAlignment-class
#'
#' @rdname SangerAlignment-class
#'
#' @exportClass SangerAlignment
#' @author Kuan-Hao Chao
#' @include ClassQualityReport.R ClassSangerRead.R
#' @examples
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' inputFilesParentDir <- file.path(rawDataDir, "Allolobophora_chlorotica")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' SangerAlignment <- new("SangerAlignment",
#'                                parentDirectory       = rawDataDir,
#'                                suffixForwardRegExp   = suffixForwardRegExp,
#'                                suffixReverseRegExp   = suffixReverseRegExp,
#'                                refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                                TrimmingMethod        = "M2",
#'                                M1TrimmingCutoff      = NULL,
#'                                M2CutoffQualityScore  = 40,
#'                                M2SlidingWindowSize   = 10,
#'                                baseNumPerRow         = 100,
#'                                heightPerRow          = 200,
#'                                signalRatioCutoff     = 0.33,
#'                                showTrimmed           = TRUE)
setClass("SangerAlignment",
         # Users need to name their ab1 files in a systematic way. Here is the
         # regulation:
         #  1. Naming: XXXXX_F[0-9]*.ab1 / XXXXX_R[0-9]*.ab1
         #        For reads in same contig, XXXXX must be same.
         #  2. Users can set
         ### -------------------------------------------------------------------
         ### Input type of each variable of 'SangerAlignment'
         ### -------------------------------------------------------------------
         representation(
             parentDirectory             = "character",
             suffixForwardRegExp         = "character",
             suffixReverseRegExp         = "character",
             minFractionCallSCSet        = "numeric",
             maxFractionLostSCSet        = "numeric",
             geneticCode                 = "character",
             refAminoAcidSeq             = "character",
             contigList                  = "list",
             contigsConsensus            = "DNAString",
             contigsAlignment            = "DNAStringSet",
             contigsTree                 = "phylo"
             ),
)

### ============================================================================
### Overwrite initialize for 'SangerContig' (New constructor)
### ============================================================================
setMethod("initialize",
          "SangerAlignment",
          function(.Object, ...,
                   parentDirectory        = parentDirectory,
                   suffixForwardRegExp    = "_[F]_[0-9]*.ab1",
                   suffixReverseRegExp    = "_[R]_[0-9]*.ab1",
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL,
                   baseNumPerRow          = 100,
                   heightPerRow           = 200,
                   signalRatioCutoff      = 0.33,
                   showTrimmed            = TRUE,
                   refAminoAcidSeq        = "",
                   minReadsNum            = 2,
                   minReadLength          = 20,
                   minFractionCall        = 0.5,
                   maxFractionLost        = 0.5,
                   geneticCode            = GENETIC_CODE,
                   acceptStopCodons       = TRUE,
                   readingFrame           = 1,
                   minFractionCallSCSet   = 0.5,
                   maxFractionLostSCSet   = 0.5,
                   processorsNum          = 1) {
    ### ------------------------------------------------------------------------
    ### Input parameter prechecking
    ### ------------------------------------------------------------------------
    errors <- character()

    ### ------------------------------------------------------------------------
    ### 'parentDirectory' prechecking
    ### ------------------------------------------------------------------------
    errors <- checkParentDirectory (parentDirectory, errors)

    ### ------------------------------------------------------------------------
    ### Input parameter prechecking for TrimmingMethod.
    ### ------------------------------------------------------------------------
    errors <- checkTrimParam(TrimmingMethod, M1TrimmingCutoff,
                             M2CutoffQualityScore, M2SlidingWindowSize, errors)

    ##### ----------------------------------------------------------------------
    ##### Input parameter prechecking for ChromatogramParam
    ##### ----------------------------------------------------------------------
    errors <- checkBaseNumPerRow (baseNumPerRow, errors)
    errors <- checkHeightPerRow (baseNumPerRow, errors)
    errors <- checkSignalRatioCutoff (signalRatioCutoff, errors)
    errors <- checkShowTrimmed (showTrimmed, errors)

    ##### ----------------------------------------------------------------------
    ##### Input parameter prechecking for SangerContig parameter
    ##### ----------------------------------------------------------------------
    errors <- checkMinReadsNum(minReadsNum, errors)
    errors <- checkMinReadLength(minReadLength, errors)
    errors <- checkMinFractionCall(minFractionCall, errors)
    errors <- checkMaxFractionLost(maxFractionLost, errors)
    errors <- checkGeneticCode(geneticCode, errors)
    errors <- checkAcceptStopCodons(acceptStopCodons, errors)
    errors <- checkReadingFrame(readingFrame, errors)

    ##### ----------------------------------------------------------------------
    ##### Input parameter prechecking for processorsNum
    ##### ----------------------------------------------------------------------
    errors <- checkProcessorsNum(processorsNum, errors)

    if (length(errors) == 0) {
        ### --------------------------------------------------------------------
        ### Automatically finding contig name by forward&reverse suffix
        ### --------------------------------------------------------------------
        parentDirFiles <- list.files(parentDirectory, recursive = TRUE)
        forwardSelectInputFiles <- parentDirFiles[grepl(suffixForwardRegExp,
                                                        parentDirFiles)]
        reverseSelectInputFiles <- parentDirFiles[grepl(suffixReverseRegExp,
                                                        parentDirFiles)]

        # Find possible consensus Name for forward and reverse reads
        forwardContigName <-
            unlist(str_split(forwardSelectInputFiles, suffixForwardRegExp,
                             n = Inf, simplify = FALSE))[c(TRUE, FALSE)]
        reverseContigName <-
            unlist(str_split(reverseSelectInputFiles, suffixReverseRegExp,
                             n = Inf, simplify = FALSE))[c(TRUE, FALSE)]

        contigName <- union(forwardContigName, reverseContigName)
        contigNumber <- length(contigName)

        # Create contig for all list of contigNumber
        ### --------------------------------------------------------------------
        ### Creating each SangerContig (store as SangerContigList)
        ### --------------------------------------------------------------------
        SangerContigList <-
            sapply(contigName,
                   function(eachConsRead) {
                       insideDirName<- dirname(eachConsRead)
                       insideContigName <- basename(eachConsRead)
                       SangerContig(
                           file.path(parentDirectory, insideDirName),
                           insideContigName, suffixForwardRegExp,
                           suffixReverseRegExp,TrimmingMethod, M1TrimmingCutoff,
                           M2CutoffQualityScore, M2SlidingWindowSize,
                           baseNumPerRow, heightPerRow, signalRatioCutoff,
                           showTrimmed, refAminoAcidSeq, minReadsNum,
                           minReadLength, minFractionCall, maxFractionLost,
                           geneticCode, acceptStopCodons,
                           readingFrame, processorsNum)
                   })

        acResult <- alignContigs(SangerContigList, geneticCode, refAminoAcidSeq,
                                 minFractionCallSCSet, maxFractionLostSCSet,
                                 processorsNum)
        consensus <- acResult[["consensus"]]
        aln <- acResult[["aln"]]
        aln.tree <- acResult[["aln.tree"]]
    } else {
        stop(errors)
    }
    callNextMethod(.Object, ...,
                   parentDirectory       = parentDirectory,
                   suffixForwardRegExp   = suffixForwardRegExp,
                   suffixReverseRegExp   = suffixReverseRegExp,
                   contigList            = SangerContigList,
                   minFractionCallSCSet  = minFractionCallSCSet,
                   maxFractionLostSCSet  = maxFractionLostSCSet,
                   geneticCode           = geneticCode,
                   contigsConsensus      = consensus,
                   refAminoAcidSeq       = refAminoAcidSeq,
                   contigsAlignment      = aln,
                   contigsTree           = aln.tree
                   )
})



