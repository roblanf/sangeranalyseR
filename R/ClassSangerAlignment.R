#' @export
setOldClass("phylo")

#' @title SangerAlignment
#'
#' @description  An S4 class containing SangerContigs lists and contigs alignment results which corresponds to a final alignment in Sanger sequencing.
#'
#' @slot parentDirectory The parent directory of all of the reads contained in ABIF format you wish to analyse. In SangerAlignment, all reads in subdirectories will be scanned recursively.
#' @slot suffixForwardRegExp The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be \code{"_[F]_[0-9]*.ab1"}.
#' @slot suffixReverseRegExp The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be \code{"_[R]_[0-9]*.ab1"}.
#' @slot trimmingMethodSA The read trimming method for all SangerRead S4 instances in SangerAlignment. The value must be \code{"M1"} (the default) or \code{'M2'}. All SangerReads must have the same trimming method.
#' @slot minFractionCallSA Minimum fraction of the sequences required to call a consensus sequence for SangerAlignment at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
#' @slot maxFractionLostSA Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerAlignment (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
#' @slot geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @slot refAminoAcidSeq An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is \code{""}.
#' @slot contigList A list storing all SangerContigs S4 instances.
#' @slot contigsConsensus The consensus read of all SangerContig S4 instances in DNAString object.
#' @slot contigsAlignment The alignment of all SangerContig S4 instances with the called consensus sequence in DNAStringSet object. Users can use BrowseSeqs() to view the alignment.
#' @slot contigsTree A phylo instance returned by bionj function in ape package. It can be used to draw the tree.
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
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' newAlignment <- new("SangerAlignment",
#'                     parentDirectory       = rawDataDir,
#'                     suffixForwardRegExp   = suffixForwardRegExp,
#'                     suffixReverseRegExp   = suffixReverseRegExp,
#'                     refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                     TrimmingMethod        = "M1",
#'                     M1TrimmingCutoff      = 0.0001,
#'                     M2CutoffQualityScore  = NULL,
#'                     M2SlidingWindowSize   = NULL,
#'                     baseNumPerRow         = 100,
#'                     heightPerRow          = 200,
#'                     signalRatioCutoff     = 0.33,
#'                     showTrimmed           = TRUE)
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
             trimmingMethodSA            = "character",
             minFractionCallSA           = "numeric",
             maxFractionLostSA           = "numeric",
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
                   minFractionCallSA      = 0.5,
                   maxFractionLostSA      = 0.5,
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
        trimmingMethodSA = TrimmingMethod
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
                                 minFractionCallSA, maxFractionLostSA,
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
                   trimmingMethodSA      = trimmingMethodSA,
                   contigList            = SangerContigList,
                   minFractionCallSA     = minFractionCallSA,
                   maxFractionLostSA     = maxFractionLostSA,
                   geneticCode           = geneticCode,
                   contigsConsensus      = consensus,
                   refAminoAcidSeq       = refAminoAcidSeq,
                   contigsAlignment      = aln,
                   contigsTree           = aln.tree
                   )
})

