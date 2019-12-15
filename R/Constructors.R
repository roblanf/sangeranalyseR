### ============================================================================
### Self-defined constructor for AlignedConsensusSet
### ============================================================================
#' @description
#'
#' @param parentDirectory The parent directory of all of the reads contained in ABIF format you wish to analyse. In SangerAlignment, all reads in subdirectories will be scanned recursively.
#' @param suffixForwardRegExp The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be \code{"_[F]_[0-9]*.ab1"}.
#' @param suffixReverseRegExp The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be \code{"_[R]_[0-9]*.ab1"}.
#' @param TrimmingMethod TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#' @param baseNumPerRow  It defines maximum base pairs in each row. The default value is \code{100}.
#' @param heightPerRow It defines the height of each row in chromatogram. The default value is \code{200}.
#' @param signalRatioCutoff The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is \code{0.33}.
#' @param showTrimmed The logical value storing whether to show trimmed base pairs in chromatogram. The default value is \code{TRUE}.
#' @param refAminoAcidSeq An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is \code{""}.
#' @param minReadsNum The minimum number of reads required to make a consensus sequence, must be 2 or more. The default value is \code{2}.
#' @param minReadLength Reads shorter than this will not be included in the readset. The default \code{20} means that all reads with length of 20 or more will be included. Note that this is the length of a read after it has been trimmed.
#' @param minFractionCall Minimum fraction of the sequences required to call a consensus sequence for SangerContig at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
#' @param maxFractionLost Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerContig (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
#' @param geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @param acceptStopCodons The logical value \code{TRUE} or \code{FALSE}. \code{TRUE} (the defualt): keep all reads, regardless of whether they have stop codons; \code{FALSE}: reject reads with stop codons. If \code{FALSE} is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
#' @param readingFrame \code{1}, \code{2}, or \code{3}. Only used if \code{accept.stop.codons == FALSE}. This specifies the reading frame that is used to determine stop codons. If you use a \code{refAminoAcidSeq}, then the frame should always be \code{1}, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame.
#' @param processorsNum The number of processors to use, or NULL (the default) for all available processors.
#'
#' @return SangerAlignment
#' @export
#' @author Kuan-Hao Chao
#' @example
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' sangerAlignment <- SangerAlignment(
#'                        parentDirectory       = parentDir,
#'                        suffixForwardRegExp   = suffixForwardRegExp,
#'                        suffixReverseRegExp   = suffixReverseRegExp,
#'                        refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                        TrimmingMethod        = "M1",
#'                        M1TrimmingCutoff      = 0.0001,
#'                        M2CutoffQualityScore  = NULL,
#'                        M2SlidingWindowSize   = NULL,
#'                        baseNumPerRow         = 100,
#'                        heightPerRow          = 200,
#'                        signalRatioCutoff     = 0.33,
#'                        showTrimmed           = TRUE)
SangerAlignment <- function(parentDirectory        = character(0),
                            suffixForwardRegExp    = character(0),
                            suffixReverseRegExp    = character(0),
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
                            processorsNum          = NULL) {
    newAlignment <- new("SangerAlignment",
                        parentDirectory        = parentDirectory,
                        suffixForwardRegExp    = suffixForwardRegExp,
                        suffixReverseRegExp    = suffixReverseRegExp,
                        TrimmingMethod         = TrimmingMethod,
                        M1TrimmingCutoff       = M1TrimmingCutoff,
                        M2CutoffQualityScore   = M2CutoffQualityScore,
                        M2SlidingWindowSize    = M2SlidingWindowSize,
                        baseNumPerRow          = baseNumPerRow,
                        heightPerRow           = heightPerRow,
                        signalRatioCutoff      = signalRatioCutoff,
                        showTrimmed            = showTrimmed,
                        refAminoAcidSeq        = refAminoAcidSeq,
                        minReadsNum            = minReadsNum,
                        minReadLength          = minReadLength,
                        minFractionCall        = minFractionCall,
                        maxFractionLost        = maxFractionLost,
                        geneticCode            = geneticCode,
                        acceptStopCodons       = acceptStopCodons,
                        readingFrame           = readingFrame,
                        processorsNum          = processorsNum)
    return(newAlignment)
}



### ============================================================================
### Self-defined constructor for SangerContig
### ============================================================================
#' @description
#'
#' @param parentDirectory The parent directory of all of the reads contained in ABIF format you wish to analyse. In SangerContig, all reads must be in the first layer in this directory.
#' @param contigName The contig name of all the reads in \code{parentDirectory}.
#' @param suffixForwardRegExp The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be \code{"_[F]_[0-9]*.ab1"}.
#' @param suffixReverseRegExp The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be \code{"_[R]_[0-9]*.ab1"}.
#' @param TrimmingMethod TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#' @param baseNumPerRow  It defines maximum base pairs in each row. The default value is \code{100}.
#' @param heightPerRow It defines the height of each row in chromatogram. The default value is \code{200}.
#' @param signalRatioCutoff The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is \code{0.33}.
#' @param showTrimmed The logical value storing whether to show trimmed base pairs in chromatogram. The default value is \code{TRUE}.
#' @param refAminoAcidSeq An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is \code{""}.
#' @param minReadsNum The minimum number of reads required to make a consensus sequence, must be 2 or more. The default value is \code{2}.
#' @param minReadLength Reads shorter than this will not be included in the readset. The default \code{20} means that all reads with length of 20 or more will be included. Note that this is the length of a read after it has been trimmed.
#' @param minFractionCall Minimum fraction of the sequences required to call a consensus sequence for SangerContig at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
#' @param maxFractionLost Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerContig (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
#' @param geneticCode Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
#' @param acceptStopCodons The logical value \code{TRUE} or \code{FALSE}. \code{TRUE} (the defualt): keep all reads, regardless of whether they have stop codons; \code{FALSE}: reject reads with stop codons. If \code{FALSE} is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
#' @param readingFrame \code{1}, \code{2}, or \code{3}. Only used if \code{accept.stop.codons == FALSE}. This specifies the reading frame that is used to determine stop codons. If you use a \code{refAminoAcidSeq}, then the frame should always be \code{1}, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame.
#' @param processorsNum The number of processors to use, or NULL (the default) for all available processors.
#'
#' @return SangerContig
#' @export
#' @author Kuan-Hao Chao
#' @example
#' rawDataDir <- system.file("extdata", package = "sangeranalyseR")
#' parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "ACHLO")
#' contigName <- "ACHLO006-09[LCO1490_t1,HCO2198_t1]"
#' suffixForwardRegExp <- "_[F]_[0-9]*.ab1"
#' suffixReverseRegExp <- "_[R]_[0-9]*.ab1"
#' sangerContig <- SangerContig(
#'                      parentDirectory       = parentDir,
#'                      contigName            = contigName,
#'                      suffixForwardRegExp   = suffixForwardRegExp,
#'                      suffixReverseRegExp   = suffixReverseRegExp,
#'                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
#'                      TrimmingMethod        = "M2",
#'                      M1TrimmingCutoff      = NULL,
#'                      M2CutoffQualityScore  = 20,
#'                      M2SlidingWindowSize   = 10,
#'                      baseNumPerRow         = 100,
#'                      heightPerRow          = 200,
#'                      signalRatioCutoff     = 0.33,
#'                      showTrimmed           = TRUE)
SangerContig <- function(parentDirectory        = character(0),
                         contigName             = character(0),
                         suffixForwardRegExp    = character(0),
                         suffixReverseRegExp    = character(0),
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
                         processorsNum          = NULL) {
    newContig <- new("SangerContig",
                     parentDirectory        = parentDirectory,
                     contigName             = contigName,
                     suffixForwardRegExp    = suffixForwardRegExp,
                     suffixReverseRegExp    = suffixReverseRegExp,
                     TrimmingMethod         = TrimmingMethod,
                     M1TrimmingCutoff       = M1TrimmingCutoff,
                     M2CutoffQualityScore   = M2CutoffQualityScore,
                     M2SlidingWindowSize    = M2SlidingWindowSize,
                     baseNumPerRow          = baseNumPerRow,
                     heightPerRow           = heightPerRow,
                     signalRatioCutoff      = signalRatioCutoff,
                     showTrimmed            = showTrimmed,
                     refAminoAcidSeq        = refAminoAcidSeq,
                     minReadsNum            = minReadsNum,
                     minReadLength          = minReadLength,
                     minFractionCall        = minFractionCall,
                     maxFractionLost        = maxFractionLost,
                     geneticCode            = geneticCode,
                     acceptStopCodons       = acceptStopCodons,
                     readingFrame           = readingFrame,
                     processorsNum          = processorsNum)
    return(newContig)
}



### ============================================================================
### Self-defined constructor for SangerRead
### ============================================================================
#' @description
#'
#' @param readFeature The direction of the Sanger read. The value must be \code{"Forward Read"} or \code{"Reverse Read"}.
#' @param readFileName The filename of the target ABIF file.
#' @param TrimmingMethod TrimmingMethod The read trimming method for this SangerRead. The value must be \code{"M1"} (the default) or \code{'M2'}.
#' @param M1TrimmingCutoff The trimming cutoff for the Method 1. If \code{TrimmingMethod} is \code{"M1"}, then the default value is \code{0.0001}. Otherwise, the value must be \code{NULL}.
#' @param M2CutoffQualityScore The trimming cutoff quality score for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{20}. Otherwise, the value must be \code{NULL}. It works with \code{M2SlidingWindowSize}.
#' @param M2SlidingWindowSize The trimming sliding window size for the Method 2. If \code{TrimmingMethod} is \code{'M2'}, then the default value is \code{10}. Otherwise, the value must be \code{NULL}. It works with \code{M2CutoffQualityScore}.
#' @param baseNumPerRow  It defines maximum base pairs in each row. The default value is \code{100}.
#' @param heightPerRow It defines the height of each row in chromatogram. The default value is \code{200}.
#' @param signalRatioCutoff The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is \code{0.33}.
#' @param showTrimmed The logical value storing whether to show trimmed base pairs in chromatogram. The default value is \code{TRUE}.
#'
#' @return SangerRead
#' @export
#' @author Kuan-Hao Chao
#' @example
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaFdFN <- file.path(inputFilesPath,
#'                               "Allolobophora_chlorotica",
#'                               "ACHLO",
#'                               "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F_1.ab1")
#' sangerRead <- SangerRead(
#'                    readFeature           = "Forward Read",
#'                    readFileName          = A_chloroticaFdFN,
#'                    geneticCode           = GENETIC_CODE,
#'                    TrimmingMethod        = "M1",
#'                    M1TrimmingCutoff      = 0.0001,
#'                    M2CutoffQualityScore  = NULL,
#'                    M2SlidingWindowSize   = NULL,
#'                    baseNumPerRow         = 100,
#'                    heightPerRow          = 200,
#'                    signalRatioCutoff     = 0.33,
#'                    showTrimmed           = TRUE)
SangerRead <- function(readFeature           = character(0),
                       readFileName          = character(0),
                       geneticCode           = GENETIC_CODE,
                       TrimmingMethod        = "M1",
                       M1TrimmingCutoff      = 0.0001,
                       M2CutoffQualityScore  = NULL,
                       M2SlidingWindowSize   = NULL,
                       baseNumPerRow         = 100,
                       heightPerRow          = 200,
                       signalRatioCutoff     = 0.33,
                       showTrimmed           = TRUE) {
    newRead <- new("SangerRead",
                   readFeature          = readFeature,
                   readFileName         = readFileName,
                   geneticCode          = geneticCode,
                   TrimmingMethod       = TrimmingMethod,
                   M1TrimmingCutoff     = M1TrimmingCutoff,
                   M2CutoffQualityScore = M2CutoffQualityScore,
                   M2SlidingWindowSize  = M2SlidingWindowSize,
                   baseNumPerRow        = baseNumPerRow,
                   heightPerRow         = heightPerRow,
                   signalRatioCutoff    = signalRatioCutoff,
                   showTrimmed          = showTrimmed)
    return(newRead)
}
