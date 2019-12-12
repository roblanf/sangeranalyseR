# # A_chloroticaRead
# #
# # traceMatrix = A_chloroticaRead@traceMatrix
# # peakPosMatrixRaw = A_chloroticaRead@peakPosMatrixRaw
# # qualityPhredScoresRaw = A_chloroticaRead@QualityReport@qualityPhredScoresRaw
# # signalRatioCutoff = A_chloroticaRead@ChromatogramParam@signalRatioCutoff
# # readFeature = A_chloroticaRead@readFeature
# #
# #
# # MakeBaseCallsInside <- function(traceMatrix, peakPosMatrixRaw,
# #                                 qualityPhredScoresRaw,
# #                                 signalRatioCutoff, readFeature) {
# #     message("     * Making basecall !!")
# #     #get peaks for each base
# #     Apeaks <- getpeaks(traceMatrix[,1])
# #     Cpeaks <- getpeaks(traceMatrix[,2])
# #     Gpeaks <- getpeaks(traceMatrix[,3])
# #     Tpeaks <- getpeaks(traceMatrix[,4])
# #
# #     #get window around primary basecall peaks
# #     primarypeaks <- peakPosMatrixRaw[,1]
# #     diffs <- diff(c(0,primarypeaks))
# #     starts <- primarypeaks - 0.5*diffs
# #     stops <- c(primarypeaks[1:(length(primarypeaks)-1)] +
# #                    0.5*diffs[2:length(diffs)],
# #                primarypeaks[length(diffs)] + 0.5*diffs[length(diffs)]
# #     )
# #     #hack for last peak. Just uses distance preceding peak
# #     #as distance after peak
# #
# #     #Now get max peak value for each channel in each peak window.
# #     #If no peak return 0
# #     primary <- NULL
# #     secondary <- NULL
# #     tempPosMatrix <- matrix(nrow=length(starts), ncol=4)
# #     tempAmpMatrix <- matrix(nrow=length(starts), ncol=4)
# #     indexBaseCall <- c()
# #     for(i in 1:length(starts)) {
# #         Apeak <- peakvalues(Apeaks, starts[i], stops[i])
# #         Cpeak <- peakvalues(Cpeaks, starts[i], stops[i])
# #         Gpeak <- peakvalues(Gpeaks, starts[i], stops[i])
# #         Tpeak <- peakvalues(Tpeaks, starts[i], stops[i])
# #         if(is.na(Apeak[2]) &
# #            is.na(Cpeak[2]) &
# #            is.na(Gpeak[2]) &
# #            is.na(Tpeak[2])) {
# #             next #rare case where no peak found
# #         }
# #         ### ----------------------------------------------------------
# #         ### My modification here: Tracking BaseCall index
# #         ###     Add qualtiy score when making basecall
# #         ### ----------------------------------------------------------
# #         indexBaseCall <- c(indexBaseCall, i)
# #         signals <- c(Apeak[1], Cpeak[1], Gpeak[1], Tpeak[1])
# #
# #         tempAmpMatrix[i,] <- signals
# #         # print(tempAmpMatrix[i,])
# #
# #         positions <- c(Apeak[2], Cpeak[2], Gpeak[2], Tpeak[2])
# #
# #         tempPosMatrix[i,] <- positions
# #         # print(tempPosMatrix[i,])
# #
# #
# #         signalratios <- signals/max(signals, na.rm=TRUE)
# #         Bases <- c("A", "C", "G", "T")
# #         Bases[signalratios < signalRatioCutoff] <- NA
# #         #sort by decreasing signal strength
# #         Bases <- Bases[order(signals, decreasing=TRUE)]
# #         positions <- positions[order(signals, decreasing=TRUE)]
# #         if(length(Bases[!is.na(Bases)]) == 4
# #            | length(Bases[!is.na(Bases)]) == 0) {
# #             primary <- c(primary, "N")
# #             secondary <- c(secondary, "N")
# #         } else if(length(Bases[!is.na(Bases)]) > 1) {
# #             primary <- c(primary, Bases[1])
# #             Bases2 <- Bases[2:4]
# #             secondaryLetter <-
# #                 mergeIUPACLetters(paste(sort(Bases2[!is.na(Bases2)]),collapse=""))
# #             secondary <- c(secondary, secondaryLetter)
# #         }
# #         else {
# #             primary <- c(primary, Bases[1])
# #             secondary <- c(secondary, Bases[1])
# #         }
# #     }
# #     if (readFeature == "Forward Read") {
# #         qualityPhredScores <- qualityPhredScoresRaw[indexBaseCall]
# #         primarySeq <- DNAString(paste(primary, collapse=""))
# #         secondarySeq <- DNAString(paste(secondary, collapse=""))
# #     } else if (readFeature == "Reverse Read") {
# #         qualityPhredScores <- rev(qualityPhredScoresRaw[indexBaseCall])
# #         primarySeq <- reverseComplement(DNAString(paste(primary, collapse="")))
# #         secondarySeq <- reverseComplement(DNAString(paste(secondary, collapse="")))
# #     }
# #     peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
# #     peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
# #     message("     * Updating slots in 'SangerRead' instance !!")
# #     return(list("qualityPhredScores" = qualityPhredScores,
# #                 "peakPosMatrix" = peakPosMatrix,
# #                 "peakAmpMatrix" = peakAmpMatrix,
# #                 "primarySeq" = primarySeq,
# #                 "secondarySeq" = secondarySeq))
# # }
#
#
#
# IUPAC_CODE_MAP <- c(
#     A="A",
#     C="C",
#     G="G",
#     T="T",
#     M="AC",
#     R="AG",
#     W="AT",
#     S="CG",
#     Y="CT",
#     K="GT",
#     V="ACG",
#     H="ACT",
#     D="AGT",
#     B="CGT",
#     N="ACGT"
# )
#
#
# yy <- unname(IUPAC_CODE_MAP[unlist("S", use.names=FALSE)])
#
#
#
# IUPAC_CODE_MAP <- c(
#     A="A",
#     C="C",
#     G="G",
#     T="T",
#     M="AC",
#     R="AG",
#     W="AT",
#     S="CG",
#     Y="CT",
#     K="GT",
#     V="ACG",
#     H="ACT",
#     D="AGT",
#     B="CGT",
#     N="ACGT"
# )
#
# mergeIUPACLetters("ATC")
#
# mergeIUPACLetters <- function(x)
# {
#     if (!is.character(x) || any(is.na(x)) || any(nchar(x) == 0))
#         stop("'x' must be a vector of non-empty character strings")
#     x <- CharacterList(strsplit(toupper(x), "", fixed=TRUE))
#     yy <- unname(IUPAC_CODE_MAP[unlist(x, use.names=FALSE)])
#     if (any(is.na(yy)))
#         stop("some strings in 'x' contain non IUPAC letters")
#     yy <- CharacterList(strsplit(yy, "", fixed=TRUE))
#     y <- unstrsplit(sort(unique(IRanges:::regroupBySupergroup(yy, x))))
#     names(IUPAC_CODE_MAP)[match(y, IUPAC_CODE_MAP)]
# }





# fileConn<-file("/Users/chaokuan-hao/Documents/ANU_2019_Semester_2/Lanfear_Lab/sangeranalyseR/vignettes/SangerRead_Report.Rmd")
# readLines(fileConn)
# writeLines(c("Hello","World"), fileConn)
# close(fileConn)



