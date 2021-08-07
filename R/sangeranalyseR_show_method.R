# the show method:
setMethod('show', 'SangerRead', function(object){
    if (object@objectResults@creationResult) {
        if (object@inputSource == "ABIF") {
            primaryDNA <- as.character(object@primarySeq)
            secondaryDNA <- as.character(object@secondarySeq)
            trimmedStartPos <- object@QualityReport@trimmedStartPos
            trimmedFinishPos <- object@QualityReport@trimmedFinishPos
            primaryDNA <- substr(primaryDNA, trimmedStartPos+1, trimmedFinishPos)
            secondaryDNA <- substr(secondaryDNA, trimmedStartPos+1,trimmedFinishPos)
            if (object@QualityReport@TrimmingMethod == "M1") {
                cat("SangerRead S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "          Read Feature : ", object@readFeature, "\n",
                    "         Read FileName : ", basename(object@readFileName), "\n",
                    "       Trimming Method : ", object@QualityReport@TrimmingMethod, "\n",
                    "    M1 Trimming Cutoff : ", object@QualityReport@M1TrimmingCutoff, "\n",
                    "      Primary Sequence : ", primaryDNA, "\n",
                    "    Secondary Sequence : ", secondaryDNA, "\n"
                )
            } else if (object@QualityReport@TrimmingMethod == "M2") {
                cat("SangerRead S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "          Read Feature : ", object@readFeature, "\n",
                    "         Read FileName : ", basename(object@readFileName), "\n",
                    "       Trimming Method : ", object@QualityReport@TrimmingMethod, "\n",
                    "  M2 Cutoff Qual Score : ", object@QualityReport@M2CutoffQualityScore, "\n",
                    "  M2 Sliding Window Sz : ", object@QualityReport@M2SlidingWindowSize, "\n",
                    "      Primary Sequence : ", primaryDNA, "\n",
                    "    Secondary Sequence : ", secondaryDNA, "\n"
                )
            }
        } else if (object@inputSource == "FASTA") {
            cat("SangerRead S4 instance\n",
                "          Input Source : ", object@inputSource, "\n",
                "          Read Feature : ", object@readFeature, "\n",
                "         Read FileName : ", basename(object@readFileName), "\n",
                "       Fasta Read Name : ", object@fastaReadName, "\n",
                "      Primary Sequence : ", as.character(object@primarySeq), "\n"
            )
        }
        sapply(paste0("'", basename(object@readFileName), "'", " is successfully created!"), 
               log_success, simplify = FALSE)
    } else {
        sapply(paste0(object@objectResults@errorTypes, '\n',object@objectResults@errorMessages, '\n') , 
               log_error, simplify = FALSE)
    }
})


setMethod('show', 'SangerContig', function(object){
    if (object@objectResults@creationResult) {
        forReadNum <- length(object@forwardReadList)
        revReadNum <- length(object@reverseReadList)
        if (object@inputSource == "ABIF") {
            if (object@processMethod == "REGEX") {
                cat("SangerContig S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "        ABIF Directory : ", object@ABIF_Directory, "\n",
                    "  REGEX Suffix Forward : ", object@REGEX_SuffixForward, "\n",
                    "  REGEX Suffix Reverse : ", object@REGEX_SuffixReverse, "\n",
                    "           Contig Name : ", object@contigName, "\n",
                    "         'minReadsNum' : ", object@minReadsNum, "\n",
                    "       'minReadLength' : ", object@minReadLength, "\n",
                    "     'minFractionCall' : ", object@minFractionCall, "\n",
                    "     'maxFractionLost' : ", object@maxFractionLost, "\n",
                    "    'acceptStopCodons' : ", object@acceptStopCodons, "\n",
                    "        'readingFrame' : ", object@readingFrame, "\n",
                    "       Contig Sequence : ", as.character(object@contigSeq), "\n",
                    " Forward reads in the contig >> ", forReadNum, "\n",
                    " Reverse reads in the contig >> ", revReadNum, "\n"
                )
            } else if (object@processMethod == "CSV") {
                cat("SangerContig S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "        ABIF Directory : ", object@ABIF_Directory, "\n",
                    "  CSV Names Conversion : ", object@CSV_NamesConversion, "\n",
                    "           Contig Name : ", object@contigName, "\n",
                    "         'minReadsNum' : ", object@minReadsNum, "\n",
                    "       'minReadLength' : ", object@minReadLength, "\n",
                    "     'minFractionCall' : ", object@minFractionCall, "\n",
                    "     'maxFractionLost' : ", object@maxFractionLost, "\n",
                    "    'acceptStopCodons' : ", object@acceptStopCodons, "\n",
                    "        'readingFrame' : ", object@readingFrame, "\n",
                    "       Contig Sequence : ", as.character(object@contigSeq), "\n",
                    " Forward reads in the contig >> ", forReadNum, "\n",
                    " Reverse reads in the contig >> ", revReadNum, "\n"
                )
            }
        } else if (object@inputSource == "FASTA") {
            if (object@processMethod == "REGEX") {
                cat("SangerContig S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "       Fasta File Name : ", object@FASTA_File, "\n",
                    "  REGEX Suffix Forward : ", object@REGEX_SuffixForward, "\n",
                    "  REGEX Suffix Reverse : ", object@REGEX_SuffixReverse, "\n",
                    "           Contig Name : ", object@contigName, "\n",
                    "         'minReadsNum' : ", object@minReadsNum, "\n",
                    "       'minReadLength' : ", object@minReadLength, "\n",
                    "     'minFractionCall' : ", object@minFractionCall, "\n",
                    "     'maxFractionLost' : ", object@maxFractionLost, "\n",
                    "    'acceptStopCodons' : ", object@acceptStopCodons, "\n",
                    "        'readingFrame' : ", object@readingFrame, "\n",
                    "       Contig Sequence : ", as.character(object@contigSeq), "\n",
                    " Forward reads in the contig >> ", forReadNum, "\n",
                    " Reverse reads in the contig >> ", revReadNum, "\n"
                )    
            } else if (object@processMethod == "CSV") {
                cat("SangerContig S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "       Fasta File Name : ", object@FASTA_File, "\n",
                    "  CSV Names Conversion : ", object@CSV_NamesConversion, "\n",
                    "           Contig Name : ", object@contigName, "\n",
                    "         'minReadsNum' : ", object@minReadsNum, "\n",
                    "       'minReadLength' : ", object@minReadLength, "\n",
                    "     'minFractionCall' : ", object@minFractionCall, "\n",
                    "     'maxFractionLost' : ", object@maxFractionLost, "\n",
                    "    'acceptStopCodons' : ", object@acceptStopCodons, "\n",
                    "        'readingFrame' : ", object@readingFrame, "\n",
                    "       Contig Sequence : ", as.character(object@contigSeq), "\n",
                    " Forward reads in the contig >> ", forReadNum, "\n",
                    " Reverse reads in the contig >> ", revReadNum, "\n"
                )    
            }
        } 
        sapply(paste0("'", object@contigName, "'", " is successfully created!"), 
               log_success, simplify = FALSE)
    } else {
        sapply(paste0(object@objectResults@errorTypes, '\n', object@objectResults@errorMessages, '\n') , 
               log_error, simplify = FALSE)
    }
    if (length(object@objectResults@warningMessages) > 0) {
        sapply(paste0(object@warningMessages, '\n') , log_warn, simplify = FALSE)
    }
})

setMethod('show', 'SangerAlignment', function(object){
    if (object@objectResults@creationResult) {
        if (object@inputSource == "ABIF") {
            if (object@processMethod == "REGEX") {
                cat("SangerAlignment S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "        ABIF Directory : ", object@ABIF_Directory, "\n",
                    "  REGEX Suffix Forward : ", object@REGEX_SuffixForward, "\n",
                    "  REGEX Suffix Reverse : ", object@REGEX_SuffixReverse, "\n",
                    "     Contigs Consensus : ", as.character(object@contigsConsensus), "\n"
                )
            } else if (object@processMethod == "CSV") {
                cat("SangerAlignment S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "        ABIF Directory : ", object@ABIF_Directory, "\n",
                    "  CSV Names Conversion : ", object@CSV_NamesConversion, "\n",
                    "     Contigs Consensus : ", as.character(object@contigsConsensus), "\n"
                )
            }
        } else if (object@inputSource == "FASTA") {
            if (object@processMethod == "REGEX") {
                cat("SangerAlignment S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "       Fasta File Name : ", object@FASTA_File, "\n",
                    "  REGEX Suffix Forward : ", object@REGEX_SuffixForward, "\n",
                    "  REGEX Suffix Reverse : ", object@REGEX_SuffixReverse, "\n",
                    "     Contigs Consensus : ", as.character(object@contigsConsensus), "\n"
                )
            } else if (object@processMethod == "CSV") {
                cat("SangerAlignment S4 instance\n",
                    "          Input Source : ", object@inputSource, "\n",
                    "        Process Method : ", object@processMethod, "\n",
                    "       Fasta File Name : ", object@FASTA_File, "\n",
                    "  CSV Names Conversion : ", object@CSV_NamesConversion, "\n",
                    "     Contigs Consensus : ", as.character(object@contigsConsensus), "\n"
                )
            }
        }    
        sapply(paste0("'SangerAlignment'", " is successfully created!"), 
               log_success, simplify = FALSE)
    } else {
        sapply(paste0(object@objectResults@errorTypes, '\n', object@objectResults@errorMessages, '\n') , 
               log_error, simplify = FALSE)
    }
    if (length(object@objectResults@warningMessages) > 0) {
        sapply(paste0(object@warningMessages, '\n') , log_warn, simplify = FALSE)
    }
})
