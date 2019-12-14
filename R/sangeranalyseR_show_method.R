# the show method:
setMethod('show', 'SangerRead', function(object){
    cat("SangerRead S4 instance\n",
        "          Read Feature : ", object@readFeature, "\n",
        "         Read FileName : ", basename(object@readFileName), "\n",
        "      Primary Sequence : ", as.character(object@primarySeq), "\n",
        "    Secondary Sequence : ", as.character(object@secondarySeq), "\n"
    )
})

setMethod('show', 'SangerContig', function(object){
    cat("SangerContig S4 instance\n",
        "      Parent Directory : ", object@parentDirectory, "\n",
        "           Contig Name : ", object@contigName, "\n",
        " Suffix Forward RegExp : ", object@suffixForwardRegExp, "\n",
        " Suffix Reverse RegExp : ", object@suffixReverseRegExp, "\n",
        "    Trimming Method SC : ", object@trimmingMethodSC, "\n",
        "         'minReadsNum' : ", object@minReadsNum, "\n",
        "       'minReadLength' : ", object@minReadLength, "\n",
        "     'minFractionCall' : ", object@minFractionCall, "\n",
        "     'maxFractionLost' : ", object@maxFractionLost, "\n",
        "    'acceptStopCodons' : ", object@acceptStopCodons, "\n",
        "        'readingFrame' : ", object@readingFrame, "\n",
        "       Contig Sequence : ", as.character(object@contigSeq), "\n"
    )
})

setMethod('show', 'SangerAlignment', function(object){
    cat("SangerAlignment S4 instance\n",
        "      Parent Directory : ", object@parentDirectory, "\n",
        " Suffix Forward RegExp : ", object@suffixForwardRegExp, "\n",
        " Suffix Reverse RegExp : ", object@suffixReverseRegExp, "\n",
        "    'trimmingMethodSA' : ", object@trimmingMethodSA, "\n",
        "   'minFractionCallSA' : ", object@minFractionCallSA, "\n",
        "   'maxFractionLostSA' : ", object@maxFractionLostSA, "\n",
        "     Contigs Consensus : ", as.character(object@contigsConsensus), "\n"
    )
})
