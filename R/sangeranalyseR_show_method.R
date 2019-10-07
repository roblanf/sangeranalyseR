# the show method:
setMethod(f='show', signature='SangerMergeReads', function(object){
    # cat("SangerMergeReads S4 object\n",
    #     "             forwardReadFileName :", object@forwardReadFileName, "\n",
    #     "             reverseReadFileName :", object@reverseReadFileName, "\n",
    # )

    cat(is(object)[[1]], "\n",
        "       forwardReadFileName: ",
        object@forwardReadSangerseq@readFileName, "\n",
        "       reverseReadFileName:  ",
        object@reverseReadSangerseq@readFileName, "\n",
        sep = ""
    )
})
