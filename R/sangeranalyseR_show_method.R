# the show method:
setMethod(f='show', signature='SangerReads', function(object){
    # cat("SangerReads S4 object\n",
    #     "             forwardReadFileName :", object@forwardReadFileName, "\n",
    #     "             reverseReadFileName :", object@reverseReadFileName, "\n",
    # )

    cat(is(object)[[1]], "\n",
        "       forwardReadFileName: ", object@forwardReadFileName, "\n",
        "       reverseReadFileName:  ", object@reverseReadFileName, "\n",
        sep = ""
    )

})
