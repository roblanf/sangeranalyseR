

get.processors <- function(processors){
    
    if(Sys.info()["sysname"] == 'Windows'){
        # mclapply is not supported on windows
        # so we give a single processor,
        # in which case mclapply calls fall back
        # on lapply
        return(1)
    }

    if(is.null(processors)){ 
        processors = detectCores(all.tests = FALSE, logical = FALSE)
    }

    return(processors)

}
