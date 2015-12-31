

get.processors <- function(processors){
    
    if(is.null(processors)){ 
        processors = detectCores(all.tests = FALSE, logical = FALSE)
    }

    return(processors)

}
