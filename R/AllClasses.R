#' @title SangerReads
#'
#' @description  An S4 class for storing reads of sanger sequencing.
#'
#' @slot forward.read .
#' @slot reverse.read .
#' @slot consensus.read .
#' @slot quality.matrix .
#'
#' @name SangerReads-class
#'
#' @rdname SangerReads-class
#'
#' @exportClass SangerReads
#' @author Kuan-Hao Chao
#' @examples
#' inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
#' A_chloroticaForwardReadFN <- file.path(inputFilesPath,
#'                                        "Allolobophora_chlorotica",
#'                                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_F.ab1")
#' A_chloroticaReverseReadFN <- file.path(inputFilesPath,
#'                                        "Allolobophora_chlorotica",
#'                                        "ACHLO006-09[LCO1490_t1,HCO2198_t1]_R.ab1")
#' A_chloroticaRead <- new("SangerReads",
#'                          forwardReadFileName = A_chloroticaForwardReadFN,
#'                          reverseReadFileName = A_chloroticaReverseReadFN)
setClass("SangerReads",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             forwardReadFileName = "character",
             reverseReadFileName = "character",
             forwardRead         = "character",
             reverseRead         = "character",
             consensusRead       = "character",
             qualityMatrix       = "data.frame"
         ),
         # ### -------------------------------------------------------------------
         # ### Initial value
         # ### -------------------------------------------------------------------
         # prototype(
         #     forwardReadFileName    = NULL,
         #     reverseReadFileName    = NULL,
         #     forwardRead            = NULL,
         #     reverseRead            = NULL,
         #     consensusRead          = "",
         #     qualityMatrix          = data.frame()),
         ### -------------------------------------------------------------------
         ### Input validity check
         ### -------------------------------------------------------------------
         validity = function(object) {
             errors <- character()
             if (!file.exists(object@forwardReadFileName)) {
                 msg <- paste("'", object@forwardReadFileName, "'",
                              " foward read file does not exist.", sep = "")
                 errors <- c(errors, msg)
             }
             if (!file.exists(object@reverseReadFileName)) {
                 msg <- paste("'", object@reverseReadFileName, "'",
                              " reverse read file does not exist.", sep = "")
                 errors <- c(errors, msg)
             }

             if (length(errors) == 0) {
                 forwardReadRF = read.abif(object@forwardReadFileName)
                 reverseReadRF = read.abif(object@reverseReadFileName)
                 seq.sanger = sangerseq(forwardRead)
                 TRUE
            } else {
                errors
            }
         }
)





### S4 class example package Try First
setClass("Person",
         slots = c(
             name = "character",
             age = "numeric"
         )
)


john <- new("Person", name = "John Smith", age = NA_real_)
setGeneric("age", function(x) standardGeneric("age"))
setGeneric("age<-", function(x, value) standardGeneric("age<-"))
setMethod("age", "Person", function(x) x@age)
setMethod("age<-", "Person", function(x, value) {
    x@age <- value
    x
})

age(john) <- 50

