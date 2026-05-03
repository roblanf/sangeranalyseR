#' @title ObjectResults
#'
#' @description  An S4 class storing results related inputs in a SangerRead, SangerContig, and SangerAlignment S4 object.
#' 
#' @slot creationResult
#' @slot errorMessages
#' @slot errorTypes
#' @slot warningMessages
#' @slot warningTypes
#' @slot readResultTable
#' @slot printLevel
#'
#' @name ObjectResults-class
#'
#' @exportClass ObjectResults
#' @author Kuan-Hao Chao
#' @examples
#' objectResults <- new("ObjectResults",
#'                      creationResult   = TRUE,
#'                      errorMessages    = character(0),
#'                      errorTypes       = character(0),
#'                      warningMessages  = character(0),
#'                      warningTypes     = character(0),
#'                      readResultTable =  data.frame(),
#'                      printLevel       = "SangerRead")
setClass("ObjectResults",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
           creationResult     = "logical",
           errorMessages      = "character",
           errorTypes         = "character",
           warningMessages    = "character",
           warningTypes       = "character",
           readResultTable    = "data.frame",
           printLevel         = "character"
         ),
)

### ============================================================================
### Overwrite initialize for QualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "ObjectResults",
          function(.Object, ...,
                   creationResult   = TRUE,
                   errorMessages    = character(0),
                   errorTypes       = character(0),
                   warningMessages  = character(0),
                   warningTypes     = character(0),
                   readResultTable  = data.frame(),
                   printLevel       = "SangerRead") {
            callNextMethod(.Object, ...,
                           creationResult   = creationResult,
                           errorMessages    = errorMessages,
                           errorTypes       = errorTypes,
                           warningMessages  = warningMessages,
                           warningTypes     = warningTypes,
                           readResultTable  = readResultTable,
                           printLevel       = printLevel)
          })

### ============================================================================
### Post-construction invariants for ObjectResults (Phase 4)
###
### Catches any code path that lands the parallel error/warning vectors in
### inconsistent states.
### ============================================================================
setValidity("ObjectResults", function(object) {
    errs <- character()
    if (length(object@creationResult) != 1L ||
        !is.logical(object@creationResult)) {
        errs <- c(errs, "creationResult must be a single logical value")
    }
    if (length(object@errorMessages) != length(object@errorTypes)) {
        errs <- c(errs, paste0("errorMessages and errorTypes must have ",
                               "equal length"))
    }
    if (length(object@warningMessages) != length(object@warningTypes)) {
        errs <- c(errs, paste0("warningMessages and warningTypes must have ",
                               "equal length"))
    }
    if (length(errs) == 0L) TRUE else errs
})

