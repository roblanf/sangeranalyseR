#' @title BaseCall
#'
#' @description  An S4 class for base call information
#'
#' @slot maxSecondaryPeaks
#' @slot secondaryPeakRatio
#'
#' @name BaseCall-class
#'
#' @rdname BaseCall-class
#'
#' @exportClass BaseCall
#' @author Kuan-Hao Chao
#' @examples
setClass("BaseCall",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             maxSecondaryPeaks         = "numeric",
             secondaryPeakRatio        = "numeric"
         ),
)

### ============================================================================
### Overwrite initialize for QualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "QualityReport",
          function(.Object, ...,
                   readFeature         = character(0)) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              if (identical(readFeature, character(0))) {
                  msg <- paste("\nYou must assign value to 'readFeature'\n")
                  errors <- c(errors, msg)
              }
              if (length(errors) == 0) {
                  ### ----------------------------------------------------------
                  ### Prechecking success.
                  ### ----------------------------------------------------------

              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature         = readFeature)
          })
