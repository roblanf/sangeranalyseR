#' @title ChromatogramParam
#'
#' @description  An S4 class for chromatogram parameter
#'
#' @slot baseNumPerRow .
#' @slot signalRatioCutoff .
#' @slot showTrimmed .
#'
#' @name ChromatogramParam-class
#'
#' @rdname ChromatogramParam-class
#'
#' @exportClass ChromatogramParam
#' @author Kuan-Hao Chao
#' @examples
#' Chromatogram <- new("ChromatogramParam",
#'                      baseNumPerRow      = 100,
#'                      signalRatioCutoff  = 0.33,
#'                      showTrimmed        = TRUE)
setClass("ChromatogramParam",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             baseNumPerRow     = "numeric",
             signalRatioCutoff = "numeric",
             showTrimmed       = "logical"
         ),
)


### ============================================================================
### Overwrite initialize for QualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "ChromatogramParam",
          function(.Object, ...,
                   baseNumPerRow     = 100,
                   signalRatioCutoff = 0.33,
                   showTrimmed       = TRUE) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              errors <- checkBaseNumPerRow (baseNumPerRow, errors)
              errors <- checkSignalRatioCutoff (signalRatioCutoff, errors)
              errors <- checkShowTrimmed (showTrimmed, errors)
              if (length(errors) == 0) {
                  ### ----------------------------------------------------------
                  ### Prechecking success.
                  ### ----------------------------------------------------------
              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             baseNumPerRow     = baseNumPerRow,
                             signalRatioCutoff = signalRatioCutoff,
                             showTrimmed       = showTrimmed)
          })

