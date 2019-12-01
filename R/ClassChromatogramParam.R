#' @title ChromatogramParam
#'
#' @description  An S4 class for chromatogram parameter
#'
#' @slot baseNumPerRow .
#' @slot heightPerRow .
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
#'                      heightPerRow       = 200,
#'                      signalRatioCutoff  = 0.33,
#'                      showTrimmed        = TRUE)
setClass("ChromatogramParam",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             baseNumPerRow     = "numeric",
             heightPerRow      = "numeric",
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
                   heightPerRow      = 200,
                   signalRatioCutoff = 0.33,
                   showTrimmed       = TRUE) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              errors <- checkBaseNumPerRow (baseNumPerRow, errors)
              errors <- checkHeightPerRow (baseNumPerRow, errors)
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
                             heightPerRow      = heightPerRow,
                             signalRatioCutoff = signalRatioCutoff,
                             showTrimmed       = showTrimmed)
          })
