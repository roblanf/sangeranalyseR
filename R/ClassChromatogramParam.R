#' @title ChromatogramParam
#'
#' @description  An S4 class storing chromatogram related inputs in a SangerRead S4 object.
#'
#' @slot baseNumPerRow  It defines maximum base pairs in each row. The default value is \code{100}.
#' @slot heightPerRow It defines the height of each row in chromatogram. The default value is \code{200}.
#' @slot signalRatioCutoff The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is \code{0.33}.
#' @slot showTrimmed The logical value storing whether to show trimmed base pairs in chromatogram. The default value is \code{TRUE}.
#'
#' @name ChromatogramParam-class
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
              errors <- list(character(0), character(0))
              errors <- checkBaseNumPerRow (baseNumPerRow, errors[[1]], errors[[2]])
              errors <- checkHeightPerRow (baseNumPerRow, errors[[1]], errors[[2]])
              errors <- checkSignalRatioCutoff (signalRatioCutoff, errors[[1]], errors[[2]])
              errors <- checkShowTrimmed (showTrimmed, errors[[1]], errors[[2]])
              if (length(errors[[1]]) == 0) {
                  ### ----------------------------------------------------------
                  ### Prechecking success.
                  ### ----------------------------------------------------------
              } else {
                log_error(paste(errors, collapse = ""))
              }
              callNextMethod(.Object, ...,
                             baseNumPerRow     = baseNumPerRow,
                             heightPerRow      = heightPerRow,
                             signalRatioCutoff = signalRatioCutoff,
                             showTrimmed       = showTrimmed)
          })

setClassUnion("ChromatogramParamORNULL", c("ChromatogramParam", "NULL"))

