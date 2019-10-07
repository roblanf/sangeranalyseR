#' @title qualityReport
#'
#' @description  An S4 class for quality report for a SangeranalyseSeq S4 object
#'
#' @slot forward.read .
#'
#' @name qualityReport-class
#'
#' @rdname qualityReport-class
#'
#' @exportClass qualityReport
#' @author Kuan-Hao Chao
#' @examples
setClass("qualityReport",
         ### -------------------------------------------------------------------
         ### Input type of each variable
         ### -------------------------------------------------------------------
         representation(
             readFeature             = "character",
             qualityScoreNumeric     = "numeric",
             qualityBaseScore        = "numeric",
             trimmingStartPos        = "integer",
             trimmingFinishPos       = "integer",
             cutoffQualityScore      = "integer",
             slidingWindowSize       = "integer"
         ),
)

### ============================================================================
### Overwrite initialize for qualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "qualityReport",
          function(.Object, ...,
                   readFeature         = character(0),
                   qualityScoreNumeric = qualityScoreNumeric,
                   qualityBaseScore    = 0,
                   trimmingStartPos    = 0L,
                   trimmingFinishPos   = 0L,
                   cutoffQualityScore  = 20L,
                   slidingWindowSize   = 5L) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              errors <- character()
              if (identical(readFeature, character(0))) {
                  msg <- paste("\nYou must assign value to 'readFeature'\n")
                  errors <- c(errors, msg)
              }
              if (length(errors) == 0) {
                  ### ----------------------------------------------------------
                  ### Prechecking success.
                  ### ----------------------------------------------------------

                  ### ----------------------------------------------------------
                  ### Quality Trimming (Using slideing window VERSION 1)
                  ### ----------------------------------------------------------
                  # calculate base score
                  # Calculate probability error per base (through column)
                  #     ==> Q = -10log10(P)
                  readLen <- length(qualityScoreNumeric)

                  qualityPbCutoff <- 10** (cutoffQualityScore / (-10.0))
                  qualityBaseScore <- 10** (qualityScoreNumeric / (-10.0))

                  remainingIndex <- c()
                  for (i in 1:(readLen-slidingWindowSize+1)) {
                      meanSLidingWindow <-
                          mean(qualityBaseScore[i:(i+slidingWindowSize-1)])
                      if (meanSLidingWindow < qualityPbCutoff) {
                          remainingIndex <- c(remainingIndex, i)
                          # or ==> i + floor(slidingWindowSize/3)
                      }
                  }
                  trimmingStartPos = remainingIndex[1]
                  trimmingFinishPos = remainingIndex[length(remainingIndex)]

              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             readFeature         = readFeature,
                             qualityScoreNumeric = qualityScoreNumeric,
                             qualityBaseScore    = qualityBaseScore,
                             trimmingStartPos    = trimmingStartPos,
                             trimmingFinishPos   = trimmingFinishPos,
                             cutoffQualityScore  = cutoffQualityScore,
                             slidingWindowSize   = slidingWindowSize)
          })
