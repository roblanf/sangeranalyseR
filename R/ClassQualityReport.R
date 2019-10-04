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
             qualityScoreNumeric     = "numeric",
             qualityBaseScore        = "numeric",
             trimmingStartPos        = "integer",
             trimmingFinishPos       = "integer"
         ),
)

### ============================================================================
### Overwrite initialize for qualityReport (New constructor)
### ============================================================================
setMethod("initialize",
          "qualityReport",
          function(.Object, ...,
                   qualityScoreNumeric = qualityScoreNumeric,
                   qualityBaseScore    = 0,
                   trimmingStartPos    = 0L,
                   trimmingFinishPos   = 0L) {
              ### --------------------------------------------------------------
              ### Input parameter prechecking
              ### --------------------------------------------------------------
              errors <- character()
              # if (typeof(qualityScoreNumeric) == "integer") {
              #     msg <- paste("\n'", qualityScoreNumeric, "'",
              #                  " data type should be integer.\n", sep = "")
              #     errors <- c(errors, msg)
              # }
              if (length(errors) == 0) {
                  ### ----------------------------------------------------------
                  ### Prechecking success.
                  ### ----------------------------------------------------------

                  ### ----------------------------------------------------------
                  ### Quality Trimming (Still need to add)
                  ### ----------------------------------------------------------
                  # calculate base score
                  # Calculate probability error per base (through column) ==> Q = -10log10(P)
                  qualityBaseScore = 10** (qualityScoreNumeric / (-10.0))
                  cumsum(qualityBaseScore)

                  cutoff = 0.0001
                  qualityBaseScoreCut = cutoff - (10 ** (qualityScoreNumeric / -10.0))
                  cumsum(qualityBaseScoreCut)

                  readLen = length(qualityScoreNumeric)

                  ### ----------------------------------------------------------
                  ### Quality Report Visualization for single read
                  ###    1. Phed score for each pair
                  ###    2. Cumulative trimming percentage
                  ### ----------------------------------------------------------
                  qualityPlotDf<- data.frame(1:length(qualityScoreNumeric),
                                             qualityScoreNumeric)
                  colnames(qualityPlotDf) <- c("Index", "Score")

                  ggplot(as.data.frame(qualityPlotDf),
                         aes(Index, Score)) +
                      geom_point()


                  trimmingStartPos = 3L
                  trimmingFinishPos = length(qualityScoreNumeric) - 5

                  stepRatio = 1 / readLen
                  trimmingStartPos / readLen
                  trimmingFinishPos / readLen

                  trimmedPer <- c()
                  remainedPer <- c()

                  for (i in 1:trimmingStartPos) {
                      if (i != trimmingStartPos) {
                          print("Start Trimming")
                          trimmedPer <- c(trimmedPer, stepRatio)
                          remainedPer <- c(remainedPer, 0)
                      }
                  }

                  for (i in trimmingStartPos:trimmingFinishPos) {
                      trimmedPer <- c(trimmedPer, 0)
                      remainedPer <- c(remainedPer, stepRatio)
                  }


                  for (i in trimmingFinishPos:readLen) {
                      if (i != trimmingFinishPos) {
                          print("End Trimming")
                          trimmedPer <- c(trimmedPer, stepRatio)
                          remainedPer <- c(remainedPer, 0)
                      }
                  }

                  trimmedPer <- cumsum(trimmedPer)
                  remainedPer <- cumsum(remainedPer)

                  trimmedPerPlot <- data.frame(1:length(trimmedPer),
                                               trimmedPer,
                                               remainedPer)
                  colnames(trimmedPerPlot) <- c("Index",
                                                "TrimmedPercent",
                                                "RemainingPercent")

                  ggplot(as.data.frame(trimmedPerPlot),
                         aes(Index, Percent)) +
                      geom_line() +
                      geom_point()



              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             qualityScoreNumeric = qualityScoreNumeric,
                             qualityBaseScore    = qualityBaseScore,
                             trimmingStartPos    = 0L,
                             trimmingFinishPos   = 0L)
          })
