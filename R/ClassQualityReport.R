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
                  ### Quality Trimming (Using slideing window VERSION 1)
                  ### ----------------------------------------------------------
                  # calculate base score
                  # Calculate probability error per base (through column) ==> Q = -10log10(P)
                  readLen <- length(qualityScoreNumeric)

                  qualityPbCutoff <- 10** (cutoffQualityScore / (-10.0))
                  qualityBaseScore <- 10** (qualityScoreNumeric / (-10.0))

                  remainingIndex <- c()
                  for (i in 1:(readLen-slidingWindowSize+1)) {
                      meanSLidingWindow <- mean(qualityBaseScore[i:(i+slidingWindowSize-1)])
                      if (meanSLidingWindow < qualityPbCutoff) {
                          remainingIndex <- c(remainingIndex, i)
                          # or ==> i + floor(slidingWindowSize/3)
                      }
                  }
                  trimmingStartPos = remainingIndex[1]
                  trimmingFinishPos = remainingIndex[length(remainingIndex)]

                  ### ----------------------------------------------------------
                  ### Quality Report Visualization for single read
                  ###    1. Phed score for each pair
                  ###    2. Cumulative trimming percentage
                  ### ----------------------------------------------------------
                  # qualityPlotDf<- data.frame(1:length(qualityScoreNumeric),
                  #                            qualityScoreNumeric)
                  # colnames(qualityPlotDf) <- c("Index", "Score")
                  #
                  # ggplot(as.data.frame(qualityPlotDf),
                  #        aes(Index, Score)) +
                  #     geom_point()
                  #
                  #
                  # trimmingStartPos = 50L
                  # trimmingFinishPos = length(qualityScoreNumeric) - 20L
                  #
                  # stepRatio = 1 / readLen
                  # trimmingStartPos / readLen
                  # trimmingFinishPos / readLen
                  #
                  # trimmedPer <- c()
                  # remainingPer <- c()
                  #
                  # for (i in 1:trimmingStartPos) {
                  #     if (i != trimmingStartPos) {
                  #         print("Start Trimming")
                  #         trimmedPer <- c(trimmedPer, stepRatio)
                  #     }
                  # }
                  #
                  # for (i in trimmingStartPos:trimmingFinishPos) {
                  #     trimmedPer <- c(trimmedPer, 0)
                  # }
                  #
                  #
                  # for (i in trimmingFinishPos:readLen) {
                  #     if (i != trimmingFinishPos) {
                  #         print("End Trimming")
                  #         trimmedPer <- c(trimmedPer, stepRatio)
                  #     }
                  # }
                  #
                  # trimmedPer <- cumsum(trimmedPer)
                  # remainingPer = 1 - trimmedPer
                  #
                  # PerData <- data.frame(1:length(trimmedPer),
                  #                       trimmedPer, remainingPer)
                  #
                  # colnames(PerData) <- c("Base",
                  #                        "Trimmed Percent",
                  #                        "Remaining Percent")
                  #
                  # PerDataPlot <- melt(PerData, id.vars = c("Base"))
                  #
                  # ggplot(as.data.frame(PerDataPlot),
                  #        aes(x=Base, y=value, colour=variable)) +
                  #     geom_line() +
                  #     geom_point() +
                  #     theme_bw() +
                  #     xlab("Base Index") + ylab("Percentage") +
                  #     ggtitle("Quality Trimming Percentage Plot") +
                  #     theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  #           plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  #           axis.title.x = element_text(size = 10),
                  #           axis.title.y = element_text(size = 10),
                  #           legend.position="top",
                  #           legend.text = element_text(size = 6),
                  #           legend.title = element_blank())
                  #

              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             qualityScoreNumeric = qualityScoreNumeric,
                             qualityBaseScore    = qualityBaseScore,
                             trimmingStartPos    = trimmingStartPos,
                             trimmingFinishPos   = trimmingFinishPos,
                             cutoffQualityScore  = cutoffQualityScore,
                             slidingWindowSize   = slidingWindowSize)
          })


### ============================================================================
###
### ============================================================================
#' @export
setGeneric("trimmingRatioPlot", function(object) {
    standardGeneric("trimmingRatioPlot")
})

### ============================================================================
###
### ============================================================================
#' @export
setGeneric("qualityBasePlot", function(object) {
    standardGeneric("qualityBasePlot")
})



### ============================================================================
###
### ============================================================================
setMethod("trimmingRatioPlot",  "qualityReport", function(object){
    trimmingStartPos = object@trimmingStartPos
    trimmingFinishPos = object@trimmingFinishPos
    readLen = length(object@qualityScoreNumeric)

    stepRatio = 1 / readLen
    trimmingStartPos / readLen
    trimmingFinishPos / readLen

    trimmedPer <- c()
    remainingPer <- c()

    for (i in 1:trimmingStartPos) {
        if (i != trimmingStartPos) {
            print("Start Trimming")
            trimmedPer <- c(trimmedPer, stepRatio)
        }
    }

    for (i in trimmingStartPos:trimmingFinishPos) {
        trimmedPer <- c(trimmedPer, 0)
    }


    for (i in trimmingFinishPos:readLen) {
        if (i != trimmingFinishPos) {
            print("End Trimming")
            trimmedPer <- c(trimmedPer, stepRatio)
        }
    }

    trimmedPer <- cumsum(trimmedPer)
    remainingPer = 1 - trimmedPer

    PerData <- data.frame(1:length(trimmedPer),
                          trimmedPer, remainingPer)

    colnames(PerData) <- c("Base",
                           "Trimmed Percent",
                           "Remaining Percent")

    PerDataPlot <- melt(PerData, id.vars = c("Base"))

    ggplot(as.data.frame(PerDataPlot),
           aes(x=Base, y=value, colour=variable)) +
        geom_line() +
        geom_point() +
        theme_bw() +
        xlab("Base Index") + ylab("Percentage") +
        ggtitle("Quality Trimming Percentage Plot") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10),
              legend.position="top",
              legend.text = element_text(size = 6),
              legend.title = element_blank())
})



### ============================================================================
###
### ============================================================================
setMethod("qualityBasePlot",  "qualityReport", function(object){
    trimmingStartPos = object@trimmingStartPos
    trimmingFinishPos = object@trimmingFinishPos
    readLen = length(object@qualityScoreNumeric)

    qualityPlotDf<- data.frame(1:length(object@qualityScoreNumeric),
                               object@qualityScoreNumeric)
    colnames(qualityPlotDf) <- c("Index", "Score")

    ggplot(as.data.frame(qualityPlotDf),
           aes(Index, Score)) +
        geom_point() + theme_bw() +
        xlab("Base Index") + ylab("Phred Quality Score") +
        geom_vline(xintercept = trimmingStartPos,
                   color = "red", size=1) +
        geom_vline(xintercept = trimmingFinishPos,
                   color = "red", size=1) +
        ggtitle("Quality Trimming Percentage Plot") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10),
              legend.position="top",
              legend.text = element_text(size = 6),
              legend.title = element_blank())
})
