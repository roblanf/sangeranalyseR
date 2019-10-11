### ============================================================================
### Plotting trimmed and remaining ratio for "QualityReport" S4 object
### ============================================================================
setMethod("preTrimmingRatioPlot",  "QualityReport", function(object){
    readFeature <- object@readFeature
    trimmingStartPos = object@trimmingStartPos
    trimmingFinishPos = object@trimmingFinishPos
    readLen = length(object@qualityPhredScores)

    stepRatio = 1 / readLen
    trimmingStartPos / readLen
    trimmingFinishPos / readLen

    trimmedPer <- c()
    remainingPer <- c()

    for (i in 1:trimmingStartPos) {
        if (i != trimmingStartPos) {
            trimmedPer <- c(trimmedPer, stepRatio)
        }
    }

    for (i in trimmingStartPos:trimmingFinishPos) {
        trimmedPer <- c(trimmedPer, 0)
    }


    for (i in trimmingFinishPos:readLen) {
        if (i != trimmingFinishPos) {
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

    p <- ggplot(as.data.frame(PerDataPlot),
               aes(x=Base, y=value, colour=variable)) +
            geom_line() +
            geom_point() +
            theme_bw() +
            xlab("Base Index") + ylab("Percentage") +
            ggtitle(paste(readFeature, " Quality Trimming Percentage Plot")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(size = 15,
                                            face = "bold",
                                            hjust = 0.5),
                  axis.title.x = element_text(size = 10),
                  axis.title.y = element_text(size = 10),
                  legend.position="top",
                  legend.text = element_text(size = 6),
                  legend.title = element_blank())
})



### ============================================================================
### Plotting quality for each base for "QualityReport" S4 object
### ============================================================================
setMethod("preQualityBasePlot",  "QualityReport", function(object){
    readFeature <- object@readFeature
    trimmingStartPos = object@trimmingStartPos
    trimmingFinishPos = object@trimmingFinishPos
    readLen = length(object@qualityPhredScores)

    qualityPlotDf<- data.frame(1:length(object@qualityPhredScores),
                               object@qualityPhredScores)
    colnames(qualityPlotDf) <- c("Index", "Score")

    p <- ggplot(as.data.frame(qualityPlotDf),
               aes(Index, Score)) +
            geom_point() + theme_bw() +
            xlab("Base Index") + ylab("Phred Quality Score") +
            geom_vline(xintercept = trimmingStartPos,
                       color = "red", size=1) +
            geom_vline(xintercept = trimmingFinishPos,
                       color = "red", size=1) +
            ggtitle(paste(readFeature, " Quality Trimming Percentage Plot")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1),
                  plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  axis.title.x = element_text(size = 10),
                  axis.title.y = element_text(size = 10),
                  legend.position="top",
                  legend.text = element_text(size = 6),
                  legend.title = element_blank())
})



### ============================================================================
### Plotting trimmed and remaining ratio for "QualityReport" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead@QualityReport)
setMethod("trimmingRatioPlot",  "QualityReport", function(object){
    plotting <- preTrimmingRatioPlot(object)
    plotting
})



### ============================================================================
### Plotting quality for each base for "QualityReport" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' qualityBasePlot(A_chloroticaSingleRead@QualityReport)
setMethod("qualityBasePlot",  "QualityReport", function(object){
    plotting <- preQualityBasePlot(object)
    plotting
})























## =============================================================================
## Updating quality parameters for QualityReport object.
## =============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' QualityReport <- A_chloroticaSingleRead@QualityReport
#' trimmingRatioPlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' QualityReport@cutoffQualityScore
#' QualityReport@slidingWindowSize
#'
#' QualityReport <- updateQualityParam(QualityReport, 20L, 5L)
#'
#' trimmingRatioPlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' QualityReport@cutoffQualityScore
#' QualityReport@slidingWindowSize
setMethod("updateQualityParam",  "QualityReport",
          function(object,
                   cutoffQualityScore = 20L,
                   slidingWindowSize  = 5L){
              ### --------------------------------------------------------------
              ### Updating QualityReport quality parameters
              ### --------------------------------------------------------------
              qualityBaseScore <- object@qualityBaseScore
              trimmingPos <- inside_calculate_trimming(qualityBaseScore,
                                                       cutoffQualityScore,
                                                       slidingWindowSize)
              object@cutoffQualityScore <- cutoffQualityScore
              object@slidingWindowSize <- slidingWindowSize
              object@trimmingStartPos <- trimmingPos[1]
              object@trimmingFinishPos <- trimmingPos[2]
              return(object)
          })
