### ============================================================================
### Plotting quality for each base for "QualityReport" S4 object
### ============================================================================
setMethod("preQualityBasePlot",  "QualityReport", function(object, readFeature){
    # trimmedStartPos = object@trimmedStartPos
    # trimmedFinishPos = object@trimmedFinishPos
    # readLen = length(object@qualityPhredScores)
    #
    # qualityPlotDf<- data.frame(1:length(object@qualityPhredScores),
    #                            object@qualityPhredScores)
    # colnames(qualityPlotDf) <- c("Index", "Score")
    #
    # p <- ggplot(as.data.frame(qualityPlotDf),
    #            aes(Index, Score)) +
    #         geom_point() + theme_bw() +
    #         xlab("Base Index") + ylab("Phred Quality Score") +
    #         geom_vline(xintercept = trimmedStartPos,
    #                    color = "red", size=1) +
    #         geom_vline(xintercept = trimmedFinishPos,
    #                    color = "red", size=1) +
    #         ggtitle(paste(readFeature, " Quality Trimming Percentage Plot")) +
    #         theme(axis.text.x = element_text(angle = 90, hjust = 1),
    #               plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    #               axis.title.x = element_text(size = 10),
    #               axis.title.y = element_text(size = 10),
    #               legend.position="top",
    #               legend.text = element_text(size = 6),
    #               legend.title = element_blank())
    trimmedStartPos = object@trimmedStartPos
    trimmedFinishPos = object@trimmedFinishPos
    qualityPhredScores = object@qualityPhredScores
    readLen = length(qualityPhredScores)
    qualityPlotDf<- data.frame(1:length(qualityPhredScores),
                               qualityPhredScores)
    colnames(qualityPlotDf) <- c("Index", "Score")
    x <- list(
        title = "Base Pair Index"
        # titlefont = f
    )
    y <- list(
        title = "Phred Quality Score"
        # titlefont = f
    )
    p <- QualityBasePlotly(trimmedStartPos, trimmedFinishPos,
                           readLen, qualityPlotDf, x,  y)
})

### ============================================================================
### Plotting trimmed and remaining ratio for "QualityReport" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead@QualityReport)
setMethod("trimmingRatioPlot",  "QualityReport", function(object){
    plotting <- preTrimmingRatioPlot(object)
    plotting
})



### ============================================================================
### Plotting quality for each base for "QualityReport" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaRead.RDdata")
#' qualityBasePlot(A_chloroticaSingleRead@QualityReport)
setMethod("qualityBasePlot",  "QualityReport", function(object){
    plotting <- preQualityBasePlot(object)
    plotting
})

## =============================================================================
## Updating quality parameters for QualityReport object.
## =============================================================================
#' @example
#' load("data/A_chloroticaRead.RDdata")
#' QualityReport <- A_chloroticaSingleRead@QualityReport
#' trimmingRatioPlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' QualityReport@TrimmingMethod
#' QualityReport@M1TrimmingCutoff
#' QualityReport@M2CutoffQualityScore
#' QualityReport@M2SlidingWindowSize
#'
#' QualityReport <- updateQualityParam(QualityReport, "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' qualityBasePlot(QualityReport)
#' QualityReport@TrimmingMethod
#' QualityReport@M1TrimmingCutoff
#' QualityReport@M2CutoffQualityScore
#' QualityReport@M2SlidingWindowSize
setMethod("updateQualityParam",  "QualityReport",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL) {
              ### --------------------------------------------------------------
              ### Updating QualityReport quality parameters
              ### --------------------------------------------------------------
              qualityBaseScores <- object@qualityBaseScores
              trimmingPos <- M2inside_calculate_trimming(qualityBaseScores,
                                                       TrimmingMethod,
                                                       M1TrimmingCutoff,
                                                       M2CutoffQualityScore,
                                                       M2SlidingWindowSize)

              object@TrimmingMethod <- TrimmingMethod
              object@M1TrimmingCutoff <- M1TrimmingCutoff
              object@M2CutoffQualityScore <- M2CutoffQualityScore
              object@M2SlidingWindowSize <- M2SlidingWindowSize

              object@rawSeqLength <- trimmingPos[1]
              object@rawMeanQualityScore <- trimmingPos[2]
              object@rawMinQualityScore <- trimmingPos[3]
              object@trimmedStartPos <- trimmingPos[4]
              object@trimmedFinishPos <- trimmingPos[5]
              object@trimmedSeqLength <- trimmingPos[6]
              object@trimmedMeanQualityScore <- trimmingPos[7]
              object@trimmedMinQualityScore <- trimmingPos[8]
              return(object)
          })
