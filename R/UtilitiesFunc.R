inside_calculate_trimming <- function(qualityBaseScore,
                                      cutoffQualityScore,
                                      slidingWindowSize) {
    readLen <- length(qualityBaseScore)
    qualityPbCutoff <- 10** (cutoffQualityScore / (-10.0))
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

    return(c(trimmingStartPos, trimmingFinishPos))
}
