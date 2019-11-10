inside_calculate_trimming <- function(qualityBaseScore,
                                      cutoffQualityScore,
                                      slidingWindowSize) {
    readLen <- length(qualityBaseScore)
    qualityPbCutoff <- 10** (cutoffQualityScore / (-10.0))
    remainingIndex <- c()
    if (slidingWindowSize > 20 || slidingWindowSize < 0 ||
        slidingWindowSize%%1!=0 ||
        cutoffQualityScore > 60 || cutoffQualityScore < 0 ||
        cutoffQualityScore%%1!=0) {
        trimmingStartPos = NULL
        trimmingFinishPos = NULL
    } else {
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
    }
    message("trimmingStartPos: ", trimmingStartPos, " trimmingFinishPos: ", trimmingFinishPos)
    return(c(trimmingStartPos, trimmingFinishPos))
}

dynamicMenuSideBar <- function(input, output, session, SangerSingleReadNum, SangerSingleReadFeature) {
    output$singleReadMenu <- renderMenu({
        menu_list <- sapply(1:SangerSingleReadNum, function(i) {
            list(menuItem(SangerSingleReadFeature[i],
                          tabName = SangerSingleReadFeature[i],
                          selected = TRUE, icon = icon("angle-right")))
        })
        sidebarMenu(.list = menu_list)
    })
    # Select consensus Read Menu first
    isolate({updateTabItems(session, "sidebar_menu", "Overview")})
}

