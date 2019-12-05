### ============================================================================
### Plotting trimmed and remaining ratio for "SangerSingleRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead)
setMethod("trimmingRatioPlot",  "SangerSingleRead", function(object){
    QualityReportObject = object@QualityReport
    plotting <- preTrimmingRatioPlot(QualityReportObject)
    plotting
})



### ============================================================================
### Plotting quality for each base for "SangerSingleRead" S4 object
### ============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' qualityBasePlot(A_chloroticaSingleRead)
setMethod("qualityBasePlot",  "SangerSingleRead", function(object){
    QualityReportObject = object@QualityReport
    plotting <- preQualityBasePlot(QualityReportObject)
    plotting
})

## =============================================================================
## Updating quality parameters for SangerSingleRead object.
## =============================================================================
#' @example
#' load("data/A_chloroticaSingleRead.RDdata")
#' trimmingRatioPlot(A_chloroticaSingleRead)
#' qualityBasePlot(A_chloroticaSingleRead)
#' A_chloroticaSingleRead@QualityReport@TrimmingMethod
#' A_chloroticaSingleRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaSingleRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@M2SlidingWindowSize
#'
#' A_chloroticaSingleRead <- updateQualityParam(A_chloroticaSingleRead,
#'                                              "M1", 0.0001, NULL, NULL)
#'
#' trimmingRatioPlot(A_chloroticaSingleRead)
#' qualityBasePlot(A_chloroticaSingleRead)
#' A_chloroticaSingleRead@QualityReport@TrimmingMethod
#' A_chloroticaSingleRead@QualityReport@M1TrimmingCutoff
#' A_chloroticaSingleRead@QualityReport@M2CutoffQualityScore
#' A_chloroticaSingleRead@QualityReport@M2SlidingWindowSize
setMethod("updateQualityParam",  "SangerSingleRead",
          function(object,
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL){
              ### --------------------------------------------------------------
              ### Updating SangerSingleRead quality parameters
              ### --------------------------------------------------------------
              object@QualityReport <- updateQualityParam(object@QualityReport,
                                                         TrimmingMethod,
                                                         M1TrimmingCutoff,
                                                         M2CutoffQualityScore,
                                                         M2SlidingWindowSize)
              return(object)
          })

setMethod("MakeBaseCalls", "SangerSingleRead",
          function(obj, ratio=.33) {
              #get peaks for each base
              Apeaks <- getpeaks(obj@traceMatrix[,1])
              Cpeaks <- getpeaks(obj@traceMatrix[,2])
              Gpeaks <- getpeaks(obj@traceMatrix[,3])
              Tpeaks <- getpeaks(obj@traceMatrix[,4])

              #get window around primary basecall peaks
              primarypeaks <- obj@peakPosMatrix[,1]
              diffs <- diff(c(0,primarypeaks))
              starts <- primarypeaks - 0.5*diffs
              stops <- c(primarypeaks[1:(length(primarypeaks)-1)] +
                             0.5*diffs[2:length(diffs)],
                         primarypeaks[length(diffs)] + 0.5*diffs[length(diffs)]
              )
              #hack for last peak. Just uses distance preceding peak
              #as distance after peak

              #Now get max peak value for each channel in each peak window.
              #If no peak return 0
              primary <- NULL
              secondary <- NULL
              tempPosMatrix <- matrix(nrow=length(starts), ncol=4)
              tempAmpMatrix <- matrix(nrow=length(starts), ncol=4)
              for(i in 1:length(starts)) {
                  Apeak <- peakvalues(Apeaks, starts[i], stops[i])
                  Cpeak <- peakvalues(Cpeaks, starts[i], stops[i])
                  Gpeak <- peakvalues(Gpeaks, starts[i], stops[i])
                  Tpeak <- peakvalues(Tpeaks, starts[i], stops[i])
                  if(is.na(Apeak[2]) &
                     is.na(Cpeak[2]) &
                     is.na(Gpeak[2]) &
                     is.na(Tpeak[2])) {
                      ### ------------------------------------------------------
                      ### My modification here: Add "N"
                      ###     Total length won't change.
                      ### ------------------------------------------------------
                      primary <- c(primary, "N")
                      secondary <- c(secondary, "N")
                      next #rare case where no peak found
                  }
                  signals <- c(Apeak[1], Cpeak[1], Gpeak[1], Tpeak[1])
                  tempAmpMatrix[i,] <- signals
                  positions <- c(Apeak[2], Cpeak[2], Gpeak[2], Tpeak[2])
                  tempPosMatrix[i,] <- positions
                  signalratios <- signals/max(signals, na.rm=TRUE)
                  Bases <- c("A", "C", "G", "T")
                  Bases[signalratios < ratio] <- NA
                  #sort by decreasing signal strength
                  Bases <- Bases[order(signals, decreasing=TRUE)]
                  positions <- positions[order(signals, decreasing=TRUE)]
                  if(length(Bases[!is.na(Bases)]) == 4
                     | length(Bases[!is.na(Bases)]) == 0) {
                      print(i)
                      primary <- c(primary, "N")
                      secondary <- c(secondary, "N")
                  }
                  else if(length(Bases[!is.na(Bases)]) > 1) {
                      print(i)
                      primary <- c(primary, Bases[1])
                      Bases2 <- Bases[2:4]
                      secondary <- c(secondary,
                                     mergeIUPACLetters(paste(sort(Bases2[!is.na(Bases2)]),
                                                             collapse="")))
                  }
                  else {
                      print(i)
                      primary <- c(primary, Bases[1])
                      secondary <- c(secondary, Bases[1])
                  }
              }
              obj@peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
              obj@peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
              obj@primarySeqID <- "sangerseq package primary basecalls"
              obj@primarySeq <- DNAString(paste(primary, collapse=""))
              obj@secondarySeqID <- "sangerseq package secondary basecalls"
              obj@secondarySeq <- DNAString(paste(secondary, collapse=""))

              return(obj)
          })
