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

                  # calculate base score
                  qualityBaseScore = - (10** (qualityScoreNumeric / (-10.0)))


                  # score_list = cutoff - (10 ** (qual / -10.0))
                  # # calculate cummulative score
                  # # if cumulative value < 0, set it to 0
                  # # the BioPython implementation always trims the first base,
                  # # this implementation does not.
                  # score = score_list[1]
                  # if(score < 0){
                  #     score = 0
                  # }else{
                  #     trim_start = 1
                  #     start = TRUE
                  # }
                  #
                  # cummul_score = c(score)
                  #
                  # for(i in 2:length(score_list)){
                  #     score = cummul_score[length(cummul_score)] + score_list[i]
                  #     if(score <= 0){
                  #         cummul_score = c(cummul_score, 0)
                  #     }else{
                  #         cummul_score = c(cummul_score, score)
                  #         if(start == FALSE){
                  #             # trim_start = value when cummulative score is first > 0
                  #             trim_start = i
                  #             start = TRUE
                  #         }
                  #     }
                  #
                  #     # trim_finish = index of highest cummulative score,
                  #     # marking the end of sequence segment with highest cummulative score
                  #     trim_finish = which.max(cummul_score)
                  #
                  # }
                  #
                  # # fix an edge case, where all scores are worse than the cutoff
                  # # in this case you wouldn't want to keep any bases at all
                  # if(sum(cummul_score)==0){trim_finish = 0}


              } else {
                  stop(errors)
              }
              callNextMethod(.Object, ...,
                             qualityScoreNumeric = qualityScoreNumeric,
                             qualityBaseScore    = qualityBaseScore,
                             trimmingStartPos    = 0L,
                             trimmingFinishPos   = 0L)
          })
