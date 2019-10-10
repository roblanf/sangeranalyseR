#' Trim a sequences with Mott's modified trimming algorithm
#' 
#' Removes low quality bases from the start and end of a sequence. This version in R is based on the Biopython implementation of the same algorithm: http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-module.html#_abi_trim, with one minor change: the BioPython implementation always trims off the first base, while this implementation does not. For more information on the alogrithm, see http://www.phrap.org/phredphrap/phred.html http://www.clcbio.com/manual/genomics/Quality_abif_trimming.html
#' Note that the resulting trim.start and trim.end values correspond to the sequence contained in the abif.seq object, which is not always the same as the sequence contained in a sangerseq object derived from the same data (e.g. if the latter has had bases recalled using sangerseqR's makeBaseCalls() method).  
#'  
#' @param abif.seq a abif.seq s4 object from the sangerseqR package
#' @param cutoff the cutoff at which you consider a base to be bad. This works on a logarithmic scale, such that if you want to consider a score of 10 as bad, you set cutoff to 0.1; for 20 set it at 0.01; for 30 set it at 0.001; for 40 set it at 0.0001; and so on. Contiguous runs of bases below this quality will be removed from the start and end of the sequence. Given the high quality reads expected of most modern ABI sequencers, the defualt is 0.0001.
#' 
#' @return a list of two integers: "start" the start position of the sequence to keep; "finish" the end position of the sequence to keep.
#'
#' @keywords quality, phred, trimming
#'
#' @export trim.mott
#'

trim.mott <- function(abif.seq, cutoff = 0.0001){

    if(class(cutoff)!='numeric' | cutoff < 0){
        stop("cutoff must be a number of at least 0")
    }

    if(class(abif.seq)!='abif'){
        stop("abif.seq must be an 'abif' object from the sangerseqR package")
    }



    abif.seq = abif.seq@data
    start = FALSE # flag for starting position of trimmed sequence
    trim_start = 0 # init start index

    seqlen = nchar(abif.seq$PBAS.2)
    qual = abif.seq$PCON.2

    # calculate base score 
    score_list = cutoff - (10 ** (qual / -10.0))

    # calculate cummulative score 
    # if cumulative value < 0, set it to 0 
    # the BioPython implementation always trims the first base, 
    # this implementation does not. 
    score = score_list[1]
    if(score < 0){ 
        score = 0 
    }else{
        trim_start = 1
        start = TRUE
    }

    cummul_score = c(score)

    for(i in 2:length(score_list)){
        score = cummul_score[length(cummul_score)] + score_list[i]
        if(score <= 0){
            cummul_score = c(cummul_score, 0)
        }else{
            cummul_score = c(cummul_score, score)
            if(start == FALSE){
                # trim_start = value when cummulative score is first > 0 
                trim_start = i
                start = TRUE
            }
        }

        # trim_finish = index of highest cummulative score, 
        # marking the end of sequence segment with highest cummulative score 
        trim_finish = which.max(cummul_score)

    }

    # fix an edge case, where all scores are worse than the cutoff
    # in this case you wouldn't want to keep any bases at all
    if(sum(cummul_score)==0){trim_finish = 0}

    return(list("start" = trim_start, "finish" = trim_finish))

}
