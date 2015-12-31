#' Trim a sequences with Mott's modified trimming algorithm
#' 
#' Removes low quality bases from the start and end of a sequence. This version in R was ported from the Biopython implementation of the same algorithm: http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-module.html#_abi_trim, with one minor change: the BioPython implementation always trims off the first base, while this implementation does not. For more information on the alogrithm, see http://www.phrap.org/phredphrap/phred.html http://www.clcbio.com/manual/genomics/Quality_abif_trimming.html
#' Note that the resulting trim.start and trim.end values correspond to the sequence contained in the abif.seq object, which is not always the same as the sequence contained in a sangerseq object derived from the same data (e.g. if the latter has had bases recalled using sangerseqR's makeBaseCalls() method).  
#'  
#' @param abif.seq a abif.seq s4 object from the sangerseqR package
#' @param segment minimum sequence length to return. If your input sequence is shorter than this, you just get your input sequence back.
#' 
#' @return a list of two integers: "trim.start" the start position of the sequence to keep; "trim.finish" the end position of the sequence to keep
#'
#' @keywords quality, phred, trimming
#'
#' @export trim.mott
#'

trim.mott <- function(abif.seq, cutoff = 0.05, segment = 20){

    abif.seq = abif.seq@data
    start = FALSE # flag for starting position of trimmed sequence
    trim_start = 0 # init start index

    seqlen = nchar(abif.seq$PBAS.2)
    qual = abif.seq$PCON.2

    if(seqlen <= segment){

        trim_start = 1
        trim_finish = seqlen

    }else{

        # calculate base score 
        score_list = cutoff - (10 ** (qual / -10.0))

        # calculate cummulative score 
        # if cummulative value < 0, set it to 0 
        # the BioPython implementation always trims the first base, 
        # this implementation does not. 

        cummul_score = c(score_list[1])

        for(i in 2:length(score_list)){
            score = cummul_score[length(cummul_score)] + score_list[i]
            if(score < 0){
                cummul_score = c(cummul_score, 0)
            }else{
                cummul_score = c(cummul_score, score)
                if(start == FALSE){
                    # trim_start = value when cummulative score is first > 0 
                    trim_start = i
                    start = TRUE
                }
            }
        }

        # trim_finish = index of highest cummulative score, 
        # marking the end of sequence segment with highest cummulative score 
        trim_finish = which.max(cummul_score)

    }

    return(list("trim.start" = trim_start, "trim.finish" = trim_finish))

}
