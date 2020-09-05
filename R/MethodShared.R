#' @title Method launchApp
#'
#' @description  A method which launches Shiny application of the SangerContig
#'  and SangerAlignment instance.
#'
#' @param object A SangerContig or SangerAlignment S4 instance.
#' @param outputDir The output directory of the saved new SangerContig or SangerAlignment S4 instance.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#'
#' @return A \code{SangerContig} or \code{SangerAlignment} object.
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' launchApp(sangerContigData)
#' launchApp(sangerAlignmentData)}
launchApp <- function(object, outputDir = NULL, colors="default") {
    if(isS4(object)) {
        if (class(object)[1] == 'SangerAlignment') {
            log_info("Your input is 'SangerAlignment' S4 instance")
            shinyApp <- launchAppSA(object, outputDir = outputDir, colors=colors)
            return(shinyApp)
        } else if (class(object)[1] == 'SangerContig') {
            log_info("Your input is 'SangerContig' S4 instance")
            shinyApp <- launchAppSC(object, outputDir = outputDir, colors=colors)
            return(shinyApp)
        } else {
            log_error("'object' must be 'SangerAlignment' or ",
                 "'SangerContig' S4 instance!")
        }
    } else {
        log_error("'object' must be S4 instance!")
    }
}

#' @title Method writeFasta
#'
#' @description  A method which writes FASTA files of the SangerRead,
#'  SangerContig, and SangerAlignment instance.
#'
#' @param object A SangerRead, SangerContig, or SangerAlignment S4 instance.
#' @param outputDir The output directory of generated FASTA files.
#' @param compress Like for the \code{save} function in base R, must be \code{TRUE} or \code{FALSE} (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip". This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param compression_level This parameter will be passed to \code{writeXStringSet} function in Biostrings package.
#' @param selection This parameter will be passed to \code{writeFastaSC} or \code{writeFastaSA}.
#'
#' @return A \code{SangerRead}, \code{SangerContig}, or \code{SangerAlignment} object.
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(sangerReadFData)
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' writeFasta(sangerReadFData)
#' writeFasta(sangerContigData)
#' writeFasta(sangerAlignmentData)}
writeFasta <- function(object, outputDir = NULL, compress  = FALSE,
                       compression_level = NA, selection = "all") {
    if(isS4(object)) {
        if (class(object)[1] == 'SangerAlignment') {
            log_info("Your input is 'SangerAlignment' S4 instance")
            writeFastaSA(object, outputDir = outputDir,
                         compress  = compress,
                         compression_level = compression_level,
                         selection = selection)
        } else if (class(object)[1] == 'SangerContig') {
            log_info("Your input is 'SangerContig' S4 instance")
            writeFastaSC(object, outputDir = outputDir,
                         compress  = compress,
                         compression_level = compression_level,
                         selection = selection)
        } else if (class(object)[1] == 'SangerRead') {
            log_info("Your input is 'SangerRead' S4 instance")
            writeFastaSR(object, outputDir = outputDir,
                         compress  = compress,
                         compression_level = compression_level)
        } else {
            log_error("'object' must be 'SangerAlignment', 'SangerContig', or ",
                 "'SangerRead' S4 instance!")
        }
    } else {
        log_error("'object' must be S4 instance!")
    }
}

#' @title Method generateReport
#'
#' @description  A method which generates final reports of the SangerRead,
#'  SangerContig, and SangerAlignment instance.
#'
#' @param object A SangerRead, SangerContig, or SangerAlignment S4 instance.
#' @param outputDir The output directory of the generated HTML report.
#' @param includeSangerContig The parameter that decides whether to include SangerContig level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param includeSangerRead The parameter that decides whether to include SangerRead level report. The value is \code{TRUE} or \code{FALSE} and the default is \code{TRUE}.
#' @param colors A vector for users to set the colors of (A, T, C, G, else). 
#'   There are three options for users to choose from. 
#'     1. "default":  (green, blue, black, red, purple). 
#'     2. "cb_friendly":  ((0, 0, 0), (199, 199, 199), (0, 114, 178), (213, 94, 0), (204, 121, 167)). 
#'     3. Users can set their own colors with a vector with five elements.
#' @param ... Further generateReportSR, generateReportSC, and generateReportSA related parameters.
#'
#' @return A \code{SangerRead}, \code{SangerContig}, or \code{SangerAlignment} object.
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(sangerReadFData)
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' generateReport(sangerReadFData)
#' generateReport(sangerContigData)
#' generateReport(sangerAlignmentData)}
generateReport <- function(object, outputDir = NULL,
                           includeSangerContig = TRUE,
                           includeSangerRead = TRUE, ...) {
    if(isS4(object)) {
        if (class(object)[1] == 'SangerAlignment') {
            log_info("Your input is 'SangerAlignment' S4 instance")
            outputHtml <-
                generateReportSA(object, outputDir = outputDir,
                                 includeSangerContig = includeSangerContig,
                                 includeSangerRead = includeSangerContig, 
                                 colors=colors, ...)
            return(outputHtml)
        } else if (class(object)[1] == 'SangerContig') {
            log_info("Your input is 'SangerContig' S4 instance")
            outputHtml <-
                generateReportSC(object, outputDir = outputDir,
                                 includeSangerRead = includeSangerRead, 
                                 colors=colors, ...)
            return(outputHtml)
        } else if (class(object)[1] == 'SangerRead') {
            log_info("Your input is 'SangerRead' S4 instance")
            outputHtml <- generateReportSR(object, outputDir = outputDir, 
                                           colors=colors, ...)
            return(outputHtml)
        } else {
            log_error("'object' must be 'SangerAlignment', 'SangerContig', or ",
                 "'SangerRead' S4 instance!")
        }
    } else {
        log_error("'object' must be S4 instance!")
    }
}
