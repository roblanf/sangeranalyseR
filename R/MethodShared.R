#' @title Method launchApp
#'
#' @description  A method which launches Shiny application of the SangerContig
#'  and SangerAlignment instance.
#'
#' @param object A SangerContig or SangerAlignment S4 instance.
#' @param outputDir The output directory of the saved new SangerContig or SangerAlignment S4 instance.
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(sangerContigData)
#' data(sangerAlignmentData)
#' \dontrun{
#' launchApp(sangerContigData)
#' launchApp(sangerAlignmentData)}
launchApp <- function(object, outputDir = NULL) {
    if(isS4(object)) {
        if (class(object)[1] == 'SangerAlignment') {
            message("Your input is 'SangerAlignment' S4 instance")
            shinyApp <- launchAppSA(object, outputDir = outputDir)
            return(shinyApp)
        } else if (class(object)[1] == 'SangerContig') {
            message("Your input is 'SangerContig' S4 instance")
            shinyApp <- launchAppSC(object, outputDir = outputDir)
            return(shinyApp)
        } else {
            stop("'object' must be 'SangerAlignment' or ",
                 "'SangerContig' S4 instance!")
        }
    } else {
        stop("'object' must be S4 instance!")
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
            message("Your input is 'SangerAlignment' S4 instance")
            writeFastaSA(object, outputDir = outputDir,
                         compress  = compress,
                         compression_level = compression_level,
                         selection = selection)
        } else if (class(object)[1] == 'SangerContig') {
            message("Your input is 'SangerContig' S4 instance")
            writeFastaSC(object, outputDir = outputDir,
                         compress  = compress,
                         compression_level = compression_level,
                         selection = selection)
        } else if (class(object)[1] == 'SangerRead') {
            message("Your input is 'SangerRead' S4 instance")
            writeFastaSR(object, outputDir = outputDir,
                         compress  = compress,
                         compression_level = compression_level)
        } else {
            stop("'object' must be 'SangerAlignment', 'SangerContig', or ",
                 "'SangerRead' S4 instance!")
        }
    } else {
        stop("'object' must be S4 instance!")
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
#' @param ... Further generateReportSR, generateReportSC, and generateReportSA related parameters.
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
            message("Your input is 'SangerAlignment' S4 instance")
            outputHtml <-
                generateReportSA(object, outputDir = outputDir,
                                 includeSangerContig = includeSangerContig,
                                 includeSangerRead = includeSangerContig, ...)
            return(outputHtml)
        } else if (class(object)[1] == 'SangerContig') {
            message("Your input is 'SangerContig' S4 instance")
            outputHtml <-
                generateReportSC(object, outputDir = outputDir,
                                 includeSangerRead = includeSangerRead, ...)
            return(outputHtml)
        } else if (class(object)[1] == 'SangerRead') {
            message("Your input is 'SangerRead' S4 instance")
            outputHtml <- generateReportSR(object, outputDir = outputDir, ...)
            return(outputHtml)
        } else {
            stop("'object' must be 'SangerAlignment', 'SangerContig', or ",
                 "'SangerRead' S4 instance!")
        }
    } else {
        stop("'object' must be S4 instance!")
    }
}
