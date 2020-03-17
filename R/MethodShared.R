#' @title Method launchApp
#'
#' @description  A method which launches Shiny application of the SangerContig
#'  and SangerAlignment instance.
#'
#' @param object object
#' @param outputDir outputDir
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(sangerContig)
#' data(sangerAlignment)
#' \dontrun{
#' launchApp(sangerContig)
#' launchApp(sangerAlignment)}
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
#' @param object object
#' @param outputDir outputDir
#' @param compress compress
#' @param compression_level compression_level
#' @param selection selection
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(sangerReadF)
#' data(sangerContig)
#' data(sangerAlignment)
#' \dontrun{
#' writeFasta(sangerReadF)
#' writeFasta(sangerContig)
#' writeFasta(sangerAlignment)}
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
#' @param object object
#' @param outputDir outputDir
#' @param includeSangerContig includeSangerContig
#' @param includeSangerRead includeSangerRead
#' @param ... ...
#'
#' @export
#' @author Kuan-Hao Chao
#' @examples
#' data(sangerReadF)
#' data(sangerContig)
#' data(sangerAlignment)
#' \dontrun{
#' generateReport(sangerReadF)
#' generateReport(sangerContig)
#' generateReport(sangerAlignment)}
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
