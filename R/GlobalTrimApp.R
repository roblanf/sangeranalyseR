# =============================================================================
# Phase 8 — Global Trimming Controls dashboard.
#
# A standalone Shiny gadget that lets a user adjust M1/M2 trimming
# parameters across an *entire* SangerAlignment in one place and instantly
# preview how the consensus changes. Built on top of the existing
# updateQualityParam,SangerAlignment-method (R/MethodSangerAlignment.R).
#
# Usage:
#     sa <- SangerAlignment(...)
#     globalTrimApp(sa)              # opens a Shiny gadget in the viewer
#
# The gadget returns the re-trimmed SangerAlignment (with the user's last
# applied parameters) when the user clicks "Done", or NULL on cancel —
# matching the convention of shiny::runGadget().
# =============================================================================

#' @title Launch a global trimming controls dashboard
#'
#' @description Opens a Shiny gadget that exposes M1/M2 trimming
#'   parameters as sliders/numeric inputs, applies any change to *every*
#'   SangerRead in the supplied SangerAlignment via
#'   `updateQualityParam(SA, ...)`, and live-previews the consensus length,
#'   contig count, and number of reads that survive the new trim.
#'
#'   Unlike `launchAppSA()` (which exposes per-read trimming) this app
#'   operates globally — useful when a whole batch needs the same
#'   re-trimming policy.
#'
#' @param SA A SangerAlignment instance.
#' @return The (last-applied) SangerAlignment when the user clicks "Done",
#'   or NULL if the user cancels.
#'
#' @examples
#' \dontrun{
#' data(sangerAlignmentData)
#' SA2 <- globalTrimApp(sangerAlignmentData)
#' }
#'
#' @export
globalTrimApp <- function(SA) {
    if (!is(SA, "SangerAlignment")) {
        stop("'SA' must be a SangerAlignment instance.")
    }
    if (SA@inputSource != "ABIF") {
        stop("globalTrimApp only supports inputSource = 'ABIF' SangerAlignments. ",
             "FASTA-derived alignments have no quality scores to retrim.")
    }

    ui <- shiny::fluidPage(
        shiny::titlePanel("Global Trimming Controls"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::radioButtons("trimMethod", "Trimming method",
                                     choices = c("M1 (modified Mott)" = "M1",
                                                 "M2 (sliding window)" = "M2"),
                                     selected = "M1"),
                shiny::conditionalPanel(
                    "input.trimMethod == 'M1'",
                    shiny::sliderInput("M1cutoff", "M1 trimming cutoff",
                                        min = 0.00001, max = 0.1,
                                        value = 0.0001, step = 0.00005)
                ),
                shiny::conditionalPanel(
                    "input.trimMethod == 'M2'",
                    shiny::sliderInput("M2score", "M2 cutoff Phred score",
                                        min = 0L, max = 60L,
                                        value = 20L, step = 1L),
                    shiny::sliderInput("M2window", "M2 sliding window size",
                                        min = 1L, max = 40L,
                                        value = 10L, step = 1L)
                ),
                shiny::actionButton("apply", "Apply",
                                     icon = shiny::icon("play"),
                                     class = "btn-primary"),
                shiny::tags$hr(),
                shiny::actionButton("done", "Done",
                                     icon = shiny::icon("check")),
                shiny::actionButton("cancel", "Cancel",
                                     icon = shiny::icon("xmark"))
            ),
            shiny::mainPanel(
                shiny::h4("Live preview"),
                shiny::tableOutput("summary"),
                shiny::tags$hr(),
                shiny::h4("Consensus (post-trim)"),
                shiny::verbatimTextOutput("consensus_preview"),
                shiny::tags$hr(),
                shiny::h4("Per-contig snapshot"),
                shiny::tableOutput("contig_table")
            )
        )
    )

    server <- function(input, output, session) {
        rv <- shiny::reactiveValues(SA = SA, applied_count = 0L)

        # --- summary panel ---
        render_summary <- function(sa) {
            data.frame(
                Metric = c("Method",
                           "Total contigs",
                           "Total reads (post-trim)",
                           "Consensus length (bp)"),
                Value  = c(input$trimMethod,
                            length(sa@contigList),
                            sum(vapply(sa@contigList, function(sc)
                                length(sc@forwardReadList) + length(sc@reverseReadList),
                                integer(1L))),
                            length(sa@contigsConsensus))
            )
        }

        output$summary <- shiny::renderTable({ render_summary(rv$SA) })

        output$consensus_preview <- shiny::renderText({
            cs <- as.character(rv$SA@contigsConsensus)
            if (nchar(cs) == 0L) return("(empty consensus)")
            # Wrap at 80 columns for readable preview.
            substr_lines <- vapply(seq.int(1L, nchar(cs), by = 80L),
                                    function(i) substr(cs, i, i + 79L),
                                    character(1L))
            paste(substr_lines, collapse = "\n")
        })

        output$contig_table <- shiny::renderTable({
            sa <- rv$SA
            data.frame(
                Contig    = vapply(sa@contigList, function(sc) sc@contigName, character(1L)),
                Forward   = vapply(sa@contigList, function(sc) length(sc@forwardReadList), integer(1L)),
                Reverse   = vapply(sa@contigList, function(sc) length(sc@reverseReadList), integer(1L)),
                ContigLen = vapply(sa@contigList, function(sc) length(sc@contigSeq), integer(1L))
            )
        })

        # --- apply button ---
        shiny::observeEvent(input$apply, {
            shiny::showNotification("Re-trimming all reads ...", type = "message",
                                     duration = 2L)
            new_sa <- tryCatch({
                if (input$trimMethod == "M1") {
                    updateQualityParam(rv$SA,
                                        TrimmingMethod        = "M1",
                                        M1TrimmingCutoff      = input$M1cutoff,
                                        M2CutoffQualityScore  = NULL,
                                        M2SlidingWindowSize   = NULL)
                } else {
                    updateQualityParam(rv$SA,
                                        TrimmingMethod        = "M2",
                                        M1TrimmingCutoff      = NULL,
                                        M2CutoffQualityScore  = input$M2score,
                                        M2SlidingWindowSize   = input$M2window)
                }
            }, error = function(e) {
                shiny::showNotification(paste("Re-trim failed:",
                                               conditionMessage(e)),
                                         type = "error")
                NULL
            })
            if (!is.null(new_sa)) {
                rv$SA <- new_sa
                rv$applied_count <- rv$applied_count + 1L
                shiny::showNotification(
                    sprintf("Applied (%d so far). New consensus length: %d bp.",
                             rv$applied_count, length(rv$SA@contigsConsensus)),
                    type = "default", duration = 3L)
            }
        })

        # --- done / cancel ---
        shiny::observeEvent(input$done, {
            shiny::stopApp(rv$SA)
        })
        shiny::observeEvent(input$cancel, {
            shiny::stopApp(NULL)
        })
    }

    shiny::runGadget(ui, server,
                      viewer = shiny::dialogViewer("Global Trimming Controls",
                                                    width = 1100L, height = 760L))
}
