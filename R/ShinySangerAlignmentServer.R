### ============================================================================
### R shiny SangerAlignment server function
### ============================================================================
SangerAlignmentServer <- function(input, output, session) {
    # Suppress Warning
    options(warn = -1)
    ### ------------------------------------------------------------------------
    ### SangerAlignment parameters initialization.
    ### ------------------------------------------------------------------------
    SangerAlignment <- getShinyOption("sangerAlignment")
    shinyDirectory <- getShinyOption("shinyDirectory")
    colors <- getShinyOption("colors")
    SangerAlignment <- SangerAlignment[[1]]
    SangerAlignmentNum <- length(SangerAlignment@contigList)
    SangerAlignmentParam <- lapply(seq_len(SangerAlignmentNum), function(i) {
        ### --------------------------------------------------------------------
        ### SangerContig-related parameters initialization.
        ### --------------------------------------------------------------------
        # readFeature
        SCName <- paste0(i, " Contig")
        # Forward & reverse reads list
        SangerReadFReadList <-
            SangerAlignment@contigList[[i]]@forwardReadList
        SangerReadRReadsList <-
            SangerAlignment@contigList[[i]]@reverseReadList
        forwardReadNum <- length(SangerReadFReadList)
        reverseReadNum <- length(SangerReadRReadsList)
        # readFeature
        forwardReadFeature <- vapply(seq_len(forwardReadNum), function(j)
            paste0(j, " ",
                   SangerReadFReadList[[j]]@readFeature), character(1))
        reverseReadFeature <- vapply(seq_len(reverseReadNum), function(j)
            paste0(j, " ",
                   SangerReadRReadsList[[j]]@readFeature), character(1))
        SangerReadFeature <- c(forwardReadFeature, reverseReadFeature)
        # readFileName (basename)
        forwardReadBFN <- vapply(seq_len(forwardReadNum), function(j)
            basename(SangerReadFReadList[[j]]@readFileName), character(1))
        reverseReadBFN <- vapply(seq_len(reverseReadNum), function(j)
            basename(SangerReadRReadsList[[j]]@readFileName), character(1))
        SangerReadBFN <- c(forwardReadBFN, reverseReadBFN)
        return(list(SCName = SCName,
                    forwardReadNum = forwardReadNum,
                    reverseReadNum = reverseReadNum,
                    forwardReadFeature = forwardReadFeature,
                    reverseReadFeature = reverseReadFeature,
                    SangerReadBFN = SangerReadBFN))
    })
    SATrimmingMethod <- SangerAlignment@trimmingMethodSA
    if (SATrimmingMethod == "M1") {
        SATrimmingMethodName = "Method 1: 'Modified Mott Trimming'"
    } else if (SATrimmingMethod == "M2") {
        SATrimmingMethodName =
            "Method 2: 'Trimmomatics Sliding Window Trimming'"
    }

    ### ------------------------------------------------------------------------
    ### SangerAlignment reactiveValue
    ### ------------------------------------------------------------------------
    sangerAlignmentParam <-
        reactiveValues(
            contigsConsensus = SangerAlignment@contigsConsensus,
            contigsAlignment = as.character(SangerAlignment@contigsAlignment),
            contigsTree = SangerAlignment@contigsTree)
    ### ------------------------------------------------------------------------
    ### SangerContig reactiveValue
    ### ------------------------------------------------------------------------
    contigParam <-
        reactiveValues(
            contigSeq   = NULL,
            differencesDF   = NULL,
            alignment       = NULL,
            distanceMatrix  = NULL,
            dendrogram      = NULL,
            indelsDF        = NULL,
            stopCodonsDF    = NULL,
            secondaryPeakDF = NULL)
    ### ------------------------------------------------------------------------
    ### SingleRead reactiveValue
    ### ------------------------------------------------------------------------
    sequenceParam <- reactiveValues(primarySeq = "",
                                    secondarySeq = "",
                                    primaryAASeqS1 = "",
                                    primaryAASeqS2 = "",
                                    primaryAASeqS3 = "")
    trimmedRV <- reactiveValues(rawSeqLength            = 0,
                                rawMeanQualityScore     = 0,
                                rawMinQualityScore      = 0,
                                trimmedStartPos         = 0,
                                trimmedFinishPos        = 0,
                                trimmedSeqLength        = 0,
                                trimmedMeanQualityScore = 0,
                                trimmedMinQualityScore  = 0,
                                remainingRatio          = 0)
    trimmedParam <- reactiveValues(M1TrimmingCutoff     = 0,
                                   M2CutoffQualityScore = 0,
                                   M2SlidingWindowSize  = 0)
    ChromatogramParam <- reactiveValues(baseNumPerRow     = 0,
                                        heightPerRow      = 0,
                                        signalRatioCutoff = 0,
                                        showTrimmed       = TRUE)

    ############################################################################
    ### Main panel dynamic UI
    ############################################################################
    output$aligned_contigSeq_content <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (sidebar_menu[[1]] == "Contigs" && 
            sidebar_menu[[2]] == "Alignment" && 
            sidebar_menu[[3]] == "Overview" && 
            sidebar_menu[[4]] == "Page" && 
            sidebar_menu[[5]] == "_" ) {
            h1(input$sidebar_menu)
            log_info(">>>>>>>> Inside 'Contigs Alignment Overview Page _'")
            shinyjs::disable("closeUI")
            shinyjs::disable("recalculateButtonSA")
            log_info("######## Reactive button clicked !!!")
            log_info("######## Start re-aligning contigs")
            CSSetResult <-
                alignContigs (SangerAlignment@contigList,
                              SangerAlignment@geneticCode,
                              SangerAlignment@refAminoAcidSeq,
                              SangerAlignment@minFractionCallSA,
                              SangerAlignment@maxFractionLostSA,
                              1)
            SangerAlignment@contigsConsensus <<- CSSetResult$consensus
            SangerAlignment@contigsAlignment <<- CSSetResult$aln
            SangerAlignment@contigsTree <<- CSSetResult$aln.tree
            sangerAlignmentParam[["contigsConsensus"]] <<- SangerAlignment@contigsConsensus
            sangerAlignmentParam[["contigsAlignment"]] <<- as.character(SangerAlignment@contigsAlignment)
            sangerAlignmentParam[["contigsTree"]] <<- SangerAlignment@contigsTree
            log_info("######## Finish contigs re-alignment")
            shinyjs::enable("recalculateButtonSA")
            shinyjs::enable("closeUI")
            fluidRow(
                useShinyjs(),
                box(title = tags$p("Input Parameters: ",
                                   style = "font-size: 26px;
                                                   font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 0.2px hidden #A9A9A9;")),
                    fluidRow(
                        column(width = 12,
                               actionBttn("recalculateButtonSA",
                                          "Re-calculate Contigs Alignment",
                                          icon = icon("calculator"),
                                          style = "simple", color = "danger",
                                          block = TRUE, size = "lg")
                        ),
                        column(12,
                               tags$hr(
                                   style = ("border-top: 2px hidden #A9A9A9;"))
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Output Directory: "),
                                             style = "font-size: 20px;
                                                   font-weight: bold;"),
                               ),
                               column(9,
                                      h4(shinyDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Raw ABI Parent Directory:"),
                                          style = "font-size: 20px;
                                                   font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerAlignment@parentDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Trimming Method: "),
                                             style = "font-size: 20px;
                                                   font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SATrimmingMethodName),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Forward Suffix RegExp: "),
                                             style = "font-size: 20px;
                                                   font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerAlignment@suffixForwardRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Reverse Suffix RegExp: "),
                                             style = "font-size: 20px;
                                                   font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerAlignment@suffixReverseRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Contigs Number: "),
                                             style = "font-size: 20px;
                                                   font-weight: bold;"),
                               ),
                               column(9,
                                      h4(length(SangerAlignment@
                                                    contigList)),
                               )
                        ),
                        tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                        box(title = tags$p("Alignment Parameters",
                                           style = "font-size: 24px;
                                       font-weight: bold;"),
                            collapsible = TRUE,
                            status = "success", width = 12,
                            column(4,
                                   uiOutput("SAMinFractionCallSA") ,
                            ),
                            column(4,
                                   uiOutput("SAMaxFractionLostSA") ,
                            ),
                        ),
                        tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                        box(title = tags$p("Genetic Code Data Frame",
                                           style = "font-size: 24px;
                                                   font-weight: bold;"),
                            collapsible = TRUE,
                            status = "success", width = 12,
                            column(width = 2,
                                   tags$p("Tri-nucleotide:",
                                          style = "font-size: 15px;
                                                   font-weight: bold;"),
                                   tags$p("Amino Acid : ",
                                          style = "font-size: 15px;
                                                   font-weight: bold;"),
                                   tags$p("('*' : stop codon) ",
                                          style = "font-size: 12px;
                                                   font-weight: italic;"),
                            ),
                            column(width = 10,
                                   excelOutput("geneticCodeDF",
                                               width = "100%", height = "50"),
                                   style = paste("height:100%; ",
                                                 "overflow-y: hidden;",
                                                 "overflow-x: scroll;")
                            ),
                        ),
                        uiOutput("SArefAminoAcidSeq") ,
                    ),
                ),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "SangerAlignment Results: "),
                                   style = "font-size: 26px;
                                                   font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Contigs Alignment",
                                       style = "font-size: 24px;
                                                   font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               htmlOutput("sangerAlignmentAlignmentHTML"),
                        ),
                    ),
                    box(title = tags$p("Contigs Tree",
                                       style = "font-size: 24px;
                                                   font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("SATreePlotUI"),
                               style = paste("height:100%; overflow-y:",
                                             "scroll;overflow-x: scroll;")
                        )
                    ),
                ),
            )
        } else if (!is.na(strtoi(sidebar_menu[[1]])) &&
                   sidebar_menu[[2]] == "Sanger" &&
                   sidebar_menu[[3]] == "Contig" &&
                   sidebar_menu[[4]] == "Overview" &&
                   sidebar_menu[[5]] == "Page") {
            contigIndex <- strtoi(sidebar_menu[[1]])
            log_info(">>>>>>>> Inside '", input$sidebar_menu, "'")
            shinyjs::disable("closeUI")
            shinyjs::disable("recalculateButton")
            log_info("@@@@@@@ 'Reactive button' has been clicked")
            log_info("In the main")
            log_info("######## Start recalculating contig")
            CSResult<-
                calculateContigSeq (
                    SangerAlignment@contigList[[contigIndex]]@inputSource,
                    SangerAlignment@contigList[[contigIndex]]@forwardReadList,
                    SangerAlignment@contigList[[contigIndex]]@reverseReadList,
                    SangerAlignment@contigList[[contigIndex]]@refAminoAcidSeq,
                    SangerAlignment@contigList[[contigIndex]]@minFractionCall,
                    SangerAlignment@contigList[[contigIndex]]@maxFractionLost,
                    SangerAlignment@contigList[[contigIndex]]@geneticCode,
                    SangerAlignment@contigList[[contigIndex]]@acceptStopCodons,
                    SangerAlignment@contigList[[contigIndex]]@readingFrame)
            SangerAlignment@contigList[[contigIndex]]@contigSeq <<- CSResult$consensusGapfree
            SangerAlignment@contigList[[contigIndex]]@differencesDF <<- CSResult$diffsDf
            SangerAlignment@contigList[[contigIndex]]@alignment <<- CSResult$aln2
            SangerAlignment@contigList[[contigIndex]]@distanceMatrix <<- CSResult$dist
            SangerAlignment@contigList[[contigIndex]]@dendrogram <<- CSResult$dend
            SangerAlignment@contigList[[contigIndex]]@indelsDF <<- CSResult$indels
            SangerAlignment@contigList[[contigIndex]]@stopCodonsDF <<- CSResult$stopsDf
            SangerAlignment@contigList[[contigIndex]]@secondaryPeakDF <<- CSResult$spDf
            contigParam[["contigSeq"]] <<- SangerAlignment@contigList[[contigIndex]]@contigSeq
            contigParam[["differencesDF"]] <<- SangerAlignment@contigList[[contigIndex]]@differencesDF
            contigParam[["alignment"]] <<- as.character(SangerAlignment@contigList[[contigIndex]]@alignment)
            contigParam[["distanceMatrix"]] <<-SangerAlignment@contigList[[contigIndex]]@distanceMatrix
            contigParam[["dendrogram"]] <<- SangerAlignment@contigList[[contigIndex]]@dendrogram
            contigParam[["indelsDF"]] <<- SangerAlignment@contigList[[contigIndex]]@indelsDF
            contigParam[["stopCodonsDF"]] <<- SangerAlignment@contigList[[contigIndex]]@stopCodonsDF
            contigParam[["secondaryPeakDF"]] <<- SangerAlignment@contigList[[contigIndex]]@secondaryPeakDF
            log_info("######## Finish recalculating contig")
            shinyjs::enable("recalculateButton")
            shinyjs::enable("closeUI")
            fluidRow(
                useShinyjs(),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Basic Information: "),
                                   style = "font-size: 26px;
                                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(
                        style = ("border-top: 0.2px hidden #A9A9A9;")),
                    fluidRow(
                        column(width = 12,
                               actionBttn("recalculateButton",
                                          "Re-calculate Contig",
                                          icon = icon("calculator"),
                                          style = "simple",
                                          color = "danger",
                                          block = TRUE, size = "lg")
                        ),
                        column(12,
                               tags$hr(
                                   style = ("border-top: 2px
                                                    hidden #A9A9A9;")),
                        ),
                        column(12,
                               column(3,
                                      tags$p(tagList(icon("caret-right"),
                                                     "Output Directory:"),
                                             style = "font-size: 20px;
                                                           font-weight: bold;"),
                               ),
                               column(9,
                                      h4(shinyDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Raw ABI Parent Directory: "),
                                          style = "font-size: 20px;
                                                           font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerAlignment@
                                             parentDirectory),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Consenesus Read Name:"),
                                          style = "font-size: 20px;
                                                           font-weight: bold;"),
                               ),
                               column(9,
                                      h4(contigParam[["contigName"]]),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Trimming Method: "),
                                          style = "font-size: 20px;
                                                           font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SATrimmingMethodName),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Forward Suffix RegExp:"),
                                          style = "font-size: 20px;
                                                           font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerAlignment@
                                             suffixForwardRegExp),
                               )
                               
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Forward Read Number:"),
                                          style = "font-size: 20px;
                                                           font-weight: bold;"),
                               ),
                               column(9,
                                      h4(length(SangerAlignment@
                                                    contigList[[contigIndex]]@
                                                    forwardReadList)),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Reverse Suffix RegExp:"),
                                          style = "font-size: 20px;
                                                           font-weight: bold;"),
                               ),
                               column(9,
                                      h4(SangerAlignment@
                                             suffixReverseRegExp),
                               )
                        ),
                        column(12,
                               column(3,
                                      tags$p(
                                          tagList(icon("caret-right"),
                                                  "Reverse Read Number:"),
                                          style = "font-size: 20px;
                                                                  font-weight: bold;"),
                               ),
                               column(9,
                                      h4(length(SangerAlignment@
                                                    contigList[[contigIndex]]@
                                                    reverseReadList)),
                               )
                        ),
                    ),
                    ################################################
                    #### Add this after having reference sample ####
                    ################################################
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Contig Parameters",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(4,
                               uiOutput("SAMinReadsNum") ,
                        ),
                        column(4,
                               uiOutput("SAMinReadLength")  ,
                        ),
                        column(4,
                               uiOutput("SAMinFractionCall") ,
                        ),
                        column(4,
                               uiOutput("SAMaxFractionLost") ,
                        ),
                        column(4,
                               uiOutput("SAAcceptStopCodons") ,
                        ),
                        column(4,
                               uiOutput("SAReadingFrame") ,
                        ),
                    ),
                    box(title = tags$p("Genetic Code Data Frame",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 2,
                               tags$p("Tri-nucleotide:",
                                      style = "font-size: 15px;
                                                           font-weight: bold;"),
                               tags$p("Amino Acid : ",
                                      style = "font-size: 15px;
                                                           font-weight: bold;"),
                               tags$p("('*' : stop codon) ",
                                      style = "font-size: 12px;
                                                           font-weight: italic;"),
                        ),
                        column(width = 10,
                               excelOutput("geneticCodeDF",
                                           width = "100%",
                                           height = "50"),
                               style = paste("height:100%; ",
                                             "overflow-y: hidden;",
                                             "overflow-x: scroll;")
                        ),
                    ),
                    uiOutput("SArefAminoAcidSeq") ,
                ),
                
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Consensus Read Results: "),
                                   style = "font-size: 26px;
                                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style =("border-top: 4px hidden #A9A9A9;")),
                    box(title = tags$p("Alignment",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               htmlOutput("contigAlignmentHTML"),
                        ),
                    ),
                    box(title = tags$p("Differences Data frame",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               # uiOutput("SADifferencesDFUI"),
                               dataTableOutput("SADifferencesDF"),
                               style =
                                   paste("height:100%; overflow-y:",
                                         "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Dendrogram",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               plotOutput("dendrogramPlot"),
                               style =
                                   paste("height:100%; overflow-y:",
                                         "scroll;overflow-x: scroll;")
                        ),
                        column(width = 12,
                               tags$hr(
                                   style =("border-top: 4px hidden #A9A9A9;")),
                        ),
                        column(width = 12,
                               dataTableOutput("dendrogramDF"),
                               style =
                                   paste("height:100%; overflow-y:",
                                         "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Samples Distance",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               # plot()
                               # uiOutput("SADistanceMatrixPlotUI"),
                               plotlyOutput("SADistanceMatrixPlot"),
                               style =
                                   paste("height:100%; overflow-y:",
                                         "scroll;overflow-x: scroll;")
                        ),
                        column(width = 12,
                               tags$hr(
                                   style =("border-top: 4px hidden #A9A9A9;")),
                        ),
                        column(width = 12,
                               # uiOutput("SADistanceMatrixUI"),
                               dataTableOutput("SADistanceMatrix"),
                               style =
                                   paste("height:100%; overflow-y:",
                                         "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Indels Data frame",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               # uiOutput("SAIndelsDFUI"),
                               dataTableOutput("SAIndelsDF"),
                               style =
                                   paste("height:100%; overflow-y:",
                                         "scroll;overflow-x: scroll;")
                        )
                    ),
                    box(title = tags$p("Stop Codons Data frame",
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               # uiOutput("SAStopCodonsDFUI"),
                               dataTableOutput("SAStopCodonsDF"),
                               style =
                                   paste("height:100%; overflow-y:",
                                         "scroll;overflow-x: scroll;")
                        )
                    )
                )
            )
            
        } else if (!is.na(strtoi(sidebar_menu[[1]])) &&
                   sidebar_menu[[2]] == "Contig" &&
                   sidebar_menu[[3]] == "-" &&
                   !is.na(strtoi(sidebar_menu[[4]])) &&
                   (sidebar_menu[[5]] == "Forward" ||
                    sidebar_menu[[5]] == "Reverse") &&
                   sidebar_menu[[6]] == "Read") {
            contigIndex <- strtoi(sidebar_menu[[1]])
            readIndex <- strtoi(sidebar_menu[[4]])
            directionParam <- sidebar_menu[[5]]
            if (directionParam == "Forward") {
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primaryAASeqS3)
                trimmedParam[["M1TrimmingCutoff"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    M1TrimmingCutoff
                trimmedParam[["M2CutoffQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    M2CutoffQualityScore
                trimmedParam[["M2SlidingWindowSize"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    M2SlidingWindowSize
                
                trimmedRV[["rawSeqLength"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerAlignment@
                              contigList[[contigIndex]]@
                              forwardReadList[[readIndex]]@QualityReport@
                              remainingRatio * 100, 2)
                
                ChromatogramParam[["baseNumPerRow"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    baseNumPerRow
                ChromatogramParam[["heightPerRow"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    heightPerRow
                ChromatogramParam[["signalRatioCutoff"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    signalRatioCutoff
                ChromatogramParam[["showTrimmed"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    showTrimmed
                
                SangerReadBFN <-
                    basename(SangerAlignment@
                                 contigList[[contigIndex]]@
                                 forwardReadList[[readIndex]]@readFileName)
                SangerReadAFN <-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@readFileName
            } else if (directionParam == "Reverse") {
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerAlignment@
                                     contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primaryAASeqS3)
                
                trimmedParam[["M1TrimmingCutoff"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    M1TrimmingCutoff
                trimmedParam[["M2CutoffQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    M2CutoffQualityScore
                trimmedParam[["M2SlidingWindowSize"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    M2SlidingWindowSize
                
                trimmedRV[["rawSeqLength"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerAlignment@
                              contigList[[contigIndex]]@
                              reverseReadList[[readIndex]]@QualityReport@
                              remainingRatio * 100, 2)
                
                ChromatogramParam[["baseNumPerRow"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    baseNumPerRow
                ChromatogramParam[["heightPerRow"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    heightPerRow
                ChromatogramParam[["signalRatioCutoff"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    signalRatioCutoff
                ChromatogramParam[["showTrimmed"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    showTrimmed
                
                SangerReadBFN <-
                    basename(SangerAlignment@
                                 contigList[[contigIndex]]@
                                 reverseReadList[[readIndex]]@readFileName)
                SangerReadAFN <-
                    SangerAlignment@
                    contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@readFileName
            }
            log_info(">>>>>>>> Inside '", input$sidebar_menu, "'")
            h1(input$sidebar_menu)
            fluidRow(
                useShinyjs(),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Raw File: "),
                                   style = "font-size: 26px;
                                                             font-weight: bold;"),
                    solidHeader = TRUE,
                    status = "success", width = 12,
                    h1(SangerReadBFN),
                    tags$h5(paste("( full path:",
                                  SangerReadAFN,
                                  ")"), style = "font-style:italic")),
                box(title = tags$p(
                    tagList(icon("dot-circle"),
                            "Primary, Secondary DNA Sequences & Amino Acid Sequence (Before Trimming):"),
                    style = "font-size: 26px;
                                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    column(width = 12,
                           tags$p(tagList(icon("bars"),
                                          "Primary Sequence"),
                                  style = "font-size: 22px;
                                                           font-weight: bold;"),
                           excelOutput("primarySeqDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Secondary Sequence"),
                                  style = "font-size: 22px;
                                                           font-weight: bold;"),
                           excelOutput("secondSeqDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Quality Phred Score"),
                                  style = "font-size: 22px;
                                                           font-weight: bold;"),
                           excelOutput("qualityScoreDF",
                                       width = "100%", height = "50"),
                           style = paste("overflow-y: hidden;",
                                         "overflow-x: scroll;")
                    ),
                ),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Quality Report: "),
                                   style = "font-size: 26px;
                                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(
                        style = ("border-top: 4px hidden #A9A9A9;")),
                    box(title =
                            tags$p(tagList(icon("arrow-circle-right"),
                                           "Trimming Parameters Input"),
                                   style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        fluidRow(
                            column(width = 12,
                                   uiOutput("TrimmingMethodSelectionOutput"),
                            ),
                        ),
                        column(width = 12,
                               uiOutput("TrimmingMethodUI"),
                        ),
                        actionBttn("startTrimmingButton",
                                   "Apply Trimming Parameters",
                                   style = "simple", color = "danger",
                                   block = TRUE, size = "lg")
                    ),
                    box(title = tags$p(tagList(icon("arrow-circle-left"),
                                               "Trimmed Result Output"),
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        fluidRow(
                            box(title = tags$p("Before Trimming",
                                               style = "font-size: 21px;
                                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       column(4,
                                              uiOutput("rawSeqLength") ,
                                       ),
                                       column(4,
                                              uiOutput("rawMeanQualityScore"),
                                       ),
                                       column(4,
                                              uiOutput("rawMinQualityScore"),
                                       ),
                                ),
                            ),
                        ),
                        fluidRow(
                            box(title = tags$p("After Trimming",
                                               style = "font-size: 21px;
                                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                column(width = 12,
                                       column(4,
                                              uiOutput("trimmedSeqLength"),
                                       ),
                                       column(4,
                                              uiOutput("trimmedMeanQualityScore"),
                                       ),
                                       column(4,
                                              uiOutput("trimmedMinQualityScore"),
                                       ),
                                ),
                                
                                column(width = 12,
                                       column(4,
                                              uiOutput("trimmedStartPos") ,
                                       ),
                                       column(4,
                                              uiOutput("trimmedFinishPos") ,
                                       ),
                                       column(4,
                                              uiOutput("remainingRatio") ,
                                       )
                                ),
                            ),
                        ),
                        tags$hr(
                            style = ("border-top: 6px double #A9A9A9;")),
                        fluidRow(
                            box(title = tags$p("Quality Trimming Plot",
                                               style = "font-size: 21px;
                                                           font-weight: bold;"),
                                collapsible = TRUE,
                                status = "success", width = 12,
                                plotlyOutput("qualityQualityBasePlot") %>%
                                    withSpinner()),
                        ),
                    ),
                ),
                box(title =
                        tags$p(tagList(icon("dot-circle"),
                                       "DNA Sequence
                                       (After Trimming):"),
                               style = "font-size: 26px;
                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    column(width = 12,
                           tags$p(tagList(icon("bars"),
                                          "Primary Sequence"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("primarySeqTrimmedDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Secondary Sequence"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("secondSeqTrimmedDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "Quality Phred Score"),
                                  style = "font-size: 22px;
                                           font-weight: bold;"),
                           excelOutput("qualityScoreTrimmedDF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "AA Sequence 1"),
                                  style = "font-size: 22px;
                                                           font-weight: bold;"),
                           excelOutput("PrimAASeqS1DF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "AA Sequence 2"),
                                  style = "font-size: 22px;
                                                           font-weight: bold;"),
                           excelOutput("PrimAASeqS2DF",
                                       width = "100%", height = "50"),
                           tags$br(),
                           tags$br(),
                           tags$p(tagList(icon("bars"),
                                          "AA Sequence 3"),
                                  style = "font-size: 22px;
                                                           font-weight: bold;"),
                           excelOutput("PrimAASeqS3DF",
                                       width = "100%", height = "50"),
                           style = paste("overflow-y: hidden;",
                                         "overflow-x: scroll;")
                    )
                ),
                box(title = tags$p(tagList(icon("dot-circle"),
                                           "Chromatogram: "),
                                   style = "font-size: 26px;
                                                           font-weight: bold;"),
                    solidHeader = TRUE, collapsible = TRUE,
                    status = "success", width = 12,
                    tags$hr(style = ("border-top: 4px hidden #A9A9A9;")),
                    
                    box(title = tags$p(tagList(icon("arrow-circle-right"),
                                               "Chromatogram Input"),
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(3,
                               sliderInput("ChromatogramBasePerRow",
                                           label = h4("Base Number Per Row"),
                                           min = 5,
                                           max = 200,
                                           value = ChromatogramParam[["baseNumPerRow"]]),
                               sliderInput("ChromatogramHeightPerRow",
                                           label = h4("Height Per Row"),
                                           min = 50,
                                           max = 600,
                                           value =ChromatogramParam[["heightPerRow"]]),
                        ),
                        column(3,
                               numericInput(
                                   "ChromatogramSignalRatioCutoff",
                                   h3("Signal Ratio Cutoff"),
                                   value = ChromatogramParam[["signalRatioCutoff"]]),
                               checkboxInput(
                                   "ChromatogramCheckShowTrimmed",
                                   "Show trimmed region",
                                   value = ChromatogramParam[["showTrimmed"]])
                        ),
                        column(3,
                               uiOutput("ChromatogramtrimmedStartPos"),
                        ),
                        column(3,
                               uiOutput("ChromatogramtrimmedFinishPos"),
                        ),
                        actionBttn("saveChromatogramParam",
                                   "Apply Chromatogram Parameters",
                                   style = "simple", color = "danger",
                                   block = TRUE, size = "lg")
                    ),
                    box(title = tags$p(tagList(icon("arrow-circle-left"),
                                               "Chromatogram Output"),
                                       style = "font-size: 24px;
                                                           font-weight: bold;"),
                        collapsible = TRUE,
                        status = "success", width = 12,
                        column(width = 12,
                               uiOutput("chromatogramUIOutput"),
                        )
                    ),
                )
            )
        }
    })
    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    observeEvent(input$sidebar_menu, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        if (sidebar_menu[[1]] == "Contigs" && 
            sidebar_menu[[2]] == "Alignment" && 
            sidebar_menu[[3]] == "Overview" && 
            sidebar_menu[[4]] == "Page" && 
            sidebar_menu[[5]] == "_") {
            shinyjs::html("rightHeader", "SangerAlignment Overview Page")
        } else if (!is.na(strtoi(sidebar_menu[[1]])) &&
                   sidebar_menu[[2]] == "Sanger" &&
                   sidebar_menu[[3]] == "Contig" &&
                   sidebar_menu[[4]] == "Overview" &&
                   sidebar_menu[[5]] == "Page") {
            contigIndex <- strtoi(sidebar_menu[[1]])
            shinyjs::html("rightHeader", paste(contigIndex, "SangerContig Overview Page"))
        } else if (!is.na(strtoi(sidebar_menu[[1]])) &&
                   sidebar_menu[[2]] == "Contig" &&
                   sidebar_menu[[3]] == "-" &&
                   !is.na(strtoi(sidebar_menu[[4]])) &&
                   (sidebar_menu[[5]] == "Forward" ||
                    sidebar_menu[[5]] == "Reverse") &&
                   sidebar_menu[[6]] == "Read") {
            contigIndex <- strtoi(sidebar_menu[[1]])
            readIndex <- strtoi(sidebar_menu[[4]])
            directionParam <- sidebar_menu[[5]]
            shinyjs::html("rightHeader",
                          paste(contigIndex, "SangerContig -",
                                readIndex, directionParam,
                                "SangerRead", "Page"))
        }
    })

    ### ------------------------------------------------------------------------
    ### observeEvent: Adding dynamic sidebar
    ### ------------------------------------------------------------------------
    # !!! Fix !!!
    dynamicMenuSideBarSA(input, output, session, SangerAlignmentParam)
    ############################################################################
    ### All other features (dynamic header / button save / button close)
    ############################################################################
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Close UI
    ### ------------------------------------------------------------------------
    observeEvent(input$closeUI, {
        log_info("@@@@@@@ 'close button' has been clicked")
        btn <- input$closeUI
        stopApp()
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button Save S4 object
    ### ------------------------------------------------------------------------
    observeEvent(input$saveS4, {
        shinyjs::disable("closeUI")
        shinyjs::disable("saveS4")
        log_info("@@@@@@@ 'save button' has been clicked")
        newS4Object <- file.path(shinyDirectory,
                                 "SangerAlignment.Rda")
        showNotification(
            ui = column(12,
                        tags$p(tagList(icon("dot-circle"),
                                       "New S4 object is store as: "),
                               style = "font-size: 28px;
                                   font-weight: bold;"),
                        tags$p(paste0("'", newS4Object, "'"),
                               style = "font-size: 26px;
                                   font-style: italic"),
                        tags$br(),
                        tags$p(">> Run",
                               style = "font-size: 18px;
                                   font-style: italic"),
                        tags$p(paste0("readRDS(\"",newS4Object,
                                      "\")"),
                               style = "font-size: 18px;
                                      font-weight: bold;
                                      font-family: Courier, Monospace;"),
                        tags$p("   to load saved S4 object into R
                                     environment",
                               style = "font-size: 18px;
                                      font-style: italic")),
            closeButton = TRUE, id = "saveNotification",
            type = "message", duration = 10)
        ### --------------------------------------------------------------------
        ### Save SangerContig quality S4 object
        ### --------------------------------------------------------------------
        saveRDS(SangerAlignment, file=newS4Object)
        log_info("New S4 object is store as: ", newS4Object)
        NEW_SANGER_ALIGNED_CONSENSUS_READ_SET <<- readRDS(file=newS4Object)
        shinyjs::enable("saveS4")
        shinyjs::enable("closeUI")
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button contigs alignment re-calculating (SA) UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButtonSA, {
        shinyjs::disable("closeUI")
        shinyjs::disable("recalculateButtonSA")
        log_info("######## Reactive button clicked !!!")
        log_info("######## Start re-aligning contigs")
        if (input$sidebar_menu == "Contigs Alignment Overview Page _") {
            CSSetResult <-
                alignContigs (SangerAlignment@contigList,
                                     SangerAlignment@geneticCode,
                                     SangerAlignment@refAminoAcidSeq,
                                     SangerAlignment@minFractionCallSA,
                                     SangerAlignment@maxFractionLostSA,
                                     1)
            SangerAlignment@contigsConsensus <<- CSSetResult$consensus
            SangerAlignment@contigsAlignment <<- CSSetResult$aln
            SangerAlignment@contigsTree <<- CSSetResult$aln.tree
            sangerAlignmentParam[["contigsConsensus"]] <<- SangerAlignment@contigsConsensus
            sangerAlignmentParam[["contigsAlignment"]] <<- as.character(SangerAlignment@contigsAlignment)
            sangerAlignmentParam[["contigsTree"]] <<- SangerAlignment@contigsTree
            log_info("######## Finish contigs re-alignment")
        }
        shinyjs::enable("recalculateButtonSA")
        shinyjs::enable("closeUI")
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button SangerContig re-calculating (SA) UI
    ### ------------------------------------------------------------------------
    observeEvent(input$recalculateButton, {
        shinyjs::disable("closeUI")
        shinyjs::disable("recalculateButton")
        log_info("@@@@@@@ 'Reactive button' has been clicked")
        log_info("######## Start recalculating contig")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            CSResult<-
                calculateContigSeq (
                    SangerAlignment@contigList[[contigIndex]]@inputSource,
                    SangerAlignment@contigList[[contigIndex]]@forwardReadList,
                    SangerAlignment@contigList[[contigIndex]]@reverseReadList,
                    SangerAlignment@contigList[[contigIndex]]@refAminoAcidSeq,
                    SangerAlignment@contigList[[contigIndex]]@minFractionCall,
                    SangerAlignment@contigList[[contigIndex]]@maxFractionLost,
                    SangerAlignment@contigList[[contigIndex]]@geneticCode,
                    SangerAlignment@contigList[[contigIndex]]@acceptStopCodons,
                    SangerAlignment@contigList[[contigIndex]]@readingFrame)
            SangerAlignment@contigList[[contigIndex]]@contigSeq <<- CSResult$consensusGapfree
            SangerAlignment@contigList[[contigIndex]]@differencesDF <<- CSResult$diffsDf
            SangerAlignment@contigList[[contigIndex]]@alignment <<- CSResult$aln2
            SangerAlignment@contigList[[contigIndex]]@distanceMatrix <<- CSResult$dist
            SangerAlignment@contigList[[contigIndex]]@dendrogram <<- CSResult$dend
            SangerAlignment@contigList[[contigIndex]]@indelsDF <<- CSResult$indels
            SangerAlignment@contigList[[contigIndex]]@stopCodonsDF <<- CSResult$stopsDf
            SangerAlignment@contigList[[contigIndex]]@secondaryPeakDF <<- CSResult$spDf
            contigParam[["contigSeq"]] <<- SangerAlignment@contigList[[contigIndex]]@contigSeq
            contigParam[["differencesDF"]] <<- SangerAlignment@contigList[[contigIndex]]@differencesDF
            contigParam[["alignment"]] <<- as.character(SangerAlignment@contigList[[contigIndex]]@alignment)
            contigParam[["distanceMatrix"]] <<-SangerAlignment@contigList[[contigIndex]]@distanceMatrix
            contigParam[["dendrogram"]] <<- SangerAlignment@contigList[[contigIndex]]@dendrogram
            contigParam[["indelsDF"]] <<- SangerAlignment@contigList[[contigIndex]]@indelsDF
            contigParam[["stopCodonsDF"]] <<- SangerAlignment@contigList[[contigIndex]]@stopCodonsDF
            contigParam[["secondaryPeakDF"]] <<- SangerAlignment@contigList[[contigIndex]]@secondaryPeakDF
        }
        log_info("######## Finish recalculating contig")
        shinyjs::enable("recalculateButton")
        shinyjs::enable("closeUI")
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button apply trimming parameters
    ### ------------------------------------------------------------------------
    observeEvent(input$startTrimmingButton, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            if (SangerAlignment@trimmingMethodSA == "M1") {
                if (!is.na(as.numeric(input$M1TrimmingCutoffText)) &&
                    as.numeric(input$M1TrimmingCutoffText) > 0 &&
                    as.numeric(input$M1TrimmingCutoffText) <= 1) {
                    inputM1TrimmingCutoffText <- input$M1TrimmingCutoffText
                } else {
                    inputM1TrimmingCutoffText <- 0.0001
                }
                if (directionParam == "Forward") {
                    trimmingPos <-
                        M1inside_calculate_trimming(
                            SangerAlignment@contigList[[contigIndex]]@
                                forwardReadList[[readIndex]]@QualityReport@
                                qualityPhredScores,
                            SangerAlignment@contigList[[contigIndex]]@
                                forwardReadList[[readIndex]]@QualityReport@
                                qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerAlignment@contigList[[contigIndex]]@
                        forwardReadList[[readIndex]]@QualityReport@
                        M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerAlignment@contigList[[contigIndex]]@
                        forwardReadList[[readIndex]]@QualityReport@M1TrimmingCutoff
                } else if (directionParam == "Reverse") {
                    trimmingPos <-
                        M1inside_calculate_trimming(
                            SangerAlignment@contigList[[contigIndex]]@
                                reverseReadList[[readIndex]]@QualityReport@
                                qualityPhredScores,
                            SangerAlignment@contigList[[contigIndex]]@
                                reverseReadList[[readIndex]]@QualityReport@
                                qualityBaseScores,
                            as.numeric(inputM1TrimmingCutoffText))
                    SangerAlignment@contigList[[contigIndex]]@
                        reverseReadList[[readIndex]]@QualityReport@
                        M1TrimmingCutoff <<- as.numeric(inputM1TrimmingCutoffText)
                    trimmedParam[["M1TrimmingCutoff"]] <<-
                        SangerAlignment@contigList[[contigIndex]]@
                        reverseReadList[[readIndex]]@QualityReport@M1TrimmingCutoff
                }
            } else if (SangerAlignment@trimmingMethodSA == "M2") {
                if (!is.na(strtoi(input$M2CutoffQualityScoreText)) &&
                    strtoi(input$M2CutoffQualityScoreText) > 0 &&
                    strtoi(input$M2CutoffQualityScoreText) <= 60 &&
                    strtoi(input$M2CutoffQualityScoreText) %% 1 ==0) {
                    inputM2CutoffQualityScoreText <- input$M2CutoffQualityScoreText
                } else {
                    inputM2CutoffQualityScoreText <- 20
                }
                if (!is.na(strtoi(input$M2SlidingWindowSizeText)) &&
                    strtoi(input$M2SlidingWindowSizeText) > 0 &&
                    strtoi(input$M2SlidingWindowSizeText) <= 40 &&
                    strtoi(input$M2SlidingWindowSizeText) %% 1 ==0) {
                    inputM2SlidingWindowSizeText <- input$M2SlidingWindowSizeText
                } else {
                    inputM2SlidingWindowSizeText <- 10
                }
                if (directionParam == "Forward") {
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerAlignment@contigList[[contigIndex]]@
                                forwardReadList[[readIndex]]@QualityReport@
                                qualityPhredScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText)
                        )
                    SangerAlignment@contigList[[contigIndex]]@
                        forwardReadList[[readIndex]]@QualityReport@
                        M2CutoffQualityScore <<- strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerAlignment@contigList[[contigIndex]]@
                        forwardReadList[[readIndex]]@QualityReport@
                        M2CutoffQualityScore
                    SangerAlignment@contigList[[contigIndex]]@
                        forwardReadList[[readIndex]]@QualityReport@
                        M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerAlignment@contigList[[contigIndex]]@
                        forwardReadList[[readIndex]]@QualityReport@
                        M2SlidingWindowSize
                } else if ("Reverse") {
                    trimmingPos <-
                        M2inside_calculate_trimming(
                            SangerAlignment@contigList[[contigIndex]]@
                                reverseReadList[[readIndex]]@QualityReport@
                                qualityPhredScores,
                            strtoi(inputM2CutoffQualityScoreText),
                            strtoi(inputM2SlidingWindowSizeText)
                        )
                    SangerAlignment@contigList[[contigIndex]]@
                        reverseReadList[[readIndex]]@QualityReport@
                        M2CutoffQualityScore <<- strtoi(inputM2CutoffQualityScoreText)
                    trimmedParam[["M2CutoffQualityScore"]] <<-
                        SangerAlignment@contigList[[contigIndex]]@
                        reverseReadList[[readIndex]]@QualityReport@
                        M2CutoffQualityScore
                    SangerAlignment@contigList[[contigIndex]]@
                        reverseReadList[[readIndex]]@QualityReport@
                        M2SlidingWindowSize <<- strtoi(inputM2SlidingWindowSizeText)
                    trimmedParam[["M2SlidingWindowSize"]] <<-
                        SangerAlignment@contigList[[contigIndex]]@
                        reverseReadList[[readIndex]]@QualityReport@
                        M2SlidingWindowSize
                }
            }
            if (directionParam == "Forward") {
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]

                trimmedRV[["rawSeqLength"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerAlignment@contigList[[contigIndex]]@
                              forwardReadList[[readIndex]]@QualityReport@
                              remainingRatio * 100, 2)
            } else if (directionParam == "Reverse") {
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    rawSeqLength <<- trimmingPos[["rawSeqLength"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore <<- trimmingPos[["rawMeanQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    rawMinQualityScore <<- trimmingPos[["rawMinQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedStartPos <<- trimmingPos[["trimmedStartPos"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedFinishPos <<- trimmingPos[["trimmedFinishPos"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedSeqLength <<- trimmingPos[["trimmedSeqLength"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore <<- trimmingPos[["trimmedMeanQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore <<- trimmingPos[["trimmedMinQualityScore"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    remainingRatio <<- trimmingPos[["remainingRatio"]]

                trimmedRV[["rawSeqLength"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@rawSeqLength
                trimmedRV[["rawMeanQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    rawMeanQualityScore
                trimmedRV[["rawMinQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@rawMinQualityScore
                trimmedRV[["trimmedStartPos"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@trimmedStartPos
                trimmedRV[["trimmedFinishPos"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@trimmedFinishPos
                trimmedRV[["trimmedSeqLength"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@trimmedSeqLength
                trimmedRV[["trimmedMeanQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedMeanQualityScore
                trimmedRV[["trimmedMinQualityScore"]] <<-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@
                    trimmedMinQualityScore
                trimmedRV[["remainingRatio"]] <<-
                    round(SangerAlignment@contigList[[contigIndex]]@
                              reverseReadList[[readIndex]]@QualityReport@
                              remainingRatio * 100, 2)
            }
        }
    })
    ### ------------------------------------------------------------------------
    ### observeEvent: Button chromatogram parameters re-calculating UI
    ### ------------------------------------------------------------------------
    observeEvent(input$saveChromatogramParam, {
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        ### ------------------------------------------------------------
        ### Update ChromatogramBasePerRow
        ### ------------------------------------------------------------
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            if (directionParam == "Forward") {
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@ChromatogramParam@
                    showTrimmed <<- input$ChromatogramCheckShowTrimmed
            } else if (directionParam == "Reverse") {
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    baseNumPerRow <<- input$ChromatogramBasePerRow
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    heightPerRow <<- input$ChromatogramHeightPerRow
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    signalRatioCutoff <<- input$ChromatogramSignalRatioCutoff
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@ChromatogramParam@
                    showTrimmed <<- input$ChromatogramCheckShowTrimmed
            }
            ### ----------------------------------------------------------------
            ### Save 'ChromatogramParam' dynamic value
            ### ----------------------------------------------------------------
            ChromatogramParam[["baseNumPerRow"]] <<-
                input$ChromatogramBasePerRow
            ChromatogramParam[["heightPerRow"]] <<-
                input$ChromatogramHeightPerRow
            ChromatogramParam[["signalRatioCutoff"]] <<-
                input$ChromatogramSignalRatioCutoff
            ChromatogramParam[["showTrimmed"]] <<-
                input$ChromatogramCheckShowTrimmed
        }
    })



    ############################################################################
    ### SangerAlignment
    ############################################################################
    ### ------------------------------------------------------------------------
    ### valuebox
    ### ------------------------------------------------------------------------
    valueBoxSAMinFractionCallSA(input, output,
                               SangerAlignment@minFractionCallSA, session)
    valueBoxSAMaxFractionLostSA(input, output,
                               SangerAlignment@maxFractionLostSA, session)
    ### ------------------------------------------------------------------------
    ### Alignment
    ### ------------------------------------------------------------------------
    output$sangerAlignmentAlignmentHTML<-renderUI({
        if (input$sidebar_menu == "Contigs Alignment Overview Page _") {
            browseSeqHTML <-
                file.path(shinyDirectory,
                          "_SangerAlignment_BrowseSeqs.html")
            if (!identical(sangerAlignmentParam[["contigsAlignment"]], character(0))) {
                BrowseSeqs(DNAStringSet(sangerAlignmentParam[["contigsAlignment"]]) ,
                           openURL=FALSE, htmlFile=browseSeqHTML)
            }
            log_info("file.exists(browseSeqHTML): ", file.exists(browseSeqHTML))
            log_info("(browseSeqHTML): ", (browseSeqHTML))
            if (file.exists(browseSeqHTML)) {
                column(width = 12,
                       includeHTML(browseSeqHTML),
                       style = paste("height:100%; ",
                                     "overflow-y: hidden;",
                                     "overflow-x: scroll;")
                )                
            } else {
                h4("*** The number of contigs is less than 2.
                   sCannot create 'BrowseSeqs' alignment. ***",
                   style="font-weight: bold; font-style: italic;")
            }
        }
    })
    ### ------------------------------------------------------------------------
    ### Contigs Tree
    ### ------------------------------------------------------------------------
    output$SATreePlotUI <- renderUI({
        if (input$sidebar_menu == "Contigs Alignment Overview Page _") {
            if (sangerAlignmentParam[["contigsTree"]]$tip.label != '' && 
                !is.null(sangerAlignmentParam[["contigsTree"]]$tip.label)) {
                plotOutput("SATreePlot")
            } else {
                h4("*** The number of contigs is less than 3 or quality of reads are too low.
                   'Contigs Tree' cannot be created. ***",
                   style="font-weight: bold; font-style: italic;")
            }
        }
    })
    output$SATreePlot <- renderPlot({
        if (input$sidebar_menu == "Contigs Alignment Overview Page _") {
            plot(sangerAlignmentParam[["contigsTree"]])
        }
    })



    ############################################################################
    ### SangerContig / SangerAlignment Share
    ############################################################################
    ### ------------------------------------------------------------------------
    ### genetic code
    ### ------------------------------------------------------------------------
    output$geneticCodeDF <- renderExcel({
        SAGeneticCode <- SangerAlignment@geneticCode
        suppressMessages(
            excelTable(data = t(data.frame(SAGeneticCode)),
                       defaultColWidth = 50, editable = FALSE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })
    ### ------------------------------------------------------------------------
    ### refAminoAcidSeq
    ### ------------------------------------------------------------------------
    output$SArefAminoAcidSeq <- renderUI({
        if (SangerAlignment@refAminoAcidSeq == "") {
            box(title = tags$p("Reference Amino Acid Sequence",
                               style = "font-size: 24px;
                                        font-weight: bold;"),
                collapsible = TRUE,
                status = "success", width = 12,
                column(width = 1),
                column(width = 11,
                       h4("Reference Amino Acid Sequence is not provided."))
            )
        } else {
            box(title = tags$p("Reference Amino Acid Sequence",
                               style = "font-size: 24px;
                                        font-weight: bold;"),
                collapsible = TRUE,
                status = "success", width = 12,
                column(width = 2,
                       tags$br(),
                       tags$p("AA Sequence:",
                              style = "font-size: 15px;
                                           font-weight: bold;"),
                ),
                column(width = 10,
                       excelOutput("SArefAminoAcidSeqDF",
                                   width = "100%", height = "50"),
                       style = paste("height:100%; ",
                                     "overflow-y: hidden;",
                                     "overflow-x: scroll;")
                ),
            )
        }
    })
    output$SArefAminoAcidSeqDF <- renderExcel({
        refAminoAcidSeqVec <-
            strsplit(SangerAlignment@refAminoAcidSeq, "")[[1]]
        names(refAminoAcidSeqVec) <- c(seq_len(length(refAminoAcidSeqVec)))
        suppressMessages(
            excelTable(data =
                           t(data.frame(refAminoAcidSeqVec)),
                       defaultColWidth = 50, editable = FALSE, rowResize = FALSE,
                       columnResize = FALSE, allowInsertRow = FALSE,
                       allowInsertColumn = FALSE, allowDeleteRow = FALSE,
                       allowDeleteColumn = FALSE, allowRenameColumn = FALSE)
        )
    })



    ############################################################################
    ### SangerContig
    ############################################################################
    valueBoxSAMinReadsNum (input, output, SangerAlignment, session)
    valueBoxSAMinReadLength (input, output, SangerAlignment, session)
    valueBoxSAMinFractionCall (input, output, SangerAlignment, session)
    valueBoxSAMaxFractionLost (input, output, SangerAlignment, session)
    valueBoxSAAcceptStopCodons (input, output, SangerAlignment, session)
    valueBoxSAReadingFrame (input, output, SangerAlignment, session)
    ### ------------------------------------------------------------------------
    ### Alignment
    ### ------------------------------------------------------------------------
    output$contigAlignmentHTML<-renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["contigName"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@contigName
            contigParam[["alignment"]] <-
                SangerAlignment@
                contigList[[contigIndex]]@alignment

            browseSeqHTML <-
                file.path(shinyDirectory, "BrowseSeqs_html",
                          paste0(sidebar_menu[[1]], "_",
                                 contigParam[["contigName"]],
                                 "_Alignment_BrowseSeqs.html"))
            if (!dir.exists(file.path(shinyDirectory, "BrowseSeqs_html"))) {
                dir.create(file.path(shinyDirectory, "BrowseSeqs_html"))
            }
            BrowseSeqs(DNAStringSet(contigParam[["alignment"]]),
                       openURL=FALSE, htmlFile=browseSeqHTML)
            column(width = 12,
                   includeHTML(browseSeqHTML),
                   style = paste("height:100%; ",
                                 "overflow-y: hidden;",
                                 "overflow-x: scroll;")
            )
        }
    })
    ### ------------------------------------------------------------------------
    ### difference
    ### ------------------------------------------------------------------------
    output$SADifferencesDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["differencesDF"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@differencesDF
            if (all(dim(contigParam[["differencesDF"]]) == c(0,0))) {
                h4("*** 'Differences' dataframe is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SADifferencesDF")
            }
        }
    })
    output$SADifferencesDF = renderDataTable({
        contigParam[["differencesDF"]]
    })
    ### ------------------------------------------------------------------------
    ### dendrogram
    ### ------------------------------------------------------------------------
    output$dendrogramPlot <- renderPlot({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["dendrogram"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@dendrogram
            plot(contigParam[["dendrogram"]][[2]])
            ggdendrogram(contigParam[["dendrogram"]][[2]], rotate = TRUE)
        }
    })
    output$dendrogramDF <- renderDataTable({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["dendrogram"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@dendrogram
            contigParam[["dendrogram"]][[1]]
        }
    })
    ### ------------------------------------------------------------------------
    ### distance
    ### ------------------------------------------------------------------------
    output$SADistanceMatrixPlotUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["distanceMatrix"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@distanceMatrix
            if (all(dim(contigParam[["distanceMatrix"]]) == c(0,0))) {
                h4("*** 'Distance' dataframe is empty. (Cannot plot)***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                plotlyOutput("SADistanceMatrixPlot")
            }
        }
    })
    output$SADistanceMatrixPlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            SangerReadBFN <-
                SangerAlignmentParam[[contigIndex]]$SangerReadBFN
            suppressPlotlyMessage(
                plot_ly(x = SangerReadBFN,
                        y = SangerReadBFN,
                        z = contigParam[["distanceMatrix"]],
                        colors = colorRamp(c("white", "#32a852")),
                        type = "heatmap")
            )
        }
    })
    output$SADistanceMatrixUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["distanceMatrix"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@distanceMatrix
            if (all(dim(contigParam[["distanceMatrix"]]) == c(0,0))) {
                h4("*** 'Distance' dataframe is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SADistanceMatrix")
            }
        }
    })
    output$SADistanceMatrix = renderDataTable({
        contigParam[["distanceMatrix"]]
    })
    ### ------------------------------------------------------------------------
    ### SAIndelsDF
    ### ------------------------------------------------------------------------
    output$SAIndelsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["indelsDF"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@indelsDF
            if (all(dim(contigParam[["indelsDF"]] ) == c(0,0))) {
                h4("*** 'Indels' data frame is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SAIndelsDF")
            }
        }
    })
    output$SAIndelsDF <- renderDataTable({
        contigParam[["indelsDF"]]
    })
    ### ------------------------------------------------------------------------
    ### SAStopCodons
    ### ------------------------------------------------------------------------
    output$SAStopCodonsDFUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        if (!is.na(contigIndex)) {
            contigParam[["stopCodonsDF"]] <<-
                SangerAlignment@
                contigList[[contigIndex]]@stopCodonsDF
            if (all(dim(contigParam[["stopCodonsDF"]]) == c(0,0))) {
                h4("*** 'Stop Codons' dataframe is empty. ***",
                   style="font-weight: bold; font-style: italic;")
            } else {
                dataTableOutput("SAStopCodonsDF")
            }
        }
    })
    output$SAStopCodonsDF <- renderDataTable({
        contigParam[["stopCodonsDF"]]
    })




    ############################################################################
    ### SangerRead (Function for singel read in contigSeq)
    ############################################################################
    output$primarySeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            primarySeqDisplay (sequenceParam, colors=colors)
        }
    })
    output$secondSeqDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            secondarySeqDisplay (sequenceParam, colors=colors)
        }
    })
    output$qualityScoreDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            if (directionParam == "Forward") {
                PhredScore <-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                PhredScore <-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityScoreDisplay (PhredScore)
        }
    })
    output$PrimAASeqS1DF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            PrimAASeqS1Display (sequenceParam)
        }
    })
    output$PrimAASeqS2DF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            PrimAASeqS2Display (sequenceParam)
        }
    })
    output$PrimAASeqS3DF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            PrimAASeqS3Display (sequenceParam)
        }
    })

    output$primarySeqTrimmedDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            primarySeqTrimmedDisplay (input, output, session,
                                      sequenceParam, trimmedRV, colors=colors)
        }
    })

    output$secondSeqTrimmedDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            secondSeqTrimmedDisplay (input, output, session,
                                     sequenceParam, trimmedRV, colors=colors)
        }
    })
    output$qualityScoreTrimmedDF <- renderExcel({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            if (directionParam == "Forward") {
                PhredScore <-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@
                    QualityReport@qualityPhredScores[
                        (trimmedRV[["trimmedStartPos"]]+1):
                            trimmedRV[["trimmedFinishPos"]]]
            } else if (directionParam == "Reverse") {
                PhredScore <-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@
                    QualityReport@qualityPhredScores[
                        (trimmedRV[["trimmedStartPos"]]+1):
                            trimmedRV[["trimmedFinishPos"]]]
            }
            qualityScoreDisplay (PhredScore)
        }
    })
    ############################################################################
    ### Trimming method selection functions
    ############################################################################
    ### ------------------------------------------------------------------------
    ### Trimming Method Selection (Just check forward reads => not dynamic)
    ### ------------------------------------------------------------------------
    output$TrimmingMethodSelectionOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            if (SangerAlignment@trimmingMethodSA == "M1") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                            'Modified Mott Trimming'")
            } else if (SangerAlignment@trimmingMethodSA == "M2") {
                tagList(icon("check-circle"),
                        "Your trimming method selection :
                            'Trimmomatics Sliding Window Trimming'")
            }
        }
    })
    output$TrimmingMethodUI <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            ## For method, everyone is same, so just pick forward one.
            if (SangerAlignment@trimmingMethodSA == "M1") {
                if (directionParam == "Forward") {
                    if (is.null(SangerAlignment@contigList[[contigIndex]]@
                                forwardReadList[[readIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerAlignment@contigList[[contigIndex]]@
                            forwardReadList[[readIndex]]@
                            QualityReport@M1TrimmingCutoff <<-  0.0001
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerAlignment@contigList[[contigIndex]]@
                                reverseReadList[[readIndex]]@
                                QualityReport@M1TrimmingCutoff)) {
                        SangerAlignment@contigList[[contigIndex]]@
                            reverseReadList[[readIndex]]@
                            QualityReport@M1TrimmingCutoff <<-  0.0001
                    }
                }
                fluidRow(
                    column(6,
                           uiOutput("M1TrimmingCutoff") ,
                           tags$ul(
                               textInput("M1TrimmingCutoffText",
                                         label = p("Change Value"),
                                         value = toString(
                                             trimmedParam[["M1TrimmingCutoff"]]),
                                         width = '70%')
                           ),
                    ),
                )
            } else if (SangerAlignment@trimmingMethodSA == "M2") {


                if (directionParam == "Forward") {
                    if (is.null(SangerAlignment@contigList[[contigIndex]]@
                                forwardReadList[[readIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerAlignment@contigList[[contigIndex]]@
                            forwardReadList[[readIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerAlignment@contigList[[contigIndex]]@
                                forwardReadList[[readIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerAlignment@contigList[[contigIndex]]@
                            forwardReadList[[readIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  10
                    }
                } else if (directionParam == "Reverse") {
                    if (is.null(SangerAlignment@contigList[[contigIndex]]@
                                reverseReadList[[readIndex]]@
                                QualityReport@M2CutoffQualityScore)) {
                        SangerAlignment@contigList[[contigIndex]]@
                            reverseReadList[[readIndex]]@
                            QualityReport@M2CutoffQualityScore <<-  20
                    }
                    if (is.null(SangerAlignment@contigList[[contigIndex]]@
                                reverseReadList[[readIndex]]@
                                QualityReport@M2SlidingWindowSize )) {
                        SangerAlignment@contigList[[contigIndex]]@
                            reverseReadList[[readIndex]]@
                            QualityReport@M2SlidingWindowSize <<-  10
                    }
                }
                fluidRow(
                    column(6,
                           uiOutput("M2CutoffQualityScore") ,
                           tags$ul(
                               textInput("M2CutoffQualityScoreText",
                                         label = p("Change Value"),
                                         value = toString(
                                             trimmedParam[["M2CutoffQualityScore"]]),
                                         width = '70%')
                           ),
                    ),
                    column(6,
                           uiOutput("M2SlidingWindowSize") ,
                           tags$ul(
                               textInput("M2SlidingWindowSizeText",
                                         label = p("Change Value"),
                                         value = toString(
                                             trimmedParam[["M2SlidingWindowSize"]]),
                                         width = '70%')
                           ),
                    ),
                )
            }
        }
    })

    valueBoxM1TrimmingCutoff (input, output, session)
    valueBoxM2CutoffQualityScore (input, output, session)
    valueBoxM2SlidingWindowSize (input, output, session)

    valueBoxRawSeqLength (input, output, session, trimmedRV)
    valueBoxRawMeanQualityScore (input, output, session, trimmedRV)
    valueBoxRawMinQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxTrimmedFinishPos (input, output, session, trimmedRV)
    valueBoxTrimmedSeqLength (input, output, session, trimmedRV)
    valueBoxTrimmedMeanQualityScore (input, output, session, trimmedRV)
    valueBoxTrimmedMinQualityScore (input, output, session, trimmedRV)
    valueBoxRemainingRatio (input, output, session, trimmedRV)

    output$qualityQualityBasePlot <- renderPlotly({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (!is.na(contigIndex) &&
            !is.na(readIndex)) {
            if (directionParam == "Forward") {
                qualityPhredScores <-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            } else if (directionParam == "Reverse") {
                qualityPhredScores <-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@
                    QualityReport@qualityPhredScores
            }
            qualityQualityBasePlotDisplay(input, output,session,
                                          trimmedRV, qualityPhredScores)
        }
    })

    valueBoxChromTrimmedStartPos (input, output, session, trimmedRV)
    valueBoxChromTrimmedFinishPos (input, output, session, trimmedRV)

    # chromatogram
    output$chromatogramUIOutput <- renderUI({
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (sidebar_menu[[2]] == "Contig" &&
            sidebar_menu[[3]] == "-" &&
            !is.na(readIndex) &&
            (directionParam == "Forward" ||
             directionParam == "Reverse") &&
            sidebar_menu[[6]] == "Read") {
            trimmedRV[["trimmedSeqLength"]]

            if (directionParam == "Forward") {
                rawSeqLength <-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                trimmedSeqLength <-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@
                    QualityReport@trimmedSeqLength
            } else if (directionParam == "Reverse") {
                rawSeqLength <-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@
                    QualityReport@rawSeqLength
                trimmedSeqLength <-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@
                    QualityReport@trimmedSeqLength
            }
            chromatogramRowNumAns <- chromatogramRowNum (
                strtoi(ChromatogramParam[["baseNumPerRow"]]),
                rawSeqLength, trimmedSeqLength,
                ChromatogramParam[["showTrimmed"]]) *
                strtoi(ChromatogramParam[["heightPerRow"]])
            plotOutput("chromatogram", height = chromatogramRowNumAns)
        }
    })
    output$chromatogram <- renderPlot({
        ## !!!!! Update !!!!
        shinyjs::disable("startTrimmingButton")
        shinyjs::disable("saveChromatogramParam")
        shinyjs::disable("ChromatogramBasePerRow")
        shinyjs::disable("ChromatogramHeightPerRow")
        shinyjs::disable("ChromatogramSignalRatioCutoff")
        shinyjs::disable("ChromatogramCheckShowTrimmed")
        shinyjs::disable("M1TrimmingCutoffText")
        shinyjs::disable("M2CutoffQualityScoreText")
        shinyjs::disable("M2SlidingWindowSizeText")
        sidebar_menu <- tstrsplit(input$sidebar_menu, " ")
        contigIndex <- strtoi(sidebar_menu[[1]])
        readIndex <- strtoi(sidebar_menu[[4]])
        directionParam <- sidebar_menu[[5]]
        if (sidebar_menu[[2]] == "Contig" &&
            sidebar_menu[[3]] == "-" &&
            !is.na(readIndex) &&
            (directionParam == "Forward" ||
             directionParam == "Reverse") &&
            sidebar_menu[[6]] == "Read") {
            removeNotification(id = "saveNotification")
            showNotification(
                ui = fluidRow(
                    column(1),
                    column(11,
                           tags$p(tagList(icon("hourglass-start"),
                                          "Plotting Chromatogram ...  "),
                                  style = "font-size: 28px;
                                   font-weight: bold;"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Base pairs per row : ",
                                          ChromatogramParam[["baseNumPerRow"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Height per row : ",
                                          ChromatogramParam[["heightPerRow"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Signal Ratio Cutoff : ",
                                          ChromatogramParam[["signalRatioCutoff"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p(tagList(icon("dot-circle"),
                                          "Show trimmed : ",
                                          ChromatogramParam[["showTrimmed"]]),
                                  style = "font-size: 20px;
                                   font-style: italic"),
                           tags$p("( If 'Signal Ratio Cutoff' is too small,
                                  it would need more time to replot
                                  the chromatogram)",
                                  style = "font-size: 16px;
                                   font-style: italic"),)
                ),
                action = NULL, duration = NULL, closeButton = FALSE,
                id = "chromatogramNotification",
                type = "message")
            if (directionParam == "Forward") {
                rawSeqLength <-
                    SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@QualityReport@rawSeqLength
                log_info(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (forward)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerAlignment@contigList[[contigIndex]]@
                                      forwardReadList[[readIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ----------------------------------------------------------------
                ### Update 'SangerContig'!
                ### ----------------------------------------------------------------
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ----------------------------------------------------------------
                ### Updating AASeqs
                ### ----------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerAlignment@contigList[[contigIndex]]@
                                                   forwardReadList[[readIndex]]@QualityReport@trimmedStartPos,
                                               SangerAlignment@contigList[[contigIndex]]@
                                                   forwardReadList[[readIndex]]@QualityReport@trimmedFinishPos,
                                               SangerAlignment@geneticCode)
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerAlignment@contigList[[contigIndex]]@
                    forwardReadList[[readIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ----------------------------------------------------------------
                ### Updating reactive values
                ### ----------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     forwardReadList[[readIndex]]@primaryAASeqS3)
            } else if (directionParam == "Reverse") {
                rawSeqLength <-
                    SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@QualityReport@rawSeqLength
                log_info(">>>>>>>>>>>> Re-running 'MakeBaseCalls' function (reverse)")
                ### ----------------------------------------------------------------
                ### Re-run 'MakeBaseCall' function
                ### ----------------------------------------------------------------
                hetcalls <-
                    MakeBaseCalls(SangerAlignment@contigList[[contigIndex]]@
                                      reverseReadList[[readIndex]],
                                  signalRatioCutoff = as.numeric(
                                      ChromatogramParam[["signalRatioCutoff"]]))
                ### ------------------------------------------------------------
                ### Update 'SangerContig'!
                ### ------------------------------------------------------------
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@peakPosMatrix <<-
                    hetcalls@peakPosMatrix
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@peakAmpMatrix <<-
                    hetcalls@peakAmpMatrix
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@primarySeq <<-
                    hetcalls@primarySeq
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@secondarySeq <<-
                    hetcalls@secondarySeq
                ### ------------------------------------------------------------
                ### Updating AASeqs
                ### ------------------------------------------------------------
                AASeqResult <- calculateAASeq (hetcalls@primarySeq,
                                               SangerAlignment@contigList[[contigIndex]]@
                                                   reverseReadList[[readIndex]]@QualityReport@trimmedStartPos,
                                               SangerAlignment@contigList[[contigIndex]]@
                                                   reverseReadList[[readIndex]]@QualityReport@trimmedFinishPos,
                                               SangerAlignment@geneticCode)
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@primaryAASeqS1 <<-
                    AASeqResult[["primaryAASeqS1"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@primaryAASeqS2 <<-
                    AASeqResult[["primaryAASeqS2"]]
                SangerAlignment@contigList[[contigIndex]]@
                    reverseReadList[[readIndex]]@primaryAASeqS3 <<-
                    AASeqResult[["primaryAASeqS3"]]
                ### ------------------------------------------------------------
                ### Updating reactive values
                ### ------------------------------------------------------------
                sequenceParam[["primarySeq"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primarySeq)
                sequenceParam[["secondarySeq"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@secondarySeq)
                sequenceParam[["primaryAASeqS1"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primaryAASeqS1)
                sequenceParam[["primaryAASeqS2"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primaryAASeqS2)
                sequenceParam[["primaryAASeqS3"]] <<-
                    as.character(SangerAlignment@contigList[[contigIndex]]@
                                     reverseReadList[[readIndex]]@primaryAASeqS3)
            }
            # log_info(">>>>>>>>>>>> 'MakeBaseCalls' finished")
            chromatogram_overwrite(hetcalls,
                                   width = strtoi(
                                       ChromatogramParam[["baseNumPerRow"]]),
                                   height = 2,
                                   trim5 = trimmedRV[["trimmedStartPos"]],
                                   trim3 = rawSeqLength -
                                       trimmedRV[["trimmedFinishPos"]],
                                   showtrim = (ChromatogramParam[["showTrimmed"]]),
                                   showcalls = "both", colors = colors)
            removeNotification(id = "chromatogramNotification")
            shinyjs::enable("ChromatogramBasePerRow")
            shinyjs::enable("ChromatogramHeightPerRow")
            shinyjs::enable("ChromatogramSignalRatioCutoff")
            shinyjs::enable("ChromatogramCheckShowTrimmed")
            shinyjs::enable("M1TrimmingCutoffText")
            shinyjs::enable("M2CutoffQualityScoreText")
            shinyjs::enable("M2SlidingWindowSizeText")
            shinyjs::enable("startTrimmingButton")
            shinyjs::enable("saveChromatogramParam")
        }
    })
}

