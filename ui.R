library(shiny)
library(shinyBS)

shinyUI(navbarPage("HybRIDS",
                   
                   tabPanel("Sequence Data",
                            h1("Input DNA sequence file."),hr(),
                            p("Upload your FASTA file containing sequences from your different populations to get started."),
                            p("Note: Duplicate Sequences will be removed."), hr(),
                            div(align='center',
                            fileInput('fastafile', 'Choose FASTA file to upload',
                                      accept = c(
                                        '.fas',
                                        '.fasta',
                                        '.FAS',
                                        '.FASTA'
                                      )
                            )),
                            hr(),
                            htmlOutput("SeqInfo")
                            ),
                   
                   tabPanel("Specify Triplets to Analyze",
                            h1("Specify sequence combinations to scan."), hr(),
                            fluidRow(column(6,
                                            h3("Specify what triplets to scan:"),
                                            "By default, every possible combination of 3 sequences will be analyzed.",
                                            "You can refine this by eliminating triplets containing sequences that are insufficiently diverged.",
                                            "",
                                            inputPanel(
                                              checkboxInput("betweengroups", "Compare between groups", value=TRUE),
                                              checkboxInput("distancebased", "Raw p-distance based"),
                                              conditionalPanel("input.betweengroups",
                                                               h4("Group Options"),
                                                               numericInput("numGroups", "How many groups?", value=3, min=1),
                                                               radioButtons("partitionStrictness", "How many sequences from the same group are allowed in a triplet?",
                                                                            c("One" = 1L,
                                                                              "Two" = 2L), selected=2L),
                                                               uiOutput("tripletgen")),
                                              conditionalPanel("input.distancebased",
                                                               h4("Distance based elimination options."),
                                                               radioButtons("distancemethod", "Use automatic threshold detection?",
                                                                            c("Yes" = "yes",
                                                                              "No, specify a raw p-distance to use as a threshold manually." = "no")),
                                              conditionalPanel("input.distancemethod == 'no'",
                                                               sliderInput("mandistthreshold", "Raw p-distance threshold.", 0.01, 1, value=0.01, step=0.01))),
                                              actionButton("comboGen", "Generate Combinations to Scan")
                                              )  
                                     ),
                                     column(6,
                                            fluidRow(column(6, strong("Number of combinations to scan: "),
                                                            textOutput("NumCombos"),
                                                            strong("Combinations to scan:"),
                                                            htmlOutput("GeneratedTriplets")))
                                            ))),
                   
                   tabPanel("Analyze Triplets",
                            h1("Analyze and Explore Sequence Triplets"),
                            hr(),
                            h2("Settings:"),
                            fluidRow(
                              column(4,
                                     inputPanel(
                                       h4("How to scan sequence similarity:"),
                                       numericInput("windowSize", "Size of sliding window (in bp):", min=1, value=100),
                                       numericInput("stepSize", "Step size of sliding window (in bp):", min=1, value=1)),
                                     inputPanel(
                                       uiOutput("TripletToAnalyze"),
                                       br(),
                                       actionButton("analysisGO", "Run Analysis")
                                       )),
                              column(4,
                                     inputPanel(
                                       h4("Block detection settings:"),
                                       radioButtons("detectionMethod", "Use a manual threshold?", c("Yes" = "yes", "No, automatically decide thresholds from data." = "no"), selected="no"),
                                       checkboxInput("fallbackManual", "Fallback to a manual threshold?", value = TRUE),
                                       conditionalPanel("input.detectionMethod == 'yes' || input.fallbackManual", 
                                                        sliderInput("manBlockDetectDist", "Sequence Similarity Threshold:", 1, 100, value=95, step=1))
                                     )),
                              column(4,
                                     inputPanel(
                                       h4("Block dating settings:"),
                                       numericInput("mu", "Mutation Rate:", 10e-8),
                                       numericInput("alpha", "Critical Value (alpha):", 0.05),
                                       checkboxInput("bonf", "Bonferoni Correct Critical Value", TRUE),
                                       checkboxInput("eliminateinsignificant", "Eliminate insignificant blocks", TRUE),
                                       selectInput("correctionModel", "Mutation correction model:", c("JC69", "K80", "F81",
                                                                                                      "K81", "F84", "BH87",
                                                                                                      "T92", "TN93", "GG95"))))),
                            bsCollapse(
                                       bsCollapsePanel("Plotting Settings",
                                                       fluidRow(
                                                         column(4,
                                                                inputPanel(
                                                                  h4("Plot title settings:"),
                                                                  br(),
                                                                  checkboxInput("plotTitle", "Include a title in the plot?", TRUE),
                                                                  br(),
                                                                  numericInput("plotTitleSize", "Font size of plot titles", 14),
                                                                  br(),
                                                                  textInput("plotTitleFace", "Font face of plot titles", "bold"),
                                                                  br(),
                                                                  textInput("plotTitleColour", "Colour of plot titles", "black"),
                                                                  br(),
                                                                  h4("Legend settings:"),
                                                                  br(),
                                                                  checkboxInput("plotLegends", "Add legends to the plots?", TRUE),
                                                                  br(),
                                                                  numericInput("plotLegendFontSize", "Font size for legend", 12)
                                                                  )),
                                                         column(4,
                                                                inputPanel(
                                                                  h4("X axis options:"),
                                                                  br(),
                                                                  checkboxInput("plotXTitle", "Include the title of the x-axis in the plots?", TRUE),
                                                                  br(),
                                                                  numericInput("plotXTitleFontSize", "Font size of x-axis title", 12),
                                                                  br(),
                                                                  textInput("plotXTitleColour", "Colour of x-axis title", "black"),
                                                                  br(),
                                                                  checkboxInput("plotXLabels", "Include the value labels of the x-axis in the plots?", TRUE),
                                                                  br(),
                                                                  numericInput("plotXLabelSize", "Font size of x-axis labels", 10),
                                                                  br(),
                                                                  textInput("plotXLabelColour", "Colour of x-axis labels", "black")
                                                                )),
                                                         column(4,
                                                                inputPanel(
                                                                  h4("Y axis options:"),
                                                                  br(),
                                                                  checkboxInput("plotYTitle", "Include the title of the y-axis in the plots?", TRUE),
                                                                  br(),
                                                                  numericInput("plotYTitleFontSize", "Font size of y-axis title", 12),
                                                                  br(),
                                                                  textInput("plotYTitleColour", "Colour of y-axis title", "black"),
                                                                  br(),
                                                                  checkboxInput("plotYLabels", "Include the value labels of the y-axis in the plots?", TRUE),
                                                                  br(),
                                                                  numericInput("plotYLabelSize", "Font size of y-axis labels", 10),
                                                                  br(),
                                                                  textInput("plotYLabelColour", "Colour of y-axis labels", "black")
                                                                ))),
                                                       fluidRow(
                                                         inputPanel(
                                                           h4("Other Settings:"),
                                                           numericInput("plotMosaicScale", "Number of segments in RGB bars", 500),
                                                           numericInput("heightSave", "Height of plots when saved to file (inches)", 10),
                                                           numericInput("widthSave", "Width of plots when saved to file (inches)", 12),
                                                           numericInput("resSave", "Resolution of plots when saved (dpi)", 300),
                                                           radioButtons("whichPlotToSave", "Which plot to save?", list(Bars = "Bars", Lines = "Lines", Both = "Both"))
                                                           ))
                                                       )),
                            hr(),
                            h2("View Results:"),
                            fluidRow(div(inputPanel(uiOutput("TripletSelector"), downloadButton("saveBarPlots", "Save Plots"), downloadButton("saveTable", "Save Table"), align="center"))),
                            br(),
                            fluidRow(column(12,
                            bsCollapse(multiple = FALSE, open = "col1", id = "collapse1",
                                       bsCollapsePanel("Sequence Similarity, RGB plot",
                                                       plotOutput("barsPlot"),
                                                       id="col1",
                                                       value="test1"),
                                       bsCollapsePanel("Sequence Similarity, lines plot",
                                                       plotOutput("linesPlot"),
                                                       id="col2",
                                                       value="test2")),
                            dataTableOutput("blocksTable"))),
                            id = "tab"
                            ),
                   
                   tabPanel("Analyze user defined regions",
                            h1("Analyze user defined regions"),
                            hr(),
                            "You can specify a recombinant block between two sequences,
                            and analyze it according to the settings.",
                            fluidRow(
                              column(4,
                                     inputPanel(
                                       h4("Add or clear user blocks:"),
                                       uiOutput("userBlocksSeqSelect"),
                                       numericInput("startPosition", "Start of block in BP", value = NULL,
                                                    min = 1),
                                       numericInput("endPosition", "End of block in BP", value = NULL,
                                                    min = 1),
                                       actionButton("addUBButton", "Add user defined block between sequences"),
                                       actionButton("clearUBButton", "Remove user defined blocks between sequences"),
                                       actionButton("dateUBButton", "Test and date user defined blocks between sequences")
                                       )),
                              column(4,
                                     inputPanel(
                                       h4("Block dating settings:"),
                                       numericInput("mu2", "Mutation Rate:", 10e-8),
                                       numericInput("alpha2", "Critical Value (alpha):", 0.05),
                                       checkboxInput("bonf2", "Bonferoni Correct Critical Value", TRUE),
                                       checkboxInput("eliminateinsignificant2", "Eliminate insignificant blocks", TRUE),
                                       selectInput("correctionModel2", "Mutation correction model:", c("JC69", "K80", "F81",
                                                                                                      "K81", "F84", "BH87",
                                                                                                      "T92", "TN93", "GG95"))))
                              
                              ),
                            dataTableOutput("userBlocksTable")
                            )
))
