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
                            fluidRow(
                              column(4,
                                     inputPanel(
                                       h4("How to scan sequence similarity:"),
                                       numericInput("windowSize", "Size of sliding window (in bp):", min=1, value=100),
                                       numericInput("stepSize", "Step size of sliding window (in bp):", min=1, value=1)),
                                     inputPanel(
                                       uiOutput("TripletSelector")
                                       )),
                              column(4,
                                     inputPanel(
                                       h4("Block detection settings:"),
                                       radioButtons("detectionMethod", "Use a manual threshold?", c("Yes" = "yes", "No, automatically decide thresholds from data." = "no"), selected="no"),
                                       checkboxInput("fallbackManual", "Fallback to a manual threshold?"),
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
                            hr(),
                            bsCollapse(multiple = FALSE, open = "col1", id = "collapse1",
                                       bsCollapsePanel("Sequence Similarity, RGB plot",
                                                       plotOutput("barsPlot"),
                                                       id="col1",
                                                       value="test1"),
                                       bsCollapsePanel("Sequence Similarity, lines plot",
                                                       plotOutput("linesPlot"),
                                                       id="col2",
                                                       value="test2")
                            ),
                            dataTableOutput("blocksTable")                    
)))
