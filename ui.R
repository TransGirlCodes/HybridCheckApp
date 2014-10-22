library(shiny)

shinyUI(navbarPage("HybRIDS",
                   
                   tabPanel("Sequence Data",
                            HTML("<h1>Input DNA sequence file.</h1>"),
                            HTML("<p>Upload your FASTA file containing sequences from your different populations to get started.</p>"),
                            HTML("<p>Note: Duplicate Sequences will be removed.</p>"),
                            HTML("<hr><div align='center'>"),
                            fileInput('fastafile', 'Choose FASTA file to upload',
                                      accept = c(
                                        '.fas',
                                        '.fasta',
                                        '.FAS',
                                        '.FASTA'
                                      )
                            ),
                            HTML("</div><hr>"),
                            htmlOutput("SeqInfo")
                            ),
                   
                   tabPanel("Selecting Comparrisons",
                            HTML("<h1>Specify sequence comparrisons to scan.</h1>"),
                            sidebarLayout(
                              sidebarPanel(radioButtons("combomethod", "Triplet Generation Method:",
                                                       c("All possible triplets" = "all",
                                                         "Between Groups" = "between")),
                                           conditionalPanel("input.combomethod == 'between'",
                                                            numericInput("numGroups", "How many groups?", value=3, min=1),
                                                            uiOutput("tripletgen"))
                                           ),
                              mainPanel("Ho!")
                              
                              )
                            
                            )
))
