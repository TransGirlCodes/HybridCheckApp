library(shiny)
library(shinyBS)

shinyServer(function(input, output) {
  # Require the HybRIDS package.
  require(HybRIDS)
  # Start up a new HybRIDS session.
  hybridsobj <- HybRIDS$new()
  # Makes sure the DNA sequences are updated.
  # Reactive so when the DNA sequence is changed, then so is everything that depends on this.
  updateSequence <- reactive({
    validate(
      need(grepl(".fas", input$fastafile$name) || grepl(".fasta", input$fastafile$name), "Provide a FASTA format sequence file.")
    )
    hybridsobj$inputDNA(input$fastafile$datapath, "fasta")
  })
  
  output$SeqInfo <- renderText({
    updateSequence()
    validate(
      need(hybridsobj$DNA$hasDNA(), "No sequence file is loaded into HybRIDS.")
    )
    hybridsobj$DNA$htmlSummary()
  })
  
  # Gets HybRIDS to update the settings for generating triplets and update the triplets that are to be generated and analyzed.
  updateTripletGenSettings <- reactive({
    input$comboGen
    newoptions <- c()
    if(isolate(input$betweengroups)){
      newoptions <- c(newoptions, 1L)
    }
    if(isolate(input$distancebased)){
      if(isolate(input$distancemethod) == "yes"){
        newoptions <- c(newoptions, 3L)
      } else {
        newoptions <- c(newoptions, 2L)
      }
    }
    validate(
      need(!is.null(newoptions), "Select either to generate triplet combinations based on group specification or based on distance, or both.")
    )
    validate(
      need(hybridsobj$DNA$hasDNA(), "No sequence file is loaded into HybRIDS.")
    )
    groups <- lapply(unlist(lapply(1:isolate(input$numGroups), function(i) paste0("group", i))), function(n) isolate(input[[n]]))
    groups[sapply(groups, is.null)] <- NULL
    hybridsobj$setParameters("TripletGeneration",
                             Method = newoptions,
                             DistanceThreshold = isolate(input$mandistthreshold),
                             PartitionStrictness = as.integer(isolate(input$partitionStrictness)),
                             Groups = groups)
  })
  
  output$NumCombos <- renderText({
    updateTripletGenSettings()
    length(hybridsobj$comparrisonSettings$AcceptedCombinations)
  })
  
  output$GeneratedTriplets <- renderText({
    updateTripletGenSettings()
    hybridsobj$comparrisonSettings$htmlCombinations()
    })
  
  output$tripletgen <- renderUI({
    updateSequence()
    lapply(1:input$numGroups, function(i) {
      selectInput(inputId=paste0("group", i), label=paste0("Sequence Group ", i), choices = hybridsobj$DNA$getSequenceNames(), multiple=TRUE)
    })
  })
  
  # Server functions for the analyze triplets page.
  output$TripletSelector <- renderUI({
    updateTripletGenSettings()
    selectInput("tripletSelection", tags$strong("Selected Triplet"), hybridsobj$triplets$tripletCombinations())
  })
  
  currentTripletSelection <- reactive({
    return(unlist(strsplit(input$tripletSelection, ", ")))
  })
  
  scanSeqSim <- reactive({
    hybridsobj$setParameters("SSAnalysis", WindowSize = as.integer(input$windowSize),
                             StepSize = as.integer(input$stepSize))
    hybridsobj$analyzeSS(currentTripletSelection())
  })
  
  findBlocks <- reactive({
    scanSeqSim()
    hybridsobj$setParameters("BlockDetection", ManualThresholds = input$manBlockDetectDist,
                             AutoThresholds = (input$detectionMethod == "no"),
                             ManualFallback = input$fallbackManual)
    hybridsobj$findBlocks(currentTripletSelection())
  })
  
  dateBlocks <- reactive({
    dateanyway <- !input$eliminateinsignificant
    hybridsobj$setParameters("BlockDating", MutationRate = input$mu, PValue = input$alpha,
                             BonfCorrection = input$bonf, DateAnyway = dateanyway,
                             MutationCorrection = input$correctionModel)
    hybridsobj$dateBlocks(currentTripletSelection())
  })
  
  output$barsPlot <- renderPlot({
    scanSeqSim()
    print(hybridsobj$plotTriplets(currentTripletSelection(), What="Bars"))
  })
  
  output$linesPlot <- renderPlot({
    scanSeqSim()
    print(hybridsobj$plotTriplets(currentTripletSelection(), What="Lines"))
  })
  
  output$blocksTable <- renderDataTable({
    findBlocks()
    dateBlocks()
    hybridsobj$tabulateDetectedBlocks(currentTripletSelection(), Neat=TRUE)
  })
  
})
