# Following code makes sure all the dependencies of HybRIDS and HybRIDSapp are loaded and installed.
# Note it only considered installed/noninstalled - it does not check version numbers.
installed <- installed.packages()
# Packages required from R's package manager:
chooseCRANmirror(ind = 83)
pkg <- c("devtools", "shiny", "ggplot2", "png", "grid", "gridExtra", "ape", "Rcpp")
new.pkg <- pkg[!(pkg %in% installed)]
if(length(new.pkg)){install.packages(new.pkg)}
if(!("Biostrings" %in% installed)){
  source("http://bioconductor.org/biocLite.R")
  biocLite()
  biocLite("Biostrings")
}
if(!("HybRIDS" %in% installed)){
  library(devtools)
  install_github("Ward9250/HybRIDS")
}
if(!("shinyBS" %in% installed)){
  library(devtools)
  install_github("ebailey78/shinyBS")
}

library(shiny)
library(shinyBS)
library(HybRIDS)

shinyServer(function(input, output, session){
  
  # Start up a new HybRIDS session.
  hybridsobj <- HybRIDS$new()
  
  # DNA sequence loading.
  updateSequence <- reactive({
    validate(
      need(grepl(".fas", input$fastafile$name) || grepl(".fasta", input$fastafile$name), "Provide a FASTA format sequence file.")
    )
    hybridsobj$inputDNA(input$fastafile$datapath)
  })
  
  updatePopulations <- reactive({
    input$setPops
    if(!isolate(input$oneSeqOnePop)){
      hybridsobj$setPopulations(NULL)
    } else {
      nameSel <- unlist(lapply(1:isolate(input$numPops), function(i) paste0("PopulationName", i)))
      popSel <- unlist(lapply(1:isolate(input$numPops), function(i) paste0("Population", i)))
      populations <- setNames(object = lapply(popSel, function(i) isolate(input[[i]])), 
                              unlist(lapply(nameSel, function(i) isolate(input[[i]]))))
      hybridsobj$setPopulations(populations)
    }
  })
  
  output$SeqInfo <- renderText({
    updateSequence()
    updatePopulations()
    validate(need(hybridsobj$DNA$hasDNA(), "Provide a FASTA format alignment to begin."))
    hybridsobj$DNA$htmlSummary()
  })
  
  output$populationsGen <- renderUI({
    updateSequence()
    validate(need(!is.na(input$numPops), "Enter a number of populations."))
    lapply(1:input$numPops, function(i) {
      fluidRow(
        column(6, textInput(inputId=paste0("PopulationName", i), label=paste0("Population Name"))),
        column(6,
               selectInput(inputId=paste0("Population", i), label=paste0("Sequences"),
                           choices = hybridsobj$DNA$getSequenceNames(), multiple = TRUE) 
        )
      )
    })
  })
                    
                    
  
  
  # ABBA-BABA section
  
  output$ABBAtree <- renderImage({
    width  <- session$clientData$output_ABBAtree_width
    height <- session$clientData$output_ABBAtree_height
    pixelratio <- session$clientData$pixelratio
    outfile <- tempfile(fileext='.png')
    png(outfile, width=width*pixelratio, height=height*pixelratio,
        res=72*pixelratio)
    plot(read.tree(text="(((P1,P2),P3),P4);"))
    nodelabels("A", c(1,4), adj = c(-2.5, 0.5), bg = "red", col="white")
    nodelabels("B", c(2,3), adj = c(-2.5, 0.5), bg = "blue", col="white")
    dev.off()
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")}, deleteFile = TRUE)
             
  output$BABAtree <- renderImage({
    width  <- session$clientData$output_BABAtree_width
    height <- session$clientData$output_BABAtree_height
    pixelratio <- session$clientData$pixelratio
    outfile <- tempfile(fileext='.png')
    png(outfile, width=width*pixelratio, height=height*pixelratio,
        res=72*pixelratio)
    plot(read.tree(text="(((P1,P2),P3),P4);"))
    nodelabels("A", c(2,4), adj = c(-2.5, 0.5), bg = "red", col="white")
    nodelabels("B", c(1,3), adj = c(-2.5, 0.5), bg = "blue", col="white")
    dev.off()
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")}, deleteFile = TRUE)
  
  output$fttGen <- renderUI({
    updateSequence()
    updatePopulations()
    validate(need(!is.na(input$fttNumCombos), "Enter a number of taxa combos to analyse."))
    validate(need(length(hybridsobj$DNA$Populations) >= 4, "You need 4 or more populations defined."))
    fluidRow(
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("P1", i), label = "P1",
                           choices = hybridsobj$DNA$namesOfPopulations())
             })),
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("P2", i), label = "P2",
                           choices = hybridsobj$DNA$namesOfPopulations())
             })),
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("P3", i), label = "P3",
                           choices = hybridsobj$DNA$namesOfPopulations())
             })),
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("A", i), label = "P4",
                           choices = hybridsobj$DNA$namesOfPopulations())
             }))
    )
  })
  
  updateTaxaCombos <- reactive({
    input$setCombinations
    if(isolate(input$fttAutoSets)){
      hybridsobj$prepareFourTaxonTests()
    } else {
      combos <- lapply(1:isolate(input$fttNumCombos), function(i){
        c(P1 = isolate(input[[paste0("P1", i)]]), P2 = isolate(input[[paste0("P2", i)]]),
          P3 = isolate(input[[paste0("P3", i)]]), A = isolate(input[[paste0("A", i)]]))
      })
      hybridsobj$prepareFourTaxonTests(combos)
    }
  })
  
  output$generatedCombos <- renderText({
    updateTaxaCombos()
    paste(hybridsobj$FTTmodule$printAllNames(), collapse = "\n\n")
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
    updateTripletGenSettings()
    lapply(1:input$numGroups, function(i) {
      selectInput(inputId=paste0("group", i), label=paste0("Sequence Group ", i),
                  choices = hybridsobj$DNA$getSequenceNames(), multiple = TRUE)
    })
  })
  
  # Server functions for the analyze triplets page.
  
  output$TripletToAnalyze <- renderUI({
    updateTripletGenSettings()
    selectInput("tripletToAnalyze", 
                tags$strong("Run / rerun analysis for triplets:"),
                c("ALL", hybridsobj$triplets$tripletCombinations()),
                selected = "ALL", multiple = TRUE)
  })
  
  output$TripletSelector <- renderUI({
    updateTripletGenSettings()
    selectInput("tripletSelection", tags$strong("View triplet:"),
                hybridsobj$triplets$tripletCombinations())
  })
  
  updatePlottingSettings <- reactive({
    hybridsobj$setParameters("Plotting", PlotTitle = input$plotTitle,
                             TitleSize = input$plotTitleSize,
                             TitleFace = input$plotTitleFace,
                             TitleColour = input$plotTitleColour,
                             XLabels = input$plotXLabels,
                             YLabels = input$plotYLabels,
                             XTitle = input$plotXTitle,
                             XTitleFontSize = input$plotXTitleFontSize,
                             XTitleColour = input$plotXTitleColour,
                             XLabelSize = input$plotXLabelSize,
                             XLabelColour = input$plotXLabelColour,
                             YTitle = input$plotYTitle,
                             YTitleFontSize = input$plotYTitleFontSize,
                             YTitleColour = input$plotYTitleColour,
                             YLabelSize = input$plotYLabelSize,
                             YLabelColour = input$plotYLabelColour,
                             Legends = input$plotLegends,
                             LegendFontSize = input$plotLegendFontSize,
                             MosaicScale = input$plotMosaicScale)
  })
  
  analysis <- reactive({
    input$analysisGO
    if(hybridsobj$DNA$hasDNA() && (input$analysisGO != 0)){
      hybridsobj$setParameters("SSAnalysis", WindowSize = as.integer(isolate(input$windowSize)),
                               StepSize = as.integer(isolate(input$stepSize)))
      hybridsobj$setParameters("BlockDetection", ManualThresholds = isolate(input$manBlockDetectDist),
                               AutoThresholds = (isolate(input$detectionMethod) == "no"),
                               ManualFallback = isolate(input$fallbackManual))
      dateanyway <- !isolate(input$eliminateinsignificant)
      hybridsobj$setParameters("BlockDating", MutationRate = isolate(input$mu), PValue = isolate(input$alpha),
                               BonfCorrection = isolate(input$bonf), DateAnyway = dateanyway,
                               MutationCorrection = isolate(input$correctionModel))
      selections <- strsplit(isolate(input$tripletToAnalyze), ", ")
      tripletsToDo <- hybridsobj$triplets$getTriplets(selections)
      numToDo <- length(tripletsToDo)
      prog <- Progress$new(session, min = 1, max = numToDo)
      prog$set(value = 0, message = "Scanning SS in triplets...")
      for(i in 1:numToDo){
        HybRIDS:::scan.similarity(hybridsobj$DNA, tripletsToDo[[i]], hybridsobj$ssAnalysisSettings)
        prog$set(value = i)
      }
      prog$close()
      prog <- Progress$new(session, min = 1, max = numToDo)
      prog$set(value = 0, message = "Finding blocks in triplets...")
      for(i in 1:numToDo){
        tripletsToDo[[i]]$putativeBlockFind(hybridsobj$blockDetectionSettings)
        prog$set(value = i)
      }
      prog$close()
      prog <- Progress$new(session, min = 1, max = numToDo)
      prog$set(value = 0, message = "Dating blocks...")
      for(i in 1:numToDo){
        tripletsToDo[[i]]$blockDate(hybridsobj$DNA, hybridsobj$blockDatingSettings)
        prog$set(value = i)
      }
      prog$close() 
    }
  })
  
  output$barsPlot <- renderPlot({
    analysis()
    updatePlottingSettings()
    validate(need(is.character(input$tripletSelection), 
                  "Triplets to choose from have not been generated yet."))
    selection <- unlist(strsplit(input$tripletSelection, ", "))
    validate(need(!hybridsobj$triplets$getTriplets(selection)[[1]]$noScanPerformed(),
             "Sequence similarity scan has not been performed for this triplet yet."))
    barsPlot <- hybridsobj$plotTriplets(unlist(strsplit(input$tripletSelection, ", ")), What="Bars")[[1]]
    print(barsPlot)
  })

  output$linesPlot <- renderPlot({
    analysis()
    updatePlottingSettings()
    validate(need(is.character(input$tripletSelection), 
                  "Triplets to choose from have not been generated yet."))
    selection <- unlist(strsplit(input$tripletSelection, ", "))
    validate(need(!hybridsobj$triplets$getTriplets(selection)[[1]]$noScanPerformed(),
             "Sequence similarity scan has not been performed for this triplet yet."))
    linesPlot <- hybridsobj$plotTriplets(selection, What="Lines")[[1]]
    print(linesPlot)
  })

  output$blocksTable <- renderDataTable({
    analysis()
    validate(need(is.character(input$tripletSelection),
                  "Triplets to choose from have not been generated yet."))
    validate(
      need(!hybridsobj$triplets$getTriplets(
        unlist(strsplit(input$tripletSelection, ", ")))[[1]]$blocksNotDated(),
           "Recombinant regions have not been found or dated yet for the selected triplet."))
    hybridsobj$tabulateDetectedBlocks(unlist(strsplit(input$tripletSelection, ", ")), Neat=TRUE)
  })
  
  output$userBlocksTable <- renderDataTable({
    updateSequence()
    clearUserBlocks()
    addUserBlocks()
    dateUserBlocks()
    hybridsobj$tabulateUserBlocks()
  })
  
  output$saveTable <- downloadHandler(filename = 
                                        function(){
                                          paste0(strsplit(input$fastafile$name, ".fas")[1], "_Triplet_",
                                                 paste(unlist(strsplit(input$tripletSelection, ", ")), collapse = ":"), ".csv")
                                        },
                                      content =
                                        function(file){
                                          write.csv(hybridsobj$tabulateDetectedBlocks(unlist(strsplit(input$tripletSelection, ", ")), Neat=TRUE), file)
                                        }
  )
  
  output$saveBarPlots <- downloadHandler(filename =
                                        function(){
                                          paste0(strsplit(input$fastafile$name, ".fas")[1], "_Triplet_",
                                                 paste(unlist(strsplit(input$tripletSelection, ", ")), collapse = ":"), "_Bars.png")
                                        },
                                        content = 
                                          function(file){
                                            selection <- unlist(strsplit(isolate(input$tripletSelection), ", "))
                                            whichToPlot <- isolate(input$whichPlotToSave)
                                            if(whichToPlot == "Both"){
                                              whichToPlot <- c("Bars", "Lines")
                                            }
                                            Plot <- hybridsobj$plotTriplets(selection, What = whichToPlot)[[1]]
                                            ggsave("plot.png", plot = Plot, width = isolate(input$widthSave), height = isolate(input$heightSave),
                                                   dpi = isolate(input$resSave), units = "in")
                                            file.copy("plot.png", file, overwrite=TRUE)
                                          })
    
  output$userBlocksSeqSelect <- renderUI({
    updateSequence()
    selectInput("seqChoice", "Select two sequences",
                hybridsobj$DNA$getSequenceNames(), multiple = TRUE)
  })
  
  addUserBlocks <- reactive({
    input$addUBButton
    validate(need(isolate(input$startPosition) < isolate(input$endPosition), "Start position must be less than the end position."))
    if(length(isolate(input$seqChoice)) == 2){
      hybridsobj$addUserBlock(isolate(input$seqChoice), isolate(input$startPosition), isolate(input$endPosition))
    }
  })
  
  clearUserBlocks <- reactive({
    input$clearUBButton
    if(length(isolate(input$seqChoice)) == 2){
      hybridsobj$clearUserBlocks(isolate(input$seqChoice))
    }
  })
  
  dateUserBlocks <- reactive({
    input$dateUBButton
    dateanyway <- !isolate(input$eliminateinsignificant2)
    hybridsobj$setParameters("BlockDating", MutationRate = isolate(input$mu2), PValue = isolate(input$alpha2),
                             BonfCorrection = isolate(input$bonf2), DateAnyway = dateanyway,
                             MutationCorrection = isolate(input$correctionModel2))
    hybridsobj$dateUserBlocks()
  })
  
  output$activeTab <- reactive({
    return(input$tab)
  })
  outputOptions(output, 'activeTab', suspendWhenHidden=FALSE)
  
})
