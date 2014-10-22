library(shiny)

shinyServer(function(input, output) {
  
  require(HybRIDS)
  
  hybridsobj <- HybRIDS$new()
  
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
  
  output$tripletgen <- renderUI({
    updateSequence()
    validate(
      need(input$numGroups >= 1, "Enter how many sequence groups there are.")
    )
    lapply(1:input$numGroups, function(i) {
      selectInput(inputId=paste0("group", i), label=paste0("Sequence Group ", i), choices = hybridsobj$DNA$getSequenceNames(), multiple=TRUE)
    })
  })

})
