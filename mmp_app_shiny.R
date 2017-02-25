library(shiny)

server <- function(input, output, session) {
    observe({
        strain_list <<- input$input_strains
    })
    observeEvent(input$start, {
        print(strain_list)
    })
        
}

ui <- basicPage(
    
    h3("Calculate gene enrichment for a list of MMP strains"),
    textInput(inputId = "input_strains", label = "Strains to analyze"),
    actionButton("start", "Analyze strains")
    
)

shinyApp(ui = ui, server = server)