library(shiny)
source("mmp_gene_enrichment.R", local = TRUE)

setwd("/home/sfrenk/Documents/Scripts/mmp_app")
mmp_dir <- "./data"
mmp <- read.table(gzfile(paste0(mmp_dir, "/mmp_mutations.txt.gz")), sep = "\t", header = TRUE)

server <- function(input, output, session) {
    observe({
        strain_list <<- input$input_strains
    })
    observeEvent(input$start, {
        output$result <- renderTable({
            mmp_gene_enrichment(strain_list, mmp = mmp)
        })
    })
        
}

ui <- basicPage(
    
    h3("Calculate gene enrichment for a list of MMP strains"),
    textInput(inputId = "input_strains", label = "Strains to analyze"),
    actionButton("start", "Analyze strains"),
    tableOutput("result")
    
)

shinyApp(ui = ui, server = server)