library(shiny)
source("mmp_gene_enrichment.R", local = TRUE)

mmp_dir <- "./data"
mmp_data <- read.table(gzfile(paste0(mmp_dir, "/mmp_mutations.txt.gz")), sep = "\t", header = TRUE)

server <- function(input, output, session) {
    observe({
        input_strain_list <<- input$input_strains
    })
    observeEvent(input$start, {
  
            output$result <- renderTable({
                mmp_gene_enrichment(strains = input_strain_list, ncrna = input$ncrna_box, censor = input$censor, mmp = mmp_data)
            }, digits = 4)
            
    })
        
}

ui <- basicPage(
    
    h3("Calculate gene enrichment for a list of MMP strains"),
    textInput(inputId = "input_strains", label = "Strains to analyze"),
    checkboxInput("ncrna_box", "ncRNA"),
    checkboxInput("censor", "censor strains with missing sequencing data"),
    actionButton("start", "Analyze strains"),
    tableOutput("result")
    
)

shinyApp(ui = ui, server = server)