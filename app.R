library(shiny)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(plotly)
source("cell_select_ggnetwork.R")

celltype <- c(unique(all_results_merged$CellType))
ui <- fluidPage(
  navbarPage(theme = shinytheme("flatly"),title = h2("Cell-Phenotype"),
             tabPanel(h3("Cell type"),fluidRow(
               column(2, wellPanel(
                 pickerInput(inputId = "cellId",label = "Select cell type",choices = celltype,selected = "Microglia",
                             options = list(style = "btn-success")),
                 numericInput("qvalue", "q-value threshold", value = 0.0005,min = 0,max =1,step=0.0005),
                 numericInput("fold", "Minimum fold change", value = 1),
                 actionBttn(
                   inputId = "do",label = "Plot",style = "jelly", color = "success"),
                 br(),
                 br(),
                 sliderInput("Cell_plot_height", "Download image height (px)",
                             min = 400, max = 4000, value = 1500, step = 100),
                 sliderInput("Cell_plot_width", "Download image width (px)",
                             min = 400, max = 4000, value = 1500, step = 100),
                 downloadLink("download", "Download figure")
               )),
               column(10, plotlyOutput("cell_select",height="700px")),
               fluidPage(DTOutput('Celltype'))
             ))
  )


)


server <- function(input, output) {

  ggnetwork_plot <- eventReactive(input$do,ggnetwork_plot_full(phenotype_to_genes = phenotype_to_genes, 
                                                               all_results_merged = all_results_merged, 
                                                               mpo = mpo, 
                                                               disease_descriptions = disease_descriptions, 
                                                               cell_type = input$cellId, 
                                                               q_threshold = input$qvalue,
                                                               fold_threshold = input$fold))
  ggnetwork_DF <- eventReactive(input$do,subset_phenos(phenotype_to_genes = phenotype_to_genes, 
                                                       all_results_merged = all_results_merged, 
                                                       mpo = mpo, 
                                                       cell_type = input$cellId, 
                                                       q_threshold = input$qvalue,
                                                       fold_threshold = input$fold))
  output$cell_select <- plotly::renderPlotly(
    plotly::ggplotly(ggnetwork_plot(), tooltip = "hover"))
  
  output$Celltype <- renderDT(ggnetwork_DF(), 
                              rownames = FALSE,extensions = 'Buttons', 
                              options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  
  output$download <- downloadHandler(
    filename = paste0("sig_pheno_cell",Sys.Date(),".png"),
    content = function(filename) {
      png(filename, width = input$Cell_plot_width, height = input$Cell_plot_height)
      print(ggnetwork_plot())
      dev.off()
    }, contentType = "image/png")
}

 
shinyApp(ui = ui, server = server)
