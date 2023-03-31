#' module_levels UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_module_levels_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    shiny::fluidRow(
    shiny::column(12,
      shiny::h1("Explore normalized gene expression"),
      shiny::hr()
    ),
    
    shinybusy::add_busy_spinner(
      spin = "self-building-square",
      position = 'top-left',
      margins = c(70, 1200)
    ),
    
    shinydashboard::tabBox(
      title = "Explore normalized data",
      width = 12,
      height = "80vh",
      shiny::tabPanel(title = "PCA",
                      shiny::fluidRow(
                      shiny::uiOutput(ns("pca_ui")),
                      shiny::column(12, shiny::includeMarkdown(system.file(
                        "extdata", "pca.md", package = "DIANE"))))
      ),
        shiny::tabPanel(title = "Visualize gene expression levels",
                        shiny::fluidRow(
                        shinydashboardPlus::box(
                          title = "Genes and conditions choice",
                          solidHeader = FALSE,
                          status = "success",
                          collapsible = TRUE,
                          closable = FALSE,
                          width = 12,
                          
                          
                          shiny::uiOutput(ns("gene_choice")),
                          
                          
                          shiny::uiOutput(ns("condition_choice"))
                          
                        ),
                        shinydashboardPlus::box(solidHeader = FALSE,
                                                status = "success",
                                                collapsible = TRUE,
                                                closable = FALSE,
                                                width = 12,
                                                shiny::plotOutput(ns("expression_plot"), height = "700px"))
                        )  
        )
  )
  ))
}
    
#' module_levels Server Function
#'
#' @noRd 
mod_module_levels_server <- function(input, output, session, r){
  ns <- session$ns
  
  ###PCA computation when page load.
  pca_raw_results <- shiny::reactive({
    shiny::req(r$normalized_counts, r$conditions)
    compute_pca(r$normalized_counts, kept_axes = 10)
  })
  
  output$condition_choice <- shiny::renderUI({
    shiny::req(r$normalized_counts, r$conditions)
    
    shinyWidgets::checkboxGroupButtons(
      inputId = ns('input_conditions'),
      label = "Conditions to include to the expression levels plot:",
      choices = unique(r$conditions),
      justified = TRUE,
      checkIcon = list(yes = shiny::icon("ok",
                                         lib = "glyphicon")),
      selected = unique(r$conditions)
    )
  })
  
  output$gene_choice <- shiny::renderUI({
    shiny::req(r$normalized_counts)
    
    shiny::textInput(ns("genes"), 
                     label = "Genes to plot, as identified in the Gene 
                     column of expression data. For several genes, they must be comma separated, without space, as in the example:", 
                     width = '100%',
                     value = paste0(sample(rownames(r$normalized_counts), 4), collapse = ','))
  })
  
  #   ____________________________________________________________________________
  #   profiles                                                                ####
  
  output$expression_plot <- shiny::renderPlot({
    
    shiny::req(r$normalized_counts, r$conditions, input$genes)
    
    genes <- unlist(strsplit(input$genes, ','))
    shiny::req(length(genes) > 0)
    shiny::req(length(genes) < 10)
    
    shiny::req(sum(genes %in% rownames(r$normalized_counts)) > 0)
    
  
  draw_expression_levels(as.data.frame(r$normalized_counts),
                           conds = input$input_conditions,
                           genes = genes, gene.name.size = 22)
  })
  
  
  #   ____________________________________________________________________________
  #   pca                                                                     ####
  
  output$pca_plot <- shiny::renderPlot({
    shiny::req(r$normalized_counts)
    # draw_PCA(r$normalized_counts)
    quick_pca(r$normalized_counts)
  })
  
  ###UI for all the pca related things. Contains 3 tabs with 3 plotoutputs.
  output$pca_ui <- shiny::renderUI({
    if(is.null(r$normalized_counts)) {
      shinydashboardPlus::descriptionBlock(
        number = "Please normalize and filter raw data in previous tab",
        numberColor = "orange",
        rightBorder = FALSE
      )
    } else {
      # else shiny::plotOutput(ns('pca_plot'), height = "800px")
      shinydashboard::tabBox(
        id = "tabset_advanced_pca",
        # height = "450px",
        width = NULL,
        tabPanel("PCA Summary",
                 shiny::plotOutput(ns('pca_plot'), height = "800px")
        ),
        tabPanel("Specific PCA plot",
                 shiny::tagList(
                   shiny::column(6, align = "center",
                                 shiny::selectInput(ns("component_1_choice"), "First component",
                                                    as.character(1:(ncol(pca_raw_results()$co)-2)),
                                                    selected = "1")
                   ),
                   shiny::column(6, align = "center",
                                 shiny::selectInput(ns("component_2_choice"), "Second component",
                                                    as.character(1:(ncol(pca_raw_results()$co)-2)),
                                                    selected = "2")
                   ),
                   shiny::column(12, align = "center",
                                 shiny::plotOutput(ns("specific_pca_plot"), height = "800px")
                   )
                 )
        ),
        tabPanel("PCA correlation plot",
                 shiny::plotOutput(outputId = ns("pca_plot_correlation"))
        )
      )
    }
    
    
    
  })
  
  ###For the second tab. Plot two specific components.
  output$specific_pca_plot <- shiny::renderPlot({
    req(pca_raw_results(), input$component_1_choice, input$component_2_choice)
    draw_specific_pca(
      pca_raw_results(),
      component_1 = input$component_1_choice,
      component_2 = input$component_2_choice,
      legend = TRUE
    )
  })
  
  ###Third tab, PCA correlation plot.
  output$pca_plot_correlation <- shiny::renderPlot({
    req(pca_raw_results())
    if (is.null(r$design)) {
      golem::print_dev("In dev\n")("no design provided")
      pca_plot_correlation(pca_raw_results())
    } else {
      golem::print_dev("A design is provided")
      pca_plot_correlation(pca_raw_results(), design = r$design)
    }
  })
  
 
}
    
## To be copied in the UI
# mod_module_levels_ui("module_levels_ui_1")
    
## To be copied in the server
# callModule(mod_module_levels_server, "module_levels_ui_1")
 
