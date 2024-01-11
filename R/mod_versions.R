#' versions UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList 
mod_versions_ui <- function(id){
  ns <- NS(id)
  tagList(
    shiny::includeMarkdown(system.file("extdata", "NEWS.md", package = "DIANE")),
    shiny::h1("Informations about loaded packages"),
    # htmlOutput(ns("sessionInfo"))
    shiny::hr(),
    shiny::textOutput(ns("r_version")),
    shiny::br(),
    DT::dataTableOutput(ns("sessionInfo"))
  )
}
    
#' versions Server Functions
#'
#' @noRd
mod_versions_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    package_list <- t(sapply(sessionInfo()[["loadedOnly"]], function(x){
      x[c("Package","Title","URL","Version")]
    }))[,-1]
    
    output$r_version <- shiny::renderText({
      sessionInfo()[1][["R.version"]][["version.string"]]
    })
    
    output$sessionInfo <- DT::renderDataTable({
      package_list
    })
 
  })
}
    
## To be copied in the UI
# mod_versions_ui("versions_ui_1")
    
## To be copied in the server
# mod_versions_server("versions_ui_1")
 
