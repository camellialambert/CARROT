library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(plotly)

options(shiny.maxRequestSize = 1000 * 1024^2) # 100 MB

source("~/Desktop/Institut_Cochin/Code/CARROT/pipeline_cytometrie.R")
source("~/Desktop/Institut_Cochin/Code/CARROT/R/module_import.R.R")
source("~/Desktop/Institut_Cochin/Code/CARROT/R/utils.R")
source("~/Desktop/Institut_Cochin/Code/CARROT/R/module_compensation.R")
source("~/Desktop/Institut_Cochin/Code/CARROT/R/module_qc.R")
source("~/Desktop/Institut_Cochin/Code/CARROT/R/module_pretraitement.R")

CANAUX_CONNUS <- c(
  "", "FITC-A", "Alexa Fluor 488-A", "Alexa Fluor 700-A",
  "PE-A", "PE-Cy5-A", "PE-Cy7-A", "PerCP-Cy5-5-A", "GFP-A",
  "APC-A", "APC-R700-A", "APC-H7-A",
  "AmCyan-A", "BV 650-A" , "BV 711-A"  ,
  "Pacific Blue-A", "BV 786-A", "DAPI-A", 
  "PE-Texas Red-A", "BV 421-A", "BV 605-A"
)

ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "CARROT"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs_sidebar",
      menuItem("Importation",    tabName = "import_tab",       icon = icon("file-import")),
      menuItem("Compensation",   tabName = "compensation_tab", icon = icon("calculator")),
      menuItem("Unmixing",       tabName = "unmixing_tab",     icon = icon("bolt")),
      menuItem("Quality Control",tabName = "nettoyage_tab",    icon = icon("broom")),
      menuItem("Prétraitement",  tabName = "pretraitement_tab",icon = icon("filter")),
      menuItem("Analyses",       tabName = "analyses_tab",     icon = icon("chart-pie"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .datatables, .dataTables_wrapper, table.dataTable tbody td { color: black !important; }
      .disabled-tab  { pointer-events: none; opacity: 0.4; }
      .btn-purple    { color: white; background-color: #605ca8; border-color: #555299; }
    "))),
    tabItems(
      tabItem(tabName = "import_tab",
              import_data_ui("mon_module_import")),
      tabItem(tabName = "compensation_tab",
              compensation_ui("mon_module_comp")),
      tabItem(tabName = "nettoyage_tab",
              qc_ui("mon_module_qc")),
      tabItem(tabName = "pretraitement_tab",
              pretraitement_ui("mon_module_pretrait")),
      tabItem(tabName = "unmixing_tab",  "Module Spectral en attente."),
      tabItem(tabName = "analyses_tab",  "Analyses statistiques en attente.")
    )
  )
)

server <- function(input, output, session) {
  
  # ── Deux reactiveVal partagés entre les modules ──
  my_pipeline         <- reactiveVal(CARROT$new())
  my_pipeline_version <- reactiveVal(0L)          # s'incrémente après charger_fcs()
  
  import_data_server("mon_module_import",
                     pipeline  = my_pipeline,
                     pipeline_version = my_pipeline_version,
                     canaux = CANAUX_CONNUS)
  
  compensation_server("mon_module_comp",
                      pipeline         = my_pipeline,
                      pipeline_version = my_pipeline_version)
  
  qc_server("mon_module_qc",
            pipeline         = my_pipeline,
            pipeline_version = my_pipeline_version)
  
  pretraitement_server("mon_module_pretrait",
                       pipeline         = my_pipeline,
                       pipeline_version = my_pipeline_version)
  
}

shinyApp(ui, server)