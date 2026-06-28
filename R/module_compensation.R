library(shiny)
library(shinydashboard)
library(shinyjs)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

compensation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Utilisation de shinyjs pour d'éventuelles interactions avancées
    useShinyjs(),
    
    tabBox(
      title = tagList(icon("sliders-h"), "Outils de Compensation"),
      id = ns("comp_steps"), width = 12,
      
      # ── Onglet Transformation ────────────────────────────────────────────
      tabPanel("Transformation",
               fluidRow(
                 column(width = 4,
                        wellPanel(
                          h4("Choix des paramètres", style = "margin-top:0;"),
                          numericInput(ns("cofacteur"), "Cofacteur Arcsinh",
                                       value = 600, min = 1, step = 50),
                          hr(),
                          uiOutput(ns("ui_trans_file")),
                          uiOutput(ns("ui_trans_cx")),
                          uiOutput(ns("ui_trans_cy")),
                          hr(),
                          actionButton(ns("btn_apply_trans"), " Appliquer la transformation",
                                       class = "btn-primary", style = "width:100%; font-weight:bold;"),
                          br(), br(),
                          uiOutput(ns("ui_trans_status"))
                        )
                 ),
                 column(width = 8,
                        # Centrage du plot pour garantir le respect du format carré
                        box(title = "Visualisation", width = NULL, status = "info", solidHeader = TRUE,
                            div(style = "display: flex; justify-content: center;",
                                plotOutput(ns("plot_transformation"), width = "400px", height = "400px")
                            )
                        )
                 )
               )
      ),
      
      # ── Onglet Gating ────────────────────────────────────────────────────
      tabPanel("Gating",
               fluidRow(
                 column(width = 4,
                        wellPanel(
                          h4("Paramètres de Gating"),
                          uiOutput(ns("ui_gate_canal")),
                          uiOutput(ns("ui_gate_use_unstained_wrapper")),
                          hr(),
                          uiOutput(ns("ui_slider_neg")),
                          uiOutput(ns("ui_slider_pos")),
                          actionButton(ns("save_gate"), " Enregistrer les gates", class = "btn-success", style = "width:100%"),
                          br(), br(),
                          uiOutput(ns("ui_gates_status"))
                        )
                 ),
                 column(width = 8,
                        box(title = "Gating interactif", width = NULL, status = "danger", solidHeader = TRUE,
                            plotOutput(ns("plot_gating"), height = "400px")
                        )
                 )
               )
      )
      
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════

compensation_server <- function(id, pipeline, pipeline_version) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Accès réactif à l'objet CARROT
    carrot_obj <- reactive({
      pipeline_version()
      pipeline()
    })
    
    # Signaux locaux pour déclencher les mises à jour UI
    trans_done <- reactiveVal(0L)
    gates_done <- reactiveVal(0L)
    
    # Fonction utilitaire : Récupérer le flowFrame via le mapping canal -> fichier
    get_fcs_by_canal <- function(p, canal) {
      nom_cle <- p$mapping_canal_fichier[[canal]]
      if (is.null(nom_cle)) return(NULL)
      
      # Priorité aux données transformées
      if (!is.null(p$monomarques_trans) && !is.null(p$monomarques_trans[[nom_cle]])) {
        return(p$monomarques_trans[[nom_cle]])
      }
      return(p$tubes_monomarques[[nom_cle]])
    }
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 1 — TRANSFORMATION
    # ════════════════════════════════════════════════════════════════════════
    
    # Liste des fichiers disponibles (clés de tubes_monomarques)
    output$ui_trans_file <- renderUI({
      p <- carrot_obj()
      selectInput(ns("trans_file_sel"), "Tube à visualiser :", choices = names(p$tubes_monomarques))
    })
    
    output$ui_trans_cx <- renderUI({
      choices <- get_canaux_filtres(carrot_obj())
      selectInput(ns("trans_cx"), "Axe X :", choices = choices)
    })
    
    output$ui_trans_cy <- renderUI({
      choices <- get_canaux_filtres(carrot_obj())
      selectInput(ns("trans_cy"), "Axe Y :", choices = choices, selected = choices[min(2, length(choices))])
    })
    
    # Application de la transformation
    observeEvent(input$btn_apply_trans, {
      req(input$cofacteur)
      p <- carrot_obj()
      
      withProgress(message = "Transformation en cours...", value = 0.5, {
        # 1. On applique la transformation sur l'objet R6
        p$transformer_fcs(cofacteur = input$cofacteur)
        
        # 2. On notifie les autres composants du module que l'état a changé
        # Cela forcera le renderPlot et les autres éléments dépendants à se mettre à jour
        trans_done(trans_done() + 1L)
      })
    })
    
    output$ui_trans_status <- renderUI({
      trans_done()
      message("Debug: Mise à jour du statut transformation")
      p <- carrot_obj()
      if (is.null(p$monomarques_trans)) {
        tags$span(icon("times-circle"), "Données brutes (Non transformées)", style="color:red;")
      } else {
        tags$span(icon("check-circle"), "Transformation Arcsinh appliquée", style="color:green;")
      }
    })
    
    output$plot_transformation <- renderPlot({
      trans_done()
      req(input$trans_file_sel, input$trans_cx, input$trans_cy)
      p <- carrot_obj()
      fcs_t <- p$monomarques_trans[[input$trans_file_sel]]
      
      validate(need(!is.null(fcs_t), "Veuillez cliquer sur 'Appliquer la transformation'."))
      
      mat <- data.frame(x = flowCore::exprs(fcs_t)[, input$trans_cx],
                        y = flowCore::exprs(fcs_t)[, input$trans_cy])
      
      ggplot2::ggplot(mat, ggplot2::aes(x = x, y = y)) +
        ggpointdensity::geom_pointdensity(size = 0.2) +
        ggplot2::theme_bw()
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 2 — GATING INTERACTIF
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_gate_canal <- renderUI({
      p <- carrot_obj()
      trans_done() # Dépend de la transformation
      noms <- names(p$mapping_canal_fichier)
      selectInput(ns("gate_canal_tube"), "Canal à gater :", choices = noms[noms != "TUBE_UNSTAINED"])
    })
    
    output$ui_gate_use_unstained_wrapper <- renderUI({
      p <- carrot_obj()
      if ("TUBE_UNSTAINED" %in% names(p$tubes_monomarques)) {
        checkboxInput(ns("gate_use_unstained"), "Utiliser Unstained (Négatif)", value = TRUE)
      }
    })
    
    # Génération dynamique des bornes et sliders
    output$ui_slider_neg <- renderUI({
      req(input$gate_canal_tube)
      p <- carrot_obj()
      fcs <- get_fcs_by_canal(p, input$gate_canal_tube)
      rng <- range(flowCore::exprs(fcs)[, input$gate_canal_tube], na.rm = TRUE)
      sliderInput(ns("gate_slider_neg"), "Bornes Négatives", min=floor(rng[1]), max=ceiling(rng[2]), value=c(0, 2))
    })
    
    output$ui_slider_pos <- renderUI({
      req(input$gate_canal_tube)
      p <- carrot_obj()
      fcs <- get_fcs_by_canal(p, input$gate_canal_tube)
      rng <- range(flowCore::exprs(fcs)[, input$gate_canal_tube], na.rm = TRUE)
      sliderInput(ns("gate_slider_pos"), "Bornes Positives", min=floor(rng[1]), max=ceiling(rng[2]), value=c(4, 6))
    })
    
    output$plot_gating <- renderPlot({
      req(input$gate_canal_tube, input$gate_slider_neg, input$gate_slider_pos)
      p <- carrot_obj()
      p$graphiques_gates(
        nom_canal = input$gate_canal_tube,
        shiny_neg = input$gate_slider_neg,
        shiny_pos = input$gate_slider_pos,
        afficher_unstained_neg = isTRUE(input$gate_use_unstained)
      )
    })
    
    observeEvent(input$save_gate, {
      p <- carrot_obj()
      p$definir_et_extraire(input$gate_canal_tube, input$gate_slider_neg, input$gate_slider_pos, isTRUE(input$gate_use_unstained))
      gates_done(gates_done() + 1L)
    })
    
 
  })
}