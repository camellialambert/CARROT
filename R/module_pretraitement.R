library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

pretraitement_ui <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    
    tags$style(HTML("
      .pretrait-instr {
        background:#eef6fb; border-left:3px solid #0077b6;
        padding:8px 12px; border-radius:4px;
        font-size:12px; color:#444; margin-top:6px;
      }
      .pretrait-instr b { color:#0077b6; }
      .gate-chip {
        display:inline-block; cursor:pointer;
        background:#e3f2fd; border:1px solid #0077b6;
        border-radius:12px; padding:3px 10px; margin:3px 2px;
        font-size:12px; color:#0077b6; font-weight:600;
      }
      .gate-chip:hover { background:#bbdefb; }
      .gate-chip.active { background:#0077b6; color:#fff; }
      .recap-row { font-size:12px; margin-bottom:6px; padding:6px 8px;
                   background:#f9f9f9; border-radius:4px; }
      .recap-row .label { font-weight:bold; color:#333; }
      .recap-row .value { float:right; }
      .source-info-box {
        background:#f3f0fb; border-left:3px solid #605ca8;
        padding:8px 12px; border-radius:4px;
        font-size:12px; color:#444; margin-bottom:8px;
      }
      .source-info-box b { color:#605ca8; }
      .canal-info-box {
        background:#fff8e1; border-left:3px solid #f9a825;
        padding:6px 10px; border-radius:4px;
        font-size:12px; color:#444; margin:8px 0;
      }
    ")),
    
    tabBox(
      title = tagList(icon("filter"), "Prétraitement"),
      id = ns("pretrait_steps"), width = 12,
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET — RETRAIT DES DÉBRIS
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("broom"), " Retrait des Débris"),
        fluidRow(
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE GAUCHE — Contrôles
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 wellPanel(
                   h4("Source des données"),
                   radioButtons(ns("source_debris"), NULL,
                                choices = c(
                                  "Compensées (brutes)" = "brutes",
                                  "Après PeacoQC"       = "peacoqc",
                                  "Après flowAI"        = "flowai"
                                ),
                                selected = "brutes"),
                   
                   hr(),
                   h4("Retrait des bordures"),
                   div(class = "pretrait-instr",
                       icon("info-circle"),
                       " Utilise automatiquement la dernière étape du pipeline (PeacoQC > flowAI > compensées)."),
                   checkboxInput(ns("activer_bordures"),
                                 "Retirer les événements de bordure (margins)", FALSE),
                   conditionalPanel(
                     condition = paste0("input['", ns("activer_bordures"), "']"),
                     uiOutput(ns("ui_canaux_bordures")),
                     br(),
                     actionButton(ns("btn_apply_bordures"),
                                  tagList(icon("play"), " Appliquer les bordures"),
                                  class = "btn-info btn-sm", style = "width:100%;")
                   ),
                   
                   hr(),
                   h4("Paramètres du gate"),
                   uiOutput(ns("ui_select_echantillon_debris")),
                   uiOutput(ns("ui_canal_x_debris")),
                   uiOutput(ns("ui_canal_y_debris")),
                   
                   hr(),
                   h4("Nommer et valider"),
                   textInput(ns("nom_gate"), "Nom du gate :", value = "Débris",
                             placeholder = "ex: Débris, Cellules vivantes..."),
                   
                   div(class = "pretrait-instr",
                       icon("hand-pointer"),
                       " Cliquez sur le graphique pour ajouter des sommets. ",
                       tags$b("Double-clic"), " pour fermer le polygone."),
                   br(),
                   
                   fluidRow(
                     column(6, actionButton(ns("btn_undo_sommet"),
                                            tagList(icon("undo"), " Annuler sommet"),
                                            class = "btn-default btn-sm", style = "width:100%;")),
                     column(6, actionButton(ns("btn_reset_gate"),
                                            tagList(icon("eraser"), " Effacer"),
                                            class = "btn-default btn-sm", style = "width:100%;"))
                   ),
                   br(),
                   actionButton(ns("btn_save_gate"),
                                tagList(icon("check-circle"), " Enregistrer le gate"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE CENTRALE — Graphiques
          # ──────────────────────────────────────────────────────────────────
          column(width = 6,
                 box(title = tagList(icon("draw-polygon"), " Dessin interactif du gate"),
                     width = NULL, status = "info", solidHeader = TRUE,
                     plotlyOutput(ns("plot_gate_dessin"), height = "480px")
                 ),
                 
                 box(title = uiOutput(ns("titre_resultat")), width = NULL,
                     status = "success", solidHeader = TRUE,
                     uiOutput(ns("ui_gate_chips")),
                     plotOutput(ns("plot_gate_resultat"), width = "100%", height = "420px")
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE DROITE — Résumé cumulatif
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 box(title = "Résumé du pipeline", width = NULL,
                     status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_recap_pipeline"))
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET — RETRAIT DES DOUBLETS
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("clone"), " Retrait des Doublets"),
        fluidRow(
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE GAUCHE — Contrôles
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 wellPanel(
                   h4("Étape de discrimination"),
                   radioButtons(ns("etape_doublet"), NULL,
                                choices = c(
                                  "Doublets FSC (Forward Scatter)" = "FSC",
                                  "Doublets SSC (Side Scatter)"    = "SSC"
                                ),
                                selected = "FSC"),
                   uiOutput(ns("ui_source_info_doublet")),
                   
                   hr(),
                   uiOutput(ns("ui_select_echantillon_doublet")),
                   selectInput(ns("axe_doublet"), "Axe de discrimination :",
                               choices = c(
                                 "Height vs Area (H_A)"  = "H_A",
                                 "Width vs Area (W_A)"   = "W_A",
                                 "Height vs Width (H_W)" = "H_W"
                               ),
                               selected = "H_A"),
                   uiOutput(ns("ui_canaux_info_doublet")),
                   
                   hr(),
                   h4("Méthode de retrait"),
                   radioButtons(ns("methode_doublet"), NULL,
                                choices = c(
                                  "Automatique (statistique, MAD)" = "auto",
                                  "Gating manuel"                  = "gate"
                                ),
                                selected = "auto"),
                   
                   # ── Méthode automatique ──────────────────────────────────
                   conditionalPanel(
                     condition = paste0("input['", ns("methode_doublet"), "'] == 'auto'"),
                     div(class = "pretrait-instr",
                         icon("info-circle"),
                         " Le seuil d'exclusion est fixé à la médiane du ratio + ",
                         tags$b("facteur de sensibilité"), " x MAD. Plus le facteur est ",
                         "élevé, plus le filtre est permissif (moins de cellules retirées)."),
                     numericInput(ns("facteur_sensibilite_doublet"),
                                  "Facteur de sensibilité :",
                                  value = 4, min = 0.5, step = 0.5),
                     checkboxInput(ns("doublet_tous_echantillons_auto"),
                                   "Appliquer à tous les échantillons", TRUE),
                     conditionalPanel(
                       condition = paste0("!input['", ns("doublet_tous_echantillons_auto"), "']"),
                       uiOutput(ns("ui_select_echantillon_cible_auto"))
                     ),
                     br(),
                     actionButton(ns("btn_appliquer_doublet_auto"),
                                  tagList(icon("play"), " Appliquer le retrait automatique"),
                                  class = "btn-success", style = "width:100%; font-weight:bold;")
                   ),
                   
                   # ── Méthode par gating manuel ────────────────────────────
                   conditionalPanel(
                     condition = paste0("input['", ns("methode_doublet"), "'] == 'gate'"),
                     numericInput(ns("doublet_max_points"), "Points affichés :",
                                  value = 10000, min = 1000, step = 1000),
                     div(class = "pretrait-instr",
                         icon("hand-pointer"),
                         " Cliquez sur le graphique pour ajouter des sommets. ",
                         tags$b("Double-clic"), " pour fermer le polygone."),
                     br(),
                     fluidRow(
                       column(6, actionButton(ns("btn_undo_sommet_doublet"),
                                              tagList(icon("undo"), " Annuler sommet"),
                                              class = "btn-default btn-sm", style = "width:100%;")),
                       column(6, actionButton(ns("btn_reset_gate_doublet"),
                                              tagList(icon("eraser"), " Effacer"),
                                              class = "btn-default btn-sm", style = "width:100%;"))
                     ),
                     br(),
                     checkboxInput(ns("doublet_tous_echantillons_gate"),
                                   "Appliquer le polygone à tous les échantillons", TRUE),
                     conditionalPanel(
                       condition = paste0("!input['", ns("doublet_tous_echantillons_gate"), "']"),
                       uiOutput(ns("ui_select_echantillon_cible_gate"))
                     ),
                     br(),
                     actionButton(ns("btn_save_gate_doublet"),
                                  tagList(icon("check-circle"), " Enregistrer le gate"),
                                  class = "btn-success", style = "width:100%; font-weight:bold;")
                   )
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE CENTRALE — Graphiques
          # ──────────────────────────────────────────────────────────────────
          column(width = 6,
                 conditionalPanel(
                   condition = paste0("input['", ns("methode_doublet"), "'] == 'gate'"),
                   box(title = tagList(icon("draw-polygon"), " Dessin interactif du gate"),
                       width = NULL, status = "info", solidHeader = TRUE,
                       plotlyOutput(ns("plot_gate_doublet_dessin"), height = "480px")
                   )
                 ),
                 box(title = uiOutput(ns("titre_resultat_doublet")), width = NULL,
                     status = "success", solidHeader = TRUE,
                     plotOutput(ns("plot_doublet_resultat"), width = "100%", height = "420px")
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE DROITE — Résumé cumulatif
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 box(title = "Résumé du pipeline", width = NULL,
                     status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_recap_pipeline_doublet"))
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET — RETRAIT DES CELLULES MORTES (VIABILITÉ)
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("heart-pulse"), " Cellules mortes"),
        fluidRow(
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE GAUCHE — Contrôles
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 wellPanel(
                   h4("Canal de viabilité"),
                   checkboxInput(ns("viabilite_afficher_tous_canaux"),
                                 "Afficher tous les canaux", FALSE),
                   uiOutput(ns("ui_canal_viabilite")),
                   
                   conditionalPanel(
                     condition = paste0("typeof input['", ns("canal_viabilite"), "'] !== 'undefined' && input['",
                                        ns("canal_viabilite"), "'] != ''"),
                     
                     hr(),
                     h4("1. Transformation Arcsinh"),
                     div(class = "pretrait-instr",
                         icon("info-circle"),
                         " La transformation Arcsinh doit être appliquée sur le canal de viabilité ",
                         "avant de pouvoir tracer le gate des cellules vivantes."),
                     numericInput(ns("cofacteur_viabilite"), "Cofacteur :",
                                  value = 400, min = 1, step = 10),
                     checkboxInput(ns("viabilite_tous_echantillons_transfo"),
                                   "Appliquer à tous les échantillons", TRUE),
                     conditionalPanel(
                       condition = paste0("!input['", ns("viabilite_tous_echantillons_transfo"), "']"),
                       uiOutput(ns("ui_select_echantillon_viabilite_transfo_cible"))
                     ),
                     actionButton(ns("btn_appliquer_transfo_viabilite"),
                                  tagList(icon("play"), " Appliquer la transformation Arcsinh"),
                                  class = "btn-info", style = "width:100%; font-weight:bold;"),
                     br(), br(),
                     uiOutput(ns("ui_transfo_viabilite_status")),
                     
                     hr(),
                     h4("2. Gating des cellules vivantes"),
                     uiOutput(ns("ui_select_echantillon_viabilite")),
                     numericInput(ns("viabilite_max_points"), "Points affichés :",
                                  value = 10000, min = 1000, step = 1000),
                     div(class = "pretrait-instr",
                         icon("hand-pointer"),
                         " Cliquez sur le graphique pour ajouter des sommets autour des cellules ",
                         tags$b("vivantes"), ". ", tags$b("Double-clic"), " pour fermer le polygone."),
                     br(),
                     fluidRow(
                       column(6, actionButton(ns("btn_undo_sommet_viabilite"),
                                              tagList(icon("undo"), " Annuler sommet"),
                                              class = "btn-default btn-sm", style = "width:100%;")),
                       column(6, actionButton(ns("btn_reset_gate_viabilite"),
                                              tagList(icon("eraser"), " Effacer"),
                                              class = "btn-default btn-sm", style = "width:100%;"))
                     ),
                     br(),
                     checkboxInput(ns("viabilite_tous_echantillons_gate"),
                                   "Appliquer le polygone à tous les échantillons", TRUE),
                     conditionalPanel(
                       condition = paste0("!input['", ns("viabilite_tous_echantillons_gate"), "']"),
                       uiOutput(ns("ui_select_echantillon_cible_gate_viabilite"))
                     ),
                     br(),
                     actionButton(ns("btn_save_gate_viabilite"),
                                  tagList(icon("check-circle"), " Enregistrer le gate"),
                                  class = "btn-success", style = "width:100%; font-weight:bold;")
                   )
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE CENTRALE — Graphiques
          # ──────────────────────────────────────────────────────────────────
          column(width = 6,
                 box(title = tagList(icon("draw-polygon"), " Dessin interactif du gate"),
                     width = NULL, status = "info", solidHeader = TRUE,
                     plotlyOutput(ns("plot_gate_viabilite_dessin"), height = "480px")
                 ),
                 box(title = "Résultat — cellules vivantes", width = NULL,
                     status = "success", solidHeader = TRUE,
                     plotOutput(ns("plot_viabilite_resultat"), width = "100%", height = "420px")
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE DROITE — Résumé cumulatif
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 box(title = "Résumé du pipeline", width = NULL,
                     status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_recap_pipeline_viabilite"))
                 )
          )
        )
      )
      
    ) # /tabBox
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════

pretraitement_server <- function(id, pipeline, pipeline_version) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ── Accès réactif au pipeline ─────────────────────────────────────────────
    carrot_obj <- reactive({
      pipeline_version()
      pipeline()
    })
    
    # ── Signaux et état local ─────────────────────────────────────────────────
    bordures_trigger  <- reactiveVal(0L)
    gates_trigger     <- reactiveVal(0L)
    doublet_trigger   <- reactiveVal(0L)
    
    # Sommets du polygone en cours de dessin (débris)
    sommets_rv        <- reactiveVal(list())
    # Nom du gate actuellement affiché dans le panneau résultat (débris)
    gate_actif_rv     <- reactiveVal(NULL)
    
    # Sommets du polygone en cours de dessin (doublets)
    sommets_doublet_rv <- reactiveVal(list())
    
    # Signaux et état local (viabilité)
    viabilite_trigger   <- reactiveVal(0L)
    sommets_viabilite_rv <- reactiveVal(list())
    # Mémorise le dernier canal effectivement transformé via cet onglet,
    # afin de n'autoriser le gating qu'une fois la transformation appliquée
    # sur le canal actuellement sélectionné.
    canal_transforme_rv <- reactiveVal(NULL)
    
    # ════════════════════════════════════════════════════════════════════════
    # RÉSOLUTION DE LA SOURCE (DÉBRIS)
    # ════════════════════════════════════════════════════════════════════════
    
    obtenir_source <- function(p, choix) {
      if (choix == "peacoqc" && length(p$post_PeacoQC) > 0) return(p$post_PeacoQC)
      if (choix == "flowai"  && length(p$post_flowAI)  > 0) return(p$post_flowAI)
      return(p$get_derniere_source())
    }
    
    # ════════════════════════════════════════════════════════════════════════
    # BORDURES (MARGINS)
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_canaux_bordures <- renderUI({
      p <- carrot_obj()
      src <- p$get_derniere_source()
      req(length(src) > 0)
      cols <- flowCore::colnames(src[[1]])
      selectInput(ns("canaux_bordures_sel"), "Canaux :", choices = cols,
                  selected = cols, multiple = TRUE)
    })
    
    observeEvent(input$btn_apply_bordures, {
      p  <- carrot_obj()
      ch <- input$canaux_bordures_sel
      req(length(ch) >= 1)
      withProgress(message = "Retrait des bordures...", value = 0.4, {
        tryCatch({
          p$retirer_les_bordures(
            canal1 = ch[1],
            canal2 = if (length(ch) >= 2) ch[2] else ch[1]
          )
          pipeline(p)
          bordures_trigger(bordures_trigger() + 1L)
          showNotification("Bordures retirées.", type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # SÉLECTEURS DYNAMIQUES (DÉBRIS)
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_select_echantillon_debris <- renderUI({
      bordures_trigger(); gates_trigger()
      p  <- carrot_obj()
      src <- obtenir_source(p, input$source_debris %||% "brutes")
      req(length(src) > 0)
      selectInput(ns("sel_ech_debris"), "Échantillon :", choices = names(src))
    })
    
    output$ui_canal_x_debris <- renderUI({
      bordures_trigger()
      p   <- carrot_obj()
      src <- obtenir_source(p, input$source_debris %||% "brutes")
      req(length(src) > 0, input$sel_ech_debris)
      fcs <- src[[input$sel_ech_debris]]
      req(!is.null(fcs))
      cols <- flowCore::colnames(fcs)
      def  <- grep("FSC", cols, value = TRUE, ignore.case = TRUE)[1]
      selectInput(ns("canal_x_debris"), "Canal X :", choices = cols,
                  selected = if (!is.na(def)) def else cols[1])
    })
    
    output$ui_canal_y_debris <- renderUI({
      bordures_trigger()
      p   <- carrot_obj()
      src <- obtenir_source(p, input$source_debris %||% "brutes")
      req(length(src) > 0, input$sel_ech_debris)
      fcs <- src[[input$sel_ech_debris]]
      req(!is.null(fcs))
      cols <- flowCore::colnames(fcs)
      def  <- grep("SSC", cols, value = TRUE, ignore.case = TRUE)[1]
      selectInput(ns("canal_y_debris"), "Canal Y :", choices = cols,
                  selected = if (!is.na(def)) def else cols[min(2, length(cols))])
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # DONNÉES POUR LE GRAPHIQUE INTERACTIF (DÉBRIS)
    # ════════════════════════════════════════════════════════════════════════
    
    donnees_plot <- reactive({
      bordures_trigger()
      p <- carrot_obj()
      req(input$sel_ech_debris, input$canal_x_debris, input$canal_y_debris)
      src <- obtenir_source(p, input$source_debris %||% "brutes")
      fcs <- src[[input$sel_ech_debris]]
      req(!is.null(fcs))
      cx <- input$canal_x_debris
      cy <- input$canal_y_debris
      req(cx %in% flowCore::colnames(fcs), cy %in% flowCore::colnames(fcs))
      mat <- flowCore::exprs(fcs)[, c(cx, cy), drop = FALSE]
      df <- as.data.frame(mat)
      colnames(df) <- c("X", "Y")
      df
    })
    
    # Réinitialise le polygone si l'utilisateur change de canal/échantillon/source
    observeEvent(list(input$sel_ech_debris, input$canal_x_debris,
                      input$canal_y_debris, input$source_debris), {
                        sommets_rv(list())
                      }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # GESTION DU POLYGONE — clic par clic + double-clic pour fermer (DÉBRIS)
    # ════════════════════════════════════════════════════════════════════════
    
    # Ajout d'un sommet via plotly_click
    observeEvent(event_data("plotly_click", source = ns("plot_gate_dessin")), {
      ev <- event_data("plotly_click", source = ns("plot_gate_dessin"))
      req(ev)
      soms <- sommets_rv()
      soms[[length(soms) + 1]] <- list(x = ev$x, y = ev$y)
      sommets_rv(soms)
    })
    
    # Fermeture du polygone via double-clic (ne rajoute pas de sommet)
    observeEvent(event_data("plotly_doubleclick", source = ns("plot_gate_dessin")), {
      soms <- sommets_rv()
      if (length(soms) >= 3) {
        showNotification(
          paste0("Polygone fermé (", length(soms), " sommets). Vérifiez et enregistrez le gate."),
          type = "message"
        )
      } else {
        showNotification("Tracez au moins 3 sommets avant de fermer le polygone.", type = "warning")
      }
    })
    
    # Annuler le dernier sommet
    observeEvent(input$btn_undo_sommet, {
      soms <- sommets_rv()
      if (length(soms) > 0) sommets_rv(soms[-length(soms)])
    })
    
    # Effacer tout le polygone
    observeEvent(input$btn_reset_gate, {
      sommets_rv(list())
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # GRAPHIQUE INTERACTIF — rendu de base (DÉBRIS)
    # Ce rendu ne dépend QUE des données (échantillon/canaux/source/N points),
    # PAS des sommets du polygone : le tracé du polygone et les poignées de
    # sommets sont ensuite injectés via plotlyProxy, sans jamais redessiner
    # le nuage de points complet. C'est ce qui rend le placement des points
    # fluide et réactif, y compris avec 10 000+ événements affichés.
    # ════════════════════════════════════════════════════════════════════════
    
    output$plot_gate_dessin <- renderPlotly({
      df <- donnees_plot()
      req(nrow(df) > 0)
      
      p_obj <- carrot_obj()
      src   <- obtenir_source(p_obj, input$source_debris %||% "brutes")
      fcs   <- src[[input$sel_ech_debris]]
      lbl_x <- p_obj$get_label(fcs, input$canal_x_debris)
      lbl_y <- p_obj$get_label(fcs, input$canal_y_debris)
      
      # Densité par binning raster (rapide, indépendant du nombre d'événements affichés),
      # rendue en heatmap plotly avec la même palette pseudo-spectrale que le reste de l'app,
      # plutôt qu'un nuage de points uni (plus lourd à afficher/interagir avec beaucoup d'événements).
      lim_x <- range(df$X, na.rm = TRUE)
      lim_y <- range(df$Y, na.rm = TRUE)
      dens  <- calculer_densite_matrice_plotly(df$X, df$Y, lim_x, lim_y)
      
      plt <- plot_ly(source = ns("plot_gate_dessin"))
      if (!is.null(dens)) {
        plt <- plt %>%
          add_trace(x = dens$x, y = dens$y, z = dens$z, type = "heatmap",
                    colorscale = COLORSCALE_DENSITE_PLOTLY, showscale = FALSE,
                    hoverinfo = "none")
      } else {
        plt <- plt %>%
          add_trace(data = df, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    name = "Cellules",
                    marker = list(size = 2, color = "#0077b6", opacity = 0.4),
                    hoverinfo = "none")
      }
      plt %>%
        # Trace vide dès le départ : elle sera mise à jour par plotlyProxy
        # à chaque ajout/déplacement de sommet, sans toucher au fond de densité.
        add_trace(x = numeric(0), y = numeric(0), type = "scatter",
                  mode = "lines+markers", name = "Gate en cours",
                  line   = list(color = "#e65100", width = 2, dash = "dot"),
                  marker = list(size = 5, color = "#e65100", symbol = "circle"),
                  hoverinfo = "x+y") %>%
        event_register("plotly_click") %>%
        event_register("plotly_doubleclick") %>%
        event_register("plotly_relayout") %>%
        layout(
          dragmode = "zoom",
          xaxis = list(title = lbl_x),
          yaxis = list(title = lbl_y),
          legend = list(orientation = "h", y = -0.15),
          margin = list(b = 60),
          shapes = list()
        ) %>%
        config(displayModeBar = TRUE, editable = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE)
    })
    
    # Pousse les sommets courants vers le graphique déjà rendu, sans redraw complet
    observeEvent(sommets_rv(), {
      soms <- sommets_rv()
      proxy <- plotlyProxy(ns("plot_gate_dessin"), session)
      
      if (length(soms) == 0) {
        plotlyProxyInvoke(proxy, "restyle", list(x = list(list()), y = list(list())), list(1))
        plotlyProxyInvoke(proxy, "relayout", list(shapes = list()))
        return(invisible(NULL))
      }
      
      xs <- c(vapply(soms, `[[`, numeric(1), "x"), soms[[1]]$x)   # ferme visuellement
      ys <- c(vapply(soms, `[[`, numeric(1), "y"), soms[[1]]$y)
      
      plotlyProxyInvoke(proxy, "restyle", list(x = list(xs), y = list(ys)), list(1))
      
      if (length(soms) >= 2) {
        df <- isolate(donnees_plot())
        rx <- diff(range(df$X, na.rm = TRUE)); if (!is.finite(rx) || rx == 0) rx <- 1
        ry <- diff(range(df$Y, na.rm = TRUE)); if (!is.finite(ry) || ry == 0) ry <- 1
        
        shapes <- lapply(seq_along(soms), function(i) {
          list(type = "circle",
               x0 = soms[[i]]$x - 0.012 * rx, x1 = soms[[i]]$x + 0.012 * rx,
               y0 = soms[[i]]$y - 0.012 * ry, y1 = soms[[i]]$y + 0.012 * ry,
               xref = "x", yref = "y",
               fillcolor = "rgba(230,101,0,0.3)",
               line = list(color = "#e65100"))
        })
        plotlyProxyInvoke(proxy, "relayout", list(shapes = shapes))
      } else {
        plotlyProxyInvoke(proxy, "relayout", list(shapes = list()))
      }
    }, ignoreNULL = FALSE)
    
    # Capture les déplacements de shapes (sommets éditables)
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_dessin")), {
      ev   <- event_data("plotly_relayout", source = ns("plot_gate_dessin"))
      req(ev)
      soms <- sommets_rv()
      
      # Les shapes[i].x0 / shapes[i].y0 sont renvoyés quand on déplace un cercle
      modifie <- FALSE
      for (nm in names(ev)) {
        m_idx <- regmatches(nm, regexpr("[0-9]+", nm))
        if (length(m_idx) == 0) next
        i <- as.integer(m_idx) + 1L   # 0-indexé en JS → 1-indexé en R
        if (i < 1 || i > length(soms)) next
        
        if (grepl("\\.x0$", nm)) {
          soms[[i]]$x <- ev[[nm]] + 0     # centrage : on récupère le bord gauche du cercle
          modifie <- TRUE
        }
        if (grepl("\\.y0$", nm)) {
          soms[[i]]$y <- ev[[nm]] + 0
          modifie <- TRUE
        }
      }
      if (modifie) sommets_rv(soms)
    }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # ENREGISTREMENT DU GATE (DÉBRIS)
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(input$btn_save_gate, {
      soms    <- sommets_rv()
      nom_g   <- trimws(input$nom_gate %||% "")
      
      if (length(soms) < 3) {
        showNotification("Tracez au moins 3 sommets pour définir un gate.", type = "error")
        return(invisible(NULL))
      }
      if (nchar(nom_g) == 0) {
        showNotification("Donnez un nom au gate.", type = "error")
        return(invisible(NULL))
      }
      
      p <- carrot_obj()
      
      # Construction de la matrice polygone à partir des sommets courants
      xs <- sapply(soms, `[[`, "x")
      ys <- sapply(soms, `[[`, "y")
      cx <- input$canal_x_debris
      cy <- input$canal_y_debris
      
      # Les shapes éditables décalent le centre d'un rayon → on utilise directement x0/y0
      # (le centre réel est x0 + rayon, mais pour polygonGate on veut les vrais sommets)
      # On calcule l'enveloppe convexe pour garantir un polygone propre
      idx_hull <- grDevices::chull(xs, ys)
      mat_poly <- cbind(xs[idx_hull], ys[idx_hull])
      colnames(mat_poly) <- c(cx, cy)
      
      withProgress(message = paste0("Application du gate '", nom_g, "'..."), value = 0.4, {
        tryCatch({
          p$appliquer_gate_nomme(
            nom_gate       = nom_g,
            matrice_points = mat_poly,
            canal_x        = cx,
            canal_y        = cy,
            source_nettoyage = input$source_debris %||% "brutes"
          )
          pipeline(p)
          gates_trigger(gates_trigger() + 1L)
          gate_actif_rv(nom_g)
          sommets_rv(list())   # remet à zéro pour le prochain gate
          showNotification(paste0("Gate '", nom_g, "' enregistré."), type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # PANNEAU RÉSULTAT — chips de navigation + graphique (DÉBRIS)
    # ════════════════════════════════════════════════════════════════════════
    
    output$titre_resultat <- renderUI({
      ga <- gate_actif_rv()
      if (is.null(ga)) "Résultat après filtration"
      else paste0("Résultat — gate : ", ga)
    })
    
    # Chips cliquables (un bouton par gate nommé enregistré)
    output$ui_gate_chips <- renderUI({
      gates_trigger()
      p  <- carrot_obj()
      nms <- names(p$gates_history)
      if (length(nms) == 0) return(NULL)
      
      ga <- gate_actif_rv()
      
      chips <- lapply(nms, function(nm) {
        cls <- paste0("gate-chip", if (!is.null(ga) && ga == nm) " active" else "")
        actionLink(ns(paste0("chip_", gsub("[^a-zA-Z0-9]", "_", nm))),
                   label = nm,
                   class = cls,
                   onclick = paste0("Shiny.setInputValue('", ns("gate_chip_clicked"), "', '",
                                    nm, "', {priority: 'event'})"))
      })
      
      tagList(
        div(style = "margin-bottom:8px; font-size:11px; color:#888;",
            "Cliquer sur un gate pour l'afficher ou le modifier :"),
        do.call(tagList, chips)
      )
    })
    
    # Réaction au clic sur un chip → charge le polygone pour réajustement
    observeEvent(input$gate_chip_clicked, {
      nm <- input$gate_chip_clicked
      req(nchar(nm) > 0)
      
      p <- carrot_obj()
      infos_gate <- p$gates_history[[nm]]
      req(!is.null(infos_gate))
      
      # Prend les coordonnées du premier échantillon (référence)
      premier_ech <- infos_gate[[input$sel_ech_debris %||% names(infos_gate)[1]]]
      if (is.null(premier_ech)) premier_ech <- infos_gate[[1]]
      req(!is.null(premier_ech))
      
      mat <- premier_ech$polygone
      soms <- lapply(seq_len(nrow(mat)), function(i) {
        list(x = mat[i, 1], y = mat[i, 2])
      })
      
      sommets_rv(soms)
      gate_actif_rv(nm)
      updateTextInput(session, "nom_gate", value = nm)
    }, ignoreInit = TRUE)
    
    # Graphique résultat (ggplot via visualiser_debris)
    output$plot_gate_resultat <- renderPlot({
      gates_trigger()
      p  <- carrot_obj()
      ga <- gate_actif_rv()
      req(input$sel_ech_debris)
      
      if (is.null(ga) || length(p$gates_history) == 0) {
        validate("Enregistrez un gate pour afficher le résultat.")
      }
      
      # Assure que post_debris correspond au gate actif avant de visualiser
      infos <- p$gates_history[[ga]]
      req(!is.null(infos))
      p$post_debris <- lapply(infos, function(r) r$post_data)
      
      res <- tryCatch(
        p$visualiser_debris(nom_echantillon = input$sel_ech_debris),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
      )
      
      if (is.null(res)) {
        plot.new(); text(0.5, 0.5, "Aucune donnée.", cex = 1.1, col = "grey40"); return()
      }
      print(res)
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # RÉSOLUTION DE LA SOURCE ET DES CANAUX (DOUBLETS)
    # ════════════════════════════════════════════════════════════════════════
    
    # Reproduit exactement la chaîne de priorité utilisée côté pipeline
    # (retirer_doublets_FSC/SSC, gate_les_doublets_FSC/SSC) : pour l'étape FSC,
    # on part du retrait des débris ou, à défaut, de la dernière étape valide
    # (QC, importation...). Pour l'étape SSC, on chaîne après le FSC si déjà
    # fait, sinon on retombe sur les débris ou la dernière étape valide.
    obtenir_source_doublet <- function(p, etape) {
      if (identical(etape, "SSC")) {
        if (!is.null(p$post_doublets_FSC) && length(p$post_doublets_FSC) > 0) return(p$post_doublets_FSC)
        if (!is.null(p$post_debris) && length(p$post_debris) > 0) return(p$post_debris)
        return(p$get_derniere_source())
      }
      if (!is.null(p$post_debris) && length(p$post_debris) > 0) return(p$post_debris)
      return(p$get_derniere_source())
    }
    
    libelle_source_doublet <- function(p, etape) {
      if (identical(etape, "SSC")) {
        if (!is.null(p$post_doublets_FSC) && length(p$post_doublets_FSC) > 0) return("post-retrait doublets FSC")
        if (!is.null(p$post_debris) && length(p$post_debris) > 0) return("post-retrait débris")
        return("dernière étape disponible (QC / importation)")
      }
      if (!is.null(p$post_debris) && length(p$post_debris) > 0) return("post-retrait débris")
      return("dernière étape disponible (QC / importation)")
    }
    
    # Détermine les deux canaux (X, Y) utilisés pour l'axe de discrimination choisi,
    # en cohérence stricte avec la logique interne du pipeline (préfixe FSC/SSC).
    resoudre_canaux_doublet <- function(cols, etape, axe) {
      prefixe <- etape  # "FSC" ou "SSC"
      if (axe == "H_A") {
        cx <- grep(paste0(prefixe, "-H"), cols, value = TRUE, ignore.case = TRUE)[1]
        cy <- grep(paste0(prefixe, "-A"), cols, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe == "W_A") {
        cx <- grep(paste0(prefixe, "-W"), cols, value = TRUE, ignore.case = TRUE)[1]
        cy <- grep(paste0(prefixe, "-A"), cols, value = TRUE, ignore.case = TRUE)[1]
      } else {
        cx <- grep(paste0(prefixe, "-H"), cols, value = TRUE, ignore.case = TRUE)[1]
        cy <- grep(paste0(prefixe, "-W"), cols, value = TRUE, ignore.case = TRUE)[1]
      }
      c(x = cx, y = cy)
    }
    
    output$ui_source_info_doublet <- renderUI({
      doublet_trigger(); gates_trigger(); bordures_trigger()
      p <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      div(class = "source-info-box",
          icon("info-circle"), " Données utilisées : ",
          tags$b(libelle_source_doublet(p, etape)), ".")
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # SÉLECTEURS DYNAMIQUES (DOUBLETS)
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_select_echantillon_doublet <- renderUI({
      doublet_trigger(); gates_trigger(); bordures_trigger()
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      src   <- obtenir_source_doublet(p, etape)
      req(length(src) > 0)
      selectInput(ns("sel_ech_doublet"), "Échantillon (référence pour le tracé) :",
                  choices = names(src))
    })
    
    output$ui_select_echantillon_cible_auto <- renderUI({
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      src   <- obtenir_source_doublet(p, etape)
      req(length(src) > 0)
      selectInput(ns("sel_ech_cible_auto"), "Échantillon ciblé :", choices = names(src))
    })
    
    output$ui_select_echantillon_cible_gate <- renderUI({
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      src   <- obtenir_source_doublet(p, etape)
      req(length(src) > 0)
      selectInput(ns("sel_ech_cible_gate"), "Échantillon ciblé :", choices = names(src))
    })
    
    output$ui_canaux_info_doublet <- renderUI({
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      axe   <- input$axe_doublet %||% "H_A"
      src   <- obtenir_source_doublet(p, etape)
      req(length(src) > 0)
      fcs <- src[[1]]
      req(!is.null(fcs))
      canaux <- resoudre_canaux_doublet(flowCore::colnames(fcs), etape, axe)
      if (is.na(canaux["x"]) || is.na(canaux["y"])) {
        return(div(class = "canal-info-box",
                   icon("exclamation-triangle"),
                   paste0(" Canaux ", etape, " introuvables pour l'axe ", axe, ".")))
      }
      div(class = "canal-info-box",
          icon("crosshairs"), " Canal X : ", tags$b(canaux["x"]),
          " — Canal Y : ", tags$b(canaux["y"]))
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # DONNÉES POUR LE GRAPHIQUE INTERACTIF (DOUBLETS — méthode gating)
    # ════════════════════════════════════════════════════════════════════════
    
    donnees_plot_doublet <- reactive({
      doublet_trigger()
      p     <- carrot_obj()
      req(input$sel_ech_doublet)
      etape <- input$etape_doublet %||% "FSC"
      axe   <- input$axe_doublet %||% "H_A"
      src   <- obtenir_source_doublet(p, etape)
      fcs   <- src[[input$sel_ech_doublet]]
      req(!is.null(fcs))
      canaux <- resoudre_canaux_doublet(flowCore::colnames(fcs), etape, axe)
      req(!is.na(canaux["x"]), !is.na(canaux["y"]))
      mat <- flowCore::exprs(fcs)[, c(canaux["x"], canaux["y"]), drop = FALSE]
      n    <- nrow(mat)
      nmax <- input$doublet_max_points %||% 10000
      if (!is.null(p$seed)) set.seed(p$seed)
      idx <- if (n > nmax) sample(seq_len(n), nmax) else seq_len(n)
      df  <- as.data.frame(mat[idx, , drop = FALSE])
      colnames(df) <- c("X", "Y")
      attr(df, "canal_x") <- unname(canaux["x"])
      attr(df, "canal_y") <- unname(canaux["y"])
      df
    })
    
    # Réinitialise le polygone si l'utilisateur change étape/axe/échantillon/méthode
    observeEvent(list(input$etape_doublet, input$axe_doublet,
                      input$sel_ech_doublet, input$methode_doublet), {
                        sommets_doublet_rv(list())
                      }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # GESTION DU POLYGONE — clic par clic + double-clic pour fermer (DOUBLETS)
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(event_data("plotly_click", source = ns("plot_gate_doublet_dessin")), {
      ev <- event_data("plotly_click", source = ns("plot_gate_doublet_dessin"))
      req(ev)
      soms <- sommets_doublet_rv()
      soms[[length(soms) + 1]] <- list(x = ev$x, y = ev$y)
      sommets_doublet_rv(soms)
    })
    
    observeEvent(event_data("plotly_doubleclick", source = ns("plot_gate_doublet_dessin")), {
      soms <- sommets_doublet_rv()
      if (length(soms) >= 3) {
        showNotification(
          paste0("Polygone fermé (", length(soms), " sommets). Vérifiez et enregistrez le gate."),
          type = "message"
        )
      } else {
        showNotification("Tracez au moins 3 sommets avant de fermer le polygone.", type = "warning")
      }
    })
    
    observeEvent(input$btn_undo_sommet_doublet, {
      soms <- sommets_doublet_rv()
      if (length(soms) > 0) sommets_doublet_rv(soms[-length(soms)])
    })
    
    observeEvent(input$btn_reset_gate_doublet, {
      sommets_doublet_rv(list())
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # GRAPHIQUE INTERACTIF — rendu de base + plotlyProxy (DOUBLETS)
    # Même principe que pour les débris : le nuage de points n'est redessiné
    # que si les données changent (échantillon/axe/étape), jamais lors de
    # l'ajout ou du déplacement d'un sommet.
    # ════════════════════════════════════════════════════════════════════════
    
    output$plot_gate_doublet_dessin <- renderPlotly({
      df <- donnees_plot_doublet()
      req(nrow(df) > 0)
      
      p_obj <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      src   <- obtenir_source_doublet(p_obj, etape)
      fcs   <- src[[input$sel_ech_doublet]]
      lbl_x <- p_obj$get_label(fcs, attr(df, "canal_x"))
      lbl_y <- p_obj$get_label(fcs, attr(df, "canal_y"))
      
      # Densité par binning raster, rendue en heatmap plotly (même palette que le reste
      # de l'app), plutôt qu'un nuage de points uni.
      lim_x <- range(df$X, na.rm = TRUE)
      lim_y <- range(df$Y, na.rm = TRUE)
      dens  <- calculer_densite_matrice_plotly(df$X, df$Y, lim_x, lim_y)
      
      plt <- plot_ly(source = ns("plot_gate_doublet_dessin"))
      if (!is.null(dens)) {
        plt <- plt %>%
          add_trace(x = dens$x, y = dens$y, z = dens$z, type = "heatmap",
                    colorscale = COLORSCALE_DENSITE_PLOTLY, showscale = FALSE,
                    hoverinfo = "none")
      } else {
        plt <- plt %>%
          add_trace(data = df, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    name = "Cellules",
                    marker = list(size = 2, color = "#0077b6", opacity = 0.4),
                    hoverinfo = "none")
      }
      plt %>%
        add_trace(x = numeric(0), y = numeric(0), type = "scatter",
                  mode = "lines+markers", name = "Gate en cours",
                  line   = list(color = "#e65100", width = 2, dash = "dot"),
                  marker = list(size = 5, color = "#e65100", symbol = "circle"),
                  hoverinfo = "x+y") %>%
        event_register("plotly_click") %>%
        event_register("plotly_doubleclick") %>%
        event_register("plotly_relayout") %>%
        layout(
          dragmode = "zoom",
          xaxis = list(title = lbl_x),
          yaxis = list(title = lbl_y),
          legend = list(orientation = "h", y = -0.15),
          margin = list(b = 60),
          shapes = list()
        ) %>%
        config(displayModeBar = TRUE, editable = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE)
    })
    
    observeEvent(sommets_doublet_rv(), {
      soms <- sommets_doublet_rv()
      proxy <- plotlyProxy(ns("plot_gate_doublet_dessin"), session)
      
      if (length(soms) == 0) {
        plotlyProxyInvoke(proxy, "restyle", list(x = list(list()), y = list(list())), list(1))
        plotlyProxyInvoke(proxy, "relayout", list(shapes = list()))
        return(invisible(NULL))
      }
      
      xs <- c(vapply(soms, `[[`, numeric(1), "x"), soms[[1]]$x)
      ys <- c(vapply(soms, `[[`, numeric(1), "y"), soms[[1]]$y)
      
      plotlyProxyInvoke(proxy, "restyle", list(x = list(xs), y = list(ys)), list(1))
      
      if (length(soms) >= 2) {
        df <- isolate(donnees_plot_doublet())
        rx <- diff(range(df$X, na.rm = TRUE)); if (!is.finite(rx) || rx == 0) rx <- 1
        ry <- diff(range(df$Y, na.rm = TRUE)); if (!is.finite(ry) || ry == 0) ry <- 1
        
        shapes <- lapply(seq_along(soms), function(i) {
          list(type = "circle",
               x0 = soms[[i]]$x - 0.012 * rx, x1 = soms[[i]]$x + 0.012 * rx,
               y0 = soms[[i]]$y - 0.012 * ry, y1 = soms[[i]]$y + 0.012 * ry,
               xref = "x", yref = "y",
               fillcolor = "rgba(230,101,0,0.3)",
               line = list(color = "#e65100"))
        })
        plotlyProxyInvoke(proxy, "relayout", list(shapes = shapes))
      } else {
        plotlyProxyInvoke(proxy, "relayout", list(shapes = list()))
      }
    }, ignoreNULL = FALSE)
    
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_doublet_dessin")), {
      ev   <- event_data("plotly_relayout", source = ns("plot_gate_doublet_dessin"))
      req(ev)
      soms <- sommets_doublet_rv()
      
      modifie <- FALSE
      for (nm in names(ev)) {
        m_idx <- regmatches(nm, regexpr("[0-9]+", nm))
        if (length(m_idx) == 0) next
        i <- as.integer(m_idx) + 1L
        if (i < 1 || i > length(soms)) next
        
        if (grepl("\\.x0$", nm)) { soms[[i]]$x <- ev[[nm]] + 0; modifie <- TRUE }
        if (grepl("\\.y0$", nm)) { soms[[i]]$y <- ev[[nm]] + 0; modifie <- TRUE }
      }
      if (modifie) sommets_doublet_rv(soms)
    }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # APPLICATION — MÉTHODE AUTOMATIQUE (STATISTIQUE)
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(input$btn_appliquer_doublet_auto, {
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      fs    <- input$facteur_sensibilite_doublet %||% 4
      axe   <- input$axe_doublet %||% "H_A"
      cible <- if (isTRUE(input$doublet_tous_echantillons_auto)) NULL else input$sel_ech_cible_auto
      
      withProgress(message = paste0("Retrait automatique des doublets ", etape, "..."), value = 0.4, {
        tryCatch({
          if (etape == "FSC") {
            p$retirer_doublets_FSC(facteur_sensibilite = fs, axe_discrimination = axe,
                                   nom_echantillon = cible)
          } else {
            p$retirer_doublets_SSC(facteur_sensibilite = fs, axe_discrimination = axe,
                                   nom_echantillon = cible)
          }
          pipeline(p)
          doublet_trigger(doublet_trigger() + 1L)
          showNotification(paste0("Doublets ", etape, " retirés (méthode statistique, facteur = ", fs, ")."),
                           type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # APPLICATION — MÉTHODE PAR GATING MANUEL
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(input$btn_save_gate_doublet, {
      soms  <- sommets_doublet_rv()
      etape <- input$etape_doublet %||% "FSC"
      axe   <- input$axe_doublet %||% "H_A"
      
      if (length(soms) < 3) {
        showNotification("Tracez au moins 3 sommets pour définir un gate.", type = "error")
        return(invisible(NULL))
      }
      
      p <- carrot_obj()
      
      xs <- sapply(soms, `[[`, "x")
      ys <- sapply(soms, `[[`, "y")
      idx_hull <- grDevices::chull(xs, ys)
      mat_poly <- cbind(xs[idx_hull], ys[idx_hull])
      
      cible <- if (isTRUE(input$doublet_tous_echantillons_gate)) NULL else input$sel_ech_cible_gate
      
      withProgress(message = paste0("Application du gate de doublets ", etape, "..."), value = 0.4, {
        tryCatch({
          if (etape == "FSC") {
            p$gate_les_doublets_FSC(points_utilisateur = mat_poly, axe_discrimination = axe,
                                    nom_echantillon = cible)
          } else {
            p$gate_les_doublets_SSC(points_utilisateur = mat_poly, axe_discrimination = axe,
                                    nom_echantillon = cible)
          }
          pipeline(p)
          doublet_trigger(doublet_trigger() + 1L)
          sommets_doublet_rv(list())
          showNotification(paste0("Gate de doublets ", etape, " enregistré."), type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # PANNEAU RÉSULTAT (DOUBLETS)
    # ════════════════════════════════════════════════════════════════════════
    
    output$titre_resultat_doublet <- renderUI({
      etape <- input$etape_doublet %||% "FSC"
      paste0("Résultat — doublets ", etape)
    })
    
    output$plot_doublet_resultat <- renderPlot({
      doublet_trigger()
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      req(input$sel_ech_doublet)
      
      res <- tryCatch(
        p$visualiser_doublets(nom_echantillon = input$sel_ech_doublet,
                              type_analyse     = etape,
                              max_points       = input$doublet_max_points %||% 10000),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
      )
      
      if (is.null(res)) {
        plot.new()
        text(0.5, 0.5, paste0("Aucun retrait de doublets ", etape, " appliqué pour cet échantillon."),
             cex = 1.05, col = "grey40")
        return()
      }
      print(res)
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # RÉSOLUTION DE LA SOURCE ET DÉTECTION DU CANAL (VIABILITÉ)
    # ════════════════════════════════════════════════════════════════════════
    
    # Reproduit exactement la chaîne de priorité interne de transformation_arcsinh :
    # doublets (final) > débris > échantillons compensés bruts.
    obtenir_source_viabilite_pretransfo <- function(p) {
      if (!is.null(p$post_doublets_final) && length(p$post_doublets_final) > 0) return(p$post_doublets_final)
      if (!is.null(p$post_debris) && length(p$post_debris) > 0) return(p$post_debris)
      return(p$echantillons_traites)
    }
    
    # Détecte les canaux susceptibles de correspondre à un marqueur de viabilité
    # (colorants Live/Dead usuels), sans présumer d'une nomenclature imposée.
    detecter_canaux_viabilite <- function(cols) {
      motif <- "viab|zombie|live.?dead|l/?d|7-?aad|fixable|ghost.?dye|amine.?react"
      grep(motif, cols, value = TRUE, ignore.case = TRUE)
    }
    
    resoudre_canal_fsc_viabilite <- function(cols) {
      cx <- grep("FSC-A", cols, value = TRUE, ignore.case = TRUE)[1]
      if (is.na(cx)) cx <- grep("FSC", cols, value = TRUE, ignore.case = TRUE)[1]
      cx
    }
    
    output$ui_canal_viabilite <- renderUI({
      p   <- carrot_obj()
      src <- obtenir_source_viabilite_pretransfo(p)
      if (length(src) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucune donnée disponible en amont."))
      }
      cols <- flowCore::colnames(src[[1]])
      candidats <- detecter_canaux_viabilite(cols)
      afficher_tout <- isTRUE(input$viabilite_afficher_tous_canaux)
      
      if (length(candidats) == 0 && !afficher_tout) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"),
                   " Aucun canal de viabilité détecté automatiquement dans ces données. ",
                   "Cette étape n'est disponible que si un tel canal est présent — ",
                   "cochez « Afficher tous les canaux » pour vérifier et sélectionner manuellement."))
      }
      
      choix <- if (afficher_tout || length(candidats) == 0) cols else candidats
      selectInput(ns("canal_viabilite"), "Canal :", choices = choix, selected = choix[1])
    })
    
    # Réinitialise l'état de transformation si l'utilisateur change de canal
    observeEvent(input$canal_viabilite, {
      canal_transforme_rv(NULL)
      sommets_viabilite_rv(list())
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    output$ui_select_echantillon_viabilite_transfo_cible <- renderUI({
      p   <- carrot_obj()
      src <- obtenir_source_viabilite_pretransfo(p)
      req(length(src) > 0)
      selectInput(ns("sel_ech_viabilite_transfo_cible"), "Échantillon ciblé :", choices = names(src))
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # ÉTAPE 1 — TRANSFORMATION ARCSINH DU CANAL DE VIABILITÉ
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(input$btn_appliquer_transfo_viabilite, {
      p     <- carrot_obj()
      canal <- input$canal_viabilite
      req(canal)
      cof   <- input$cofacteur_viabilite %||% 150
      cible <- if (isTRUE(input$viabilite_tous_echantillons_transfo)) NULL else input$sel_ech_viabilite_transfo_cible
      
      withProgress(message = "Transformation Arcsinh en cours...", value = 0.4, {
        tryCatch({
          p$transformation_arcsinh(canaux = canal, echantillon = cible, cofactor = cof)
          pipeline(p)
          canal_transforme_rv(canal)
          viabilite_trigger(viabilite_trigger() + 1L)
          showNotification(
            paste0("Transformation Arcsinh appliquée sur '", canal, "' (cofacteur = ", cof, ")."),
            type = "message"
          )
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    output$ui_transfo_viabilite_status <- renderUI({
      viabilite_trigger()
      canal   <- input$canal_viabilite
      req(canal)
      transforme <- identical(canal_transforme_rv(), canal)
      if (transforme) {
        div(class = "source-info-box",
            icon("check-circle"),
            " Transformation Arcsinh appliquée sur ce canal — le gating est disponible ci-dessous.")
      } else {
        div(class = "canal-info-box",
            icon("info-circle"),
            " Appliquez la transformation Arcsinh sur ce canal pour activer le gating.")
      }
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # SÉLECTEURS DYNAMIQUES (VIABILITÉ) — actifs seulement après transformation
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_select_echantillon_viabilite <- renderUI({
      viabilite_trigger()
      p     <- carrot_obj()
      canal <- input$canal_viabilite
      req(canal, identical(canal_transforme_rv(), canal))
      req(!is.null(p$post_transformation), length(p$post_transformation) > 0)
      selectInput(ns("sel_ech_viabilite"), "Échantillon (référence pour le tracé) :",
                  choices = names(p$post_transformation))
    })
    
    output$ui_select_echantillon_cible_gate_viabilite <- renderUI({
      p     <- carrot_obj()
      canal <- input$canal_viabilite
      req(canal, identical(canal_transforme_rv(), canal))
      req(!is.null(p$post_transformation), length(p$post_transformation) > 0)
      selectInput(ns("sel_ech_cible_gate_viabilite"), "Échantillon ciblé :",
                  choices = names(p$post_transformation))
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # DONNÉES POUR LE GRAPHIQUE INTERACTIF (VIABILITÉ)
    # ════════════════════════════════════════════════════════════════════════
    
    donnees_plot_viabilite <- reactive({
      viabilite_trigger()
      p     <- carrot_obj()
      canal <- input$canal_viabilite
      req(canal, identical(canal_transforme_rv(), canal))
      req(input$sel_ech_viabilite)
      fcs <- p$post_transformation[[input$sel_ech_viabilite]]
      req(!is.null(fcs))
      cols <- flowCore::colnames(fcs)
      cx   <- resoudre_canal_fsc_viabilite(cols)
      req(!is.na(cx), canal %in% cols)
      mat <- flowCore::exprs(fcs)[, c(cx, canal), drop = FALSE]
      n    <- nrow(mat)
      nmax <- input$viabilite_max_points %||% 10000
      if (!is.null(p$seed)) set.seed(p$seed)
      idx <- if (n > nmax) sample(seq_len(n), nmax) else seq_len(n)
      df  <- as.data.frame(mat[idx, , drop = FALSE])
      colnames(df) <- c("X", "Y")
      attr(df, "canal_x") <- cx
      attr(df, "canal_y") <- canal
      df
    })
    
    observeEvent(list(input$sel_ech_viabilite, input$canal_viabilite), {
      sommets_viabilite_rv(list())
    }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # GESTION DU POLYGONE — clic par clic + double-clic pour fermer (VIABILITÉ)
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(event_data("plotly_click", source = ns("plot_gate_viabilite_dessin")), {
      ev <- event_data("plotly_click", source = ns("plot_gate_viabilite_dessin"))
      req(ev)
      soms <- sommets_viabilite_rv()
      soms[[length(soms) + 1]] <- list(x = ev$x, y = ev$y)
      sommets_viabilite_rv(soms)
    })
    
    observeEvent(event_data("plotly_doubleclick", source = ns("plot_gate_viabilite_dessin")), {
      soms <- sommets_viabilite_rv()
      if (length(soms) >= 3) {
        showNotification(
          paste0("Polygone fermé (", length(soms), " sommets). Vérifiez et enregistrez le gate."),
          type = "message"
        )
      } else {
        showNotification("Tracez au moins 3 sommets avant de fermer le polygone.", type = "warning")
      }
    })
    
    observeEvent(input$btn_undo_sommet_viabilite, {
      soms <- sommets_viabilite_rv()
      if (length(soms) > 0) sommets_viabilite_rv(soms[-length(soms)])
    })
    
    observeEvent(input$btn_reset_gate_viabilite, {
      sommets_viabilite_rv(list())
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # GRAPHIQUE INTERACTIF — rendu de base + plotlyProxy (VIABILITÉ)
    # ════════════════════════════════════════════════════════════════════════
    
    output$plot_gate_viabilite_dessin <- renderPlotly({
      df <- donnees_plot_viabilite()
      req(nrow(df) > 0)
      
      p_obj <- carrot_obj()
      fcs   <- p_obj$post_transformation[[input$sel_ech_viabilite]]
      lbl_x <- p_obj$get_label(fcs, attr(df, "canal_x"))
      lbl_y <- p_obj$get_label(fcs, attr(df, "canal_y"))
      
      # Densité par binning raster, rendue en heatmap plotly (même palette que le reste
      # de l'app), plutôt qu'un nuage de points uni.
      lim_x <- range(df$X, na.rm = TRUE)
      lim_y <- range(df$Y, na.rm = TRUE)
      dens  <- calculer_densite_matrice_plotly(df$X, df$Y, lim_x, lim_y)
      
      plt <- plot_ly(source = ns("plot_gate_viabilite_dessin"))
      if (!is.null(dens)) {
        plt <- plt %>%
          add_trace(x = dens$x, y = dens$y, z = dens$z, type = "heatmap",
                    colorscale = COLORSCALE_DENSITE_PLOTLY, showscale = FALSE,
                    hoverinfo = "none")
      } else {
        plt <- plt %>%
          add_trace(data = df, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    name = "Cellules",
                    marker = list(size = 2, color = "#0077b6", opacity = 0.4),
                    hoverinfo = "none")
      }
      plt %>%
        add_trace(x = numeric(0), y = numeric(0), type = "scatter",
                  mode = "lines+markers", name = "Gate en cours",
                  line   = list(color = "#e65100", width = 2, dash = "dot"),
                  marker = list(size = 5, color = "#e65100", symbol = "circle"),
                  hoverinfo = "x+y") %>%
        event_register("plotly_click") %>%
        event_register("plotly_doubleclick") %>%
        event_register("plotly_relayout") %>%
        layout(
          dragmode = "zoom",
          xaxis = list(title = lbl_x),
          yaxis = list(title = lbl_y),
          legend = list(orientation = "h", y = -0.15),
          margin = list(b = 60),
          shapes = list()
        ) %>%
        config(displayModeBar = TRUE, editable = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE)
    })
    
    observeEvent(sommets_viabilite_rv(), {
      soms <- sommets_viabilite_rv()
      proxy <- plotlyProxy(ns("plot_gate_viabilite_dessin"), session)
      
      if (length(soms) == 0) {
        plotlyProxyInvoke(proxy, "restyle", list(x = list(list()), y = list(list())), list(1))
        plotlyProxyInvoke(proxy, "relayout", list(shapes = list()))
        return(invisible(NULL))
      }
      
      xs <- c(vapply(soms, `[[`, numeric(1), "x"), soms[[1]]$x)
      ys <- c(vapply(soms, `[[`, numeric(1), "y"), soms[[1]]$y)
      
      plotlyProxyInvoke(proxy, "restyle", list(x = list(xs), y = list(ys)), list(1))
      
      if (length(soms) >= 2) {
        df <- isolate(donnees_plot_viabilite())
        rx <- diff(range(df$X, na.rm = TRUE)); if (!is.finite(rx) || rx == 0) rx <- 1
        ry <- diff(range(df$Y, na.rm = TRUE)); if (!is.finite(ry) || ry == 0) ry <- 1
        
        shapes <- lapply(seq_along(soms), function(i) {
          list(type = "circle",
               x0 = soms[[i]]$x - 0.012 * rx, x1 = soms[[i]]$x + 0.012 * rx,
               y0 = soms[[i]]$y - 0.012 * ry, y1 = soms[[i]]$y + 0.012 * ry,
               xref = "x", yref = "y",
               fillcolor = "rgba(230,101,0,0.3)",
               line = list(color = "#e65100"))
        })
        plotlyProxyInvoke(proxy, "relayout", list(shapes = shapes))
      } else {
        plotlyProxyInvoke(proxy, "relayout", list(shapes = list()))
      }
    }, ignoreNULL = FALSE)
    
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_viabilite_dessin")), {
      ev   <- event_data("plotly_relayout", source = ns("plot_gate_viabilite_dessin"))
      req(ev)
      soms <- sommets_viabilite_rv()
      
      modifie <- FALSE
      for (nm in names(ev)) {
        m_idx <- regmatches(nm, regexpr("[0-9]+", nm))
        if (length(m_idx) == 0) next
        i <- as.integer(m_idx) + 1L
        if (i < 1 || i > length(soms)) next
        
        if (grepl("\\.x0$", nm)) { soms[[i]]$x <- ev[[nm]] + 0; modifie <- TRUE }
        if (grepl("\\.y0$", nm)) { soms[[i]]$y <- ev[[nm]] + 0; modifie <- TRUE }
      }
      if (modifie) sommets_viabilite_rv(soms)
    }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # ÉTAPE 2 — ENREGISTREMENT DU GATE DE VIABILITÉ
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(input$btn_save_gate_viabilite, {
      soms  <- sommets_viabilite_rv()
      canal <- input$canal_viabilite
      
      if (!identical(canal_transforme_rv(), canal)) {
        showNotification("Appliquez d'abord la transformation Arcsinh sur ce canal.", type = "error")
        return(invisible(NULL))
      }
      if (length(soms) < 3) {
        showNotification("Tracez au moins 3 sommets pour définir un gate.", type = "error")
        return(invisible(NULL))
      }
      
      p   <- carrot_obj()
      fcs <- p$post_transformation[[input$sel_ech_viabilite]]
      req(!is.null(fcs))
      canal_fsc <- resoudre_canal_fsc_viabilite(flowCore::colnames(fcs))
      req(!is.na(canal_fsc))
      
      xs <- sapply(soms, `[[`, "x")
      ys <- sapply(soms, `[[`, "y")
      idx_hull <- grDevices::chull(xs, ys)
      mat_poly <- cbind(xs[idx_hull], ys[idx_hull])
      
      cible <- if (isTRUE(input$viabilite_tous_echantillons_gate)) NULL else input$sel_ech_cible_gate_viabilite
      
      withProgress(message = "Application du gate de viabilité...", value = 0.4, {
        tryCatch({
          p$retirer_les_cellules_mortes(
            canal_fsc          = canal_fsc,
            marqueur_viabilite  = canal,
            points_utilisateur  = mat_poly,
            nom_echantillon     = cible
          )
          pipeline(p)
          viabilite_trigger(viabilite_trigger() + 1L)
          sommets_viabilite_rv(list())
          showNotification("Gate de cellules vivantes enregistré.", type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # PANNEAU RÉSULTAT (VIABILITÉ)
    # ════════════════════════════════════════════════════════════════════════
    
    output$plot_viabilite_resultat <- renderPlot({
      viabilite_trigger()
      p <- carrot_obj()
      req(input$sel_ech_viabilite)
      
      res <- tryCatch(
        p$visualiser_viabilite(nom_echantillon = input$sel_ech_viabilite,
                               max_points       = input$viabilite_max_points %||% 10000),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
      )
      
      if (is.null(res)) {
        plot.new()
        text(0.5, 0.5, "Enregistrez un gate de viabilité pour afficher le résultat.",
             cex = 1.05, col = "grey40")
        return()
      }
      print(res)
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # RÉSUMÉ CUMULATIF — construction commune (DÉBRIS + DOUBLETS)
    # ════════════════════════════════════════════════════════════════════════
    
    construire_recap <- function(p, nom) {
      etapes <- list()
      
      if (!is.null(p$echantillons_traites[[nom]])) {
        n <- nrow(flowCore::exprs(p$echantillons_traites[[nom]]))
        etapes[["Compensées"]] <- list(n = n, ref = n)
      }
      ref <- etapes[["Compensées"]]$n %||% NA
      
      if (!is.null(p$post_PeacoQC[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_PeacoQC[[nom]]))
        etapes[["PeacoQC"]] <- list(n = n, ref = ref)
      }
      if (!is.null(p$post_flowAI[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_flowAI[[nom]]))
        etapes[["flowAI"]] <- list(n = n, ref = ref)
      }
      if (!is.null(p$post_retrait_bordures[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_retrait_bordures[[nom]]))
        etapes[["Bordures"]] <- list(n = n, ref = ref)
      }
      
      # Gates de débris nommés (empilés)
      ref_gate <- if (length(etapes) > 0) etapes[[length(etapes)]]$n else ref
      for (nm_gate in names(p$gates_history)) {
        infos <- p$gates_history[[nm_gate]][[nom]]
        if (!is.null(infos)) {
          etapes[[nm_gate]] <- list(n = infos$n_apres, ref = infos$n_avant)
          ref_gate <- infos$n_apres
        }
      }
      
      # Doublets FSC puis SSC, chaînés à la suite des débris
      if (!is.null(p$post_doublets_FSC[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_doublets_FSC[[nom]]))
        etapes[["Doublets FSC"]] <- list(n = n, ref = ref_gate)
        ref_gate <- n
      }
      if (!is.null(p$post_doublets_SSC[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_doublets_SSC[[nom]]))
        etapes[["Doublets SSC"]] <- list(n = n, ref = ref_gate)
        ref_gate <- n
      }
      
      # Viabilité, chaînée à la suite des doublets/débris
      if (!is.null(p$post_viabilite[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_viabilite[[nom]]))
        etapes[["Viabilité"]] <- list(n = n, ref = ref_gate)
        ref_gate <- n
      }
      
      etapes
    }
    
    rendre_recap <- function(etapes, nom) {
      if (length(etapes) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucune étape réalisée."))
      }
      
      lignes <- lapply(names(etapes), function(step) {
        info <- etapes[[step]]
        n    <- info$n
        ref_n <- info$ref
        pct  <- if (!is.null(ref_n) && !is.na(ref_n) && ref_n > 0)
          round(n / ref_n * 100, 1) else NA
        
        couleur <- if (is.na(pct)) "#666"
        else if (pct >= 90) "#2e7d32"
        else if (pct >= 70) "#e65100"
        else "#c62828"
        
        div(class = "recap-row",
            span(class = "label", step), " : ",
            span(class = "value",
                 style = paste0("color:", couleur, "; font-weight:bold;"),
                 format(n, big.mark = " "),
                 if (!is.na(pct)) paste0(" (", pct, "%)") else ""
            )
        )
      })
      
      tagList(
        div(style = "font-size:11px; color:#888; margin-bottom:8px;",
            icon("user"), " ", nom),
        do.call(tagList, lignes)
      )
    }
    
    output$ui_recap_pipeline <- renderUI({
      bordures_trigger(); gates_trigger(); doublet_trigger(); viabilite_trigger()
      p <- carrot_obj()
      req(input$sel_ech_debris)
      rendre_recap(construire_recap(p, input$sel_ech_debris), input$sel_ech_debris)
    })
    
    output$ui_recap_pipeline_doublet <- renderUI({
      bordures_trigger(); gates_trigger(); doublet_trigger(); viabilite_trigger()
      p <- carrot_obj()
      req(input$sel_ech_doublet)
      rendre_recap(construire_recap(p, input$sel_ech_doublet), input$sel_ech_doublet)
    })
    
    output$ui_recap_pipeline_viabilite <- renderUI({
      bordures_trigger(); gates_trigger(); doublet_trigger(); viabilite_trigger()
      p <- carrot_obj()
      req(input$sel_ech_viabilite)
      rendre_recap(construire_recap(p, input$sel_ech_viabilite), input$sel_ech_viabilite)
    })
    
  })
}

# Opérateur null-coalescing utilitaire
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}