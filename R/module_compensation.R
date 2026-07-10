library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(ggplot2)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

compensation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    
    # jQuery UI (module "sortable") : requis par le plugin selectize "drag_drop" pour pouvoir
    # glisser-déposer les canaux et réordonner la matrice de spillover. Shiny ne charge pas
    # cette dépendance par défaut, il faut l'ajouter explicitement.
    tags$head(
      tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/jqueryui/1.13.2/jquery-ui.min.js")
    ),
    
    # CSS pour le curseur drag sur les lignes de gate
    tags$style(HTML("
      .gate-instructions {
        background:#eef6fb; border-left:3px solid #0077b6;
        padding:8px 12px; border-radius:4px;
        font-size:12px; color:#444; margin-top:8px;
      }
      .gate-instructions b { color:#0077b6; }
    ")),
    
    tabBox(
      title = tagList(icon("sliders-h"), "Outils de Compensation"),
      id = ns("comp_steps"), width = 12,
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 1 — TRANSFORMATION (en premier)
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("wave-square"), " Transformation"),
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
                   actionButton(ns("btn_apply_trans"), tagList(icon("play"), " Appliquer la transformation"),
                                class = "btn-primary", style = "width:100%; font-weight:bold;"),
                   br(), br(),
                   uiOutput(ns("ui_trans_status"))
                 )
          ),
          column(width = 8,
                 box(title = "Visualisation", width = NULL, status = "info", solidHeader = TRUE,
                     div(style = "display:flex; justify-content:center;",
                         plotOutput(ns("plot_transformation"), width = "400px", height = "400px")
                     )
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 2 — DÉFINITION DES GATES (en second)
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("crosshairs"), " Définition des Gates"),
        value = "tab_gates",
        
        fluidRow(
          
          # ───────────────────────────────────────────────
          # COLONNE GAUCHE — Contrôles
          # ───────────────────────────────────────────────
          column(
            width = 3,
            wellPanel(
              style = "background:#f8f9fa; border:1px solid #dee2e6;",
              
              h4(tagList(icon("vial"), " Tube à gater"),
                 style = "margin-top:0; color:#605ca8;"),
              
              uiOutput(ns("ui_gate_canal_select")),
              uiOutput(ns("ui_gate_unstained_info")),
              uiOutput(ns("ui_unstained_option")),
              
              hr(style = "border-color:#dee2e6;"),
              
              # Négatif
              h5(tagList(icon("minus-circle"), " Population Négative"),
                 style = "color:#0077b6; margin-bottom:4px;"),
              fluidRow(
                column(6, numericInput(ns("gate_neg_min"), "Min", value = 0, step = 0.05)),
                column(6, numericInput(ns("gate_neg_max"), "Max", value = 0.5, step = 0.05))
              ),
              
              # Positif
              h5(tagList(icon("plus-circle"), " Population Positive"),
                 style = "color:#d90429; margin-bottom:4px;"),
              fluidRow(
                column(6, numericInput(ns("gate_pos_min"), "Min", value = 1, step = 0.05)),
                column(6, numericInput(ns("gate_pos_max"), "Max", value = 1.5, step = 0.05))
              ),
              
              hr(style = "border-color:#dee2e6;"),
              
              actionButton(
                ns("btn_apply_gate"),
                tagList(icon("check-circle"), " Valider ce gate"),
                class = "btn-success",
                style = "width:100%; font-weight:bold;"
              ),
              
              br(), br(),
              uiOutput(ns("ui_gate_status_canal"))
            )
          ),
          
          # ───────────────────────────────────────────────
          # COLONNE CENTRALE — Graphique interactif
          # ───────────────────────────────────────────────
          column(
            width = 6,
            box(
              title = tagList(icon("chart-area"), " Ajustement interactif des bornes"),
              width = NULL, status = "info", solidHeader = TRUE,
              
              uiOutput(ns("ui_gate_plot_placeholder")),
              
              # Graphique Plotly avec shapes déplaçables
              plotlyOutput(ns("plot_gate_interactive"), height = "420px"),
              
              div(
                class = "gate-instructions",
                icon("hand-pointer"),
                " Faites glisser les lignes verticales directement sur le graphique pour ajuster les bornes. ",
                tags$b("Bleu = négatif / Unstained"), " — ",
                tags$b("Rouge = monomarqué (positif)"), ".",
                br(),
                icon("info-circle"),
                " Les champs numériques à gauche se synchronisent automatiquement."
              )
            )
          ),
          
          # ───────────────────────────────────────────────
          # COLONNE DROITE — Récapitulatif
          # ───────────────────────────────────────────────
          column(
            width = 3,
            box(
              title = tagList(icon("clipboard-list"), " Récapitulatif"),
              width = NULL, status = "warning", solidHeader = TRUE,
              
              uiOutput(ns("ui_gates_recap")),
              hr(),
              
              div(style = "text-align:center;",
                  uiOutput(ns("ui_btn_calculer_spillover")))
            )
          )
        )
      ),
      
      # ───────────────────────────────────────────────
      # MATRICE DE SPILLOVER
      # ───────────────────────────────────────────────
      
      tabPanel(
        title = "Matrice de Spillover",
        fluidRow(
          # --- Colonne de Gauche : Flux logique ---
          column(width = 4,
                 wellPanel(
                   h4("Paramètres de Compensation"),
                   
                   # 1. Bouton de calcul en premier
                   uiOutput(ns("ui_btn_calculer_spillover")),
                   
                   hr(),
                   
                   # 2. Sélecteur d'échantillon ensuite
                   uiOutput(ns("ui_select_echantillon")),
                   
                   hr(),
                   
                   # 2 bis. Personnalisation de l'affichage (ordre des canaux / transposition)
                   h5(tagList(icon("table"), " Personnaliser l'affichage"),
                      style = "color:#605ca8;"),
                   h5(tagList(icon("table"), " Personnaliser l'affichage"), style = "color:#605ca8;"),
                   selectizeInput(
                     ns("select_ordre_canaux"),
                     "Glissez les canaux pour réordonner lignes/colonnes :",
                     choices  = NULL, # Sera mis à jour via le serveur
                     multiple = TRUE,
                     options  = list(plugins = list("drag_drop"), placeholder = "En attente du calcul...")
                   ),
                   fluidRow(
                     column(6, actionButton(ns("btn_inverser_ordre"), tagList(icon("exchange-alt"), " Inverser l'ordre"),
                                            class = "btn-default", style = "width:100%; font-size:11px;")),
                     column(6, actionButton(ns("btn_transposer_matrice"), tagList(icon("sync"), " Transposer L/C"),
                                            class = "btn-default", style = "width:100%; font-size:11px;"))
                   ),
                   
                   hr(),
                   
                   # 3. Boutons undo / redo
                   fluidRow(
                     column(6,
                            actionButton(ns("btn_undo_spillover"), icon("undo"), " Annuler",
                                         class = "btn-default", style = "width:100%;")
                     ),
                     column(6,
                            actionButton(ns("btn_redo_spillover"), icon("redo"), " Rétablir",
                                         class = "btn-default", style = "width:100%;")
                     )
                   ),
                   br(),
                   # 4. Bouton valider à la fin
                   actionButton(ns("btn_save_spillover"), 
                                "Valider et Enregistrer", 
                                class = "btn-success", 
                                icon = icon("save"),
                                style = "font-weight:bold; width:100%;")
                 )
          ),
          # --- Colonne de Droite : Édition ---
          column(width = 8,
                 box(title = "Matrice de Spillover éditée (%)", 
                     width = NULL, status = "info", solidHeader = TRUE,
                     DTOutput(ns("dt_spillover_matrix")),
                     hr(),
                     div(style = "font-size:12px; color:#666; margin-bottom:6px;",
                         icon("info-circle"),
                         " Téléchargez la matrice telle qu'affichée (ordre et orientation courants) :"),
                     fluidRow(
                       column(6, downloadButton(ns("dl_matrix_png"), "Télécharger en PNG",
                                                class = "btn-default", style = "width:100%;")),
                       column(6, downloadButton(ns("dl_matrix_pdf"), "Télécharger en PDF",
                                                class = "btn-default", style = "width:100%;"))
                     )
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # COMPARAISON DES MÉDIANES (AVANT / APRÈS COMPENSATION)
      # ══════════════════════════════════════════════════════════════════════
      
      tabPanel(
        title = tagList(icon("balance-scale"), " Comparaison des Médianes"),
        fluidRow(
          column(width = 4,
                 wellPanel(
                   h4("Comparaison Avant / Après Compensation"),
                   p("Pour chaque fluorochrome (canal principal) et chaque canal, cette analyse affiche l'écart entre la médiane de la population positive et celle de la population négative, avant et après application de la matrice de compensation.",
                     style = "font-size:12px; color:#555;"),
                   hr(),
                   uiOutput(ns("ui_select_echantillon_medianes")),
                   hr(),
                   uiOutput(ns("ui_btn_comparer_medianes")),
                   hr(),
                   div(class = "gate-instructions",
                       icon("info-circle"),
                       " ", tags$b("Critère de bonne compensation :"),
                       " la compensation est correcte quand la médiane de la population positive (compensée) devient égale à la médiane de la population négative, dans tous les canaux hors diagonale — c'est-à-dire quand l'écart affiché tend vers 0.",
                       br(),
                       "La couleur reflète l'ampleur de l'écart (vert = proche de 0, rouge = écart important). La diagonale (canal principal) n'est pas concernée par ce test, car elle correspond au signal réel du marqueur — elle est affichée en gris neutre."
                   )
                 )
          ),
          column(width = 8,
                 box(title = tagList(icon("exclamation-triangle"), " Écart de médianes AVANT compensation (Positif − Négatif)"),
                     width = NULL, status = "danger", solidHeader = TRUE,
                     DTOutput(ns("dt_medianes_avant"))
                 ),
                 box(title = tagList(icon("check-circle"), " Écart de médianes APRÈS compensation (Positif − Négatif)"),
                     width = NULL, status = "success", solidHeader = TRUE,
                     DTOutput(ns("dt_medianes_apres"))
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # BIPLOTS
      # ══════════════════════════════════════════════════════════════════════
      
      tabPanel(
        title = "Biplots de Compensation",
        fluidRow(
          # --- Colonne de Gauche : Contrôles ---
          column(width = 3,
                 wellPanel(
                   h4("Paramètres de Visualisation"),
                   actionButton(ns("btn_apply_compensation"), "Appliquer la compensation",
                                class = "btn-warning", icon = icon("gears"), style = "width:100%"),
                   hr(),
                   uiOutput(ns("ui_select_echantillon_plot")),
                   checkboxInput(ns("mode_ensemble"), "Vue d'ensemble (toutes combinaisons)", FALSE),
                   conditionalPanel(
                     condition = paste0("!input['", ns("mode_ensemble"), "']"),
                     selectInput(ns("canal_x"), "Canal X :", choices = NULL),
                     selectInput(ns("canal_y"), "Canal Y :", choices = NULL)
                   ),
                   radioButtons(ns("affichage_type"), "Type d'affichage :",
                                choices = c("Both", "Before compensation only", "After compensation only"))
                 )
          ),
          # --- Colonne de Droite : Visualisation ---
          column(width = 9,
                 box(title = "Contrôle Qualité Compensation", width = NULL, status = "primary",
                     div(
                       style = "max-height: 80vh; overflow-y: auto; overflow-x: hidden;",
                       uiOutput(ns("ui_biplots_container"))
                     )
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # EXPORT
      # ══════════════════════════════════════════════════════════════════════
      
      tabPanel(
        title = tagList(icon("download"), " Export"),
        fluidRow(
          
          # --- Colonne Gauche : FCS compensés ---
          column(width = 6,
                 box(title = tagList(icon("file-medical"), " Fichiers FCS compensés"),
                     width = NULL, status = "primary", solidHeader = TRUE,
                     
                     p("Sélectionnez les échantillons à télécharger. Les fichiers seront exportés au format FCS avec la matrice de spillover intégrée dans les métadonnées.",
                       style = "font-size:12px; color:#555;"),
                     hr(),
                     
                     uiOutput(ns("ui_select_export_fcs")),
                     
                     fluidRow(
                       column(6,
                              actionButton(ns("btn_select_all_fcs"), "Tout sélectionner",
                                           class = "btn-default btn-sm", style = "width:100%;")
                       ),
                       column(6,
                              actionButton(ns("btn_deselect_all_fcs"), "Tout désélectionner",
                                           class = "btn-default btn-sm", style = "width:100%;")
                       )
                     ),
                     br(),
                     
                     uiOutput(ns("ui_download_fcs"))
                 )
          ),
          
          # --- Colonne Droite : Session RDS ---
          column(width = 6,
                 box(title = tagList(icon("box-archive"), " Session de compensation (RDS)"),
                     width = NULL, status = "success", solidHeader = TRUE,
                     
                     p("Exportez l'intégralité des paramètres ayant servi à la compensation dans un fichier RDS : transformation Arcsinh, matrice de spillover, gates, et figures de contrôle qualité.",
                       style = "font-size:12px; color:#555;"),
                     hr(),
                     
                     textInput(ns("rds_filename"), "Nom du fichier :",
                               value = "Compensation_Session_Complete.rds",
                               placeholder = "nom_session.rds"),
                     br(),
                     
                     div(style = "text-align:center;",
                         uiOutput(ns("ui_download_rds"))
                     )
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

compensation_server <- function(id, pipeline, pipeline_version) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ── Accès à l'objet CARROT ────────────────────────────────────────────────
    carrot_obj <- reactive({
      pipeline_version()
      pipeline()
    })
    
    # ── Signaux locaux ────────────────────────────────────────────────────────
    trans_done <- reactiveVal(0L)   # incrémenté après transformer_fcs()
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 1 — TRANSFORMATION
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_trans_file <- renderUI({
      p <- carrot_obj()
      req(!is.null(p$tubes_monomarques), length(p$tubes_monomarques) > 0)
      selectInput(ns("trans_file_sel"), "Tube à visualiser :",
                  choices = names(p$tubes_monomarques))
    })
    
    output$ui_trans_cx <- renderUI({
      p <- carrot_obj()
      req(!is.null(p$tubes_monomarques), length(p$tubes_monomarques) > 0)
      choices <- get_canaux_filtres(p)
      selectInput(ns("trans_cx"), "Axe X :", choices = choices)
    })
    
    output$ui_trans_cy <- renderUI({
      p <- carrot_obj()
      req(!is.null(p$tubes_monomarques), length(p$tubes_monomarques) > 0)
      choices <- get_canaux_filtres(p)
      selectInput(ns("trans_cy"), "Axe Y :", choices = choices,
                  selected = choices[min(2, length(choices))])
    })
    
    observeEvent(input$btn_apply_trans, {
      req(input$cofacteur)
      p <- carrot_obj()
      withProgress(message = "Transformation en cours...", value = 0.5, {
        p$transformer_fcs(cofacteur = input$cofacteur)
        trans_done(trans_done() + 1L)
      })
    })
    
    output$ui_trans_status <- renderUI({
      trans_done()
      p <- pipeline()
      if (is.null(p$monomarques_trans)) {
        tags$span(icon("times-circle"), " Données brutes (non transformées)",
                  style = "color:red;")
      } else {
        tags$span(icon("check-circle"), " Transformation Arcsinh appliquée",
                  style = "color:green;")
      }
    })
    
    output$plot_transformation <- renderPlot({
      trans_done()
      req(input$trans_file_sel, input$trans_cx, input$trans_cy)
      p <- pipeline()
      fcs_t <- p$monomarques_trans[[input$trans_file_sel]]
      validate(need(!is.null(fcs_t),
                    "Cliquez sur 'Appliquer la transformation'."))
      mat <- data.frame(
        x = flowCore::exprs(fcs_t)[, input$trans_cx],
        y = flowCore::exprs(fcs_t)[, input$trans_cy]
      )
      ggplot2::ggplot(mat, ggplot2::aes(x = x, y = y)) +
        ggpointdensity::geom_pointdensity(size = 0.2) +
        ggplot2::theme_bw()
    })
    
    # ========================================================================
    # SECTION 2 — GATING INTERACTIF
    # ========================================================================
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 2 — GATING INTERACTIF (VERSION SHAPES DÉPLAÇABLES)
    # ════════════════════════════════════════════════════════════════════════
    
    # ── Bornes réactives (source de vérité) ───────────────────────────────────
    rv_neg_min <- reactiveVal(0)
    rv_neg_max <- reactiveVal(2)
    rv_pos_min <- reactiveVal(4)
    rv_pos_max <- reactiveVal(8)
    
    # ── Gates validés ─────────────────────────────────────────────────────────
    gates_valides <- reactiveVal(list())
    
    # Flag anti-boucle pour updateNumericInput
    .prog_update <- reactiveVal(FALSE)
    
    # ── Helpers réactifs ──────────────────────────────────────────────────────
    
    trans_disponibles <- reactive({
      trans_done(); pipeline_version()
      p <- pipeline()
      !is.null(p$monomarques_trans) && length(p$monomarques_trans) > 0
    })
    
    canaux_monomarques <- reactive({
      pipeline_version()
      p <- pipeline()
      if (is.null(p$chemins_monomarques)) return(character(0))
      df <- p$chemins_monomarques
      cx <- df$canal[df$type == "Monomarque"]
      cx[!is.na(cx) & nchar(trimws(cx)) > 0]
    })
    
    has_unstained <- reactive({
      trans_done(); pipeline_version()
      "TUBE_UNSTAINED" %in% names(pipeline()$monomarques_trans)
    })
    
    utiliser_unstained_rv <- reactive({
      if (has_unstained()) isTRUE(input$use_unstained) else FALSE
    })
    
    # Synchronisation rv + numericInput
    set_bornes <- function(n1, n2, p1, p2) {
      .prog_update(TRUE)
      rv_neg_min(n1); rv_neg_max(n2)
      rv_pos_min(p1); rv_pos_max(p2)
      updateNumericInput(session, "gate_neg_min", value = n1)
      updateNumericInput(session, "gate_neg_max", value = n2)
      updateNumericInput(session, "gate_pos_min", value = p1)
      updateNumericInput(session, "gate_pos_max", value = p2)
      .prog_update(FALSE)
    }
    
    # ── UI : sélecteur de canal ───────────────────────────────────────────────
    output$ui_gate_canal_select <- renderUI({
      trans_done()
      canaux <- canaux_monomarques()
      if (length(canaux) == 0)
        return(div(class="alert alert-warning", style="font-size:12px; padding:7px;",
                   icon("exclamation-triangle"),
                   " Importez et transformez les données d'abord."))
      
      gv <- gates_valides()
      labels <- sapply(canaux, function(cx)
        if (cx %in% names(gv)) paste0("✔ ", cx) else cx)
      names(labels) <- canaux
      
      selectInput(ns("gate_canal_sel"),
                  label = tagList(icon("tag"), " Canal :"),
                  choices = setNames(canaux, labels))
    })
    
    output$ui_gate_unstained_info <- renderUI({
      trans_done()
      if (has_unstained())
        div(class="alert alert-info", style="padding:5px 9px; font-size:12px;",
            icon("check-circle"), " Tube Unstained disponible.")
      else
        div(class="alert alert-secondary", style="padding:5px 9px; font-size:12px;",
            icon("info-circle"), " Pas d'Unstained — référence négative interne.")
    })
    
    output$ui_unstained_option <- renderUI({
      req(has_unstained())
      checkboxInput(ns("use_unstained"),
                    label = tagList(icon("vial"),
                                    " Utiliser l'Unstained comme référence négative"),
                    value = TRUE)
    })
    
    output$ui_gate_plot_placeholder <- renderUI({
      if (!trans_disponibles())
        div(class="alert alert-warning", style="margin:6px 0; font-size:12px;",
            icon("exclamation-triangle"),
            " Appliquez d'abord la transformation.")
    })
    
    # ── Changement de canal → recharge bornes sauvegardées ───────────────────
    observeEvent(input$gate_canal_sel, {
      req(trans_disponibles())
      cx <- input$gate_canal_sel
      gv <- gates_valides()
      
      if (!is.null(cx) && cx %in% names(gv)) {
        s <- gv[[cx]]
        set_bornes(s$neg[1], s$neg[2], s$pos[1], s$pos[2])
        if (has_unstained())
          updateCheckboxInput(session, "use_unstained", value = s$use_unstained)
      } else {
        set_bornes(0, 2, 4, 8)
      }
    }, ignoreInit = TRUE)
    
    # ── Saisie clavier → rv ───────────────────────────────────────────────────
    observeEvent(input$gate_neg_min, {
      if (.prog_update()) return()
      rv_neg_min(input$gate_neg_min)
    }, ignoreInit = TRUE)
    
    observeEvent(input$gate_neg_max, {
      if (.prog_update()) return()
      rv_neg_max(input$gate_neg_max)
    }, ignoreInit = TRUE)
    
    observeEvent(input$gate_pos_min, {
      if (.prog_update()) return()
      rv_pos_min(input$gate_pos_min)
    }, ignoreInit = TRUE)
    
    observeEvent(input$gate_pos_max, {
      if (.prog_update()) return()
      rv_pos_max(input$gate_pos_max)
    }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # GRAPHIQUE PLOTLY — SHAPES DÉPLAÇABLES
    # ════════════════════════════════════════════════════════════════════════
    
    output$plot_gate_interactive <- renderPlotly({
      req(trans_disponibles(), input$gate_canal_sel)
      trans_done()
      
      n_min <- rv_neg_min()
      n_max <- rv_neg_max()
      p_min <- rv_pos_min()
      p_max <- rv_pos_max()
      
      p     <- pipeline()
      canal <- input$gate_canal_sel
      
      use_u    <- utiliser_unstained_rv()
      tube_neg <- if (use_u && has_unstained()) "TUBE_UNSTAINED" else canal
      
      fcs_pos  <- p$monomarques_trans[[canal]]
      fcs_neg  <- p$monomarques_trans[[tube_neg]]
      
      vals_pos <- flowCore::exprs(fcs_pos)[, canal]
      vals_neg <- flowCore::exprs(fcs_neg)[, canal]
      
      dens_pos <- density(vals_pos)
      dens_neg <- density(vals_neg)
      y_max    <- max(c(dens_pos$y, dens_neg$y)) * 1.15
      
      # SHAPES DÉPLAÇABLES — VERSION PRO
      shapes <- list(
        list(type="line", x0=n_min, x1=n_min, y0=0, y1=1, yref="paper",
             line=list(color="#0077b6", width=3)),
        list(type="line", x0=n_max, x1=n_max, y0=0, y1=1, yref="paper",
             line=list(color="#0077b6", width=3)),
        list(type="line", x0=p_min, x1=p_min, y0=0, y1=1, yref="paper",
             line=list(color="#d90429", width=3)),
        list(type="line", x0=p_max, x1=p_max, y0=0, y1=1, yref="paper",
             line=list(color="#d90429", width=3))
      )
      
      p <- plot_ly(source = ns("plot_gate_interactive")) %>%
        add_lines(x = dens_neg$x, y = dens_neg$y, name = "Négatif") %>%
        add_lines(x = dens_pos$x, y = dens_pos$y, name = "Positif") %>%
        event_register("plotly_relayout") %>%
        layout(
          shapes = shapes, # Vos shapes dynamiques
          dragmode = "pan"
        ) %>%
        config(editable = TRUE)
      
      # OBLIGATOIRE : enregistrer l’événement
      event_register(p, "plotly_relayout")
    })
    
    
    # ── Capture du déplacement des shapes ─────────────────────────────────────
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_interactive")), {
      
      ev <- event_data("plotly_relayout", source = ns("plot_gate_interactive"))
      
      print("Event reçu :")
      print(ev)
      
      req(ev)
      
      changed <- FALSE
      for (nm in names(ev)) {
        if (grepl("^shapes\\[([0-9]+)\\]\\.x0$", nm)) {
          idx <- as.integer(gsub("^shapes\\[([0-9]+)\\]\\.x0$", "\\1", nm))
          new_val <- ev[[nm]]
          
          # Mise à jour des valeurs réactives selon l'index
          switch(as.character(idx),
                 "0" = if(!isTRUE(all.equal(rv_neg_min(), new_val))) { rv_neg_min(new_val); changed <<- TRUE },
                 "1" = if(!isTRUE(all.equal(rv_neg_max(), new_val))) { rv_neg_max(new_val); changed <<- TRUE },
                 "2" = if(!isTRUE(all.equal(rv_pos_min(), new_val))) { rv_pos_min(new_val); changed <<- TRUE },
                 "3" = if(!isTRUE(all.equal(rv_pos_max(), new_val))) { rv_pos_max(new_val); changed <<- TRUE }
          )
        }
      }
      
      if (changed) {
        updateNumericInput(session, "gate_neg_min", value = rv_neg_min())
        updateNumericInput(session, "gate_neg_max", value = rv_neg_max())
        updateNumericInput(session, "gate_pos_min", value = rv_pos_min())
        updateNumericInput(session, "gate_pos_max", value = rv_pos_max())
      }
    })
    
    
    # ════════════════════════════════════════════════════════════════════════
    # VALIDATION DU GATE
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(input$btn_apply_gate, {
      req(trans_disponibles(), input$gate_canal_sel)
      
      p     <- pipeline()
      canal <- input$gate_canal_sel
      
      lim_n <- c(rv_neg_min(), rv_neg_max())
      lim_p <- c(rv_pos_min(), rv_pos_max())
      use_u <- utiliser_unstained_rv()
      
      if (lim_n[1] >= lim_n[2]) {
        showNotification("⚠ Borne négative : min doit être < max.", type="warning")
        return()
      }
      if (lim_p[1] >= lim_p[2]) {
        showNotification("⚠ Borne positive : min doit être < max.", type="warning")
        return()
      }
      
      withProgress(message = paste("Extraction gate :", canal), value = 0.5, {
        tryCatch({
          canaux_presents <- names(p$monomarques_trans)
          canaux_fluo     <- canaux_presents[canaux_presents != "TUBE_UNSTAINED"]
          if (is.null(p$canaux) || !canal %in% p$canaux) p$canaux <- canaux_fluo
          
          p$definir_et_extraire(
            nom_canal               = canal,
            intervalle_gate_negatif = lim_n,
            intervalle_gate_positif = lim_p,
            utiliser_unstained      = use_u
          )
          
          gv <- gates_valides()
          gv[[canal]] <- list(neg = lim_n, pos = lim_p, use_unstained = use_u)
          gates_valides(gv)
          
          showNotification(paste0("✔ Gate validé : ", canal),
                           type="message", duration=3)
        }, error = function(e) {
          showNotification(paste("Erreur :", conditionMessage(e)), type="error")
        })
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # RÉCAPITULATIF
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_gates_recap <- renderUI({
      gv     <- gates_valides()
      canaux <- canaux_monomarques()
      if (length(canaux) == 0)
        return(p("Aucun canal.", style="color:gray; font-size:12px;"))
      
      tagList(lapply(canaux, function(cx) {
        if (cx %in% names(gv)) {
          s   <- gv[[cx]]
          src <- if (s$use_unstained && has_unstained()) "Unstained" else "Interne"
          
          div(style="border:1px solid #c3e6cb; border-radius:5px;
                 padding:7px 9px; margin-bottom:5px; background:#f0fff4;",
              div(style="font-weight:bold; color:#155724; font-size:12px;",
                  icon("check-circle"), " ", cx),
              div(style="font-size:11px; color:#555; margin-top:2px;",
                  tags$span(style="color:#0077b6;",
                            icon("minus-circle"),
                            sprintf(" Nég [%.2f ; %.2f]", s$neg[1], s$neg[2])),
                  " | ", tags$em(src)),
              div(style="font-size:11px; color:#555;",
                  tags$span(style="color:#d90429;",
                            icon("plus-circle"),
                            sprintf(" Pos [%.2f ; %.2f]", s$pos[1], s$pos[2]))),
              div(style="margin-top:4px;",
                  actionButton(
                    ns(paste0("edit_gate_", make.names(cx))),
                    label = tagList(icon("edit"), " Modifier"),
                    class = "btn-xs btn-default",
                    style = "font-size:11px; padding:2px 6px;"
                  ))
          )
        } else {
          div(style="border:1px solid #f5c6cb; border-radius:5px;
                 padding:7px 9px; margin-bottom:5px; background:#fff5f5;",
              div(style="font-weight:bold; color:#721c24; font-size:12px;",
                  icon("times-circle"), " ", cx),
              div(style="font-size:11px; color:#888;", "Gate non défini")
          )
        }
      }))
    })
    
    # Boutons "Modifier"
    observe({
      canaux <- canaux_monomarques()
      for (cx in canaux) {
        local({
          cl <- cx
          observeEvent(input[[paste0("edit_gate_", make.names(cl))]], {
            updateSelectInput(session, "gate_canal_sel", selected = cl)
          }, ignoreInit = TRUE, ignoreNULL = TRUE)
        })
      }
    })
    
    # Statut
    output$ui_gate_status_canal <- renderUI({
      gv      <- gates_valides()
      canaux  <- canaux_monomarques()
      valides <- length(gv)
      total   <- length(canaux)
      
      div(style="font-size:12px; color:#555; text-align:center;",
          tags$b(valides), " / ", tags$b(total), " gate(s) validé(s)",
          if (total > 0 && valides == total)
            div(style="color:green; margin-top:4px;",
                icon("check-circle"), " Tous les gates sont définis !")
      )
    })
    
    
    # ════════════════════════════════════════════════════════════════════════
    # MATRICE DE SPILLOVER
    # ════════════════════════════════════════════════════════════════════════
    spillover_trigger <- reactiveVal(0)
    
    # Le bouton de calcul se base sur les gates validés
    output$ui_btn_calculer_spillover <- renderUI({
      gv     <- gates_valides() # Assurez-vous que gates_valides() est accessible ici
      canaux <- canaux_monomarques()
      ok     <- length(canaux) > 0 && all(canaux %in% names(gv))
      
      if (ok) {
        tagList(
          div(class="alert alert-success", style="font-size:11px; padding:5px 9px; margin-bottom:7px;",
              icon("check-circle"), " Gates validés. Prêt pour le calcul."),
          actionButton(ns("btn_calc_spillover"), "Calculer la matrice", 
                       class="btn-primary btn-block", icon=icon("calculator"))
        )
      } else {
        div(class="alert alert-warning", style="font-size:11px; padding:5px 9px;",
            icon("lock"), " Veuillez définir tous les gates dans l'onglet précédent.")
      }
    })
    
    # L'observeEvent du calcul
    observeEvent(input$btn_calc_spillover, {
      p <- pipeline()
      withProgress(message="Calcul en cours...", value=0.5, {
        tryCatch({
          p$calculer_spillover()
          spillover_trigger(spillover_trigger() + 1)
          showNotification("✔ Matrice de spillover calculée !", type="message")
        }, error=function(e) {
          showNotification(paste("Erreur :", conditionMessage(e)), type="error")
        })
      })
    })
    
    # 1. Sélecteur d'échantillon dynamique
    output$ui_select_echantillon <- renderUI({
      p <- pipeline()
      # On liste les échantillons disponibles dans l'objet R6
      selectInput(ns("sel_echantillon"), "Sélectionner l'échantillon à ajuster :", 
                  choices = names(p$echantillons))
    })
    
    # ── Personnalisation dynamique de l'affichage (ordre des canaux + transposition) ──
    # ── Personnalisation dynamique de l'affichage (ordre des canaux + transposition) ──
    ordre_canaux_rv <- reactiveVal(NULL)    # ordre courant d'affichage des canaux
    transpose_rv    <- reactiveVal(FALSE)   # inversion lignes / colonnes
    
    # Observateur pour initialiser/mettre à jour les choix du selectize quand la matrice change
    observe({
      spillover_trigger()
      p <- pipeline()
      req(!is.null(p$S_matrix))
      canaux <- rownames(p$S_matrix)
      
      ordre_actuel <- ordre_canaux_rv()
      # Si l'ordre interne est vide ou obsolète par rapport aux canaux réels, on le réinitialise
      if (is.null(ordre_actuel) || !setequal(ordre_actuel, canaux)) {
        ordre_canaux_rv(canaux)
        ordre_actuel <- canaux
      }
      
      updateSelectizeInput(
        session, 
        "select_ordre_canaux",
        choices = canaux,
        selected = ordre_actuel,
        server = FALSE
      )
    })
    
    # Capture le changement d'ordre suite au Drag & Drop de l'utilisateur
    observeEvent(input$select_ordre_canaux, {
      p <- pipeline()
      req(!is.null(p$S_matrix))
      canaux_reels <- rownames(p$S_matrix)
      
      # Sécurité : On s'assure qu'on ne traite pas un vecteur vide ou incomplet pendant un drag
      if (length(input$select_ordre_canaux) == length(canaux_reels) && 
          setequal(input$select_ordre_canaux, canaux_reels)) {
        ordre_canaux_rv(input$select_ordre_canaux)
      }
    }, ignoreNULL = TRUE)
    
    # Inversion de l'ordre
    observeEvent(input$btn_inverser_ordre, {
      ord <- ordre_canaux_rv()
      req(!is.null(ord))
      nouvel_ordre <- rev(ord)
      ordre_canaux_rv(nouvel_ordre)
      updateSelectizeInput(session, "select_ordre_canaux", selected = nouvel_ordre)
    })
    
    # Transposition
    observeEvent(input$btn_transposer_matrice, {
      transpose_rv(!transpose_rv())
    })
    
    # Matrice réellement affichée à l'écran
    mat_affichee <- reactive({
      spillover_trigger()
      req(input$sel_echantillon)
      
      p <- pipeline()
      validate(need(!is.null(p$S_matrix), "Veuillez d'abord calculer la matrice de spillover."))
      
      mat_brute <- if (!is.null(p$S_matrices_par_echantillon[[input$sel_echantillon]])) {
        p$S_matrices_par_echantillon[[input$sel_echantillon]]
      } else {
        p$S_matrix
      }
      
      canaux <- rownames(mat_brute)
      ordre  <- ordre_canaux_rv()
      
      # Fallback de sécurité au cas où l'ordre réactif n'est pas encore synchronisé
      if (is.null(ordre) || !setequal(ordre, canaux)) ordre <- canaux
      
      m <- round(mat_brute[ordre, ordre, drop = FALSE] * 100, 2)
      if (isTRUE(transpose_rv())) m <- t(m)
      m
    })
    
    # 2. RenderDT intelligent (lit la matrice spécifique ou la générale, avec ordre/orientation personnalisés)
    output$dt_spillover_matrix <- renderDT({
      mat_display <- mat_affichee()
      
      datatable(
        mat_display, 
        editable = TRUE,
        rownames = TRUE,
        selection = "none",
        options = list(
          dom = 't',           # 't' = Tableau seul (supprime pagination, recherche, info)
          ordering = FALSE,    # Empêche le tri pour conserver l'ordre des canaux
          scrollX = TRUE,      # Permet le défilement horizontal si besoin
          scrollY = "500px",   # Fixe la hauteur totale
          scrollCollapse = TRUE,
          paging = FALSE       # Désactive explicitement la pagination
        )
      ) %>% formatStyle(
        columns = colnames(mat_display),
        backgroundColor = styleInterval(c(5, 10, 20), c('#d1e7dd', '#fff3cd', '#f8d7da', '#e2c6c6'))
      )
    })
    
    # ── Export de la matrice affichée (PNG / PDF) ───────────────────────────
    generer_plot_matrice <- function(mat, titre = "Matrice de Spillover (%)") {
      df <- as.data.frame(as.table(mat))
      colnames(df) <- c("Ligne", "Colonne", "Valeur")
      df$Ligne     <- factor(df$Ligne,   levels = rev(rownames(mat)))
      df$Colonne   <- factor(df$Colonne, levels = colnames(mat))
      df$categorie <- cut(df$Valeur, breaks = c(-Inf, 5, 10, 20, Inf),
                          labels = c('#d1e7dd', '#fff3cd', '#f8d7da', '#e2c6c6'))
      
      ggplot(df, aes(x = Colonne, y = Ligne)) +
        geom_tile(aes(fill = categorie), color = "white") +
        scale_fill_identity() +
        geom_text(aes(label = sprintf("%.2f", Valeur)), size = 3.2) +
        labs(title = titre, x = "Canal secondaire (colonne)", y = "Canal principal (ligne)") +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold")
        )
    }
    
    output$dl_matrix_png <- downloadHandler(
      filename = function() {
        paste0("matrice_spillover_", gsub("[^a-zA-Z0-9_]", "_", input$sel_echantillon), "_",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        ggsave(file, plot = generer_plot_matrice(mat_affichee()), width = 8, height = 6, dpi = 300, bg = "white")
      }
    )
    
    output$dl_matrix_pdf <- downloadHandler(
      filename = function() {
        paste0("matrice_spillover_", gsub("[^a-zA-Z0-9_]", "_", input$sel_echantillon), "_",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        ggsave(file, plot = generer_plot_matrice(mat_affichee()), width = 8, height = 6)
      }
    )
    
    # 3. Mettez à jour le trigger lors du calcul
    observeEvent(input$btn_calc_spillover, {
      p <- pipeline()
      p$calculer_spillover()
      spillover_trigger(spillover_trigger() + 1) # Incrémente pour rafraîchir
      showNotification("✔ Matrice de spillover calculée !")
    })
    
    # ── Historique undo/redo par échantillon ────────────────────────────────
    # Chaque entrée est une copie de la matrice (en fraction, pas en %)
    historique      <- reactiveVal(list())   # pile des états passés
    historique_redo <- reactiveVal(list())   # pile des états annulés (pour redo)
    
    # Capture l'état courant de la matrice de l'échantillon sélectionné
    capturer_etat <- function(p, nom_ech) {
      mat <- if (!is.null(p$S_matrices_par_echantillon[[nom_ech]])) {
        p$S_matrices_par_echantillon[[nom_ech]]
      } else {
        p$S_matrix
      }
      mat
    }
    
    # 4. Édition manuelle d'une cellule
    observeEvent(input$dt_spillover_matrix_cell_edit, {
      info <- input$dt_spillover_matrix_cell_edit
      
      # Récupérer l'ordre actuel affiché à l'écran
      ordre_actuel <- ordre_canaux_rv()
      p <- pipeline()
      req(!is.null(p$S_matrix), input$sel_echantillon)
      
      # Identifier les canaux modifiés grâce au nom (et non à l'index brut)
      # info$row et info$col sont basés sur les lignes/colonnes affichées
      nom_ligne <- ordre_actuel[info$row]
      
      # Attention : si la matrice est transposée, les lignes deviennent les colonnes !
      if (isTRUE(transpose_rv())) {
        nom_colonne <- ordre_actuel[info$row]
        nom_ligne   <- ordre_actuel[info$col]
      } else {
        nom_colonne <- ordre_actuel[info$col]
      }
      
      # Nouvelle valeur entrée par l'utilisateur (divisée par 100 car affichée en %)
      nouvelle_valeur <- as.numeric(info$value) / 100
      if (is.na(nouvelle_valeur)) return()
      
      # Déterminer quelle matrice modifier (spécifique à l'échantillon ou globale)
      if (!is.null(p$S_matrices_par_echantillon[[input$sel_echantillon]])) {
        p$S_matrices_par_echantillon[[input$sel_echantillon]][nom_ligne, nom_colonne] <- nouvelle_valeur
      } else {
        p$S_matrix[nom_ligne, nom_colonne] <- nouvelle_valeur
      }
      
      # Déclencher la mise à jour du pipeline
      pipeline(p)
      spillover_trigger(spillover_trigger() + 1)
    })
    
    # Annuler (undo)
    observeEvent(input$btn_undo_spillover, {
      hist <- historique()
      if (length(hist) == 0) {
        showNotification("Aucune action à annuler.", type = "warning")
        return()
      }
      
      dernier        <- hist[[length(hist)]]
      historique(hist[-length(hist)])
      
      p   <- pipeline()
      nom <- dernier$echantillon
      
      # Sauvegarder l'état actuel dans redo avant d'écraser
      etat_actuel <- capturer_etat(p, nom)
      historique_redo(c(historique_redo(), list(list(echantillon = nom, matrice = etat_actuel))))
      
      # Restaurer la matrice
      p$S_matrices_par_echantillon[[nom]] <- dernier$matrice
      pipeline(p)
      spillover_trigger(spillover_trigger() + 1)
      showNotification("✔ Modification annulée.", type = "message")
    })
    
    # Rétablir (redo)
    observeEvent(input$btn_redo_spillover, {
      redo <- historique_redo()
      if (length(redo) == 0) {
        showNotification("Aucune action à rétablir.", type = "warning")
        return()
      }
      
      prochain          <- redo[[length(redo)]]
      historique_redo(redo[-length(redo)])
      
      p   <- pipeline()
      nom <- prochain$echantillon
      
      # Sauvegarder l'état actuel dans undo avant d'écraser
      etat_actuel <- capturer_etat(p, nom)
      historique(c(historique(), list(list(echantillon = nom, matrice = etat_actuel))))
      
      p$S_matrices_par_echantillon[[nom]] <- prochain$matrice
      pipeline(p)
      spillover_trigger(spillover_trigger() + 1)
      showNotification("✔ Modification rétablie.", type = "message")
    })
    
    # Bouton de validation finale
    observeEvent(input$btn_save_spillover, {
      showModal(modalDialog(
        title = "Succès",
        "La matrice de spillover pour l'échantillon a été validée et enregistrée avec succès.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # COMPARAISON DES MÉDIANES (AVANT / APRÈS COMPENSATION)
    # ════════════════════════════════════════════════════════════════════════
    medianes_trigger <- reactiveVal(0)
    
    output$ui_select_echantillon_medianes <- renderUI({
      p <- pipeline()
      selectInput(ns("sel_echantillon_medianes"), "Échantillon (matrice de compensation utilisée) :",
                  choices = names(p$echantillons))
    })
    
    output$ui_btn_comparer_medianes <- renderUI({
      spillover_trigger()
      p <- pipeline()
      if (is.null(p$S_matrix)) {
        return(div(class = "alert alert-warning", style = "font-size:11px; padding:5px 9px;",
                   icon("lock"), " Calculez d'abord la matrice de spillover dans l'onglet précédent."))
      }
      actionButton(ns("btn_comparer_medianes"), tagList(icon("balance-scale"), " Comparer les médianes"),
                   class = "btn-primary btn-block")
    })
    
    observeEvent(input$btn_comparer_medianes, {
      p <- pipeline()
      req(input$sel_echantillon_medianes)
      withProgress(message = "Comparaison des médianes en cours...", value = 0.5, {
        tryCatch({
          p$comparer_medianes_spillover(nom_echantillon = input$sel_echantillon_medianes)
          pipeline(p)
          medianes_trigger(medianes_trigger() + 1)
          showNotification("✔ Comparaison des médianes calculée !", type = "message")
        }, error = function(e) {
          showNotification(paste("Erreur :", conditionMessage(e)), type = "error")
        })
      })
    })
    
    # Résultat courant (avant / après / delta) pour l'échantillon sélectionné
    resultat_medianes <- reactive({
      medianes_trigger()
      req(input$sel_echantillon_medianes)
      p   <- pipeline()
      res <- p$comparaison_medianes[[input$sel_echantillon_medianes]]
      validate(need(!is.null(res), "Cliquez sur 'Comparer les médianes' pour lancer le calcul."))
      res
    })
    
    # Construit un tableau canal x canal affichant l'écart (médiane positive - médiane négative).
    # Plus la valeur est proche de 0, meilleure est la compensation. La diagonale (canal principal,
    # signal réel du marqueur) est affichée en gris neutre avec un tiret, car elle n'est pas concernée
    # par le critère d'égalité des médianes.
    afficher_matrice_medianes <- function(mat_ecart) {
      canaux <- colnames(mat_ecart)
      
      df_affichage <- as.data.frame(mat_ecart)
      colnames(df_affichage) <- canaux
      rownames(df_affichage) <- rownames(mat_ecart)
      
      datatable(
        df_affichage, rownames = TRUE, selection = "none",
        options = list(dom = 't', ordering = FALSE, scrollX = TRUE, scrollY = "300px",
                       scrollCollapse = TRUE, paging = FALSE)
      ) %>%
        formatStyle(
          columns = canaux,
          backgroundColor = styleInterval(c(-20, -5, 5, 20), c('#f8d7da', '#fff3cd', '#d1e7dd', '#fff3cd', '#f8d7da'))
        )
    }
    
    output$dt_medianes_avant <- renderDT({
      afficher_matrice_medianes(resultat_medianes()$avant$ecart)
    })
    
    output$dt_medianes_apres <- renderDT({
      afficher_matrice_medianes(resultat_medianes()$apres$ecart)
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # Compensation et biplots
    # ════════════════════════════════════════════════════════════════════════
    
    # Trigger incrémenté après chaque compensation pour forcer le re-rendu
    comp_trigger <- reactiveVal(0L)
    
    # 1. Mise à jour des sélecteurs de canaux
    observe({
      req(input$sel_echantillon_plot)
      p <- pipeline()
      fcs <- p$echantillons[[input$sel_echantillon_plot]]
      req(!is.null(fcs))
      
      all_cols    <- flowCore::colnames(fcs)
      canaux_fluo <- all_cols[!grepl("FSC|SSC|Time", all_cols, ignore.case = TRUE)]
      req(length(canaux_fluo) > 0)
      
      choices_list <- setNames(canaux_fluo, sapply(canaux_fluo, function(c) p$get_label(fcs, c)))
      
      updateSelectInput(session, "canal_x", choices = choices_list)
      updateSelectInput(session, "canal_y", choices = choices_list,
                        selected = choices_list[min(2, length(choices_list))])
    })
    
    output$ui_select_echantillon_plot <- renderUI({
      p <- pipeline()
      req(length(p$echantillons) > 0)
      selectInput(ns("sel_echantillon_plot"), "Sélectionner l'échantillon :", 
                  choices = names(p$echantillons))
    })
    
    # 2. Application de la compensation — utilise la matrice propre à chaque échantillon
    observeEvent(input$btn_apply_compensation, {
      p <- pipeline()
      withProgress(message = "Compensation en cours...", value = 0.5, {
        tryCatch({
          p$compenser()   # compenser() utilise déjà S_matrices_par_echantillon si disponible
          pipeline(p)
          comp_trigger(comp_trigger() + 1L)
          showNotification("✔ Compensation appliquée à la cohorte.", type = "message")
        }, error = function(e) {
          showNotification(paste("Erreur :", conditionMessage(e)), type = "error")
        })
      })
    })
    
    # 3. Calcul réactif du résultat de visualisation
    res_biplots <- reactive({
      comp_trigger()
      p <- pipeline()
      req(input$sel_echantillon_plot, input$affichage_type)
      
      canal_x <- if (isTRUE(input$mode_ensemble)) "ALL" else {
        req(input$canal_x); input$canal_x
      }
      canal_y <- if (isTRUE(input$mode_ensemble)) "ALL" else {
        req(input$canal_y); input$canal_y
      }
      
      need_apres <- input$affichage_type %in% c("Both", "After compensation only")
      if (need_apres && length(p$echantillons_traites) == 0) {
        validate("Cliquez sur 'Appliquer la compensation' pour afficher les données compensées.")
      }
      
      tryCatch(
        p$visualiser_compensation(
          nom_echantillon = input$sel_echantillon_plot,
          canal_x         = canal_x,
          canal_y         = canal_y,
          affichage       = input$affichage_type
        ),
        error = function(e) {
          showNotification(paste("Erreur visualisation :", conditionMessage(e)), type = "error")
          NULL
        }
      )
    })
    
    # 4. Conteneur UI : mode ALL → N plotOutputs dans un div scrollable
    #                   mode simple → un seul plotOutput
    output$ui_biplots_container <- renderUI({
      res <- res_biplots()
      
      if (is.null(res)) {
        return(div(style = "padding:40px; text-align:center; color:grey;",
                   "Aucune donnée à afficher. Vérifiez la sélection ou appliquez la compensation."))
      }
      
      if (inherits(res, "carrot_plots_list")) {
        # Mode ALL : on crée une ligne de plotOutputs par paire de plots
        n <- length(res)
        # On groupe par 2 (2 plots par ligne)
        n_lignes <- ceiling(n / 2)
        hauteur_plot <- 400   # px par plot individuel
        
        lignes <- lapply(seq_len(n_lignes), function(i) {
          idx_gauche <- (i - 1) * 2 + 1
          idx_droite <- idx_gauche + 1
          
          id_g <- paste0("biplot_", idx_gauche)
          id_d <- if (idx_droite <= n) paste0("biplot_", idx_droite) else NULL
          
          fluidRow(
            column(6, plotOutput(ns(id_g), width = "100%", height = paste0(hauteur_plot, "px"))),
            if (!is.null(id_d))
              column(6, plotOutput(ns(id_d), width = "100%", height = paste0(hauteur_plot, "px")))
          )
        })
        
        do.call(tagList, lignes)
        
      } else {
        # Mode simple : un seul plotOutput
        plotOutput(ns("biplot_single"), width = "100%", height = "600px")
      }
    })
    
    # 5. Observateur : branche les renderPlot sur les outputs créés dynamiquement
    observe({
      res <- res_biplots()
      req(!is.null(res))
      
      if (inherits(res, "carrot_plots_list")) {
        for (i in seq_along(res)) {
          local({
            idx   <- i
            plot_i <- res[[idx]]
            id    <- paste0("biplot_", idx)
            output[[id]] <- renderPlot({
              if (inherits(plot_i, c("gtable", "grob"))) {
                grid::grid.newpage(); grid::grid.draw(plot_i)
              } else {
                print(plot_i)
              }
            })
          })
        }
      } else {
        output$biplot_single <- renderPlot({
          if (is.null(res)) {
            plot.new()
            text(0.5, 0.5, "Aucune donnée.", cex = 1.1, col = "grey40")
            return()
          }
          if (inherits(res, c("gtable", "grob"))) {
            grid::grid.newpage(); grid::grid.draw(res)
          } else {
            print(res)
          }
        })
      }
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION EXPORT
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(comp_trigger(), {
      p <- pipeline()
      noms <- names(p$echantillons_traites)
      
      # On met à jour les choix tout en préservant la sélection précédente si possible
      updateCheckboxGroupInput(session, "sel_export_fcs",
                               choices = noms,
                               selected = if (is.null(input$sel_export_fcs)) noms else input$sel_export_fcs)
    }, ignoreInit = FALSE)
    
    # ── Sélecteur multi-échantillons ─────────────────────────────────────────
    output$ui_select_export_fcs <- renderUI({
      p <- pipeline()
      noms <- names(p$echantillons_traites)
      
      if (length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"),
                   " Appliquez d'abord la compensation pour pouvoir exporter des fichiers FCS."))
      }
      
      checkboxGroupInput(ns("sel_export_fcs"), "Échantillons à exporter :",
                         choices  = noms,
                         selected = noms)
    })
    
    # Tout sélectionner / désélectionner
    observeEvent(input$btn_select_all_fcs, {
      p <- pipeline()
      updateCheckboxGroupInput(session, "sel_export_fcs",
                               selected = names(p$echantillons_traites))
    })
    observeEvent(input$btn_deselect_all_fcs, {
      updateCheckboxGroupInput(session, "sel_export_fcs", selected = character(0))
    })
    
    # ── Bouton téléchargement FCS (affiché seulement si sélection non vide) ──
    output$ui_download_fcs <- renderUI({
      req(length(input$sel_export_fcs) > 0)
      n <- length(input$sel_export_fcs)
      downloadButton(ns("dl_fcs"),
                     label = paste0("Télécharger ", n, " fichier(s) FCS"),
                     class = "btn-primary",
                     style = "width:100%; font-weight:bold;")
    })
    
    output$dl_fcs <- downloadHandler(
      filename = function() {
        if (length(input$sel_export_fcs) == 1) {
          paste0(gsub("[^a-zA-Z0-9_]", "_", input$sel_export_fcs), "_compense.fcs")
        } else {
          paste0("FCS_Compenses_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
        }
      },
      content = function(file) {
        p <- pipeline()
        noms_sel <- input$sel_export_fcs
        req(length(noms_sel) > 0)
        
        if (length(noms_sel) == 1) {
          # Export direct d'un seul FCS
          nom        <- noms_sel
          fcs_obj    <- p$echantillons_traites[[nom]]
          if (!is.null(p$S_matrix)) {
            flowCore::keyword(fcs_obj)[["$SPILLOVER"]] <- p$S_matrix
          }
          flowCore::write.FCS(fcs_obj, filename = file)
          
        } else {
          # Export multiple → ZIP
          tmp_dir <- tempfile()
          dir.create(tmp_dir)
          on.exit(unlink(tmp_dir, recursive = TRUE))
          
          withProgress(message = "Préparation des fichiers FCS...", value = 0, {
            for (i in seq_along(noms_sel)) {
              nom     <- noms_sel[i]
              fcs_obj <- p$echantillons_traites[[nom]]
              if (is.null(fcs_obj)) next
              
              # Matrice spécifique à l'échantillon ou globale
              mat <- if (!is.null(p$S_matrices_par_echantillon[[nom]])) {
                p$S_matrices_par_echantillon[[nom]]
              } else {
                p$S_matrix
              }
              if (!is.null(mat)) flowCore::keyword(fcs_obj)[["$SPILLOVER"]] <- mat
              
              nom_fichier <- paste0(gsub("[^a-zA-Z0-9_]", "_", nom), "_compense.fcs")
              flowCore::write.FCS(fcs_obj, filename = file.path(tmp_dir, nom_fichier))
              incProgress(1 / length(noms_sel), detail = nom)
            }
          })
          
          fichiers <- list.files(tmp_dir, full.names = TRUE)
          zip::zip(zipfile = file, files = fichiers, mode = "cherry-pick")
        }
      },
      contentType = "application/octet-stream"
    )
    
    # ── Bouton téléchargement RDS ────────────────────────────────────────────
    output$ui_download_rds <- renderUI({
      p <- pipeline()
      # On désactive si rien n'a encore été fait
      pret <- !is.null(p$S_matrix) || !is.null(p$trans_list)
      if (!pret) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"),
                   " Aucun paramètre de compensation disponible à exporter."))
      }
      downloadButton(ns("dl_rds"),
                     label = "Télécharger la session RDS",
                     class = "btn-success",
                     style = "font-weight:bold; min-width:220px;")
    })
    
    output$dl_rds <- downloadHandler(
      filename = function() {
        nm <- trimws(input$rds_filename)
        if (nchar(nm) == 0) nm <- "Compensation_Session_Complete.rds"
        if (!grepl("\\.rds$", nm, ignore.case = TRUE)) nm <- paste0(nm, ".rds")
        nm
      },
      content = function(file) {
        p <- pipeline()
        withProgress(message = "Sérialisation de la session...", value = 0.5, {
          sauvegarde <- list(
            meta = list(
              date_export     = Sys.time(),
              canaux_utilises = p$canaux
            ),
            configuration_technique = list(
              trans_list        = p$trans_list,
              matrice_spillover = p$S_matrix,
              matrices_par_echantillon = p$S_matrices_par_echantillon
            ),
            gating = list(
              gates_positifs = p$gates_positifs,
              gates_negatifs = p$gates_negatifs
            ),
            visualisations = list(
              plots_gates        = p$plots_gates,
              plots_compensation = p$plots_compensation
            )
          )
          saveRDS(sauvegarde, file = file)
        })
      },
      contentType = "application/octet-stream"
    )
    
  })
}

# Opérateur null-coalescing utilitaire (évite d'importer rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b