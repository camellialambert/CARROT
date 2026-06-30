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
      .pretrait-instructions {
        background:#eef6fb; border-left:3px solid #0077b6;
        padding:8px 12px; border-radius:4px;
        font-size:12px; color:#444; margin-top:8px;
      }
      .pretrait-instructions b { color:#0077b6; }
    ")),
    
    tabBox(
      title = tagList(icon("filter"), "Prétraitement"),
      id = ns("pretrait_steps"), width = 12,
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 1 — RETRAIT DES DÉBRIS (bordures + gate polygonal)
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("broom"), " Retrait des Débris"),
        fluidRow(
          # ───────────────────────────────────────────────
          # COLONNE GAUCHE — Source + Paramètres
          # ───────────────────────────────────────────────
          column(width = 3,
                 wellPanel(
                   h4("Source des données"),
                   
                   div(class = "pretrait-instructions",
                       icon("info-circle"),
                       " Choisissez la source à partir de laquelle effectuer le retrait des débris."),
                   hr(),
                   
                   radioButtons(ns("source_debris"), "Données de départ :",
                                choices = c(
                                  "Données compensées (brutes)" = "brutes",
                                  "Après PeacoQC"               = "peacoqc",
                                  "Après flowAI"                = "flowai"
                                ),
                                selected = "brutes"),
                   hr(),
                   
                   h4("Retrait des bordures (Margins)"),
                   div(class = "pretrait-instructions",
                       icon("info-circle"),
                       " Cette étape s'applique automatiquement sur la ",
                       tags$b("dernière source disponible"),
                       " du pipeline (PeacoQC > flowAI > compensées), indépendamment ",
                       "du choix de source fait ci-dessous pour le gate de débris."),
                   checkboxInput(ns("activer_bordures"), "Retirer les événements de bordure (margins)",
                                 value = FALSE),
                   conditionalPanel(
                     condition = paste0("input['", ns("activer_bordures"), "']"),
                     uiOutput(ns("ui_canaux_bordures")),
                     actionButton(ns("btn_apply_bordures"), tagList(icon("play"), " Appliquer le retrait des bordures"),
                                  class = "btn-info", style = "width:100%; font-weight:bold;")
                   ),
                   hr(),
                   
                   h4("Gate de débris (polygone)"),
                   uiOutput(ns("ui_select_echantillon_debris")),
                   uiOutput(ns("ui_canal_x_debris")),
                   uiOutput(ns("ui_canal_y_debris")),
                   numericInput(ns("debris_max_points"), "Points affichés :",
                                value = 10000, min = 1000, step = 1000),
                   hr(),
                   
                   actionButton(ns("btn_reset_gate_debris"), tagList(icon("eraser"), " Effacer la sélection"),
                                class = "btn-default", style = "width:100%;"),
                   br(), br(),
                   actionButton(ns("btn_apply_debris"), tagList(icon("check-circle"), " Valider le gate de débris"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          
          # ───────────────────────────────────────────────
          # COLONNE CENTRALE — Graphique interactif
          # ───────────────────────────────────────────────
          column(width = 6,
                 box(title = tagList(icon("draw-polygon"), " Sélection interactive du gate"),
                     width = NULL, status = "info", solidHeader = TRUE,
                     
                     plotlyOutput(ns("plot_debris_interactif"), height = "500px"),
                     
                     div(class = "pretrait-instructions",
                         icon("hand-pointer"),
                         " Utilisez l'outil de sélection ", tags$b("Lasso"),
                         " (icône dans la barre d'outils en haut à droite du graphique) pour entourer ",
                         "la population de cellules à ", tags$b("conserver"),
                         ". Cliquez ensuite sur \"Valider le gate de débris\".")
                 ),
                 
                 box(title = "Résultat après filtration", width = NULL,
                     status = "success", solidHeader = TRUE,
                     plotOutput(ns("plot_debris_resultat"), width = "100%", height = "450px")
                 )
          ),
          
          # ───────────────────────────────────────────────
          # COLONNE DROITE — Récapitulatif
          # ───────────────────────────────────────────────
          column(width = 3,
                 box(title = "Résumé", width = NULL, status = "success", solidHeader = TRUE,
                     uiOutput(ns("ui_recap_debris"))
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
    
    # ── Accès à l'objet CARROT ────────────────────────────────────────────────
    carrot_obj <- reactive({
      pipeline_version()
      pipeline()
    })
    
    # ── Signaux locaux ────────────────────────────────────────────────────────
    bordures_trigger <- reactiveVal(0L)   # incrémenté après retirer_les_bordures()
    debris_trigger    <- reactiveVal(0L)  # incrémenté après retirer_les_debris()
    selection_lasso   <- reactiveVal(NULL) # stocke le data.frame des points sélectionnés au lasso
    
    # ════════════════════════════════════════════════════════════════════════
    # SOURCE DES DONNÉES — résolution selon le choix utilisateur
    # ════════════════════════════════════════════════════════════════════════
    
    # Renvoie la liste nommée des flowFrame correspondant à la source choisie par l'utilisateur.
    # Pour "peacoqc"/"flowai", utilise directement la liste demandée (avec repli sur la logique
    # pyramidale du pipeline si elle est vide) ; pour "brutes", suit la même logique pyramidale
    # que get_derniere_source() afin de refléter le retrait des bordures s'il a déjà été appliqué.
    obtenir_source_active <- function(p, choix) {
      if (choix == "peacoqc" && length(p$post_PeacoQC) > 0) return(p$post_PeacoQC)
      if (choix == "flowai"  && length(p$post_flowAI)  > 0) return(p$post_flowAI)
      return(p$get_derniere_source())
    }
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 1 — RETRAIT DES BORDURES (Margins)
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_canaux_bordures <- renderUI({
      p <- carrot_obj()
      source_liste <- p$get_derniere_source()
      req(length(source_liste) > 0)
      
      premier_fcs  <- source_liste[[1]]
      tous_canaux  <- flowCore::colnames(premier_fcs)
      
      selectInput(ns("canaux_bordures_sel"), "Canaux à vérifier (bordures) :",
                  choices  = tous_canaux,
                  selected = tous_canaux,
                  multiple = TRUE)
    })
    
    observeEvent(input$btn_apply_bordures, {
      p <- carrot_obj()
      req(input$canaux_bordures_sel)
      
      if (length(input$canaux_bordures_sel) < 1) {
        showNotification("Veuillez sélectionner au moins un canal.", type = "error")
        return(invisible(NULL))
      }
      
      withProgress(message = "Retrait des événements de bordure...", value = 0.3, {
        tryCatch({
          # retirer_les_bordures() s'appuie en interne sur get_derniere_source() pour
          # déterminer automatiquement la dernière étape valide du pipeline (PeacoQC > flowAI > brut)
          if (length(input$canaux_bordures_sel) >= 2) {
            p$retirer_les_bordures(
              canal1 = input$canaux_bordures_sel[1],
              canal2 = input$canaux_bordures_sel[2]
            )
          } else {
            p$retirer_les_bordures(
              canal1 = input$canaux_bordures_sel[1],
              canal2 = input$canaux_bordures_sel[1]
            )
          }
          pipeline(p)
          bordures_trigger(bordures_trigger() + 1L)
          showNotification("✔ Retrait des bordures appliqué.", type = "message")
        }, error = function(e) {
          showNotification(paste("Erreur retrait des bordures :", conditionMessage(e)), type = "error")
        })
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 2 — GATE DE DÉBRIS (sélection polygonale interactive)
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_select_echantillon_debris <- renderUI({
      p <- carrot_obj()
      source_liste <- obtenir_source_active(p, input$source_debris %||% "brutes")
      req(length(source_liste) > 0)
      selectInput(ns("sel_echantillon_debris"), "Échantillon :",
                  choices = names(source_liste))
    })
    
    output$ui_canal_x_debris <- renderUI({
      p <- carrot_obj()
      source_liste <- obtenir_source_active(p, input$source_debris %||% "brutes")
      req(length(source_liste) > 0, input$sel_echantillon_debris)
      fcs <- source_liste[[input$sel_echantillon_debris]]
      req(!is.null(fcs))
      
      tous_canaux <- flowCore::colnames(fcs)
      defaut_fsc  <- grep("FSC", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      
      selectInput(ns("canal_x_debris"), "Canal X :",
                  choices  = tous_canaux,
                  selected = if (!is.na(defaut_fsc)) defaut_fsc else tous_canaux[1])
    })
    
    output$ui_canal_y_debris <- renderUI({
      p <- carrot_obj()
      source_liste <- obtenir_source_active(p, input$source_debris %||% "brutes")
      req(length(source_liste) > 0, input$sel_echantillon_debris)
      fcs <- source_liste[[input$sel_echantillon_debris]]
      req(!is.null(fcs))
      
      tous_canaux <- flowCore::colnames(fcs)
      defaut_ssc  <- grep("SSC", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      
      selectInput(ns("canal_y_debris"), "Canal Y :",
                  choices  = tous_canaux,
                  selected = if (!is.na(defaut_ssc)) defaut_ssc else tous_canaux[min(2, length(tous_canaux))])
    })
    
    # Données sous-échantillonnées du graphique interactif (réactif aux choix utilisateur)
    donnees_debris_interactif <- reactive({
      bordures_trigger()
      p <- carrot_obj()
      req(input$sel_echantillon_debris, input$canal_x_debris, input$canal_y_debris)
      
      source_liste <- obtenir_source_active(p, input$source_debris %||% "brutes")
      fcs <- source_liste[[input$sel_echantillon_debris]]
      req(!is.null(fcs))
      
      cx <- input$canal_x_debris
      cy <- input$canal_y_debris
      req(cx %in% flowCore::colnames(fcs), cy %in% flowCore::colnames(fcs))
      
      mat   <- flowCore::exprs(fcs)[, c(cx, cy), drop = FALSE]
      n_tot <- nrow(mat)
      n_max <- input$debris_max_points %||% 10000
      
      if (!is.null(p$seed)) set.seed(p$seed)
      idx <- if (n_tot > n_max) sample(seq_len(n_tot), n_max) else seq_len(n_tot)
      
      df <- as.data.frame(mat[idx, , drop = FALSE])
      colnames(df) <- c("X", "Y")
      df$key <- idx   # clé stable pour retrouver les points sélectionnés
      df
    })
    
    # Réinitialise la sélection lasso quand l'utilisateur change d'échantillon/canaux/source
    observeEvent(list(input$sel_echantillon_debris, input$canal_x_debris,
                      input$canal_y_debris, input$source_debris), {
                        selection_lasso(NULL)
                      }, ignoreInit = TRUE)
    
    observeEvent(input$btn_reset_gate_debris, {
      selection_lasso(NULL)
      showNotification("Sélection effacée.", type = "message")
    })
    
    # Graphique interactif plotly — mode lasso pour entourer la population à conserver
    output$plot_debris_interactif <- renderPlotly({
      df <- donnees_debris_interactif()
      req(nrow(df) > 0)
      
      p <- carrot_obj()
      source_liste <- obtenir_source_active(p, input$source_debris %||% "brutes")
      fcs <- source_liste[[input$sel_echantillon_debris]]
      lbl_x <- p$get_label(fcs, input$canal_x_debris)
      lbl_y <- p$get_label(fcs, input$canal_y_debris)
      
      plot_ly(df, x = ~X, y = ~Y, key = ~key, type = "scatter", mode = "markers",
              source = ns("plot_debris_interactif"),
              marker = list(size = 3, color = "#0077b6", opacity = 0.45)) %>%
        layout(
          dragmode = "lasso",
          xaxis = list(title = lbl_x),
          yaxis = list(title = lbl_y)
        ) %>%
        event_register("plotly_selected")
    })
    
    # Capture la sélection lasso et la conserve en mémoire (persiste tant que non réinitialisée)
    observeEvent(event_data("plotly_selected", source = ns("plot_debris_interactif")), {
      ev <- event_data("plotly_selected", source = ns("plot_debris_interactif"))
      req(ev)
      if (nrow(ev) < 3) {
        showNotification("Sélectionnez au moins 3 points pour former un polygone.", type = "warning")
        return(invisible(NULL))
      }
      selection_lasso(ev)
    })
    
    # Validation du gate : calcule l'enveloppe convexe des points sélectionnés
    # et l'envoie à retirer_les_debris()
    observeEvent(input$btn_apply_debris, {
      sel <- selection_lasso()
      
      if (is.null(sel) || nrow(sel) < 3) {
        showNotification("Veuillez d'abord sélectionner une zone au lasso sur le graphique.", type = "error")
        return(invisible(NULL))
      }
      
      p <- carrot_obj()
      
      # Enveloppe convexe des points sélectionnés → polygone exploitable par flowCore::polygonGate
      idx_hull   <- grDevices::chull(sel$x, sel$y)
      polygone   <- as.matrix(sel[idx_hull, c("x", "y")])
      colnames(polygone) <- c(input$canal_x_debris, input$canal_y_debris)
      
      withProgress(message = "Application du gate de débris...", value = 0.3, {
        tryCatch({
          p$retirer_les_debris(
            matrice_points    = polygone,
            canal_x           = input$canal_x_debris,
            canal_y           = input$canal_y_debris,
            source_nettoyage  = input$source_debris %||% "brutes"
          )
          pipeline(p)
          debris_trigger(debris_trigger() + 1L)
          showNotification("✔ Gate de débris appliqué à la cohorte.", type = "message")
        }, error = function(e) {
          showNotification(paste("Erreur retrait des débris :", conditionMessage(e)), type = "error")
        })
      })
    })
    
    # ── Rendu du résultat après filtration (biplot avec contour du polygone) ──
    output$plot_debris_resultat <- renderPlot({
      debris_trigger()
      p <- carrot_obj()
      req(input$sel_echantillon_debris)
      
      if (length(p$post_debris) == 0 || is.null(p$post_debris[[input$sel_echantillon_debris]])) {
        validate("Cliquez sur 'Valider le gate de débris' pour afficher le résultat.")
      }
      
      res <- tryCatch(
        p$visualiser_debris(
          nom_echantillon = input$sel_echantillon_debris,
          max_points       = input$debris_max_points %||% 10000
        ),
        error = function(e) {
          showNotification(paste("Erreur visualisation débris :", conditionMessage(e)), type = "error")
          NULL
        }
      )
      
      if (is.null(res)) {
        plot.new()
        text(0.5, 0.5, "Aucune donnée à afficher.", cex = 1.1, col = "grey40")
        return()
      }
      print(res)
    })
    
    # ── Récapitulatif des taux de conservation pour tous les échantillons traités ──
    output$ui_recap_debris <- renderUI({
      debris_trigger()
      p <- carrot_obj()
      
      if (length(p$post_debris) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucun résultat disponible."))
      }
      
      source_liste <- obtenir_source_active(p, input$source_debris %||% "brutes")
      
      lignes <- lapply(names(p$post_debris), function(nom) {
        avant <- source_liste[[nom]]
        apres <- p$post_debris[[nom]]
        if (is.null(avant) || is.null(apres)) return(NULL)
        
        n_avant <- nrow(flowCore::exprs(avant))
        n_apres <- nrow(flowCore::exprs(apres))
        pct     <- if (n_avant > 0) round(n_apres / n_avant * 100, 1) else 0
        
        couleur <- if (pct >= 90) "#2e7d32" else if (pct >= 70) "#e65100" else "#c62828"
        
        div(style = "margin-bottom:8px; font-size:12px;",
            strong(nom), br(),
            span(style = paste0("color:", couleur, ";"),
                 format(n_apres, big.mark = " "), " / ", format(n_avant, big.mark = " "),
                 " (", pct, "%)")
        )
      })
      
      tagList(Filter(Negate(is.null), lignes))
    })
    
  })
}

# Opérateur null-coalescing utilitaire (présent aussi dans module_compensation.R et module_qc.R)
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}