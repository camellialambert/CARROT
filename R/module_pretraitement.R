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
                   numericInput(ns("debris_max_points"), "Points affichés :",
                                value = 10000, min = 1000, step = 1000),
                   
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
    
    # Sommets du polygone en cours de dessin (liste de listes x/y)
    sommets_rv        <- reactiveVal(list())
    # Nom du gate actuellement affiché dans le panneau résultat
    gate_actif_rv     <- reactiveVal(NULL)
    
    # ════════════════════════════════════════════════════════════════════════
    # RÉSOLUTION DE LA SOURCE
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
    # SÉLECTEURS DYNAMIQUES
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
    # DONNÉES POUR LE GRAPHIQUE INTERACTIF
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
      n   <- nrow(mat)
      nmax <- input$debris_max_points %||% 10000
      if (!is.null(p$seed)) set.seed(p$seed)
      idx <- if (n > nmax) sample(seq_len(n), nmax) else seq_len(n)
      df  <- as.data.frame(mat[idx, , drop = FALSE])
      colnames(df) <- c("X", "Y")
      df
    })
    
    # Réinitialise le polygone si l'utilisateur change de canal/échantillon/source
    observeEvent(list(input$sel_ech_debris, input$canal_x_debris,
                      input$canal_y_debris, input$source_debris), {
                        sommets_rv(list())
                      }, ignoreInit = TRUE)
    
    # ════════════════════════════════════════════════════════════════════════
    # GESTION DU POLYGONE — clic par clic + double-clic pour fermer
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
    # GRAPHIQUE INTERACTIF — scatter + polygone en cours
    # ════════════════════════════════════════════════════════════════════════
    
    output$plot_gate_dessin <- renderPlotly({
      df   <- donnees_plot()
      soms <- sommets_rv()
      req(nrow(df) > 0)
      
      p_obj <- carrot_obj()
      src   <- obtenir_source(p_obj, input$source_debris %||% "brutes")
      fcs   <- src[[input$sel_ech_debris]]
      lbl_x <- p_obj$get_label(fcs, input$canal_x_debris)
      lbl_y <- p_obj$get_label(fcs, input$canal_y_debris)
      
      # Scatter de base
      fig <- plot_ly(source = ns("plot_gate_dessin"),
                     type   = "scatter", mode = "markers") %>%
        add_trace(data = df, x = ~X, y = ~Y, name = "Cellules",
                  marker = list(size = 3, color = "#0077b6", opacity = 0.4),
                  hoverinfo = "none") %>%
        event_register("plotly_click") %>%
        event_register("plotly_doubleclick")
      
      # Superposition du polygone si au moins 2 sommets
      if (length(soms) >= 2) {
        xs <- c(sapply(soms, `[[`, "x"), soms[[1]]$x)   # ferme visuellement
        ys <- c(sapply(soms, `[[`, "y"), soms[[1]]$y)
        
        fig <- fig %>%
          add_trace(x = xs, y = ys, type = "scatter", mode = "lines+markers",
                    name = "Gate en cours",
                    line   = list(color = "#e65100", width = 2, dash = "dot"),
                    marker = list(size = 8, color = "#e65100", symbol = "circle"),
                    hoverinfo = "x+y") %>%
          # Shapes plotly éditables pour déplacer chaque sommet
          layout(
            shapes = lapply(seq_along(soms), function(i) {
              list(type = "circle",
                   x0 = soms[[i]]$x - 0.02 * diff(range(df$X)),
                   x1 = soms[[i]]$x + 0.02 * diff(range(df$X)),
                   y0 = soms[[i]]$y - 0.02 * diff(range(df$Y)),
                   y1 = soms[[i]]$y + 0.02 * diff(range(df$Y)),
                   xref = "x", yref = "y",
                   fillcolor = "rgba(230,101,0,0.3)",
                   line = list(color = "#e65100"))
            })
          ) %>%
          config(editable = TRUE, displayModeBar = TRUE)
      }
      
      fig %>% layout(
        dragmode = "zoom",
        xaxis = list(title = lbl_x),
        yaxis = list(title = lbl_y),
        legend = list(orientation = "h", y = -0.15),
        margin = list(b = 60)
      ) %>%
        config(displayModeBar = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE)
    })
    
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
    # ENREGISTREMENT DU GATE
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
    # PANNEAU RÉSULTAT — chips de navigation + graphique
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
        p$visualiser_debris(nom_echantillon = input$sel_ech_debris,
                            max_points       = input$debris_max_points %||% 10000),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
      )
      
      if (is.null(res)) {
        plot.new(); text(0.5, 0.5, "Aucune donnée.", cex = 1.1, col = "grey40"); return()
      }
      print(res)
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # RÉSUMÉ CUMULATIF — toutes les étapes du pipeline
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_recap_pipeline <- renderUI({
      bordures_trigger(); gates_trigger()
      p <- carrot_obj()
      req(input$sel_ech_debris)
      nom <- input$sel_ech_debris
      
      # Résolution du nombre de cellules à chaque étape
      etapes <- list()
      
      # Compensées brutes
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
      
      # Gates nommés (empilés)
      ref_gate <- if (length(etapes) > 0) etapes[[length(etapes)]]$n else ref
      for (nm_gate in names(p$gates_history)) {
        infos <- p$gates_history[[nm_gate]][[nom]]
        if (!is.null(infos)) {
          etapes[[nm_gate]] <- list(n = infos$n_apres, ref = infos$n_avant)
          ref_gate <- infos$n_apres
        }
      }
      
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
    })
    
  })
}

# Opérateur null-coalescing utilitaire
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}