library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

compensation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    
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
      
      tabPanel(title = "Matrice de Spillover",
               DTOutput(ns("dt_spillover_matrix"))
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
    
    # Bouton spillover
    output$ui_btn_calculer_spillover <- renderUI({
      gv     <- gates_valides()
      canaux <- canaux_monomarques()
      ok     <- length(canaux) > 0 && all(canaux %in% names(gv))
      
      if (ok) {
        tagList(
          div(class="alert alert-success",
              style="font-size:11px; padding:5px 9px; margin-bottom:7px;",
              icon("check-circle"), " Tous les gates sont définis."),
          actionButton(
            ns("btn_calc_spillover"),
            tagList(icon("calculator"), " Calculer la matrice de spillover"),
            class="btn-primary btn-block",
            style="font-weight:bold;"
          )
        )
      } else {
        div(style="font-size:11px; color:#888; padding:3px;",
            icon("lock"), " Définissez tous les gates pour débloquer.")
      }
    })
    
    observeEvent(input$btn_calc_spillover, {
      p <- pipeline()
      withProgress(message="Calcul de la matrice de spillover...", value=0.5, {
        tryCatch({
          p$calculer_spillover()
          showNotification("✔ Matrice de spillover calculée !",
                           type="message", duration=4)
        }, error=function(e) {
          showNotification(paste("Erreur :", conditionMessage(e)),
                           type="error")
        })
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # MATRICE DE SPILLOVER
    # ════════════════════════════════════════════════════════════════════════
    
    # Dans le module Server
    output$dt_spillover_matrix <- renderDT({
      req(pipeline()$S_matrix)
      # On récupère la matrice (ou celle spécifique à l'échantillon)
      mat <- pipeline()$S_matrix
      
      datatable(mat, editable = TRUE, options = list(pageLength = 20))
    })
    
    # Capture de l'édition
    observeEvent(input$dt_spillover_matrix_cell_edit, {
      info <- input$dt_spillover_matrix_cell_edit
      # info contient row, col, value
      
      canal1 <- rownames(pipeline()$S_matrix)[info$row]
      canal2 <- colnames(pipeline()$S_matrix)[info$col]
      valeur <- as.numeric(info$value)
      
      # Appel de votre méthode R6
      # Note : Vous devrez passer le nom de l'échantillon actif ici
      tryCatch({
        pipeline()$modifier_spillover(nom_echantillon = "Echantillon_A", # À dynamiser
                                      canal1 = canal1, 
                                      canal2 = canal2, 
                                      valeur = valeur)
        showNotification("Matrice mise à jour", type = "message")
      }, error = function(e) {
        showNotification(e$message, type = "error")
      })
    })
    
  })
}

# Opérateur null-coalescing utilitaire (évite d'importer rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b