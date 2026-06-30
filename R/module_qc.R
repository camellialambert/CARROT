library(shiny)
library(shinydashboard)
library(shinyjs)
library(base64enc)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

qc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    
    tags$style(HTML("
      .qc-instructions {
        background:#eef6fb; border-left:3px solid #0077b6;
        padding:8px 12px; border-radius:4px;
        font-size:12px; color:#444; margin-top:8px;
      }
      .qc-instructions b { color:#0077b6; }
    ")),
    
    tabBox(
      title = tagList(icon("broom"), "Contrôle Qualité"),
      id = ns("qc_steps"), width = 12,
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 1 — PeacoQC
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("wave-square"), " PeacoQC"),
        fluidRow(
          # --- Colonne de Gauche : Paramètres ---
          column(width = 3,
                 wellPanel(
                   h4("Paramètres PeacoQC"),
                   
                   div(class = "qc-instructions",
                       icon("info-circle"),
                       " PeacoQC détecte les ", tags$b("instabilités temporelles"),
                       " du signal (flux, débit) sur les données compensées."),
                   hr(),
                   
                   selectInput(ns("peaco_determine_good_cells"), "Cellules évaluées :",
                               choices = c("all", "channels"), selected = "all"),
                   
                   numericInput(ns("peaco_min_cells"), "Cellules minimales par bin :",
                                value = 150, min = 1, step = 10),
                   numericInput(ns("peaco_max_bins"), "Nombre maximal de bins :",
                                value = 100, min = 1, step = 5),
                   numericInput(ns("peaco_step"), "Pas glissant (step) :",
                                value = 500, min = 1, step = 50),
                   numericInput(ns("peaco_MAD"), "Multiplicateur MAD :",
                                value = 6, min = 0, step = 0.5),
                   numericInput(ns("peaco_IT_limit"), "Limite IT (Intervalle de Temps) :",
                                value = 0.6, min = 0, max = 1, step = 0.05),
                   numericInput(ns("peaco_consecutive_bins"), "Bins consécutifs tolérés :",
                                value = 5, min = 1, step = 1),
                   checkboxInput(ns("peaco_remove_zeros"), "Retirer les valeurs nulles/négatives",
                                 value = FALSE),
                   numericInput(ns("peaco_force_IT"), "Forçage IT (force_IT) :",
                                value = 150, min = 0, step = 10),
                   numericInput(ns("peaco_peak_removal"), "Fraction de retrait des pics :",
                                value = round(1/3, 4), min = 0, max = 1, step = 0.01),
                   numericInput(ns("peaco_min_nr_bins_peakdetection"), "Bins min. détection de pics :",
                                value = 10, min = 1, step = 1),
                   hr(),
                   
                   actionButton(ns("btn_reset_peacoqc"), "Réinitialiser les valeurs par défaut",
                                class = "btn-default btn-sm", style = "width:100%;"),
                   br(), br(),
                   
                   actionButton(ns("btn_apply_peacoqc"), tagList(icon("play"), " Appliquer PeacoQC"),
                                class = "btn-warning", style = "width:100%; font-weight:bold;")
                 )
          ),
          
          # --- Colonne Centrale : Sélection échantillon + Visualisation ---
          column(width = 6,
                 wellPanel(
                   style = "padding:10px 15px;",
                   fluidRow(
                     column(width = 8, uiOutput(ns("ui_select_echantillon_peacoqc"))),
                     column(width = 4,
                            numericInput(ns("peaco_max_points"), "Points affichés :",
                                         value = 10000, min = 1000, step = 1000))
                   )
                 ),
                 box(title = "Visualisation cinétique (Temps vs FSC)", width = NULL,
                     status = "warning", solidHeader = TRUE,
                     plotOutput(ns("plot_peacoqc"), width = "100%", height = "550px")
                 ),
                 box(title = tagList(icon("image"), " Rapport diagnostique natif PeacoQC"), width = NULL,
                     status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                     uiOutput(ns("ui_plot_peacoqc_natif"))
                 )
          ),
          
          # --- Colonne de Droite : Récapitulatif ---
          column(width = 3,
                 box(title = "Résumé", width = NULL, status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_recap_peacoqc"))
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 2 — flowAI
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("filter"), " flowAI"),
        fluidRow(
          # --- Colonne de Gauche : Paramètres ---
          column(width = 3,
                 wellPanel(
                   h4("Paramètres flowAI"),
                   
                   div(class = "qc-instructions",
                       icon("info-circle"),
                       " flowAI évalue trois critères : ", tags$b("débit (Flow Rate)"),
                       ", ", tags$b("signal (Flow Signal)"), " et ",
                       tags$b("temps de vol (Flow Margin)"), "."),
                   hr(),
                   
                   selectInput(ns("flowai_remove_from"), "Critères évalués (remove_from) :",
                               choices = c("all", "FR", "FS", "FM", "FR_FS", "FR_FM", "FS_FM"),
                               selected = "all"),
                   
                   h5(tags$b("Flow Rate (débit)")),
                   numericInput(ns("flowai_second_fractionFR"), "Fraction de seconde :",
                                value = 0.1, min = 0.01, step = 0.01),
                   numericInput(ns("flowai_alphaFR"), "Seuil alpha :",
                                value = 0.01, min = 0, max = 1, step = 0.01),
                   
                   h5(tags$b("Flow Signal")),
                   numericInput(ns("flowai_max_cptFS"), "Nombre max. de points de changement :",
                                value = 3, min = 1, step = 1),
                   numericInput(ns("flowai_pen_valueFS"), "Valeur de pénalité :",
                                value = 500, min = 1, step = 50),
                   textInput(ns("flowai_ChExcludeFS"), "Canaux exclus (séparés par virgule) :",
                             value = "FSC,SSC"),
                   
                   h5(tags$b("Flow Margin")),
                   selectInput(ns("flowai_neg_valuesFM"), "Gestion valeurs négatives :",
                               choices = c("1" = 1, "2" = 2), selected = 1),
                   hr(),
                   
                   actionButton(ns("btn_reset_flowai"), "Réinitialiser les valeurs par défaut",
                                class = "btn-default btn-sm", style = "width:100%;"),
                   br(), br(),
                   
                   actionButton(ns("btn_apply_flowai"), tagList(icon("play"), " Appliquer flowAI"),
                                class = "btn-warning", style = "width:100%; font-weight:bold;")
                 )
          ),
          
          # --- Colonne Centrale : Sélection échantillon + Visualisation ---
          column(width = 6,
                 wellPanel(
                   style = "padding:10px 15px;",
                   fluidRow(
                     column(width = 8, uiOutput(ns("ui_select_echantillon_flowai"))),
                     column(width = 4,
                            numericInput(ns("flowai_max_points"), "Points affichés :",
                                         value = 10000, min = 1000, step = 1000))
                   )
                 ),
                 box(title = "Visualisation cinétique (Temps vs FSC)", width = NULL,
                     status = "warning", solidHeader = TRUE,
                     plotOutput(ns("plot_flowai"), width = "100%", height = "550px")
                 )
          ),
          
          # --- Colonne de Droite : Récapitulatif ---
          column(width = 3,
                 box(title = "Résumé", width = NULL, status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_recap_flowai"))
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

qc_server <- function(id, pipeline, pipeline_version) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ── Accès à l'objet CARROT ────────────────────────────────────────────────
    carrot_obj <- reactive({
      pipeline_version()
      pipeline()
    })
    
    # ── Signaux locaux : incrémentés après chaque exécution pour forcer le re-rendu ──
    peacoqc_trigger <- reactiveVal(0L)
    flowai_trigger  <- reactiveVal(0L)
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 1 — PeacoQC
    # ════════════════════════════════════════════════════════════════════════
    
    # Valeurs par défaut (synchronisées avec parametres_par_defaut de appliquer_peacoqc)
    defaults_peacoqc <- list(
      determine_good_cells = "all", min_cells = 150, max_bins = 100, step = 500,
      MAD = 6, IT_limit = 0.6, consecutive_bins = 5, remove_zeros = FALSE,
      force_IT = 150, peak_removal = 1/3, min_nr_bins_peakdetection = 10
    )
    
    observeEvent(input$btn_reset_peacoqc, {
      updateSelectInput(session, "peaco_determine_good_cells", selected = defaults_peacoqc$determine_good_cells)
      updateNumericInput(session, "peaco_min_cells", value = defaults_peacoqc$min_cells)
      updateNumericInput(session, "peaco_max_bins", value = defaults_peacoqc$max_bins)
      updateNumericInput(session, "peaco_step", value = defaults_peacoqc$step)
      updateNumericInput(session, "peaco_MAD", value = defaults_peacoqc$MAD)
      updateNumericInput(session, "peaco_IT_limit", value = defaults_peacoqc$IT_limit)
      updateNumericInput(session, "peaco_consecutive_bins", value = defaults_peacoqc$consecutive_bins)
      updateCheckboxInput(session, "peaco_remove_zeros", value = defaults_peacoqc$remove_zeros)
      updateNumericInput(session, "peaco_force_IT", value = defaults_peacoqc$force_IT)
      updateNumericInput(session, "peaco_peak_removal", value = round(defaults_peacoqc$peak_removal, 4))
      updateNumericInput(session, "peaco_min_nr_bins_peakdetection", value = defaults_peacoqc$min_nr_bins_peakdetection)
      showNotification("Paramètres PeacoQC réinitialisés.", type = "message")
    })
    
    # Sélecteur d'échantillon (alimenté par les échantillons compensés disponibles)
    output$ui_select_echantillon_peacoqc <- renderUI({
      p <- carrot_obj()
      req(length(p$echantillons_traites) > 0)
      selectInput(ns("sel_echantillon_peacoqc"), "Échantillon à visualiser :",
                  choices = names(p$echantillons_traites))
    })
    
    # Application de PeacoQC sur l'ensemble de la cohorte avec les paramètres utilisateur
    observeEvent(input$btn_apply_peacoqc, {
      p <- carrot_obj()
      
      if (length(p$echantillons_traites) == 0) {
        showNotification("Veuillez d'abord appliquer la compensation (onglet Compensation).",
                         type = "error")
        return(invisible(NULL))
      }
      
      if (is.null(p$canaux) || length(p$canaux) == 0) {
        showNotification("Aucun canal détecté (étape d'importation incomplète).", type = "error")
        return(invisible(NULL))
      }
      
      reglages <- list(
        determine_good_cells      = input$peaco_determine_good_cells,
        min_cells                 = input$peaco_min_cells,
        max_bins                  = input$peaco_max_bins,
        step                      = input$peaco_step,
        MAD                       = input$peaco_MAD,
        IT_limit                  = input$peaco_IT_limit,
        consecutive_bins          = input$peaco_consecutive_bins,
        remove_zeros              = input$peaco_remove_zeros,
        force_IT                  = input$peaco_force_IT,
        peak_removal              = input$peaco_peak_removal,
        min_nr_bins_peakdetection = input$peaco_min_nr_bins_peakdetection
      )
      
      withProgress(message = "Exécution de PeacoQC en cours...", value = 0.3, {
        tryCatch({
          p$appliquer_peacoqc(reglages_specifiques = reglages)
          pipeline(p)
          peacoqc_trigger(peacoqc_trigger() + 1L)
          showNotification("✔ PeacoQC appliqué à la cohorte.", type = "message")
        }, error = function(e) {
          showNotification(paste("Erreur PeacoQC :", conditionMessage(e)), type = "error")
        })
      })
    })
    
    # Rendu du graphique cinétique (réactif au changement d'échantillon ou de trigger)
    output$plot_peacoqc <- renderPlot({
      peacoqc_trigger()
      p <- carrot_obj()
      req(input$sel_echantillon_peacoqc)
      
      if (length(p$post_PeacoQC) == 0 || is.null(p$post_PeacoQC[[input$sel_echantillon_peacoqc]])) {
        validate("Cliquez sur 'Appliquer PeacoQC' pour afficher le résultat du contrôle qualité.")
      }
      
      res <- tryCatch(
        p$visualiser_peacoqc(
          nom_echantillon = input$sel_echantillon_peacoqc,
          max_points       = input$peaco_max_points %||% 10000
        ),
        error = function(e) {
          showNotification(paste("Erreur visualisation PeacoQC :", conditionMessage(e)), type = "error")
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
    
    # Rapport diagnostique natif PeacoQC (PNG généré par PeacoQC::PeacoQC, encodé en base64 pour affichage inline)
    output$ui_plot_peacoqc_natif <- renderUI({
      peacoqc_trigger()
      p <- carrot_obj()
      req(input$sel_echantillon_peacoqc)
      
      chemin_png <- p$plots_peacoqc_natif[[input$sel_echantillon_peacoqc]]
      
      if (is.null(chemin_png) || !file.exists(chemin_png)) {
        return(div(style = "padding:30px; text-align:center; color:grey; font-size:12px;",
                   icon("image"),
                   " Le rapport natif PeacoQC apparaîtra ici après l'exécution de l'algorithme."))
      }
      
      donnees_brutes <- readBin(chemin_png, "raw", file.info(chemin_png)$size)
      encode_b64     <- base64enc::base64encode(donnees_brutes)
      
      tags$img(src = paste0("data:image/png;base64,", encode_b64),
               style = "width:100%; height:auto; border:1px solid #eee; border-radius:4px;")
    })
    
    # Récapitulatif des taux de conservation pour tous les échantillons traités
    output$ui_recap_peacoqc <- renderUI({
      peacoqc_trigger()
      p <- carrot_obj()
      
      if (length(p$post_PeacoQC) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucun résultat disponible."))
      }
      
      lignes <- lapply(names(p$post_PeacoQC), function(nom) {
        avant <- p$echantillons_traites[[nom]]
        apres <- p$post_PeacoQC[[nom]]
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
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION 2 — flowAI
    # ════════════════════════════════════════════════════════════════════════
    
    defaults_flowai <- list(
      remove_from = "all", second_fractionFR = 0.1, alphaFR = 0.01,
      max_cptFS = 3, pen_valueFS = 500, neg_valuesFM = 1, ChExcludeFS = c("FSC", "SSC")
    )
    
    observeEvent(input$btn_reset_flowai, {
      updateSelectInput(session, "flowai_remove_from", selected = defaults_flowai$remove_from)
      updateNumericInput(session, "flowai_second_fractionFR", value = defaults_flowai$second_fractionFR)
      updateNumericInput(session, "flowai_alphaFR", value = defaults_flowai$alphaFR)
      updateNumericInput(session, "flowai_max_cptFS", value = defaults_flowai$max_cptFS)
      updateNumericInput(session, "flowai_pen_valueFS", value = defaults_flowai$pen_valueFS)
      updateTextInput(session, "flowai_ChExcludeFS", value = paste(defaults_flowai$ChExcludeFS, collapse = ","))
      updateSelectInput(session, "flowai_neg_valuesFM", selected = defaults_flowai$neg_valuesFM)
      showNotification("Paramètres flowAI réinitialisés.", type = "message")
    })
    
    output$ui_select_echantillon_flowai <- renderUI({
      p <- carrot_obj()
      req(length(p$echantillons_traites) > 0)
      selectInput(ns("sel_echantillon_flowai"), "Échantillon à visualiser :",
                  choices = names(p$echantillons_traites))
    })
    
    observeEvent(input$btn_apply_flowai, {
      p <- carrot_obj()
      
      if (length(p$echantillons_traites) == 0) {
        showNotification("Veuillez d'abord appliquer la compensation (onglet Compensation).",
                         type = "error")
        return(invisible(NULL))
      }
      
      ch_exclude <- trimws(strsplit(input$flowai_ChExcludeFS, ",")[[1]])
      ch_exclude <- ch_exclude[nchar(ch_exclude) > 0]
      
      reglages <- list(
        second_fractionFR = input$flowai_second_fractionFR,
        alphaFR           = input$flowai_alphaFR,
        max_cptFS         = input$flowai_max_cptFS,
        pen_valueFS       = input$flowai_pen_valueFS,
        neg_valuesFM      = as.numeric(input$flowai_neg_valuesFM),
        ChExcludeFS       = ch_exclude
      )
      
      withProgress(message = "Exécution de flowAI en cours...", value = 0.3, {
        tryCatch({
          p$appliquer_flowai(reglages_specifiques = reglages)
          pipeline(p)
          flowai_trigger(flowai_trigger() + 1L)
          showNotification("✔ flowAI appliqué à la cohorte.", type = "message")
        }, error = function(e) {
          showNotification(paste("Erreur flowAI :", conditionMessage(e)), type = "error")
        })
      })
    })
    
    output$plot_flowai <- renderPlot({
      flowai_trigger()
      p <- carrot_obj()
      req(input$sel_echantillon_flowai)
      
      if (length(p$post_flowAI) == 0 || is.null(p$post_flowAI[[input$sel_echantillon_flowai]])) {
        validate("Cliquez sur 'Appliquer flowAI' pour afficher le résultat du contrôle qualité.")
      }
      
      res <- tryCatch(
        p$visualiser_flowai(
          nom_echantillon = input$sel_echantillon_flowai,
          max_points       = input$flowai_max_points %||% 10000
        ),
        error = function(e) {
          showNotification(paste("Erreur visualisation flowAI :", conditionMessage(e)), type = "error")
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
    
    output$ui_recap_flowai <- renderUI({
      flowai_trigger()
      p <- carrot_obj()
      
      if (length(p$post_flowAI) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucun résultat disponible."))
      }
      
      lignes <- lapply(names(p$post_flowAI), function(nom) {
        avant <- p$echantillons_traites[[nom]]
        apres <- p$post_flowAI[[nom]]
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

# Opérateur null-coalescing utilitaire (présent aussi dans module_compensation.R)
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}