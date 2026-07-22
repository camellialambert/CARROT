library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(shinyFiles)
library(DT)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

demixage_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    
    tabBox(
      title = tagList(icon("dna"), " AutoSpectral - Démixage Spectral"),
      id = ns("tabs_demix"), width = 12,
      
      # ───────────────────────────────────────────────────────────────
      # ONGLET 1 — PARAMÈTRES
      # ───────────────────────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("folder-open"), " Paramètres"),
        
        # ── Bulle 1 : Dossiers, pleine largeur ─────────────────────────────
        fluidRow(
          column(
            width = 12,
            box(title = tagList(icon("folder-open"), " Dossiers"), width = NULL, status = "primary", solidHeader = TRUE,
                fluidRow(
                  column(4,
                         h5("Dossier racine (sorties AutoSpectral)", style = "margin-top:0;"),
                         shinyDirButton(ns("dir_root"), "Parcourir…", "Sélectionner le dossier racine",
                                        icon = icon("folder-open"), class = "btn-primary", style = "width:100%;"),
                         br(), uiOutput(ns("root_display"))
                  ),
                  column(4,
                         h5("Contrôles monomarqués / unstained", style = "margin-top:0;"),
                         shinyDirButton(ns("dir_monomarques"), "Parcourir…", "Sélectionner le dossier des contrôles",
                                        icon = icon("folder-open"), class = "btn-primary", style = "width:100%;"),
                         br(), uiOutput(ns("monomarques_display")),
                         uiOutput(ns("unstained_picker"))
                  ),
                  column(4,
                         h5("Échantillons biologiques", style = "margin-top:0;"),
                         shinyDirButton(ns("dir_echantillons"), "Parcourir…", "Sélectionner le dossier des échantillons",
                                        icon = icon("folder-open"), class = "btn-primary", style = "width:100%;"),
                         br(), uiOutput(ns("echantillons_display"))
                  )
                ),
                hr(),
                fluidRow(
                  column(4,
                         h5(tagList(icon("microscope"), " Type de cytomètre spectral")),
                         selectInput(ns("type_cytometre"), NULL,
                                     choices = c("aurora", "a8", "s8", "a5se", "id7000", "mosaic", "opteon", "xenith", "minimal", "auroraNL"),
                                     selected = "aurora")
                  ),
                  column(4,
                         br(),
                         actionButton(ns("btn_lancer_asp"), tagList(icon("play"), " Lancer lancer_asp()"),
                                      class = "btn-success", style = "width:100%; font-weight:bold; margin-top:5px;")
                  ),
                  column(4, uiOutput(ns("asp_status")))
                ))
          )
        ),
        
        # ── Bulle 2 : Aperçu, sur toute la largeur de la page ──────────────
        fluidRow(
          column(
            width = 12,
            box(title = tagList(icon("table"), " Aperçu et édition du fichier de contrôle (fcs_control_file.csv)"),
                width = NULL, status = "primary", solidHeader = TRUE,
                fluidRow(
                  column(6, actionButton(ns("btn_recharger_csv"), tagList(icon("rotate"), " Recharger depuis le disque"),
                                         class = "btn-default")),
                  column(6, actionButton(ns("btn_save_csv"), tagList(icon("save"), " Enregistrer les modifications"),
                                         class = "btn-success", style = "float:right;"))
                ),
                hr(),
                fluidRow(
                  column(3, actionButton(ns("btn_ajouter_ligne"), tagList(icon("plus"), " Ajouter une ligne"),
                                         class = "btn-default", style = "width:100%;")),
                  column(3, actionButton(ns("btn_supprimer_ligne"), tagList(icon("minus"), " Supprimer la/les ligne(s) sélectionnée(s)"),
                                         class = "btn-danger", style = "width:100%;")),
                  column(3, textInput(ns("nouvelle_colonne_nom"), NULL, placeholder = "Nom de la nouvelle colonne")),
                  column(3, actionButton(ns("btn_ajouter_colonne"), tagList(icon("plus"), " Ajouter une colonne"),
                                         class = "btn-default", style = "width:100%;"))
                ),
                fluidRow(
                  column(6, uiOutput(ns("colonne_a_supprimer_ui"))),
                  column(6, br(), actionButton(ns("btn_supprimer_colonne"), tagList(icon("minus"), " Supprimer cette colonne"),
                                               class = "btn-danger"))
                ),
                p("Astuce : cliquez sur une ou plusieurs lignes du tableau (Ctrl/Cmd + clic pour une sélection multiple) avant de les supprimer.",
                  style = "font-size:11px; color:#999;"),
                br(),
                div(style = "width:100%; overflow-x:auto;", DTOutput(ns("csv_preview"))),
                br(),
                uiOutput(ns("csv_save_status")))
          )
        )
      ),
      
      # ───────────────────────────────────────────────────────────────
      # ONGLET 2 — VÉRIFICATION DES CONTRÔLES
      # ───────────────────────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("check-double"), " Vérification des Contrôles"),
        
        fluidRow(
          column(
            width = 4,
            wellPanel(
              h4(tagList(icon("sliders-h"), " Seuils"), style = "margin-top:0; color:#605ca8;"),
              numericInput(ns("seuil_error"), "Seuil critique d'erreur (évènements min.)", value = 1000, min = 0, step = 100),
              numericInput(ns("seuil_warning"), "Seuil d'avertissement (évènements min.)", value = 5000, min = 0, step = 500),
              actionButton(ns("btn_verifier_asp"), tagList(icon("search"), " Lancer verifier_asp()"),
                           class = "btn-primary", style = "width:100%; font-weight:bold;")
            )
          ),
          column(
            width = 8,
            box(title = tagList(icon("clipboard-list"), " Résultat"), width = NULL, status = "warning", solidHeader = TRUE,
                verbatimTextOutput(ns("verif_result")))
          )
        )
      ),
      
      # ───────────────────────────────────────────────────────────────
      # ONGLET 3 — DÉFINITION DES GATES
      # ───────────────────────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("crosshairs"), " Définition des Gates"),
        
        fluidRow(
          column(
            width = 4,
            wellPanel(
              h4(tagList(icon("crosshairs"), " Nouvelle gate"), style = "margin-top:0; color:#605ca8;"),
              
              selectInput(ns("gate_methode"), "Méthode",
                          choices = c("Tune gate" = "tune", "Landmarks" = "landmarks", "Densité" = "density")),
              
              textInput(ns("gate_name"), "Nom de la gate (gate.name / control_name)", value = ""),
              numericInput(ns("gate_n_cells"), "n.cells", value = 2000, min = 1, step = 100),
              
              conditionalPanel(
                condition = sprintf("input['%s'] == 'tune' || input['%s'] == 'landmarks'", ns("gate_methode"), ns("gate_methode")),
                numericInput(ns("gate_percentile"), "percentile", value = 70, min = 1, max = 100, step = 1)
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'tune'", ns("gate_methode")),
                numericInput(ns("gate_bandwidth"), "bandwidth", value = 1, min = 0.1, step = 0.1)
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'landmarks' || input['%s'] == 'density'", ns("gate_methode"), ns("gate_methode")),
                numericInput(ns("gate_grid_n"), "grid.n", value = 100, min = 10, step = 10),
                numericInput(ns("gate_bandwidth_factor"), "bandwidth.factor", value = 1, min = 0.1, step = 0.1),
                textInput(ns("gate_fsc_channel"), "fsc.channel (optionnel)", value = ""),
                textInput(ns("gate_ssc_channel"), "ssc.channel (optionnel)", value = "")
              ),
              
              actionButton(ns("btn_definir_gate"), tagList(icon("play"), " Créer la gate"),
                           class = "btn-success", style = "width:100%; font-weight:bold;"),
              br(), br(), uiOutput(ns("gate_status")),
              
              hr(),
              h4(tagList(icon("broom"), " Nettoyage des contrôles"), style = "color:#605ca8;"),
              p("Appelle charger_et_nettoyer() une fois toutes les gates nécessaires définies.", style = "font-size:12px; color:#666;"),
              actionButton(ns("btn_charger_nettoyer"), tagList(icon("play"), " Lancer charger_et_nettoyer()"),
                           class = "btn-primary", style = "width:100%; font-weight:bold;"),
              br(), br(), uiOutput(ns("nettoyage_status"))
            )
          ),
          column(
            width = 8,
            box(title = tagList(icon("image"), " Aperçu des figures de gating"), width = NULL, status = "primary", solidHeader = TRUE,
                fluidRow(
                  column(6, selectInput(ns("gate_figure_dossier"), "Dossier",
                                        choices = c("figure_gate", "figure_gate_tuning"))),
                  column(6, uiOutput(ns("gate_figure_fichier_ui")))
                ),
                actionButton(ns("btn_rafraichir_figures"), tagList(icon("rotate"), " Rafraîchir la liste"), class = "btn-default"),
                br(), br(),
                uiOutput(ns("gate_figure_zone")),
                
                hr(),
                h4(tagList(icon("trash-can"), " Réinitialisation des gates"), style = "color:#c9302c;"),
                
                fluidRow(
                  column(6,
                         h5(tagList(icon("list-check"), " Suppression sélective")),
                         p("Coche une ou plusieurs gates ci-dessous puis supprime uniquement celles-ci (les autres restent intactes).",
                           style = "font-size:12px; color:#666;"),
                         uiOutput(ns("gates_liste_ui")),
                         actionButton(ns("btn_supprimer_gates_selection"), tagList(icon("trash"), " Supprimer la/les gate(s) cochée(s)"),
                                      class = "btn-danger", style = "width:100%;"),
                         br(), br(), uiOutput(ns("suppr_gates_status"))
                  ),
                  column(6,
                         h5(tagList(icon("triangle-exclamation"), " Réinitialisation complète")),
                         p("Supprime les figures de gating générées et/ou vide entièrement la liste des gates enregistrées dans le pipeline.",
                           style = "font-size:12px; color:#666;"),
                         checkboxGroupInput(ns("reset_dossiers"), "Dossiers de figures à supprimer",
                                            choices = c("figure_gate" = "figure_gate",
                                                        "figure_gate_tuning" = "figure_gate_tuning"),
                                            selected = c("figure_gate", "figure_gate_tuning")),
                         checkboxInput(ns("reset_vider_gates"), "Vider aussi la liste complète des gates (self$gates)", value = TRUE),
                         actionButton(ns("btn_reset_gates"), tagList(icon("trash-can"), " Tout réinitialiser"),
                                      class = "btn-danger", style = "width:100%; font-weight:bold;"),
                         br(), br(), uiOutput(ns("reset_gates_status"))
                  )
                ))
          )
        )
      ),
      
      # ───────────────────────────────────────────────────────────────
      # ONGLET 4 — EXTRACTION DES SPECTRES
      # ───────────────────────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("wave-square"), " Extraction des Spectres"),
        
        fluidRow(
          column(
            width = 4,
            wellPanel(
              h4(tagList(icon("wave-square"), " 1. Spectres de fluorophores"), style = "margin-top:0; color:#605ca8;"),
              actionButton(ns("btn_extract_spectra"), tagList(icon("play"), " Lancer extraire_fluorophore_spectre()"),
                           class = "btn-success", style = "width:100%; font-weight:bold;"),
              br(), br(), uiOutput(ns("spectra_status")),
              
              hr(),
              h4(tagList(icon("adjust"), " 2. Autofluorescence (optionnel)"), style = "color:#605ca8;"),
              uiOutput(ns("unstained_path_ui")),
              textInput(ns("af_tissue_name"), "tissue_name", value = "Cells"),
              checkboxInput(ns("af_refine"), "refine", value = TRUE),
              actionButton(ns("btn_extract_af"), tagList(icon("play"), " Lancer extraire_spectre_af()"),
                           class = "btn-primary", style = "width:100%; font-weight:bold;"),
              br(), br(), uiOutput(ns("af_status")),
              
              hr(),
              h4(tagList(icon("layer-group"), " 3. Variants spectraux (optionnel)"), style = "color:#605ca8;"),
              uiOutput(ns("variants_tissue_ui")),
              checkboxInput(ns("variants_refine"), "refine", value = TRUE),
              actionButton(ns("btn_variants"), tagList(icon("play"), " Lancer preparer_variants_spectraux()"),
                           class = "btn-primary", style = "width:100%; font-weight:bold;"),
              br(), br(), uiOutput(ns("variants_status"))
            )
          ),
          column(
            width = 8,
            box(title = tagList(icon("image"), " Aperçu des figures de spectres de fluorophores"), width = NULL, status = "primary", solidHeader = TRUE,
                fluidRow(
                  column(6, selectInput(ns("spectra_figure_dossier"), "Dossier",
                                        choices = c("figure_spectra", "figure_spectral_variants"))),
                  column(6, uiOutput(ns("spectra_figure_fichier_ui")))
                ),
                actionButton(ns("btn_rafraichir_figures_spectra"), tagList(icon("rotate"), " Rafraîchir la liste"), class = "btn-default"),
                br(), br(),
                uiOutput(ns("spectra_figure_zone")))
          )
        )
      ),
      
      # ───────────────────────────────────────────────────────────────
      # ONGLET 5 — DÉMIXAGE
      # ───────────────────────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("layer-group"), " Démixage"),
        
        fluidRow(
          column(
            width = 4,
            wellPanel(
              h4(tagList(icon("layer-group"), " Lancer le démixage"), style = "margin-top:0; color:#605ca8;"),
              
              radioButtons(ns("unmix_cible"), "Cible",
                           choices = c("Un seul fichier (unmix_fcs)" = "fichier",
                                       "Tout le dossier (unmix_folder)" = "dossier")),
              
              conditionalPanel(
                condition = sprintf("input['%s'] == 'fichier'", ns("unmix_cible")),
                uiOutput(ns("unmix_fichier_picker"))
              ),
              
              uiOutput(ns("unmix_tissue_ui")),
              selectInput(ns("unmix_method"), "method", choices = c("AutoSpectral", "WLS", "OLS")),
              selectInput(ns("unmix_speed"), "speed", choices = c("slow", "fast")),
              
              actionButton(ns("btn_unmix"), tagList(icon("play"), " Démixer"),
                           class = "btn-success", style = "width:100%; font-weight:bold;"),
              br(), br(), uiOutput(ns("unmix_status")),
              
              hr(),
              h4(tagList(icon("folder-open"), " Charger les résultats en mémoire"), style = "color:#605ca8;"),
              textInput(ns("unmix_dossier_resultats"), "Sous-dossier des fichiers unmixés", value = "AutoSpectral_unmixed"),
              actionButton(ns("btn_charger_unmixes"), tagList(icon("play"), " Lancer charger_fcs_unmixes()"),
                           class = "btn-primary", style = "width:100%; font-weight:bold;"),
              br(), br(), uiOutput(ns("charger_unmixes_status"))
            )
          ),
          column(
            width = 8,
            box(title = tagList(icon("chart-scatter"), " visualiser_unmixing()"), width = NULL, status = "info", solidHeader = TRUE,
                uiOutput(ns("visu_fichier_ui")),
                fluidRow(
                  column(6, uiOutput(ns("visu_canal_x_ui"))),
                  column(6, uiOutput(ns("visu_canal_y_ui")))
                ),
                fluidRow(
                  column(6, numericInput(ns("visu_cofacteur"), "cofacteur", value = 150, min = 1, step = 10)),
                  column(6, numericInput(ns("visu_max_points"), "max_points", value = 10000, min = 100, step = 1000))
                ),
                actionButton(ns("btn_visualiser"), tagList(icon("eye"), " Afficher"), class = "btn-primary"),
                br(), br(),
                plotlyOutput(ns("plot_unmixing"), height = "450px"))
          )
        )
      ),
      
      # ───────────────────────────────────────────────────────────────
      # ONGLET 6 — CONTRÔLE QUALITÉ
      # ───────────────────────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("chart-area"), " Contrôle Qualité"),
        
        fluidRow(
          column(
            width = 4,
            wellPanel(
              h4(tagList(icon("chart-area"), " verifier_qualite_unmix()"), style = "margin-top:0; color:#605ca8;"),
              
              uiOutput(ns("qc_fluorophore_ui")),
              uiOutput(ns("qc_single_stained_ui")),
              uiOutput(ns("qc_unstained_ui")),
              selectInput(ns("qc_cytometer"), "cytometer",
                          choices = c("aurora", "a8", "s8", "a5se", "id7000", "mosaic", "opteon", "xenith", "minimal", "auroraNL")),
              checkboxInput(ns("qc_gate"), "gate", value = TRUE),
              
              actionButton(ns("btn_qc"), tagList(icon("play"), " Lancer verifier_qualite_unmix()"),
                           class = "btn-success", style = "width:100%; font-weight:bold;")
            )
          ),
          column(
            width = 8,
            box(title = tagList(icon("clipboard-list"), " Résultat"), width = NULL, status = "warning", solidHeader = TRUE,
                verbatimTextOutput(ns("qc_result")))
          )
        )
      )
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SERVEUR
# ══════════════════════════════════════════════════════════════════════════════

demixage_server <- function(id, pipeline, pipeline_version) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    volumes <- c(Home = path.expand("~"), shinyFiles::getVolumes()())
    
    shinyFiles::shinyDirChoose(input, "dir_root", roots = volumes, session = session)
    shinyFiles::shinyDirChoose(input, "dir_monomarques", roots = volumes, session = session)
    shinyFiles::shinyDirChoose(input, "dir_echantillons", roots = volumes, session = session)
    
    # Chemin du dossier des échantillons (nécessaire pour unmix_folder / unmix_fcs)
    dossier_echantillons_path <- reactiveVal(NULL)
    # Fichiers .fcs détectés dans le dossier des contrôles monomarqués/unstained
    fichiers_monomarques <- reactiveVal(character(0))
    # Contenu éditable du fichier fcs_control_file.csv
    donnees_csv <- reactiveVal(NULL)
    
    # Chemin complet du fichier de contrôle généré par lancer_asp()
    chemin_csv <- reactive({
      req(pipeline())
      req(pipeline()$dossier_racine)
      file.path(pipeline()$dossier_racine, "fcs_control_file.csv")
    })
    
    # Lecture robuste au séparateur : AutoSpectral écrit toujours en virgule,
    # mais si ce fichier a été ouvert/réenregistré dans Excel (notamment en
    # français), il peut revenir en point-virgule et casser la lecture
    # standard. On essaie les deux séparateurs et on garde celui qui donne
    # le plus de colonnes.
    lire_csv_robuste <- function(chemin) {
      if (!file.exists(chemin)) {
        stop("Fichier introuvable : ", chemin)
      }
      d_virgule      <- tryCatch(read.csv(chemin, sep = ",", stringsAsFactors = FALSE, check.names = FALSE),
                                 error = function(e) NULL)
      d_pointvirgule <- tryCatch(read.csv(chemin, sep = ";", stringsAsFactors = FALSE, check.names = FALSE),
                                 error = function(e) NULL)
      
      n_virgule      <- if (is.null(d_virgule)) 0 else ncol(d_virgule)
      n_pointvirgule <- if (is.null(d_pointvirgule)) 0 else ncol(d_pointvirgule)
      
      if (n_virgule <= 1 && n_pointvirgule <= 1) {
        stop("Impossible de déterminer le séparateur du fichier CSV (aucune colonne détectée avec ',' ou ';').")
      }
      
      if (n_pointvirgule > n_virgule) d_pointvirgule else d_virgule
    }
    
    # Affiche une barre de progression ET désactive/anime le bouton pendant
    # l'exécution d'un traitement long. Les fonctions AutoSpectral appelées ici
    # sont synchrones et ne renvoient aucune information d'avancement réelle :
    # l'indicateur signale donc qu'un traitement est en cours (et son message),
    # pas un pourcentage exact. Le changement du bouton (désactivé + spinner)
    # est le signal le plus visible et le plus fiable, la barre de progression
    # (style "old", pleine largeur en haut de page) vient en complément.
    executer_avec_progression <- function(message, fn, bouton_id = NULL, label_repos = NULL, icone_repos = "play") {
      if (!is.null(bouton_id)) {
        shinyjs::disable(bouton_id)
        shinyjs::html(bouton_id, as.character(tagList(icon("spinner", class = "fa-spin"), paste0(" ", message))))
      }
      on.exit({
        if (!is.null(bouton_id)) {
          shinyjs::enable(bouton_id)
          if (!is.null(label_repos)) {
            shinyjs::html(bouton_id, as.character(tagList(icon(icone_repos), paste0(" ", label_repos))))
          }
        }
      }, add = TRUE)
      
      withProgress(message = message, value = 0.2, style = "notification", {
        resultat <- fn()
        setProgress(1)
        resultat
      })
    }
    
    # ── Dossier racine ──────────────────────────────────────────────────────
    observeEvent(input$dir_root, {
      req(pipeline())
      if (is.integer(input$dir_root)) return(invisible(NULL))
      chemin_choisi <- shinyFiles::parseDirPath(volumes, input$dir_root)
      req(length(chemin_choisi) > 0)
      if (!dir.exists(chemin_choisi)) {
        showNotification("Le dossier racine sélectionné n'existe pas.", type = "error")
        return(invisible(NULL))
      }
      p <- pipeline()
      p$dossier_racine <- chemin_choisi
      pipeline_version(pipeline_version() + 1L)
    })
    
    output$root_display <- renderUI({
      pipeline_version()
      chemin_actuel <- tryCatch(pipeline()$dossier_racine, error = function(e) NULL)
      if (is.null(chemin_actuel)) {
        div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
            icon("exclamation-triangle"), " Aucun dossier racine sélectionné.")
      } else {
        div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
            icon("check-circle"), " ", tags$b(chemin_actuel))
      }
    })
    
    # ── Dossier des contrôles monomarqués / unstained ──────────────────────
    observeEvent(input$dir_monomarques, {
      req(pipeline())
      if (is.integer(input$dir_monomarques)) return(invisible(NULL))
      chemin_choisi <- shinyFiles::parseDirPath(volumes, input$dir_monomarques)
      req(length(chemin_choisi) > 0)
      if (!dir.exists(chemin_choisi)) {
        showNotification("Le dossier des contrôles sélectionné n'existe pas.", type = "error")
        return(invisible(NULL))
      }
      
      fichiers <- list.files(chemin_choisi, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
      fichiers_monomarques(fichiers)
      
      # Par défaut : aucun fichier "Unstained" tant que l'utilisateur ne l'a pas désigné
      p <- pipeline()
      p$chemins_monomarques <- data.frame(
        chemin = fichiers,
        type   = "Monomarque",
        canal  = NA_character_,
        stringsAsFactors = FALSE
      )
      pipeline_version(pipeline_version() + 1L)
    })
    
    output$monomarques_display <- renderUI({
      n <- length(fichiers_monomarques())
      if (n == 0) {
        div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
            icon("exclamation-triangle"), " Aucun dossier de contrôles sélectionné.")
      } else {
        div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
            icon("check-circle"), " ", n, " fichier(s) .fcs détecté(s).")
      }
    })
    
    # Désignation du fichier "Unstained" parmi les fichiers détectés
    output$unstained_picker <- renderUI({
      fichiers <- fichiers_monomarques()
      req(length(fichiers) > 0)
      selectInput(ns("fichier_unstained"), "Fichier Unstained / négatif",
                  choices = basename(fichiers), selected = basename(fichiers)[1])
    })
    
    observeEvent(input$fichier_unstained, {
      req(pipeline())
      fichiers <- fichiers_monomarques()
      req(length(fichiers) > 0)
      p <- pipeline()
      p$chemins_monomarques <- data.frame(
        chemin = fichiers,
        type   = ifelse(basename(fichiers) == input$fichier_unstained, "Unstained", "Monomarque"),
        canal  = NA_character_,
        stringsAsFactors = FALSE
      )
      pipeline_version(pipeline_version() + 1L)
    })
    
    # Chemin absolu du fichier Unstained désigné (utilisé onglet Extraction)
    chemin_unstained <- reactive({
      fichiers <- fichiers_monomarques()
      req(length(fichiers) > 0, input$fichier_unstained)
      fichiers[basename(fichiers) == input$fichier_unstained][1]
    })
    
    # ── Dossier des échantillons biologiques ───────────────────────────────
    observeEvent(input$dir_echantillons, {
      req(pipeline())
      if (is.integer(input$dir_echantillons)) return(invisible(NULL))
      chemin_choisi <- shinyFiles::parseDirPath(volumes, input$dir_echantillons)
      req(length(chemin_choisi) > 0)
      if (!dir.exists(chemin_choisi)) {
        showNotification("Le dossier des échantillons sélectionné n'existe pas.", type = "error")
        return(invisible(NULL))
      }
      dossier_echantillons_path(chemin_choisi)
      pipeline_version(pipeline_version() + 1L)
    })
    
    output$echantillons_display <- renderUI({
      chemin_actuel <- dossier_echantillons_path()
      if (is.null(chemin_actuel)) {
        div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
            icon("exclamation-triangle"), " Aucun dossier d'échantillons sélectionné.")
      } else {
        div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
            icon("check-circle"), " ", tags$b(chemin_actuel))
      }
    })
    
    # ── 1. lancer_asp() ─────────────────────────────────────────────────────
    observeEvent(input$btn_lancer_asp, {
      req(pipeline())
      p <- pipeline()
      tryCatch({
        executer_avec_progression("Génération du fichier de contrôle AutoSpectral…", function() {
          p$lancer_asp(type_cytometre = input$type_cytometre)
        }, bouton_id = ns("btn_lancer_asp"), label_repos = "Lancer lancer_asp()")
        pipeline_version(pipeline_version() + 1L)
        
        # Chargement immédiat du fichier généré dans le tableau éditable
        donnees_csv(lire_csv_robuste(chemin_csv()))
        
        output$asp_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Fichier de contrôle AutoSpectral généré avec succès.")
        })
      }, error = function(e) {
        output$asp_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # ── Aperçu / édition / sauvegarde du fichier de contrôle ─────────────────
    
    observeEvent(input$btn_recharger_csv, {
      tryCatch({
        donnees_csv(lire_csv_robuste(chemin_csv()))
        output$csv_save_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Fichier rechargé depuis le disque.")
        })
      }, error = function(e) {
        output$csv_save_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    output$csv_preview <- renderDT({
      req(donnees_csv())
      datatable(
        donnees_csv(),
        editable = TRUE,
        rownames = FALSE,
        selection = "multiple",
        width = "100%",
        options = list(pageLength = 15, scrollX = TRUE, autoWidth = FALSE,
                       columnDefs = list(list(width = "auto", targets = "_all")))
      )
    })
    
    observeEvent(input$csv_preview_cell_edit, {
      info <- input$csv_preview_cell_edit
      df <- donnees_csv()
      req(df)
      df[info$row, info$col + 1] <- DT::coerceValue(info$value, df[info$row, info$col + 1])
      donnees_csv(df)
    })
    
    # ── Ajout / suppression de lignes ────────────────────────────────────────
    observeEvent(input$btn_ajouter_ligne, {
      df <- donnees_csv()
      req(df)
      nouvelle_ligne <- df[1, , drop = FALSE]
      nouvelle_ligne[1, ] <- NA
      donnees_csv(rbind(df, nouvelle_ligne))
      output$csv_save_status <- renderUI({
        div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
            icon("check-circle"), " Ligne vide ajoutée en bas du tableau.")
      })
    })
    
    observeEvent(input$btn_supprimer_ligne, {
      df <- donnees_csv()
      req(df)
      if (is.null(input$csv_preview_rows_selected) || length(input$csv_preview_rows_selected) == 0) {
        output$csv_save_status <- renderUI({
          div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
              icon("exclamation-triangle"), " Sélectionnez d'abord une ou plusieurs lignes dans le tableau (clic sur la/les ligne(s), touche Ctrl/Cmd pour en sélectionner plusieurs).")
        })
        return(invisible(NULL))
      }
      n_supprimees <- length(input$csv_preview_rows_selected)
      donnees_csv(df[-input$csv_preview_rows_selected, , drop = FALSE])
      output$csv_save_status <- renderUI({
        div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
            icon("check-circle"), " ", n_supprimees, " ligne(s) supprimée(s).")
      })
    })
    
    # ── Ajout / suppression de colonnes ──────────────────────────────────────
    output$colonne_a_supprimer_ui <- renderUI({
      df <- donnees_csv()
      req(df)
      selectInput(ns("colonne_a_supprimer"), "Colonne à supprimer", choices = colnames(df))
    })
    
    observeEvent(input$btn_ajouter_colonne, {
      df <- donnees_csv()
      req(df)
      nom <- trimws(input$nouvelle_colonne_nom)
      if (!nzchar(nom)) {
        output$csv_save_status <- renderUI({
          div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
              icon("exclamation-triangle"), " Indiquez un nom de colonne avant de l'ajouter.")
        })
        return(invisible(NULL))
      }
      if (nom %in% colnames(df)) {
        output$csv_save_status <- renderUI({
          div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
              icon("exclamation-triangle"), " Une colonne \"", nom, "\" existe déjà.")
        })
        return(invisible(NULL))
      }
      df[[nom]] <- NA_character_
      donnees_csv(df)
      updateTextInput(session, "nouvelle_colonne_nom", value = "")
      output$csv_save_status <- renderUI({
        div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
            icon("check-circle"), " Colonne \"", nom, "\" ajoutée.")
      })
    })
    
    observeEvent(input$btn_supprimer_colonne, {
      df <- donnees_csv()
      req(df, input$colonne_a_supprimer)
      if (ncol(df) <= 1) {
        output$csv_save_status <- renderUI({
          div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
              icon("exclamation-triangle"), " Impossible de supprimer la dernière colonne restante.")
        })
        return(invisible(NULL))
      }
      df[[input$colonne_a_supprimer]] <- NULL
      donnees_csv(df)
      output$csv_save_status <- renderUI({
        div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
            icon("check-circle"), " Colonne \"", input$colonne_a_supprimer, "\" supprimée.")
      })
    })
    
    observeEvent(input$btn_save_csv, {
      req(donnees_csv())
      tryCatch({
        # write.csv() fige toujours sep="," et dec="." quelle que soit la
        # locale du système : le fichier réécrit ici sera donc toujours dans
        # le format attendu par AutoSpectral, sans repasser par Excel.
        executer_avec_progression("Enregistrement du fichier de contrôle…", function() {
          write.csv(donnees_csv(), file = chemin_csv(), row.names = FALSE, quote = TRUE)
        }, bouton_id = ns("btn_save_csv"), label_repos = "Enregistrer les modifications", icone_repos = "save")
        output$csv_save_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Modifications enregistrées dans fcs_control_file.csv.")
        })
      }, error = function(e) {
        output$csv_save_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur lors de l'enregistrement : ", e$message)
        })
      })
    })
    
    # ── 2. verifier_asp() ────────────────────────────────────────────────────
    observeEvent(input$btn_verifier_asp, {
      req(pipeline())
      p <- pipeline()
      output$verif_result <- renderPrint({
        tryCatch({
          executer_avec_progression("Vérification des contrôles…", function() {
            p$verifier_asp(warning = input$seuil_warning, error = input$seuil_error)
          }, bouton_id = ns("btn_verifier_asp"), label_repos = "Lancer verifier_asp()", icone_repos = "search")
        }, error = function(e) {
          cat("Erreur :", e$message, "\n")
        })
      })
    })
    
    # ── 3. Définition des gates ──────────────────────────────────────────────
    observeEvent(input$btn_definir_gate, {
      req(pipeline())
      p <- pipeline()
      req(nzchar(input$gate_name))
      
      tryCatch({
        fsc <- if (nzchar(input$gate_fsc_channel)) input$gate_fsc_channel else NULL
        ssc <- if (nzchar(input$gate_ssc_channel)) input$gate_ssc_channel else NULL
        
        executer_avec_progression(paste0("Création de la gate \"", input$gate_name, "\"…"), function() {
          if (input$gate_methode == "tune") {
            p$definir_tune_gates(
              gate.name  = input$gate_name,
              n_cells    = input$gate_n_cells,
              percentile = input$gate_percentile,
              bandwidth  = input$gate_bandwidth
            )
          } else if (input$gate_methode == "landmarks") {
            p$definir_gates_landmarks(
              control_name     = input$gate_name,
              n.cells          = input$gate_n_cells,
              percentile       = input$gate_percentile,
              grid.n           = input$gate_grid_n,
              bandwidth.factor = input$gate_bandwidth_factor,
              fsc.channel      = fsc,
              ssc.channel      = ssc
            )
          } else if (input$gate_methode == "density") {
            p$definir_gates_density(
              control_name     = input$gate_name,
              n.cells          = input$gate_n_cells,
              grid.n           = input$gate_grid_n,
              bandwidth.factor = input$gate_bandwidth_factor,
              fsc.channel      = fsc,
              ssc.channel      = ssc
            )
          }
        }, bouton_id = ns("btn_definir_gate"), label_repos = "Créer la gate")
        
        pipeline_version(pipeline_version() + 1L)
        output$gate_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Gate \"", input$gate_name, "\" créée avec succès.")
        })
      }, error = function(e) {
        output$gate_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    output$gates_liste_ui <- renderUI({
      pipeline_version()
      noms <- tryCatch(names(pipeline()$gates), error = function(e) NULL)
      if (is.null(noms) || length(noms) == 0) {
        return(p("Aucune gate définie pour le moment.", style = "font-size:12px; color:#999;"))
      }
      checkboxGroupInput(ns("gates_selectionnees"), NULL, choices = noms)
    })
    
    observeEvent(input$btn_supprimer_gates_selection, {
      req(pipeline())
      p <- pipeline()
      if (is.null(input$gates_selectionnees) || length(input$gates_selectionnees) == 0) {
        output$suppr_gates_status <- renderUI({
          div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
              icon("exclamation-triangle"), " Cochez au moins une gate à supprimer.")
        })
        return(invisible(NULL))
      }
      tryCatch({
        for (nom in input$gates_selectionnees) {
          p$gates[[nom]] <- NULL
        }
        pipeline_version(pipeline_version() + 1L)
        output$suppr_gates_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Gate(s) supprimée(s) : ", paste(input$gates_selectionnees, collapse = ", "), ".")
        })
      }, error = function(e) {
        output$suppr_gates_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # ── charger_et_nettoyer() ────────────────────────────────────────────────
    observeEvent(input$btn_charger_nettoyer, {
      req(pipeline())
      p <- pipeline()
      tryCatch({
        executer_avec_progression("Chargement et nettoyage des contrôles…", function() {
          p$charger_et_nettoyer()
        }, bouton_id = ns("btn_charger_nettoyer"), label_repos = "Lancer charger_et_nettoyer()")
        pipeline_version(pipeline_version() + 1L)
        output$nettoyage_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Contrôles chargés et nettoyés avec succès.")
        })
      }, error = function(e) {
        output$nettoyage_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # ── Réinitialisation des gates ───────────────────────────────────────────
    observeEvent(input$btn_reset_gates, {
      req(length(input$reset_dossiers) > 0 || isTRUE(input$reset_vider_gates))
      showModal(modalDialog(
        title = tagList(icon("triangle-exclamation"), " Confirmer la réinitialisation"),
        "Cette action va supprimer définitivement les dossiers de figures sélectionnés",
        if (isTRUE(input$reset_vider_gates)) " et vider la liste des gates enregistrées" else NULL,
        ". Cette action est irréversible. Continuer ?",
        footer = tagList(
          modalButton("Annuler"),
          actionButton(ns("btn_confirm_reset_gates"), "Confirmer la suppression", class = "btn-danger")
        )
      ))
    })
    
    observeEvent(input$btn_confirm_reset_gates, {
      removeModal()
      req(pipeline())
      p <- pipeline()
      tryCatch({
        req(p$dossier_racine)
        actions <- character(0)
        
        for (nom_dossier in input$reset_dossiers) {
          chemin_dossier <- file.path(p$dossier_racine, nom_dossier)
          if (dir.exists(chemin_dossier)) {
            unlink(chemin_dossier, recursive = TRUE, force = TRUE)
            actions <- c(actions, paste0("dossier \"", nom_dossier, "\" supprimé"))
          } else {
            actions <- c(actions, paste0("dossier \"", nom_dossier, "\" introuvable (rien à supprimer)"))
          }
        }
        
        if (isTRUE(input$reset_vider_gates)) {
          p$gates <- list()
          actions <- c(actions, "liste des gates (self$gates) vidée")
        }
        
        pipeline_version(pipeline_version() + 1L)
        output$reset_gates_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " ", paste(actions, collapse = " ; "), ".")
        })
      }, error = function(e) {
        output$reset_gates_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # ── Aperçu des figures de gating (.jpg) ──────────────────────────────────
    
    # Se rafraîchit : au changement de dossier choisi, après une réinitialisation,
    # après la création d'une nouvelle gate, ou via le bouton "Rafraîchir la liste"
    fichiers_figures_gate <- reactive({
      input$btn_rafraichir_figures
      pipeline_version()
      if (is.null(pipeline()) || is.null(pipeline()$dossier_racine) || is.null(input$gate_figure_dossier)) {
        return(character(0))
      }
      chemin_dossier <- file.path(pipeline()$dossier_racine, input$gate_figure_dossier)
      if (!dir.exists(chemin_dossier)) return(character(0))
      list.files(chemin_dossier, pattern = "\\.(jpg|jpeg)$", full.names = TRUE, ignore.case = TRUE)
    })
    
    output$gate_figure_fichier_ui <- renderUI({
      fichiers <- fichiers_figures_gate()
      if (length(fichiers) == 0) {
        return(p("Aucune figure .jpg trouvée dans ce dossier.", style = "font-size:12px; color:#999; margin-top:25px;"))
      }
      selectInput(ns("gate_figure_fichier"), "Figure", choices = setNames(fichiers, basename(fichiers)))
    })
    
    output$gate_figure_zone <- renderUI({
      fichiers <- fichiers_figures_gate()
      if (length(fichiers) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucune figure disponible pour le moment dans \"",
                   input$gate_figure_dossier, "\"."))
      }
      req(input$gate_figure_fichier)
      imageOutput(ns("gate_figure_img"), height = "auto")
    })
    
    output$gate_figure_img <- renderImage({
      req(input$gate_figure_fichier)
      list(src = input$gate_figure_fichier, contentType = "image/jpeg", width = "100%", alt = basename(input$gate_figure_fichier))
    }, deleteFile = FALSE)
    
    # ── Aperçu des figures d'extraction des spectres (.jpg / .pdf) ───────────
    fichiers_figures_spectra <- reactive({
      input$btn_rafraichir_figures_spectra
      pipeline_version()
      if (is.null(pipeline()) || is.null(pipeline()$dossier_racine) || is.null(input$spectra_figure_dossier)) {
        return(character(0))
      }
      chemin_dossier <- file.path(pipeline()$dossier_racine, input$spectra_figure_dossier)
      if (!dir.exists(chemin_dossier)) return(character(0))
      list.files(chemin_dossier, pattern = "\\.(jpg|jpeg|pdf)$", full.names = TRUE, ignore.case = TRUE)
    })
    
    output$spectra_figure_fichier_ui <- renderUI({
      fichiers <- fichiers_figures_spectra()
      if (length(fichiers) == 0) {
        return(p("Aucune figure trouvée dans ce dossier.", style = "font-size:12px; color:#999; margin-top:25px;"))
      }
      selectInput(ns("spectra_figure_fichier"), "Figure", choices = setNames(fichiers, basename(fichiers)))
    })
    
    output$spectra_figure_zone <- renderUI({
      fichiers <- fichiers_figures_spectra()
      if (length(fichiers) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucune figure disponible pour le moment dans \"",
                   input$spectra_figure_dossier, "\"."))
      }
      req(input$spectra_figure_fichier)
      if (grepl("\\.pdf$", input$spectra_figure_fichier, ignore.case = TRUE)) {
        # Les PDF ne peuvent pas être prévisualisés inline simplement : on propose un téléchargement.
        div(
          div(class = "alert alert-info", style = "font-size:12px; padding:8px;",
              icon("file-pdf"), " Fichier PDF — la prévisualisation inline n'est pas disponible, téléchargez-le pour l'ouvrir."),
          downloadButton(ns("dl_spectra_figure_pdf"), "Télécharger le PDF", class = "btn-primary")
        )
      } else {
        imageOutput(ns("spectra_figure_img"), height = "auto")
      }
    })
    
    output$spectra_figure_img <- renderImage({
      req(input$spectra_figure_fichier)
      req(!grepl("\\.pdf$", input$spectra_figure_fichier, ignore.case = TRUE))
      list(src = input$spectra_figure_fichier, contentType = "image/jpeg", width = "100%", alt = basename(input$spectra_figure_fichier))
    }, deleteFile = FALSE)
    
    output$dl_spectra_figure_pdf <- downloadHandler(
      filename = function() basename(req(input$spectra_figure_fichier)),
      content = function(file) file.copy(input$spectra_figure_fichier, file)
    )
    
    # ── 4. extraire_fluorophore_spectre() ────────────────────────────────────
    observeEvent(input$btn_extract_spectra, {
      req(pipeline())
      p <- pipeline()
      tryCatch({
        executer_avec_progression("Extraction des spectres de fluorophores…", function() {
          p$extraire_fluorophore_spectre()
        }, bouton_id = ns("btn_extract_spectra"), label_repos = "Lancer extraire_fluorophore_spectre()")
        pipeline_version(pipeline_version() + 1L)
        output$spectra_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Spectres de fluorophores extraits avec succès.")
        })
      }, error = function(e) {
        output$spectra_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # Chemin du fichier Unstained, pré-rempli mais modifiable
    output$unstained_path_ui <- renderUI({
      valeur <- tryCatch(chemin_unstained(), error = function(e) "")
      textInput(ns("af_unstained_path"), "unstained_fcs_path", value = valeur)
    })
    
    # ── extraire_spectre_af() ────────────────────────────────────────────────
    observeEvent(input$btn_extract_af, {
      req(pipeline())
      p <- pipeline()
      req(nzchar(input$af_unstained_path), nzchar(input$af_tissue_name))
      tryCatch({
        executer_avec_progression("Extraction de l'autofluorescence…", function() {
          p$extraire_spectre_af(
            unstained_fcs_path = input$af_unstained_path,
            tissue_name         = input$af_tissue_name,
            refine               = input$af_refine
          )
        }, bouton_id = ns("btn_extract_af"), label_repos = "Lancer extraire_spectre_af()")
        pipeline_version(pipeline_version() + 1L)
        output$af_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Autofluorescence extraite pour \"", input$af_tissue_name, "\".")
        })
      }, error = function(e) {
        output$af_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # Sélection du tissu pour les variants spectraux (basé sur les AF déjà extraites)
    output$variants_tissue_ui <- renderUI({
      pipeline_version()
      noms <- tryCatch(names(pipeline()$af_spectra), error = function(e) NULL)
      req(length(noms) > 0)
      selectInput(ns("variants_tissue_name"), "tissue_af_name", choices = noms)
    })
    
    # ── preparer_variants_spectraux() ────────────────────────────────────────
    observeEvent(input$btn_variants, {
      req(pipeline())
      p <- pipeline()
      tryCatch({
        executer_avec_progression("Préparation des variants spectraux…", function() {
          p$preparer_variants_spectraux(
            tissue_af_name = input$variants_tissue_name,
            refine          = input$variants_refine
          )
        }, bouton_id = ns("btn_variants"), label_repos = "Lancer preparer_variants_spectraux()")
        pipeline_version(pipeline_version() + 1L)
        output$variants_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Variants spectraux préparés avec succès.")
        })
      }, error = function(e) {
        output$variants_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # ── 5. Démixage ───────────────────────────────────────────────────────────
    
    # tissue_name pour unmix_fcs/unmix_folder : basé sur les AF déjà extraites
    output$unmix_tissue_ui <- renderUI({
      pipeline_version()
      noms <- tryCatch(names(pipeline()$af_spectra), error = function(e) NULL)
      if (is.null(noms) || length(noms) == 0) {
        textInput(ns("unmix_tissue_name"), "tissue_name (optionnel)", value = "")
      } else {
        selectInput(ns("unmix_tissue_name"), "tissue_name", choices = c("(aucun)", noms))
      }
    })
    
    output$unmix_fichier_picker <- renderUI({
      dossier <- dossier_echantillons_path()
      req(!is.null(dossier))
      fichiers <- list.files(dossier, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
      req(length(fichiers) > 0)
      selectInput(ns("unmix_fichier_choisi"), "Fichier à démixer", choices = setNames(fichiers, basename(fichiers)))
    })
    
    observeEvent(input$btn_unmix, {
      req(pipeline())
      p <- pipeline()
      
      tissue <- if (is.null(input$unmix_tissue_name) || input$unmix_tissue_name %in% c("", "(aucun)")) {
        NULL
      } else {
        input$unmix_tissue_name
      }
      
      tryCatch({
        executer_avec_progression(
          if (input$unmix_cible == "fichier") "Démixage du fichier en cours…" else "Démixage du dossier en cours (peut prendre du temps)…",
          function() {
            if (input$unmix_cible == "fichier") {
              req(input$unmix_fichier_choisi)
              p$unmix_fcs(
                fcs_file_path = input$unmix_fichier_choisi,
                tissue_name    = tissue,
                method          = input$unmix_method,
                speed           = input$unmix_speed
              )
            } else {
              dossier <- dossier_echantillons_path()
              req(!is.null(dossier))
              p$unmix_folder(
                folder_path = dossier,
                tissue_name  = tissue,
                method        = input$unmix_method,
                speed         = input$unmix_speed
              )
            }
          },
          bouton_id = ns("btn_unmix"), label_repos = "Démixer"
        )
        
        pipeline_version(pipeline_version() + 1L)
        output$unmix_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Démixage (", input$unmix_method, ") terminé avec succès.")
        })
      }, error = function(e) {
        output$unmix_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    observeEvent(input$btn_charger_unmixes, {
      req(pipeline())
      p <- pipeline()
      req(nzchar(input$unmix_dossier_resultats))
      tryCatch({
        executer_avec_progression("Chargement des fichiers unmixés en mémoire…", function() {
          p$charger_fcs_unmixes(dossier = input$unmix_dossier_resultats)
        }, bouton_id = ns("btn_charger_unmixes"), label_repos = "Lancer charger_fcs_unmixes()")
        pipeline_version(pipeline_version() + 1L)
        output$charger_unmixes_status <- renderUI({
          div(class = "alert alert-success", style = "font-size:12px; padding:8px;",
              icon("check-circle"), " Fichiers unmixés chargés en mémoire.")
        })
      }, error = function(e) {
        output$charger_unmixes_status <- renderUI({
          div(class = "alert alert-danger", style = "font-size:12px; padding:8px;",
              icon("times-circle"), " Erreur : ", e$message)
        })
      })
    })
    
    # ── visualiser_unmixing() ────────────────────────────────────────────────
    output$visu_fichier_ui <- renderUI({
      pipeline_version()
      noms <- tryCatch(names(pipeline()$echantillons_traites), error = function(e) NULL)
      if (is.null(noms) || length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"),
                   " Aucun échantillon en mémoire. Démixez un fichier/dossier puis cliquez sur \"charger_fcs_unmixes()\" ci-contre."))
      }
      selectInput(ns("visu_fichier"), "nom_fichier_fcs", choices = noms)
    })
    
    canaux_disponibles <- reactive({
      req(input$visu_fichier)
      fcs <- tryCatch(pipeline()$echantillons_traites[[input$visu_fichier]], error = function(e) NULL)
      req(!is.null(fcs))
      colnames(fcs)
    })
    
    output$visu_canal_x_ui <- renderUI({
      ch <- tryCatch(canaux_disponibles(), error = function(e) NULL)
      if (is.null(ch)) {
        return(p("canal_x : sélectionnez d'abord un fichier ci-dessus.", style = "font-size:12px; color:#999;"))
      }
      selectInput(ns("visu_canal_x"), "canal_x", choices = ch)
    })
    
    output$visu_canal_y_ui <- renderUI({
      ch <- tryCatch(canaux_disponibles(), error = function(e) NULL)
      if (is.null(ch)) {
        return(p("canal_y : sélectionnez d'abord un fichier ci-dessus.", style = "font-size:12px; color:#999;"))
      }
      selectInput(ns("visu_canal_y"), "canal_y", choices = ch, selected = ch[min(2, length(ch))])
    })
    
    observeEvent(input$btn_visualiser, {
      req(pipeline(), input$visu_fichier, input$visu_canal_x, input$visu_canal_y)
      output$plot_unmixing <- renderPlotly({
        p_gg <- tryCatch({
          executer_avec_progression("Génération du graphique…", function() {
            pipeline()$visualiser_unmixing(
              nom_fichier_fcs = input$visu_fichier,
              canal_x          = input$visu_canal_x,
              canal_y          = input$visu_canal_y,
              cofacteur         = input$visu_cofacteur,
              max_points        = input$visu_max_points
            )
          }, bouton_id = ns("btn_visualiser"), label_repos = "Afficher", icone_repos = "eye")
        }, error = function(e) {
          showNotification(paste("Erreur :", e$message), type = "error")
          NULL
        })
        req(p_gg)
        plotly::ggplotly(p_gg)
      })
    })
    
    # ── 6. verifier_qualite_unmix() ──────────────────────────────────────────
    output$qc_fluorophore_ui <- renderUI({
      pipeline_version()
      noms <- tryCatch(rownames(pipeline()$spectra), error = function(e) NULL)
      req(length(noms) > 0)
      selectInput(ns("qc_fluorophore"), "fluorophore", choices = noms)
    })
    
    output$qc_single_stained_ui <- renderUI({
      fichiers <- fichiers_monomarques()
      req(length(fichiers) > 0)
      selectInput(ns("qc_single_stained"), "single_stained_fcs",
                  choices = setNames(fichiers, basename(fichiers)))
    })
    
    output$qc_unstained_ui <- renderUI({
      valeur <- tryCatch(chemin_unstained(), error = function(e) "")
      textInput(ns("qc_unstained"), "unstained_fcs", value = valeur)
    })
    
    observeEvent(input$btn_qc, {
      req(pipeline())
      p <- pipeline()
      output$qc_result <- renderPrint({
        tryCatch({
          executer_avec_progression("Contrôle qualité de l'unmixing…", function() {
            p$verifier_qualite_unmix(
              fluorophore         = input$qc_fluorophore,
              single_stained_fcs  = input$qc_single_stained,
              unstained_fcs        = input$qc_unstained,
              cytometer             = input$qc_cytometer,
              gate                   = input$qc_gate
            )
          }, bouton_id = ns("btn_qc"), label_repos = "Lancer verifier_qualite_unmix()")
        }, error = function(e) {
          cat("Erreur :", e$message, "\n")
        })
      })
    })
  })
}