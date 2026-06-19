library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(shinyjs)
library(grid)
library(gridExtra)

source("pipeline_cytometrie.R")
options(shiny.maxRequestSize = 1024 * 1024 * 1024)

# Base initiale du dictionnaire de canaux connus
if (!exists("global_canaux")) {
  global_canaux <- c("", "FITC-A", "Alexa Fluor 488-A", "Alexa Fluor 700-A", "PE-A", "PE-Cy5-A","PE-Cy7-A","PerCP-Cy5-5-A", "GFP-A", "PE-Cy7-A", "APC-A", "APC-R700-A", "APC-H7-A", 
                     "BV421-A", "AmCyan-A","BV510-A", "BV605-A", "BV650-A","Pacific Blue-A", "V450-A","BV711-A", "BV786-A")
}

####################### UI #########################

library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)

ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Cytométrie"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs_sidebar",
      menuItem("Importation", tabName = "import_tab", icon = icon("file-import")),
      menuItem("Compensation", tabName = "compensation_tab", icon = icon("calculator")),
      menuItem("Unmixing", tabName = "unmixing_tab", icon = icon("bolt")),
      menuItem("Quality Control", tabName = "nettoyage_tab", icon = icon("broom")),
      menuItem("Analyses", tabName = "analyses_tab", icon = icon("chart-pie"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .datatables, .dataTables_wrapper, table.dataTable tbody td { color: black !important; }
        .disabled-tab { pointer-events: none; opacity: 0.4; }
        .nav-tabs-custom > .nav-tabs > li.active { border-top-color: #605ca8 !important; }
        /* Style pour harmoniser les selects et inputs modifiables */
        .dt-select, .dt-datalist-input { width: 100%; padding: 4px; background: white; color: black; border: 1px solid #ccc; border-radius: 4px;}
      "))
    ),
    
    tabItems(
      
      # ============================================================
      # 📥 ONGLET IMPORTATION (NOUVELLE ARCHITECTURE À ONGLETS)
      # ============================================================
      tabItem(
        tabName = "import_tab",
        fluidRow(
          column(
            width = 12,
            tabBox(
              title = tagList(icon("file-import"), "Configuration de l'Importation"),
              id = "import_steps", width = 12,
              
              # ---- SOUS-ONGLET 1 : FICHIERS  ----
              tabPanel(
                title = "1. Fichiers",
                fluidRow(
                  column(
                    width = 4,
                    box(
                      title = "Chargement des fichiers FCS", width = 12, status = "primary", solidHeader = TRUE,
                      fileInput("load_mono", "1. Tubes Monomarqués & Unstained", multiple = TRUE, accept = ".fcs"),
                      fileInput("load_ech", "2. Échantillons Biologiques", multiple = TRUE, accept = ".fcs"),
                      hr(),
                      uiOutput("ui_cyto_select")
                    ),
                    # Bouton d'initialisation global déplacé ici pour verrouiller toute l'étape 1
                    actionButton(
                      "init_r6", "Valider la configuration & Initialiser", 
                      class = "btn-purple", style = "width:100%; color:white; background-color:#605ca8;"
                    )
                  ),
                  
                  column(
                    width = 8,
                    box(
                      title = "Annotation & Configuration des Tubes", width = 12, status = "warning", solidHeader = TRUE,
                      h4("Configuration des Monomarqués"),
                      DTOutput("table_mono"),
                      br(),
                      h4("Liste des Échantillons Biologiques"),
                      p(style = "color:gray; font-size:11px;", "Double-cliquez sur le nom pour le renommer."),
                      DTOutput("table_ech")
                    )
                  )
                )
              ),
              
              # ---- SOUS-ONGLET 2 : CONFIGURATION DES MARQUEURS ----
              tabPanel(
                title = "2. Configuration des Marqueurs",
                fluidRow(
                  box(
                    title = "Annotation", 
                    width = 12, status = "primary", solidHeader = TRUE,
                    p(style = "color:gray; font-size:11px;", 
                      "Les lignes représentent vos échantillons et les colonnes vos canaux détectés. Les cases sont pré-remplies automatiquement avec le nom du marqueur extrait du fichier FCS. Double-cliquez sur une case pour corriger manuellement un marqueur si nécessaire."),
                    br(),
                    div(style = "overflow-x: auto;",
                        DTOutput("table_matrice_marqueurs")
                    )
                  )
                )
              ),
              
              # ---- ✨ NOUVEAU : SOUS-ONGLET 3 : GROUPES D'ÉCHANTILLONS ----
              tabPanel(
                title = "3. Groupes d'Échantillons",
                fluidRow(
                  column(
                    width = 4,
                    box(
                      title = "Gestion des Groupes Biologiques", width = 12, status = "info", solidHeader = TRUE,
                      p(style = "color:gray; font-size:11px;", "Créez vos groupes (ex: Contrôle, Choc Septique) puis assignez-y vos tubes."),
                      fluidRow(
                        column(8, textInput("group_name", NULL, placeholder = "Nom du groupe...")),
                        column(4, actionButton("add_group", "Créer", class = "btn-info", style = "width:100%;"))
                      ),
                      hr(),
                      h5(tags$strong("Assignation des Échantillons :")),
                      div(style = "max-height: 300px; overflow-y: auto;",
                          uiOutput("ui_assignation_groupes")
                      )
                    )
                  ),
                  column(
                    width = 8,
                    box(
                      title = "Résumé de la Cohorte", width = 12, status = "primary", solidHeader = TRUE,
                      p(style = "color:gray; font-size:11px;", "Récapitulatif des groupes créés et des échantillons associés."),
                      DTOutput("table_resume_groupes"),
                      br(),
                      actionButton(
                        "save_cohort", "💾 Enregistrer la Cohorte", 
                        class = "btn-success", style = "width:100%; font-weight:bold; font-size:14px;"
                      )
                    )
                  )
                )
              )
            ) # Fin tabBox Importation
          ) # Fin column
        ) # Fin fluidRow
      ),
      
      # ============================================================
      # 🎛️ ONGLET COMPENSATION (INCHANGÉ)
      # ============================================================
      tabItem(
        tabName = "compensation_tab",
        fluidRow(
          column(
            width = 12,
            tabBox(
              title = tagList(icon("sliders-h"), "Outils de Compensation"),
              id = "comp_steps", width = 12,
              
              # ÉTAPE 1 : TRANSFORMATION
              tabPanel(
                title = "Transformation",
                fluidRow(
                  box(
                    title = "Choix des paramètres", width = 4, status = "primary", solidHeader = TRUE,
                    column(
                      width = 12,
                      numericInput("cofacteur", "Valeur du Cofacteur Arcsinh", value = 600, min = 1),
                      uiOutput("ui_trans_file"),
                      uiOutput("ui_trans_cx"),
                      uiOutput("ui_trans_cy"),
                      hr(),
                      actionButton(
                        "btn_apply_trans", "Appliquer la transformation", 
                        class = "btn-purple", style = "width:100%; color:white; background-color:#605ca8;"
                      )
                    )
                  ),
                  column(width = 8, plotOutput("plot_transformation", height = "500px"))
                )
              ),
              
              # ÉTAPE 2 : GATING
              tabPanel(
                title = "Gating (Sélection des populations)",
                fluidRow(
                  column(
                    width = 4,
                    uiOutput("ui_gate_canal"),
                    uiOutput("ui_gate_use_unstained_wrapper"),
                    hr(),
                    h5(style = "color:#0077b6; font-weight:bold;", "Intervalle Population NÉGATIVE"),
                    uiOutput("ui_slider_neg"),
                    br(),
                    h5(style = "color:#d90429; font-weight:bold;", "Intervalle Population POSITIVE"),
                    uiOutput("ui_slider_pos"),
                    hr(),
                    actionButton("save_gate", "Enregistrer les gates", class = "btn-success", style = "width:100%")
                  ),
                  column(width = 8, plotOutput("plot_gating", height = "520px"))
                )
              ),
              
              # ÉTAPE 3 : MATRICE
              tabPanel(
                title = "Matrice de spillover",
                value = "tab_matrice_spillover",
                
                fluidRow(
                  box(
                    title = "Consultation de la matrice de spillover", 
                    width = 12, 
                    status = "primary", 
                    solidHeader = TRUE,
                    
                    helpText("Sélectionnez un échantillon"),
                    br(),
                    
                    # Sélecteur d'échantillon dédié à cet onglet
                    uiOutput("ui_select_ech_matrice_origine"),
                    br(),
                    
                    # Tableau d'affichage de la matrice
                    div(
                      style = "overflow-x: auto;",
                      tableOutput("matrice_origine_affichage")
                    )
                  )
                )
              ),
              
              # ÉTAPE 4 : FIGURES
              tabPanel(
                title = "Figures & Validations",
                value = "tab_figures_validations",
                
                # PARTIE SUPÉRIEURE : Configuration des graphiques à gauche et affichage à droite
                fluidRow(
                  column(
                    width = 3,
                    uiOutput("ui_fig_ech"),
                    radioButtons(
                      inputId = "fig_mode", "Mode d'affichage :", 
                      choices = c("Paires de Canaux" = "pairs", "Sélection Libre" = "free")
                    ),
                    uiOutput("ui_conditionnel_axes"),
                    selectInput(
                      inputId = "fig_view_type", "Type de vue :", 
                      choices = c("Avant et après compensation", "Avant compensation uniquement", "Après compensation uniquement")
                    ),
                    numericInput("fig_max_points", "Nombre de points max :", value = 5000, min = 500, step = 500),
                    br(),
                    actionButton("run_compensation", "⚡ Lancer Compensation Globale", class = "btn-success", style = "width:100%")
                  ),
                  
                  column(
                    width = 9, 
                    div(
                      style = "overflow-y: auto; max-height: 800px; padding-right: 10px;",
                      plotOutput("plot_validation_comp", height = "auto")
                    )
                  )
                ),
                
                br(), 
                
                # PARTIE INFÉRIEURE : Édition de la matrice par Axes (Ligne/Colonne) pour l'échantillon sélectionné
                fluidRow(
                  column(
                    width = 12,
                    tabBox(
                      title = tagList(shiny::icon("sliders-h"), "Édition de la Compensation par Échantillon"), 
                      width = 12, 
                      id = "box_edition_comp",
                      
                      tabPanel(
                        title = "Ajustement ciblé de la matrice",
                        fluidRow(
                          # Bloc de modification par Axes
                          column(
                            width = 6,
                            div(
                              style = "background: #fdfdfd; padding: 15px; border: 1px solid #e3e3e3; border-radius: 4px;",
                              h4(style = "margin-top:0; color: #333;", "Modifier une valeur de Spillover"),
                              fluidRow(
                                column(4, uiOutput("ui_edit_axes_c1")), # Choix Ligne (Source)
                                column(4, uiOutput("ui_edit_axes_c2")), # Choix Colonne (Spillover vers)
                                column(4, numericInput("edit_axes_val", "Nouvelle valeur :", value = 0, min = 0, max = 1, step = 0.001))
                              ),
                              br(),
                              actionButton("valider_changement_axe", "Appliquer le changement", class = "btn-primary", style = "width: 100%;")
                            )
                          ),
                          
                          # Bloc d'Historique et de Sauvegarde
                          column(
                            width = 6,
                            div(
                              style = "padding: 15px;",
                              h4("Historique de ce tube"),
                              div(
                                style = "margin-bottom: 15px;",
                                actionButton("btn_undo_matrice", "", icon = shiny::icon("arrow-left"), class = "btn-default"),
                                actionButton("btn_redo_matrice", "", icon = shiny::icon("arrow-right"), class = "btn-default"),
                                span(textOutput("texte_version_matrice"), style = "margin-left: 10px; font-weight: bold; vertical-align: middle;")
                              ),
                              hr(),
                              p(style = "color: #666; font-size: 13px;", "Une fois vos modifications par axes appliquées, enregistrez définitivement cette version pour recalculer vos graphiques :"),
                              actionButton("sauvegarder_matrice_echantillon", "💾 Sauvegarder et Appliquer à cet échantillon", class = "btn-success", style = "width:100%; font-weight: bold;")
                            )
                          )
                        ),
                        br(),
                        # Visualisation en lecture seule de la matrice en cours pour cet échantillon
                        h4("Aperçu de la matrice actuelle pour ce tube :"),
                        div(style = "overflow-x: auto;", tableOutput("table_matrice_courante_ech"))
                      )
                    )
                  )
                )
              ),
              
              # --- Dans votre UI (juste après l'onglet Figures & Validations) ---
              tabPanel(
                title = "Sauvegarde & Export",
                icon = icon("download"),
                fluidRow(
                  column(
                    width = 6,
                    box(
                      title = "Exportation des Fichiers FCS Compensés", width = 12, status = "success", solidHeader = TRUE,
                      p(style = "color:gray; font-size:12px;", 
                        "Sélectionnez un ou plusieurs échantillons. L'application exportera un fichier .fcs unique s'il y a un seul choix, ou une archive .zip s'il y en a plusieurs."),
                      br(),
                      uiOutput("ui_select_fcs_export"),
                      br(),
                      
                      div(id = "wrapper_single_fcs", style = "display: none;",
                          downloadButton("btn_fcs_unique", "Télécharger l'échantillon (.fcs)", 
                                         class = "btn-success btn-lg", style = "width:100%; font-weight:bold;")
                      ),
                      div(id = "wrapper_zip_fcs", style = "display: none;",
                          downloadButton("btn_fcs_zip", "Télécharger les fichiers FCS (.zip)", 
                                         class = "btn-success btn-lg", style = "width:100%; font-weight:bold;")
                      )
                    )
                  ),
                  
                  column(
                    width = 6,
                    box(
                      title = "Exportation de la Session Technique", width = 12, status = "primary", solidHeader = TRUE,
                      p(style = "color:gray; font-size:12px;", 
                        "Téléchargez l'intégralité des configurations de la session (matrice de compensation, transformations arsinh, configurations de gates et graphiques associés) dans un fichier unique RDS."),
                      br(),
                      downloadButton("download_session_rds", "💾 Sauvegarder la Session Complète (.rds)", 
                                     class = "btn-primary btn-lg", style = "width:100%; font-weight:bold;")
                    )
                  )
                )
              )
            ) # Fin tabBox Compensation
          ) # Fin column
        ) # Fin fluidRow
      ),
      
      # ============================================================
      # ⚠️ ONGLETS EN ATTENTE (PLACEHOLDERS)
      # ============================================================
      tabItem(
        tabName = "unmixing_tab", 
        box(title = "Unmixing Spectral", status = "info", solidHeader = TRUE, width = 12, "Module en attente d'intégration (Débloqué uniquement en mode Spectral)")
      ),
      
      # ============================================================
      # ONGLET NETTOYAGE 
      # ============================================================
      
      tabItem(
        tabName = "nettoyage_tab", 
        fluidRow(
          tabBox(
            title = tagList(shiny::icon("broom"), "Pipeline de Nettoyage & Contrôle Qualité"),
            id = "preprocessing_tabs",
            width = 12,
            side = "left",
            
            # --------------------------------------------------------
            # SOUS-ONGLET 1 : Contrôle Qualité (PeacoQC ou flowAI)
            # --------------------------------------------------------
            tabPanel(
              title = "1. PeacoQC ou flowAI", 
              value = "tab_qc",
              fluidRow(
                # PANNEL DE GAUCHE : Configuration et Lancement GLOBAL
                box(title = "Configuration du CQ (Tous les échantillons)", width = 4, status = "primary", solidHeader = TRUE,
                    radioButtons("choix_qc", "Méthode de contrôle qualité :",
                                 choices = c("PeacoQC" = "peacoqc", "flowAI" = "flowai", "Aucune" = "none"), selected = "none"),
                    hr(),
                    uiOutput("ui_parametres_peacoqc"),
                    uiOutput("ui_parametres_flowai"),
                    hr(),
                    actionButton("valider_qc", "Lancer le Contrôle Qualité Global", class = "btn-success", width = "100%")
                ),
                
                # PANNEL DE DROITE : Inspection VISUELLE par échantillon
                box(title = "Visualisation du Contrôle Qualité par tube", width = 8, status = "info",
                    # Le sélecteur d'échantillon est maintenant judicieusement placé ici, dédié à la vue
                    selectInput("echantillon_visu_qc", "Choisir l'échantillon à inspecter :", choices = NULL, width = "50%"),
                    hr(),
                    plotOutput("plot_qc")
                )
              )
            ),
            
            # --------------------------------------------------------
            # SOUS-ONGLET 2 : Bordures & Débris
            # --------------------------------------------------------
            tabPanel(
              title = "2. Bordures & Débris", 
              value = "tab_debris",
              fluidRow(
                # PANNEAU DE GAUCHE : Configuration de la source
                box(title = "Configuration de Base", width = 3, status = "primary", solidHeader = TRUE,
                    radioButtons("source_donnees_debris", "Données de départ :",
                                 choices = c("Données compensées brutes" = "brutes",
                                             "Nettoyées par PeacoQC" = "peacoqc",
                                             "Nettoyées par flowAI" = "flowai"),
                                 selected = "brutes"),
                    hr(),
                    # Bouton dynamique pour activer/désactiver le filtre des bordures
                    uiOutput("ui_bouton_bordures"),
                    helpText("Vous pouvez activer ou désactiver ce filtre à tout moment, avant ou après avoir dessiné votre polygone.")
                ),
                
                # PANNEAU DE DROITE : Graphique, Boutons d'action et Tableau
                box(title = "Visualisation & Gating Morphologique", width = 9, status = "info",
                    fluidRow(
                      column(4, selectInput("echantillon_visu_debris", "Échantillon à inspecter :", choices = NULL, width = "100%")),
                      column(8, class = "text-right", style = "margin-top: 25px;",
                             actionButton("reset_debris_gate", "Effacer le tracé", class = "btn-warning"),
                             actionButton("appliquer_tous_debris", "Appliquer ce format à tous", class = "btn-primary"),
                             actionButton("valider_debris", "Valider la sélection pour ce tube", class = "btn-success")
                      )
                    ),
                    hr(),
                    # Zone du graphique interactif
                    plotOutput("plot_debris", click = "click_debris"),
                    hr(),
                    h4(shiny::icon("table"), "Statistiques de population pour ce tube"),
                    tableOutput("table_stats_debris")
                )
              )
            ),

            # --------------------------------------------------------
            # SOUS-ONGLET 3 : Doublets
            # --------------------------------------------------------
            tabPanel(
              title = "3. Doublets", 
              value = "tab_doublets",
              fluidRow(
                box(title = "Configuration Singlets", width = 4, status = "primary", solidHeader = TRUE,
                    checkboxGroupInput("doublets_axes", "Filtres à appliquer :", 
                                       choices = c("Selon FSC (FSC-H vs FSC-A)" = "FSC", "Selon SSC (SSC-H vs SSC-A)" = "SSC"), 
                                       selected = "FSC"),
                    hr(),
                    radioButtons("type_doublets", "Type de discrimination :",
                                 choices = c("Automatique (Seuil MAD)" = "auto", "Manuelle (Polygone)" = "manuel")),
                    conditionalPanel(
                      condition = "input.type_doublets == 'auto'",
                      numericInput("sensibilite_doublets", "Facteur de sensibilité (MAD) :", value = 4, min = 1, max = 10, step = 0.5)
                    ),
                    conditionalPanel(
                      condition = "input.type_doublets == 'manuel'",
                      helpText("Cliquez sur le graphique pour entourer la population de singulets."),
                      actionButton("reset_doublets_gate", "Effacer le tracé", class = "btn-warning btn-sm")
                    ),
                    hr(),
                    actionButton("valider_doublets", "Retirer les Doublets", class = "btn-success", width = "100%")
                ),
                box(title = "Visualisation des Événements Uniques", width = 8, status = "info",
                    radioButtons("visu_doublets_axis", "Axe d'analyse affiché :", choices = c("FSC", "SSC"), inline = TRUE),
                    plotOutput("plot_doublets", click = "click_doublets")
                )
              )
            ),
            
            # --------------------------------------------------------
            # SOUS-ONGLET 4 : Viabilité
            # --------------------------------------------------------
            tabPanel(
              title = "4. Viabilité", 
              value = "tab_viabilite",
              fluidRow(
                box(title = "Configuration Viabilité", width = 4, status = "primary", solidHeader = TRUE,
                    selectInput("canal_fsc_viab", "Axe X (Taille) :", choices = NULL),
                    selectInput("canal_fluo_viab", "Axe Y (Marqueur de Mort / Viabilité) :", choices = NULL),
                    hr(),
                    helpText("Entourez la population négative (cellules vivantes n'ayant pas capté le marqueur)."),
                    actionButton("reset_viab_gate", "Effacer le tracé", class = "btn-warning btn-sm"),
                    hr(),
                    actionButton("valider_viabilite", "Exclure les cellules mortes", class = "btn-success", width = "100%")
                ),
                box(title = "Discrimination Vivantes / Mortes", width = 8, status = "info",
                    plotOutput("plot_viabilite", click = "click_viabilite")
                )
              )
            )
          )
        )
      ),
      
      # ============================================================
      # ⚠️ ONGLETS EN ATTENTE (PLACEHOLDERS)
      # ============================================================
      
      tabItem(
        tabName = "analyses_tab", 
        box(title = "Analyses Statistiques", status = "info", solidHeader = TRUE, width = 12, "Bientôt disponible.")
      )
    )
  )
)

#########################
####### SERVER #########
######################

server <- function(input, output, session) {
  
  # ============================================================
  # 🟣 VARIABLES RÉACTIVES GLOBALES
  # ============================================================
  r_dictionnaire      <- reactiveVal(global_canaux)
  r_df_mono           <- reactiveVal(NULL)
  r_df_ech            <- reactiveVal(NULL)
  r_obj_R6            <- reactiveVal(NULL)
  r_matrice_marqueurs <- reactiveVal(NULL)   # Stocke le dataframe croisé Tubes x Canaux
  r_liste_groupes     <- reactiveVal(list())
  matrice_en_cours <- reactiveVal(NULL)
  compensation_lancee <- reactiveVal(FALSE)
  refresh_plot_trans   <- reactiveVal(0)
  refresh_matrix_table <- reactiveVal(0)
  refresh_comp_plot    <- reactiveVal(0)
  coords_polygone_debris <- reactiveVal(data.frame(x = numeric(), y = numeric()))
  status_bordures <- reactiveVal(list())
  gabarit_global_debris <- reactiveVal(NULL)
  
  # ============================================================
  # 🔒 VERROUS ONGLETS (TECHNOLOGIE CYTOMÈTRE)
  # ============================================================
  observe({
    req(r_df_mono())
    type_cytometre <- input$cyto_type_global
    if (is.null(type_cytometre)) return()
    
    if (type_cytometre == "Conventionnel") {
      shinyjs::removeClass(selector = "a[data-value='compensation_tab']", class = "disabled-tab")
      shinyjs::addClass(selector = "a[data-value='unmixing_tab']", class = "disabled-tab")
    } else {
      shinyjs::addClass(selector = "a[data-value='compensation_tab']", class = "disabled-tab")
      shinyjs::removeClass(selector = "a[data-value='unmixing_tab']", class = "disabled-tab")
    }
  })
  
  # ============================================================
  # 📥 LOGIQUE SOUS-ONGLET 1 : FICHIERS
  # ============================================================
  
  # --- 1.1 Chargement des Monomarqués ---
  observeEvent(input$load_mono, {
    f <- input$load_mono
    noms_proposes <- gsub(".fcs", "", f$name, ignore.case = TRUE)
    
    canaux_choisis <- sapply(noms_proposes, function(np) {
      match <- r_dictionnaire()[sapply(r_dictionnaire(), function(dc) grepl(dc, np, ignore.case = TRUE) || grepl(np, dc, ignore.case = TRUE))]
      if(length(match) > 0) return(match[1]) else return(r_dictionnaire()[1])
    })
    
    df <- data.frame(
      fichier = f$name,
      canal = canaux_choisis,
      type = ifelse(grepl("unstained", noms_proposes, ignore.case = TRUE), "Unstained", "Monomarque"),
      chemin = f$datapath,
      stringsAsFactors = FALSE
    )
    r_df_mono(df)
  })
  
  # --- 1.2 Chargement Flash des Échantillons (Métadonnées à la volée) ---
  observeEvent(input$load_ech, {
    f <- input$load_ech
    withProgress(message = "Lecture des en-têtes FCS...", value = 0, {
      liste_meta <- lapply(seq_len(nrow(f)), function(i) {
        incProgress(1 / nrow(f), detail = paste("Fichier", i, "/", nrow(f)))
        chemin_fichier <- f$datapath[i]
        tryCatch({
          fcs_header <- flowCore::read.FCS(chemin_fichier, dataset = 1, which.lines = 1, transformation = FALSE, truncate_max_range = FALSE)
          kw <- flowCore::keyword(fcs_header)
          nb_events <- ifelse(!is.null(kw[["$TOT"]]), as.numeric(kw[["$TOT"]]), 0)
          tube_name <- ifelse(!is.null(kw[["TUBE NAME"]]), kw[["TUBE NAME"]], "Inconnu")
          exp_name  <- ifelse(!is.null(kw[["EXPERIMENT NAME"]]), kw[["EXPERIMENT NAME"]], "Inconnu")
          cytometre <- ifelse(!is.null(kw[["$CYT"]]), kw[["$CYT"]], "Inconnu")
          date_acq  <- ifelse(!is.null(kw[["$DATE"]]), kw[["$DATE"]], "Inconnu")
        }, error = function(e) {
          nb_events <<- 0; tube_name <<- "Erreur"; exp_name <<- "Erreur"; cytometre <<- "Erreur"; date_acq <<- "Erreur"
        })
        
        taille_octets <- file.info(chemin_fichier)$size
        if (is.na(taille_octets)) {
          poids_formate <- "Inconnu"
        } else if (taille_octets >= 10^9) {
          poids_formate <- paste0(round(taille_octets / 10^9, 2), " Go")
        } else {
          poids_formate <- paste0(round(taille_octets / 10^6, 1), " Mo")
        }
        
        return(data.frame(
          nb_events = format(nb_events, big.mark = " "),
          tube_name = tube_name,
          exp_name = exp_name,
          cytometre = cytometre,
          poids = poids_formate,
          date = date_acq,
          stringsAsFactors = FALSE
        ))
      })
      
      df_meta <- do.call(rbind, liste_meta)
      df_final <- data.frame(
        tube_name = df_meta$tube_name, 
        fichier = f$name,
        chemin = f$datapath,
        stringsAsFactors = FALSE
      )
      
      df_final <- cbind(df_final, df_meta[, colnames(df_meta) != "tube_name", drop = FALSE])
      r_df_ech(df_final)
    })
  })
  
  # --- 1.3 Rendu Table Monomarqués ---
  output$table_mono <- renderDT({
    req(r_df_mono())
    df <- r_df_mono()
    dict <- r_dictionnaire()
    
    datalist_options <- paste0(sapply(dict, function(opt) paste0("<option value='", opt, "'>")), collapse = "")
    datalist_html <- paste0("<datalist id='canaux_list'>", datalist_options, "</datalist>")
    
    select_canal <- sapply(1:nrow(df), function(i) {
      paste0("<input type='text' class='dt-datalist-input' list='canaux_list' id='select_canal_", i, "' value='", df$canal[i], "' onchange='Shiny.setInputValue(\"change_table_canal\", {row: ", i, ", val: this.value}, {priority: \"event\"})'>", datalist_html)
    })
    
    select_type <- sapply(1:nrow(df), function(i) {
      opt_m <- if(df$type[i] == "Monomarque") " selected" else ""
      opt_u <- if(df$type[i] == "Unstained") " selected" else ""
      paste0("<select class='dt-select' id='select_type_", i, "' onchange='Shiny.setInputValue(\"change_table_type\", {row: ", i, ", val: this.value}, {priority: \"event\"})'>",
             "<option value='Monomarque'", opt_m, ">Monomarque</option>",
             "<option value='Unstained'", opt_u, ">Unstained</option>",
             "</select>")
    })
    
    display_df <- data.frame(fichier = df$fichier, canal = select_canal, type = select_type, stringsAsFactors = FALSE)
    
    datatable(
      display_df, colnames = c("Fichier FCS", "Canal Affecté", "Type de Tube"), escape = FALSE,
      options = list(dom = 't', paging = FALSE, ordering = FALSE,
                     preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
                     drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); }')),
      selection = 'none'
    )
  })
  
  observeEvent(input$change_table_canal, {
    info <- input$change_table_canal
    df <- r_df_mono()
    df[info$row, "canal"] <- info$val
    r_df_mono(df)
    if (!(info$val %in% r_dictionnaire()) && info$val != "") {
      r_dictionnaire(c(r_dictionnaire(), info$val))
    }
  })
  
  observeEvent(input$change_table_type, {
    info <- input$change_table_type
    df <- r_df_mono()
    df[info$row, "type"] <- info$val
    r_df_mono(df)
  })
  
  # --- 1.4 Rendu Table Échantillons ---
  output$table_ech <- renderDT({
    req(r_df_ech())
    df <- r_df_ech()
    
    datatable(
      df[, c("tube_name", "fichier", "nb_events", "exp_name", "cytometre", "poids", "date")], 
      colnames = c("Nom du tube (Modifiable)", "Fichier d'origine", "Nombre d'évènements", "Expérience", "Cytomètre", "Volume", "Date d'acquisition"),
      editable = list(target = "cell", disable = list(columns = 2:7)),
      options = list(dom = 't', paging = FALSE, scrollX = TRUE), 
      selection = 'none'
    )
  })
  
  observeEvent(input$table_ech_cell_edit, {
    info <- input$table_ech_cell_edit
    df <- r_df_ech()
    df[info$row, "tube_name"] <- info$value
    r_df_ech(df)
  }) 
  
  
  # --- 1.6 Initialisation R6 ---
  output$ui_cyto_select <- renderUI({
    selectInput("cyto_type_global", "Technologie Cytomètre :", choices = c("Conventionnel", "Spectral"))
  })
  
  observeEvent(input$init_r6, {
    req(r_df_mono(), r_df_ech())
    withProgress(message = "Instanciation de la classe R6 & Lecture FCS...", value = 0.5, {
      tryCatch({
        obj <- CARROT$new(r_df_mono(), r_df_ech())
        obj$charger_fcs()
        r_df_ech(obj$chemins_echantillons)
        r_obj_R6(obj)
        showNotification("Objet R6 créé et fichiers FCS chargés avec succès !", type = "message")
      }, error = function(e) {
        showNotification(paste("Erreur d'initialisation :", e$message), type = "error", duration = NULL)
      })
    })
  })
  
  # --- ✨ LOGIQUE SOUS-ONGLET 3 : GESTION DES GROUPES ---
  
  # ============================================================
  # 🏷️ LOGIQUE SOUS-ONGLET 3 : GESTION DES GROUPES (V2 INTUITIVE)
  # ============================================================
  
  # --- 3.1 Ajout d'un groupe dans la liste ---
  observeEvent(input$add_group, {
    req(input$group_name)
    nom_g <- trimws(input$group_name)
    groupes_actuels <- r_liste_groupes()
    
    if (nom_g != "" && !(nom_g %in% groupes_actuels)) {
      r_liste_groupes(c(groupes_actuels, nom_g))
      updateTextInput(session, "group_name", value = "") # Vide le champ de texte
    }
  })
  
  # --- 3.2 Affichage de la liste des groupes créés à gauche ---
  output$ui_assignation_groupes <- renderUI({
    groupes <- r_liste_groupes()
    if (length(groupes) == 0) {
      return(p(style = "color:gray; font-style:italic;", "Aucun groupe créé pour le moment."))
    }
    
    # Simple liste textuelle sous forme de badges pour voir ce qui est disponible
    tags$div(
      lapply(groupes, function(g) {
        span(g, class = "badge", style = "background-color: #605ca8; color: white; margin: 3px; padding: 5px 10px; display: inline-block; border-radius: 4px;")
      })
    )
  })
  
  # Variable interne pour mémoriser l'état des choix de groupes par échantillon
  # Clé = Nom du tube, Valeur = Groupe sélectionné
  r_choix_groupes_save <- reactiveVal(list())
  
  # --- 3.3 Rendu du Tableau de Résumé avec Selecteurs Intégrés ---
  output$table_resume_groupes <- renderDT({
    req(r_df_ech())
    df <- r_df_ech()
    groupes <- r_liste_groupes()
    choix_sauvegardes <- r_choix_groupes_save()
    
    # Génération des balises HTML <select> pour chaque ligne
    select_groupes_html <- sapply(seq_len(nrow(df)), function(i) {
      nom_tube <- df$tube_name[i]
      
      # Récupère la valeur précédemment choisie ou met "Non assigné" par défaut
      valeur_actuelle <- ifelse(!is.null(choix_sauvegardes[[nom_tube]]), choix_sauvegardes[[nom_tube]], "")
      
      # Construction des options du menu déroulant
      opt_default <- paste0("<option value=''", ifelse(valeur_actuelle == "", " selected", ""), ">Non assigné</option>")
      opt_groupes <- sapply(groupes, function(g) {
        selected_attr <- ifelse(valeur_actuelle == g, " selected", "")
        paste0("<option value='", g, "'", selected_attr, ">", g, "</option>")
      })
      
      options_totale <- paste0(c(opt_default, opt_groupes), collapse = "")
      
      # Code HTML injecté dans la cellule du tableau avec un trigger Shinjs / Input
      paste0(
        "<select class='dt-select' style='width:100%;' onchange='Shiny.setInputValue(\"change_table_groupe\", {tube: \"", nom_tube, "\", val: this.value}, {priority: \"event\"})'>",
        options_totale,
        "</select>"
      )
    })
    
    # Dataframe final pour l'affichage DT
    df_res <- data.frame(
      "Échantillon" = df$tube_name,
      "Groupe Assigné" = select_groupes_html,
      "Fichier FCS" = df$fichier,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    
    datatable(
      df_res, escape = FALSE, # Très important pour interpréter le code HTML du <select>
      options = list(
        dom = 't', paging = FALSE, ordering = FALSE,
        preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
        drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); }')
      ), 
      selection = 'none'
    )
  })
  
  # --- 3.4 Intercepter le changement du menu déroulant dans le tableau ---
  observeEvent(input$change_table_groupe, {
    info <- input$change_table_groupe
    choix_actuels <- r_choix_groupes_save()
    
    # On met à jour l'état en mémoire
    choix_actuels[[info$tube]] <- info$val
    r_choix_groupes_save(choix_actuels)
  })
  
  # --- 3.5 Enregistrement final de la cohorte ---
  observeEvent(input$save_cohort, {
    req(r_df_ech())
    df <- r_df_ech()
    choix_actuels <- r_choix_groupes_save()
    
    # On reconstruit proprement une liste nommée propre pour validation
    # (Sécurité : si un tube n'a pas été touché, il sera explicitement "Non assigné")
    liste_validation <- lapply(df$tube_name, function(tube) {
      val <- choix_actuels[[tube]]
      if (is.null(val) || val == "") return("Non assigné")
      return(val)
    })
    names(liste_validation) <- df$tube_name
    
    # Optionnel : Vous pouvez ici injecter cette liste directement dans votre objet R6 
    # si votre classe possède une méthode dédiée, ex: r_obj_R6()$definir_groupes(liste_validation)
    
    # Notification visuelle de succès
    showNotification(
      "Structure de la cohorte enregistrée avec succès !", 
      type = "message", 
      duration = 5
    )
  })
  
  # ============================================================
  # 🏷️ ✨ NOUVEAU : SOUS-ONGLET 2 : CONFIGURATION DES MARQUEURS
  # ============================================================
  
  # Générer dynamiquement la matrice Tubes x Canaux dès que les échantillons sont chargés
  observeEvent(input$load_ech, {
    req(r_df_ech())
    df_ech <- r_df_ech()
    
    withProgress(message = "Génération de la matrice de tous les marqueurs...", value = 0.5, {
      premier_fichier <- df_ech$chemin[1]
      fcs_header <- flowCore::read.FCS(premier_fichier, dataset = 1, which.lines = 1)
      pdata <- flowCore::pData(flowCore::parameters(fcs_header))
      tous_les_canaux <- pdata$name
      matrice <- as.data.frame(matrix("", nrow = nrow(df_ech), ncol = length(tous_les_canaux)))
      colnames(matrice) <- tous_les_canaux
      rownames(matrice) <- df_ech$tube_name
      
      for (i in seq_len(nrow(df_ech))) {
        try({
          f_head <- flowCore::read.FCS(df_ech$chemin[i], dataset = 1, which.lines = 1)
          pd <- flowCore::pData(flowCore::parameters(f_head))
          
          for (canal in tous_les_canaux) {
            idx <- which(pd$name == canal)
            desc_val <- pd$desc[idx]
            is_non_bio <- grepl("FSC|SSC|Time", canal, ignore.case = TRUE)
            
            if (is_non_bio) {
              matrice[i, canal] <- canal
            } else {
              if (length(desc_val) > 0 && !is.na(desc_val) && trimws(desc_val) != "") {
                matrice[i, canal] <- trimws(desc_val)
              } else {
                matrice[i, canal] <- "-" # Si la description $PnS est vide
              }
            }
          }
        }, silent = TRUE)
      }
      
      matrice_finale <- cbind("Nom du tube" = df_ech$tube_name, matrice)
      r_matrice_marqueurs(matrice_finale)
    })
  })
  
  # Rendu de la table éditable
  output$table_matrice_marqueurs <- renderDT({
    req(r_matrice_marqueurs())
    datatable(
      r_matrice_marqueurs(),
      editable = list(target = "cell", disable = list(columns = c(1))), # Colonne "Nom du tube" non modifiable ici
      options = list(dom = 'ltip', pageLength = 10, scrollX = TRUE),
      selection = 'none'
    )
  })
  
  # Sauvegarde en temps réel de la modification manuelle d'une cellule (Marqueur)
  observeEvent(input$table_matrice_marqueurs_cell_edit, {
    info <- input$table_matrice_marqueurs_cell_edit
    mat <- r_matrice_marqueurs()
    mat[info$row, info$col] <- info$value
    r_matrice_marqueurs(mat)
  })
  
  # Transmettre les modifications de marqueurs à l'objet R6 en temps réel
  observeEvent(r_matrice_marqueurs(), {
    req(r_obj_R6(), r_matrice_marqueurs())
    mat <- r_matrice_marqueurs()
    canaux_tech <- colnames(mat)[-1]
    mapping <- sapply(canaux_tech, function(canal) {
      valeurs <- mat[[canal]]
      valeurs_propres <- valeurs[valeurs != "" & valeurs != "-"]
      if (length(valeurs_propres) > 0) return(valeurs_propres[1]) else return(canal)
    })
    r_obj_R6()$mettre_a_jour_marqueurs(mapping)
    refresh_plot_trans(refresh_plot_trans() + 1)
  })
  
  # ============================================================
  # 🟣 3. WORKFLOW DE COMPENSATION INTERACTIVE (DÉCLENCHÉ PAR BOUTON)
  # ============================================================

  # 1. INITIALISATION IMMÉDIATE AU CHARGEMENT DES DONNÉES BRUTES
  observe({
    req(r_obj_R6())
    obj <- r_obj_R6()
    noms_ech <- names(obj$echantillons) 
    req(noms_ech)
    
    for (ech in noms_ech) {
      if (is.null(obj$matrices_spillover[[ech]])) {
        ff <- obj$echantillons[[ech]]
        req(ff)
        meta <- flowCore::keyword(ff)
        mat_fcs <- NULL
        
        if (!is.null(meta$SPILL)) mat_fcs <- meta$SPILL
        else if (!is.null(meta$`$SPILLOVER`)) mat_fcs <- meta$`$SPILLOVER`
        else if (!is.null(meta$SPILLOVER)) mat_fcs <- meta$SPILLOVER
        else if (!is.null(meta$COMP)) mat_fcs <- meta$COMP
        
        if (is.null(mat_fcs)) {
          canaux_fluo <- obj$canaux 
          if (is.null(canaux_fluo)) {
            tous_canaux <- colnames(flowCore::exprs(ff))
            canaux_fluo <- grep("FSC|SSC|Time", tous_canaux, invert = TRUE, value = TRUE, ignore.case = TRUE)
          }
          if (length(canaux_fluo) > 0) {
            mat_fcs <- diag(length(canaux_fluo))
            colnames(mat_fcs) <- rownames(mat_fcs) <- canaux_fluo
          }
        }
        if (!is.null(mat_fcs)) obj$enregistrer_matrice_echantillon(ech, mat_fcs)
      }
    }
  })
  
  # 2. CHARGEMENT DE LA MATRICE AU CHANGEMENT D'ÉCHANTILLON
  observeEvent(input$sel_fig_ech, {
    req(r_obj_R6(), input$sel_fig_ech)
    obj <- r_obj_R6()
    ech <- input$sel_fig_ech
    if (!is.null(obj$matrices_spillover[[ech]])) {
      matrice_en_cours(obj$matrices_spillover[[ech]])
    }
  })
  
  # 3. ACTION : CLIC SUR LE BOUTON "LANCER COMPENSATION GLOBALE"
  observeEvent(input$run_compensation, {
    req(r_obj_R6())
    obj <- r_obj_R6()
    
    withProgress(message = "Initialisation de la vue de compensation...", {
      # Pour chaque échantillon, si aucun calcul lourd n'a été fait, 
      # on pré-compense le FCS à l'aide de sa matrice de métadonnées d'origine
      for (ech in names(obj$echantillons)) {
        if (is.null(obj$echantillons_compenses[[ech]])) {
          base_fcs <- if (!is.null(obj$fcs_compenses[[ech]])) obj$fcs_compenses[[ech]] else obj$echantillons[[ech]]
          mat_actuelle <- obj$matrices_spillover[[ech]]
          if (is.null(mat_actuelle)) mat_actuelle <- obj$S_matrix
          
          if (!is.null(base_fcs) && !is.null(mat_actuelle)) {
            obj$echantillons_compenses[[ech]] <- flowCore::compensate(base_fcs, mat_actuelle)
          }
        }
      }
      r_obj_R6(obj)
      compensation_lancee(TRUE) # On lève le verrou, le graphique peut maintenant s'afficher
      refresh_comp_plot(refresh_comp_plot() + 1) # Force le rafraîchissement
    })
    showNotification("⚡ Visualisation globale initialisée avec succès !", type = "message")
  })
  
  # 4. RENDU DU GRAPHIQUE (VERROUILLÉ TANT QU'ON A PAS CLIQUÉ SUR LE BOUTON)
  output$plot_validation_comp <- renderPlot({
    # Si le bouton n'a pas été cliqué, on affiche un message blanc/vide ou d'instruction
    if (!compensation_lancee()) {
      plot.new()
      text(0.5, 0.5, "Veuillez cliquer sur le bouton '⚡ Lancer Compensation Globale' \npour afficher les graphiques.", 
           cex = 1.2, font = 2, col = "gray40")
      return(NULL)
    }
    
    # Déclencheurs réactifs normaux
    refresh_comp_plot() 
    req(r_obj_R6(), input$sel_fig_ech, input$fig_mode, input$fig_view_type)
    
    obj <- r_obj_R6()
    ech <- input$sel_fig_ech
    
    fcs_avant <- if (!is.null(obj$fcs_compenses[[ech]])) obj$fcs_compenses[[ech]] else obj$echantillons[[ech]]
    fcs_apres <- obj$echantillons_compenses[[ech]]
    
    req(fcs_avant, fcs_apres)
    
    # Sélection des axes X et Y
    if (input$fig_mode == "free") {
      req(input$sel_fig_x, input$sel_fig_y)
      axe_x <- input$sel_fig_x
      axe_y <- input$sel_fig_y
    } else {
      req(input$sel_axe_ligne, input$sel_axe_colonne)
      axe_x <- input$sel_axe_ligne
      axe_y <- input$sel_axe_colonne
    }
    
    # Sous-échantillonnage des points
    max_pts <- input$fig_max_points
    exprs_av <- flowCore::exprs(fcs_avant)
    exprs_ap <- flowCore::exprs(fcs_apres)
    
    idx_av <- sample(1:nrow(exprs_av), min(max_pts, nrow(exprs_av)))
    idx_ap <- sample(1:nrow(exprs_ap), min(max_pts, nrow(exprs_ap)))
    
    # Dessin des graphiques
    if (input$fig_view_type == "Avant et après compensation") {
      par(mfrow = c(1, 2))
      plot(exprs_av[idx_av, axe_x], exprs_av[idx_av, axe_y], pch = 20, col = rgb(0, 0, 1, 0.3), xlab = axe_x, ylab = axe_y, main = paste("Avant -", ech))
      plot(exprs_ap[idx_ap, axe_x], exprs_ap[idx_ap, axe_y], pch = 20, col = rgb(1, 0, 0, 0.3), xlab = axe_x, ylab = axe_y, main = paste("Après -", ech))
    } else if (input$fig_view_type == "Avant compensation uniquement") {
      par(mfrow = c(1, 1))
      plot(exprs_av[idx_av, axe_x], exprs_av[idx_av, axe_y], pch = 20, col = rgb(0, 0, 1, 0.4), xlab = axe_x, ylab = axe_y, main = paste("Avant (Brut) -", ech))
    } else if (input$fig_view_type == "Après compensation uniquement") {
      par(mfrow = c(1, 1))
      plot(exprs_ap[idx_ap, axe_x], exprs_ap[idx_ap, axe_y], pch = 20, col = rgb(1, 0, 0, 0.4), xlab = axe_x, ylab = axe_y, main = paste("Après (Compensé) -", ech))
    }
  }, height = function() {
    if (input$fig_view_type == "Avant et après compensation") 400 else 500
  })
  
  # 5. ACTION DU BOUTON EN BAS : SAUVEGARDER ET METTRE À JOUR LA MODIFICATION PAR AXE
  observeEvent(input$sauvegarder_matrice_echantillon, {
    req(r_obj_R6(), input$sel_fig_ech, matrice_en_cours())
    
    # Sécurité additionnelle : Si l'utilisateur clique en bas sans avoir lancé le bouton du haut
    if (!compensation_lancee()) {
      showNotification("⚠️ Veuillez d'abord initialiser l'affichage avec le bouton 'Lancer Compensation Globale' en haut.", type = "error")
      return(NULL)
    }
    
    obj <- r_obj_R6()
    ech <- input$sel_fig_ech
    
    withProgress(message = "Mise à jour du graphique en cours...", {
      obj$enregistrer_matrice_echantillon(ech, matrice_en_cours())
      
      base_fcs <- if (!is.null(obj$fcs_compenses[[ech]])) obj$fcs_compenses[[ech]] else obj$echantillons[[ech]]
      
      if (!is.null(base_fcs)) {
        obj$echantillons_compenses[[ech]] <- flowCore::compensate(base_fcs, matrice_en_cours())
      }
      r_obj_R6(obj)
      refresh_comp_plot(refresh_comp_plot() + 1)
    })
    showNotification(paste("🎯 Modification enregistrée et appliquée pour :", ech), type = "message")
  })
  
  # --- Les fonctions suivantes (Générateurs UI des axes, Tables, Undo/Redo) restent exactement les mêmes ---
  
  # 6. AFFICHAGE DE LA TABLE EN LECTURE SEULE (APERÇU DE LA MATRICE DU TUBE)
  output$table_matrice_courante_ech <- renderTable({
    if (is.null(matrice_en_cours())) {
      obj <- r_obj_R6()
      req(obj)
      ech_initial <- input$sel_fig_ech
      if (!is.null(ech_initial) && !is.null(obj$matrices_spillover[[ech_initial]])) {
        matrice_en_cours(obj$matrices_spillover[[ech_initial]])
      }
    }
    req(matrice_en_cours())
    return(as.data.frame(round(matrice_en_cours(), 4)))
  }, rownames = TRUE, bordered = TRUE, striped = TRUE, spacing = "s")
  
  # 7. TEXTE DE LA VERSION ACTUELLE (TIMELINE UNDO/REDO)
  output$texte_version_matrice <- renderText({
    req(r_obj_R6(), input$sel_fig_ech)
    obj <- r_obj_R6()
    ech <- input$sel_fig_ech
    
    idx <- obj$index_historique[[ech]]
    total <- length(obj$historique_matrices[[ech]])
    if (is.null(idx) || total == 0) return("Version : 1/1 (Origine)")
    paste0("Version : ", idx, " / ", total)
  })
  
  # 8. UNDO
  observeEvent(input$btn_undo_matrice, {
    req(r_obj_R6(), input$sel_fig_ech)
    obj <- r_obj_R6()
    ech <- input$sel_fig_ech
    
    idx_actuel <- obj$index_historique[[ech]]
    if (!is.null(idx_actuel) && idx_actuel > 1) {
      nouvel_idx <- idx_actuel - 1
      obj$index_historique[[ech]] <- nouvel_idx
      mat_restauree <- obj$historique_matrices[[ech]][[nouvel_idx]]
      obj$matrices_spillover[[ech]] <- mat_restauree
      
      # Choix intelligent du FCS de base pour la restauration
      base_fcs <- if (!is.null(obj$fcs_compenses[[ech]])) obj$fcs_compenses[[ech]] else obj$echantillons[[ech]]
      
      if (!is.null(base_fcs)) {
        obj$echantillons_compenses[[ech]] <- flowCore::compensate(base_fcs, mat_restauree)
      }
      matrice_en_cours(mat_restauree)
      r_obj_R6(obj)
      refresh_comp_plot(refresh_comp_plot() + 1)
      showNotification("↩️ Version précédente restaurée", type = "default")
    }
  })
  
  # 9. REDO
  observeEvent(input$btn_redo_matrice, {
    req(r_obj_R6(), input$sel_fig_ech)
    obj <- r_obj_R6()
    ech <- input$sel_fig_ech
    
    idx_actuel <- obj$index_historique[[ech]]
    total_versions <- length(obj$historique_matrices[[ech]])
    
    if (!is.null(idx_actuel) && idx_actuel < total_versions) {
      nouvel_idx <- idx_actuel + 1
      obj$index_historique[[ech]] <- nouvel_idx
      mat_restauree <- obj$historique_matrices[[ech]][[nouvel_idx]]
      obj$matrices_spillover[[ech]] <- mat_restauree
      
      # Choix intelligent du FCS de base pour la restauration
      base_fcs <- if (!is.null(obj$fcs_compenses[[ech]])) obj$fcs_compenses[[ech]] else obj$echantillons[[ech]]
      
      if (!is.null(base_fcs)) {
        obj$echantillons_compenses[[ech]] <- flowCore::compensate(base_fcs, mat_restauree)
      }
      matrice_en_cours(mat_restauree)
      r_obj_R6(obj)
      refresh_comp_plot(refresh_comp_plot() + 1)
      showNotification("⏭️ Passage à la version suivante", type = "default")
    }
  })
  
  # 10. AFFICHAGE DE LA MATRICE D'ORIGINE PAR ÉCHANTILLON (ONGLET DÉDIÉ)
  output$matrice_origine_affichage <- renderTable({
    req(r_obj_R6(), input$sel_ech_matrice_origine)
    obj <- r_obj_R6()
    ech <- input$sel_ech_matrice_origine
    
    if (!is.null(obj$historique_matrices[[ech]]) && length(obj$historique_matrices[[ech]]) > 0) {
      return(as.data.frame(obj$historique_matrices[[ech]][[1]]))
    }
    return(NULL)
  }, rownames = TRUE, bordered = TRUE, striped = TRUE)
  
  # --- Les étapes suivantes (Transformation, Gating, etc.) restent identiques à votre code initial ---
  
  # --- ÉTAPE 1 : TRANSFORMATION ---
  observeEvent(input$btn_apply_trans, {
    req(r_obj_R6(), input$cofacteur)
    tryCatch({
      r_obj_R6()$transformer_fcs(input$cofacteur)
      refresh_plot_trans(refresh_plot_trans() + 1)
      showNotification("Transformation Arcsinh appliquée", type = "message")
    }, error = function(e) {
      showNotification(paste("Erreur lors de la transformation :", e$message), type = "error")
    })
  })
  
  output$ui_trans_file <- renderUI({ req(r_obj_R6()); selectInput("sel_trans_file", "Échantillon Monomarqué :", choices = names(r_obj_R6()$tubes_monomarques)) })
  output$ui_trans_cx   <- renderUI({ req(r_obj_R6()); selectInput("sel_trans_cx", "Axe X :", choices = r_obj_R6()$canaux, selected = r_obj_R6()$canaux[1]) })
  output$ui_trans_cy   <- renderUI({ choix <- r_obj_R6()$canaux; selectInput("sel_trans_cy", "Axe Y :", choices = choix, selected = ifelse(length(choix)>1, choix[2], choix[1])) })
  
  output$plot_transformation <- renderPlot({
    req(r_obj_R6(), input$sel_trans_file, input$sel_trans_cx, input$sel_trans_cy)
    refresh_plot_trans()
    if (length(r_obj_R6()$monomarques_trans) == 0) return(NULL)
    
    fcs_trans <- r_obj_R6()$monomarques_trans[[input$sel_trans_file]]
    req(fcs_trans)
    
    mat_data <- as.data.frame(flowCore::exprs(fcs_trans))
    ggplot(mat_data, aes(x = .data[[input$sel_trans_cx]], y = .data[[input$sel_trans_cy]])) +
      ggpointdensity::geom_pointdensity(size = 0.3, alpha = 0.5) +
      scale_color_gradientn(colours = c("darkblue", "blue", "cyan", "greenyellow", "red")) +
      theme_bw() + labs(title = paste("Contrôle transformation Arcsinh sur", input$sel_trans_file))
  })
  
  # --- ÉTAPE 2 : GATING ---
  output$ui_gate_canal <- renderUI({ req(r_obj_R6()); selectInput("sel_gate_canal", "Canal à régler :", choices = r_obj_R6()$canaux) })
  output$ui_gate_use_unstained_wrapper <- renderUI({
    req(r_obj_R6())
    if ("TUBE_UNSTAINED" %in% names(r_obj_R6()$monomarques_trans)) {
      checkboxInput("gate_use_unstained", "Utiliser le tube Unstained global", value = TRUE)
    } else NULL
  })
  
  output$ui_slider_neg <- renderUI({
    req(r_obj_R6(), input$sel_gate_canal, length(r_obj_R6()$monomarques_trans) > 0)
    nom_tube_source <- ifelse(isTRUE(input$gate_use_unstained) && "TUBE_UNSTAINED" %in% names(r_obj_R6()$monomarques_trans), "TUBE_UNSTAINED", input$sel_gate_canal)
    vals <- flowCore::exprs(r_obj_R6()$monomarques_trans[[nom_tube_source]])[, input$sel_gate_canal]
    sliderInput("slide_gate_neg", NULL, min = round(min(vals)-1,1), max = round(max(vals)+1,1), value = c(round(min(vals),1), round(min(vals)+1.5,1)), step = 0.1)
  })
  
  output$ui_slider_pos <- renderUI({
    req(r_obj_R6(), input$sel_gate_canal, length(r_obj_R6()$monomarques_trans) > 0)
    vals <- flowCore::exprs(r_obj_R6()$monomarques_trans[[input$sel_gate_canal]])[, input$sel_gate_canal]
    sliderInput("slide_gate_pos", NULL, min = round(min(vals)-1,1), max = round(max(vals)+1,1), value = c(round(max(vals)-2,1), round(max(vals),1)), step = 0.1)
  })
  
  output$plot_gating <- renderPlot({
    req(r_obj_R6(), input$sel_gate_canal, input$slide_gate_neg, input$slide_gate_pos)
    if (isTRUE(input$gate_use_unstained) && !("TUBE_UNSTAINED" %in% names(r_obj_R6()$monomarques_trans))) return(NULL)
    r_obj_R6()$graphiques_gates(nom_canal = input$sel_gate_canal, shiny_neg = input$slide_gate_neg, shiny_pos = input$slide_gate_pos, afficher_unstained_neg = isTRUE(input$gate_use_unstained))
  })
  
  observeEvent(input$save_gate, {
    req(r_obj_R6(), input$sel_gate_canal, input$slide_gate_neg, input$slide_gate_pos)
    r_obj_R6()$definir_et_extraire(nom_canal = input$sel_gate_canal, intervalle_gate_negatif = input$slide_gate_neg, intervalle_gate_positif = input$slide_gate_pos, utiliser_unstained = isTRUE(input$gate_use_unstained))
    showNotification(paste("🎯 Gate enregistré pour le canal", input$sel_gate_canal), type = "message")
    tryCatch({ r_obj_R6()$calculer_spillover(); refresh_matrix_table(refresh_matrix_table() + 1) }, error = function(e){})
  })
  
  # --- ÉTAPE 3 : MATRICE DE SPILLOVER DU PANNEAU GLOBAL ---
  observeEvent(input$reset_mat, { 
    req(r_obj_R6())
    tryCatch({ r_obj_R6()$calculer_spillover(); refresh_matrix_table(refresh_matrix_table() + 1); showNotification("Matrice brute réinitialisée !", type = "message") }, error = function(e) { showNotification(paste("Impossible de réinitialiser :", e$message), type = "error") })
  })
  
  output$table_spillover <- renderDT({
    req(r_obj_R6()); refresh_matrix_table() 
    if (is.null(r_obj_R6()$S_matrix)) return(datatable(data.frame(Message = "Veuillez configurer les Gates pour générer la matrice."), options = list(dom = 't')))
    datatable(round(as.data.frame(r_obj_R6()$S_matrix), 4), options = list(dom = 't', paging = FALSE, scrollX = TRUE))
  })
  
  output$ui_edit_mat_c1 <- renderUI({ req(r_obj_R6()); selectInput("sel_mat_c1", "Ligne (Canal principal) :", choices = r_obj_R6()$canaux) })
  output$ui_edit_mat_c2 <- renderUI({ req(r_obj_R6()); selectInput("sel_mat_c2", "Colonne (Spillover vers) :", choices = r_obj_R6()$canaux) })
  
  observeEvent(input$apply_edit_mat, {
    req(r_obj_R6(), r_obj_R6()$S_matrix, input$sel_mat_c1, input$sel_mat_c2)
    tryCatch({ r_obj_R6()$modifier_spillover(input$sel_mat_c1, input$sel_mat_c2, input$edit_mat_val); refresh_matrix_table(refresh_matrix_table() + 1); showNotification("Matrice modifiée", type = "warning") }, error = function(e) { showNotification(e$message, type = "error") })
  })
  
  # --- ÉTAPE 4 : FIGURES ET DISPOSITION DE LA MOSAÏQUE ---
  output$ui_fig_ech <- renderUI({ req(r_obj_R6()); selectInput("sel_fig_ech", "Sélectionner l'échantillon :", choices = names(r_obj_R6()$echantillons)) })
  
  observeEvent(input$run_compensation, {
    req(r_obj_R6())
    withProgress(message = "Lancement de la compensation globale brute...", {
      r_obj_R6()$compenser()
      # Par défaut, on synchronise le plot avec les matrices brutes fraîchement calculées
      for(ech in names(r_obj_R6()$echantillons_compenses)) {
        if(!is.null(r_obj_R6()$S_matrix)) {
          r_obj_R6()$enregistrer_matrice_echantillon(ech, r_obj_R6()$S_matrix)
        }
      }
      refresh_comp_plot(refresh_comp_plot() + 1)
    })
    showNotification("Tous les échantillons ont été initialisés avec la compensation globale !", type = "default")
  })
  
  output$ui_conditionnel_axes <- renderUI({
    req(input$fig_mode == "free", r_obj_R6())
    tagList(
      selectInput("sel_fig_cx", "Canal Axe X :", choices = r_obj_R6()$canaux),
      selectInput("sel_fig_cy", "Canal Axe Y :", choices = r_obj_R6()$canaux, selected = rev(r_obj_R6()$canaux)[1])
    )
  })
  
  hauteur_dynamique <- reactive({
    req(r_obj_R6(), input$fig_mode, input$fig_view_type)
    if (input$fig_mode == "free") return(600) 
    
    canaux <- r_obj_R6()$canaux
    if (length(canaux) < 2) return(400)
    
    n_combinaisons <- ncol(combn(canaux, 2))
    mapping_affichage <- c(
      "Avant et après compensation" = "Both",
      "Avant compensation uniquement" = "Before compensation only",
      "Après compensation uniquement" = "After compensation only"
    )
    valeur_affichage_r6 <- mapping_affichage[[input$fig_view_type]]
    
    if (!is.null(valeur_affichage_r6) && valeur_affichage_r6 == "Both") {
      n_rows <- n_combinaisons
    } else {
      cols_grille <- ceiling(sqrt(n_combinaisons))
      n_rows <- ceiling(n_combinaisons / cols_grille)
    }
    return(max((n_rows * 350) + 50, 400))
  })
  
  output$plot_validation_comp <- renderPlot({
    # 1. Récupération des déclencheurs et données
    refresh_comp_plot() # Notre trigger réactif pour forcer la mise à jour
    req(r_obj_R6(), input$sel_fig_ech, input$fig_mode, input$fig_view_type)
    
    obj <- r_obj_R6()
    ech <- input$sel_fig_ech
    
    # 2. Récupération du fichier "Avant" compensation
    # Si la transformation/gating globale n'a pas été faite, on prend le fichier brut d'origine
    fcs_avant <- if (!is.null(obj$fcs_compenses[[ech]])) obj$fcs_compenses[[ech]] else obj$echantillons[[ech]]
    req(fcs_avant)
    
    # 3. Récupération du fichier "Après" compensation
    # Si l'utilisateur n'a pas encore cliqué sur "Sauvegarder et Appliquer" pour ce tube,
    # on génère à la volée la compensation basée sur la matrice actuelle (métadonnées ou modifiée)
    fcs_apres <- obj$echantillons_compenses[[ech]]
    
    if (is.null(fcs_apres)) {
      mat_actuelle <- obj$matrices_spillover[[ech]]
      if (is.null(mat_actuelle)) mat_actuelle <- obj$S_matrix
      
      if (!is.null(mat_actuelle)) {
        fcs_apres <- flowCore::compensate(fcs_avant, mat_actuelle)
      } else {
        fcs_apres <- fcs_avant # Repli par défaut si aucune matrice n'existe du tout
      }
    }
    
    # 4. Extraction des axes à tracer selon le mode de l'UI
    if (input$fig_mode == "free") {
      req(input$sel_fig_x, input$sel_fig_y)
      axe_x <- input$sel_fig_x
      axe_y <- input$sel_fig_y
    } else {
      # Mode paires automatiques (Exemple: si l'utilisateur modifie une ligne/colonne spécifique)
      req(input$sel_axe_ligne, input$sel_axe_colonne)
      axe_x <- input$sel_axe_ligne
      axe_y <- input$sel_axe_colonne
    }
    
    # 5. Échantillonnage des points pour éviter les lenteurs d'affichage
    max_pts <- input$fig_max_points
    exprs_av <- flowCore::exprs(fcs_avant)
    exprs_ap <- flowCore::exprs(fcs_apres)
    
    idx_av <- sample(1:nrow(exprs_av), min(max_pts, nrow(exprs_av)))
    idx_ap <- sample(1:nrow(exprs_ap), min(max_pts, nrow(exprs_ap)))
    
    # 6. Génération dynamique du graphique selon le type de vue choisi
    if (input$fig_view_type == "Avant et après compensation") {
      par(mfrow = c(1, 2))
      # Graphique Avant
      plot(exprs_av[idx_av, axe_x], exprs_av[idx_av, axe_y], 
           pch = 20, col = rgb(0, 0, 1, 0.3), xlab = axe_x, ylab = axe_y, main = paste("Avant -", ech))
      # Graphique Après
      plot(exprs_ap[idx_ap, axe_x], exprs_ap[idx_ap, axe_y], 
           pch = 20, col = rgb(1, 0, 0, 0.3), xlab = axe_x, ylab = axe_y, main = paste("Après -", ech))
      
    } else if (input$fig_view_type == "Avant compensation uniquement") {
      par(mfrow = c(1, 1))
      plot(exprs_av[idx_av, axe_x], exprs_av[idx_av, axe_y], 
           pch = 20, col = rgb(0, 0, 1, 0.4), xlab = axe_x, ylab = axe_y, main = paste("Avant (Brut) -", ech))
      
    } else if (input$fig_view_type == "Après compensation uniquement") {
      par(mfrow = c(1, 1))
      plot(exprs_ap[idx_ap, axe_x], exprs_ap[idx_ap, axe_y], 
           pch = 20, col = rgb(1, 0, 0, 0.4), xlab = axe_x, ylab = axe_y, main = paste("Après (Compensé) -", ech))
    }
  }, height = function() {
    # Ajuste la hauteur du plot dynamiquement selon le mode d'affichage
    if (input$fig_view_type == "Avant et après compensation") 400 else 500
  })
  
  # ============================================================
  # 💾 LOGIQUE DE L'ONGLET : SAUVEGARDE & EXPORT
  # ============================================================
  
  output$ui_select_fcs_export <- renderUI({
    req(r_obj_R6(), r_df_ech())
    noms_personnalises <- r_df_ech()$tube_name
    
    if (length(noms_personnalises) == 0) {
      return(p(style = "color:orange; font-style:italic;", "⚠️ Aucun échantillon trouvé."))
    }
    
    tagList(
      selectizeInput(
        inputId = "fcs_export_choices", 
        label = "Sélectionner les échantillons à exporter :", 
        choices = noms_personnalises, 
        selected = noms_personnalises,
        multiple = TRUE,
        options = list(placeholder = 'Sélectionnez vos tubes...', plugins = list('remove_button'))
      ),
      actionButton("select_all_export", "Tout sélectionner", class = "btn-xs btn-default"),
      actionButton("clear_all_export", "Tout désélectionner", class = "btn-xs btn-default")
    )
  })
  
  observeEvent(input$select_all_export, { req(r_df_ech()); updateSelectizeInput(session, "fcs_export_choices", selected = r_df_ech()$tube_name) })
  observeEvent(input$clear_all_export, { updateSelectizeInput(session, "fcs_export_choices", selected = character(0)) })
  
  observe({
    choix <- input$fcs_export_choices
    if (length(choix) == 0) {
      shinyjs::hide("wrapper_single_fcs")
      shinyjs::hide("wrapper_zip_fcs")
    } else if (length(choix) == 1) {
      txt <- paste0("🔬 Télécharger ", choix, " (.fcs)")
      shinyjs::html(id = "btn_fcs_unique", html = paste0("<i class='fa fa-download'></i> ", txt))
      shinyjs::show("wrapper_single_fcs")
      shinyjs::hide("wrapper_zip_fcs")
    } else {
      shinyjs::show("wrapper_zip_fcs")
      shinyjs::hide("wrapper_single_fcs")
    }
  })
  
  # 1. EXPORTATION : CAS UNIQUE AVEC INJECTION DE LA MATRICE DU TUBE PROPRE
  output$btn_fcs_unique <- downloadHandler(
    filename = function() {
      req(input$fcs_export_choices)
      # Nettoie le nom du fichier en remplaçant les caractères spéciaux par des "_"
      paste0(gsub("[^a-zA-Z0-9_]", "_", input$fcs_export_choices[1]), "_compense.fcs")
    },
    content = function(file) {
      # 1. Vérification des inputs requis
      req(r_obj_R6(), r_df_ech(), input$fcs_export_choices)
      obj <- r_obj_R6()
      
      # 2. SÉCURITÉ : On bloque l'export si la liste des compensés est totalement vide
      if (is.null(obj$echantillons_compenses) || length(obj$echantillons_compenses) == 0) {
        showNotification("⚠️ Veuillez d'abord lancer la compensation globale avant d'exporter.", type = "error")
        return(NULL) # Arrête proprement la fonction de téléchargement
      }
      
      tube_nom <- input$fcs_export_choices[1]
      idx <- which(r_df_ech()$tube_name == tube_nom)
      
      if (length(idx) > 0) {
        fcs_final <- NULL
        
        # 3. SÉCURITÉ : On vérifie si le nom du tube existe dans la liste pour éviter l'indice hors limites
        if (tube_nom %in% names(obj$echantillons_compenses)) {
          fcs_final <- obj$echantillons_compenses[[tube_nom]]
        } else if (idx <= length(obj$echantillons_compenses)) {
          # Si le nom n'est pas trouvé mais que l'index numérique reste valide
          fcs_final <- obj$echantillons_compenses[[idx]]
        }
        
        # 4. Écriture du fichier si l'échantillon a bien été récupéré
        if (!is.null(fcs_final)) {
          
          # Extraction de sa matrice spécifique (ou la globale par défaut)
          mat_specifique <- obj$matrices_spillover[[tube_nom]]
          if (is.null(mat_specifique)) {
            mat_specifique <- obj$S_matrix
          }
          
          # Injection de la matrice dans les métadonnées avec la méthode moderne keyword()
          if (!is.null(mat_specifique)) {
            flowCore::keyword(fcs_final)[["$SPILLOVER"]] <- mat_specifique
          }
          
          # Génération physique du fichier FCS
          flowCore::write.FCS(fcs_final, filename = file)
          
        } else {
          # Si malgré tout l'échantillon est introuvable (non compensé par exemple)
          showNotification("⚠️ Cet échantillon n'a pas encore été traité ou est introuvable.", type = "error")
        }
      } else {
        showNotification("⚠️ Échantillon introuvable dans la liste de données.", type = "error")
      }
    }
  )
  
  # 2. EXPORTATION : CAS MULTIPLE ZIP AVEC ENREGISTREMENT DES MATRICES INDIVIDUELLES
  output$btn_fcs_zip <- downloadHandler(
    filename = function() {
      paste0("FCS_Compenses_Selection_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(r_obj_R6(), r_df_ech(), input$fcs_export_choices)
      obj <- r_obj_R6()
      tubes_choisis <- input$fcs_export_choices
      
      tmp_dir <- file.path(tempdir(), paste0("export_zip_", sample(1000:9999, 1)))
      dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
      fichiers_complets <- c()
      
      withProgress(message = "Génération de l'archive...", value = 0, {
        for (i in seq_along(tubes_choisis)) {
          nom_tube <- tubes_choisis[i]
          idx <- which(r_df_ech()$tube_name == nom_tube)
          
          if (length(idx) > 0) {
            fcs_final <- obj$echantillons_compenses[[nom_tube]]
            if(is.null(fcs_final)) fcs_final <- obj$echantillons_compenses[[idx]]
            
            if (!is.null(fcs_final)) {
              # Injection de sa propre matrice customisée avant compression
              mat_specifique <- obj$matrices_spillover[[nom_tube]]
              if (is.null(mat_specifique)) mat_specifique <- obj$S_matrix
              
              if (!is.null(mat_specifique)) {
                flowCore::keyword(fcs_final)[["$SPILLOVER"]] <- mat_specifique
              }
              
              nom_fichier_propre <- paste0(gsub("[^a-zA-Z0-9_]", "_", nom_tube), "_compense.fcs")
              chemin_complet <- file.path(tmp_dir, nom_fichier_propre)
              
              flowCore::write.FCS(fcs_final, filename = chemin_complet)
              fichiers_complets <- c(fichiers_complets, chemin_complet)
            }
          }
          incProgress(1/length(tubes_choisis))
        }
        setProgress(0.9, detail = "Compression ZIP...")
        utils::zip(zipfile = file, files = fichiers_complets, flags = "-r9Xj")
      })
      unlink(tmp_dir, recursive = TRUE)
    }
  )
  
  # 3. CAS ARCHIVE COMPLÈTE RDS
  output$download_session_rds <- downloadHandler(
    filename = function() {
      paste0("Session_Compensation_Complete_", format(Sys.time(), "%Y%m%d"), ".rds")
    },
    content = function(file) {
      req(r_obj_R6())
      withProgress(message = "Génération du fichier RDS...", {
        sauvegarde <- list(
          meta = list(
            date_export     = Sys.time(),
            canaux_utilises = r_obj_R6()$canaux
          ),
          configuration_technique = list(
            trans_list        = r_obj_R6()$trans_list,
            matrices_par_tube = r_obj_R6()$matrices_spillover,
            matrice_globale   = r_obj_R6()$S_matrix
          ),
          gating = list(
            gates_positifs = r_obj_R6()$gates_positifs,
            gates_negatifs = r_obj_R6()$gates_negatifs
          )
        )
        saveRDS(sauvegarde, file = file)
      })
    }
  )
  
  
  # ============================================================
  # 🟣 4. WORKFLOW DE NETTOYAGE DES DONNEES (MIS À JOUR R6 CENTRALISÉ)
  # ============================================================
  
  # ============================================================
  # 1. LIEN INTER-MODULES : ENCHAINEMENT DEPUIS LA COMPENSATION
  # ============================================================

  observeEvent(input$choix_qc, {
    req(input$choix_qc)
    if (input$choix_qc == "none") return(NULL)
    if (is.null(r_obj_R6())) {
      showNotification("⚠️ Aucun pipeline trouvé. Veuillez charger vos données à l'étape 1.", type = "error")
      updateRadioButtons(session, "choix_qc", selected = "none")
      return(NULL)
    }
    
    obj <- r_obj_R6()
    
    liste_fcs <- NULL
    if (!is.null(obj$echantillons_compenses) && length(obj$echantillons_compenses) > 0) {
      liste_fcs <- obj$echantillons_compenses
    } else if (!is.null(obj$fcs_compenses) && length(obj$fcs_compenses) > 0) {
      liste_fcs <- obj$fcs_compenses
    }
    
    if (is.null(liste_fcs) || length(liste_fcs) == 0) {
      showNotification("⚠️ Les données compensées ne sont pas prêtes. Avez-vous cliqué sur 'Lancer la Compensation Globale' ?", type = "warning", duration = 10)
      updateRadioButtons(session, "choix_qc", selected = "none")
      return(NULL)
    }
    
    updateSelectInput(session, "echantillon_visu_qc", choices = names(liste_fcs))
    updateSelectInput(session, "echantillon_visu_debris", choices = names(liste_fcs))
    
    if (!is.null(obj$canaux)) {
      updateSelectInput(session, "canal_fluo_viab", choices = obj$canaux)
    }
    updateSelectInput(session, "canal_fsc_viab", choices = c("FSC-A", "FSC-H", "FSC-W"))
  })
  
  # ============================================================
  # 2. DEBRIS
  # ============================================================
  
  # 1. AIGUILLAGE DYNAMIQUE DU FLOWFRAME EN COURS D'INSPECTION
  fcs_courant_debris <- reactive({
    req(r_obj_R6(), input$echantillon_visu_debris, input$source_donnees_debris)
    obj <- r_obj_R6()
    nom_ech <- input$echantillon_visu_debris
    
    # Si le bouton "Retirer les Bordures" a été activé pour ce tube
    if (!is.null(bordures_cliquees()[[nom_ech]]) && bordures_cliquees()[[nom_ech]] == TRUE) {
      return(obj$post_bordures[[nom_ech]])
    }
    
    # Sinon, on suit le bouton radio de l'UI
    if (input$source_donnees_debris == "peacoqc" && !is.null(obj$post_PeacoQC[[nom_ech]])) {
      return(obj$post_PeacoQC[[nom_ech]])
    } else if (input$source_donnees_debris == "flowai" && !is.null(obj$post_flowAI[[nom_ech]])) {
      return(obj$post_flowAI[[nom_ech]])
    } else {
      # CORRECTION ICI : Remplacement de self$ par obj$
      return(if(!is.null(obj$echantillons_compenses) && length(obj$echantillons_compenses) > 0) obj$echantillons_compenses[[nom_ech]] else obj$fcs_compenses[[nom_ech]])
    }
  })
  
  
  # ============================================================
  # 2. GENERATION DYNAMIQUE DES PARAMÈTRES (UI CONDITIONNELLE)
  # ============================================================
  
  # [Cette partie reste inchangée et est très bien écrite]
  output$ui_parametres_peacoqc <- renderUI({
    req(input$choix_qc == "peacoqc")
    tagList(
      h4("Paramètres de configuration PeacoQC", 
         style = "color: #2c3e50; font-weight: bold; margin-bottom: 15px; border-bottom: 1px solid #eee; padding-bottom: 5px;"),
      fluidRow(
        column(6, numericInput("p_min_cells", "min_cells :", value = 150, min = 10, step = 10)),
        column(6, numericInput("p_max_bins", "max_bins :", value = 100, min = 10, step = 5))
      ),
      fluidRow(
        column(6, numericInput("p_step", "step :", value = 500, min = 50, step = 50)),
        column(6, numericInput("p_mad", "MAD :", value = 6, min = 1, max = 20, step = 0.5))
      ),
      fluidRow(
        column(6, numericInput("p_it_limit", "IT_limit :", value = 0.6, min = 0, max = 1, step = 0.05)),
        column(6, numericInput("p_consecutive_bins", "consecutive_bins :", value = 5, min = 1, step = 1))
      ),
      fluidRow(
        column(6, numericInput("p_force_it", "force_IT :", value = 150, min = 0, step = 10)),
        column(6, numericInput("p_peak_removal", "peak_removal :", value = 0.33, min = 0, max = 1, step = 0.05))
      ),
      fluidRow(
        column(6, numericInput("p_min_nr_bins", "min_nr_bins_peakdetection :", value = 10, min = 2, step = 1)),
        column(6, selectInput("p_remove_zeros", "remove_zeros :", choices = c("Faux" = FALSE, "Vrai" = TRUE), selected = FALSE))
      )
    )
  })
  
  output$ui_parametres_flowai <- renderUI({
    req(input$choix_qc == "flowai")
    tagList(
      h4("Paramètres flowAI", style = "color: #2c3e50; font-weight: bold;"),
      fluidRow(
        column(6, numericInput("f_second_fraction", "second_fractionFR :", value = 0.1, min = 0.01, max = 0.5, step = 0.05)),
        column(6, numericInput("f_alpha", "alphaFR :", value = 0.01, min = 0.001, max = 0.1, step = 0.005))
      ),
      fluidRow(
        column(6, numericInput("f_max_cpt", "max_cptFS :", value = 3, min = 1)),
        column(6, numericInput("f_pen_value", "pen_valueFS :", value = 500, min = 10))
      ),
      fluidRow(
        column(6, numericInput("f_neg_values", "neg_valuesFM :", value = 1, min = 1, max = 3)),
        column(6, selectizeInput("f_ch_exclude", "ChExcludeFS :", 
                                 choices = c("FSC", "SSC"), selected = c("FSC", "SSC"), multiple = TRUE))
      )
    )
  })
  
  
  # ============================================================
  # 3. EXÉCUTION DU CALCUL (BOUTON VALIDER)
  # ============================================================
  
  observeEvent(input$valider_qc, {
    req(r_obj_R6(), input$choix_qc)
    
    obj <- r_obj_R6()
    choix_clean <- tolower(trimws(input$choix_qc))
    
    if (choix_clean == "peacoqc") {
      if (is.null(input$p_min_cells) || is.null(input$p_max_bins)) {
        showNotification("Configuration en cours, veuillez recliquer sur le bouton.", type = "warning")
        return(NULL)
      }
      
      withProgress(message = 'Calcul PeacoQC en cours...', value = 0, {
        incProgress(0.1, detail = "Configuration des paramètres...")
        
        reglages_peacoqc <- list(
          min_cells                 = input$p_min_cells,
          max_bins                  = input$p_max_bins,
          step                      = input$p_step,
          MAD                       = input$p_mad,
          IT_limit                  = input$p_it_limit,
          consecutive_bins          = input$p_consecutive_bins,
          remove_zeros              = as.logical(input$p_remove_zeros),
          force_IT                  = input$p_force_it,
          peak_removal              = input$p_peak_removal,
          min_nr_bins_peakdetection = input$p_min_nr_bins
        )
        
        incProgress(0.3, detail = "Exécution de PeacoQC sur l'ensemble des fichiers...")
        
        calcul_reussi <- FALSE
        tryCatch({
          message("🧬 Appel de obj$appliquer_peacoqc()...")
          # Utilisation de la méthode de ta classe unifiée pipeline_cytometry.R
          obj$appliquer_peacoqc(reglages_specifiques = reglages_peacoqc)
          calcul_reussi <- TRUE
        }, error = function(e) {
          message("❌ CRASH DANS LA MÉTHODE R6 : ", e$message)
          showNotification(paste("Erreur PeacoQC (R6) :", e$message), type = "error", duration = NULL)
        })
        
        if (!calcul_reussi) return(NULL)
        
        incProgress(0.8, detail = "Finalisation...")
        r_obj_R6(obj)
      })
      
      showNotification("PeacoQC appliqué avec succès sur l'ensemble des tubes", type = "message")
      
    } else if (choix_clean == "flowai") {
      if (is.null(input$f_second_fraction) || is.null(input$f_alpha)) {
        showNotification("Configuration en cours, veuillez recliquer sur le bouton.", type = "warning")
        return(NULL)
      }
      
      withProgress(message = 'Calcul flowAI en cours...', value = 0, {
        incProgress(0.2, detail = "Capture des paramètres de l'UI...")
        
        reglages_flowai <- list(
          second_fractionFR = input$f_second_fraction,
          alphaFR           = input$f_alpha,
          max_cptFS         = input$f_max_cpt,
          pen_valueFS       = input$f_pen_value,
          neg_valuesFM      = input$f_neg_values,
          ChExcludeFS       = input$f_ch_exclude
        )
        
        incProgress(0.4, detail = "Analyse du débit, des signaux et de la marge...")
        
        tryCatch({
          obj$appliquer_flowai(reglages_specifiques = reglages_flowai)
          r_obj_R6(obj)
          showNotification("✨ flowAI appliqué avec succès sur l'ensemble des tubes !", type = "message")
          
        }, error = function(e) {
          showNotification(paste("Erreur flowAI (R6) :", e$message), type = "error", duration = NULL)
        })
      })
    }
  })
  
  
  # ============================================================
  # 4. RENDU DU GRAPHIQUE DE CONTRÔLE QUALITÉ (ONGLET 1)
  # ============================================================
  
  output$plot_qc <- renderPlot({
    # Dépendances réactives : on écoute l'échantillon choisi ET l'objet global
    req(r_obj_R6(), input$echantillon_visu_qc)
    
    # On force le graphique à dépendre du clic sur le bouton de calcul
    input$valider_qc 
    
    obj <- r_obj_R6()
    ech_courant <- input$echantillon_visu_qc
    methode <- input$choix_qc
    
    # --------------------------------------------------------
    # CAS 1 : L'UTILISATEUR A SÉLECTIONNÉ PEACOQC
    # --------------------------------------------------------
    if (methode == "peacoqc") {
      if (is.null(obj$post_PeacoQC) || is.null(obj$post_PeacoQC[[ech_courant]])) {
        grid::grid.newpage()
        grid::grid.text(
          "📊 L'analyse PeacoQC n'a pas encore été exécutée.\nConfigurez vos paramètres à gauche puis cliquez sur\n'Lancer le Contrôle Qualité Global'.", 
          gp = grid::gpar(col = "gray50", fontsize = 14, fontface = "italic")
        )
        return(NULL)
      }
      
      withProgress(message = 'Génération du graphique PeacoQC...', value = 0.5, {
        p <- obj$visualiser_peacoqc(nom_echantillon = ech_courant, max_points = 10000)
      })
      return(p)
    }
    
    # --------------------------------------------------------
    # CAS 2 : L'UTILISATEUR A SÉLECTIONNÉ FLOWAI
    # --------------------------------------------------------
    if (methode == "flowai") {
      if (is.null(obj$post_flowAI) || is.null(obj$post_flowAI[[ech_courant]])) {
        grid::grid.newpage()
        grid::grid.text(
          "📊 L'analyse flowAI n'a pas encore été exécutée.\nConfigurez vos paramètres à gauche puis cliquez sur\n'Lancer le Contrôle Qualité Global'.", 
          gp = grid::gpar(col = "gray50", fontsize = 14, fontface = "italic")
        )
        return(NULL)
      }
      
      withProgress(message = 'Génération du graphique flowAI...', value = 0.5, {
        p <- obj$visualiser_flowai(nom_echantillon = ech_courant, max_points = 10000)
      })
      return(p)
    }
    
    # --------------------------------------------------------
    # CAS 3 : AUCUNE MÉTHODE SÉLECTIONNÉE ("none")
    # --------------------------------------------------------
    grid::grid.newpage()
    grid::grid.text(
      "ℹ️ Veuillez sélectionner une méthode de contrôle qualité (PeacoQC ou flowAI)\ndans le panneau de gauche pour commencer.", 
      gp = grid::gpar(col = "gray40", fontsize = 13)
    )
    return(NULL)
  })
  
  
  # ============================================================
  # 5. RENDU DU GRAPHIQUE ET LOGIQUE DÉBRIS (ONGLET 2)
  # ============================================================
  
  # Structure réactive pour enregistrer pas à pas les coordonnées cliquées par la souris
  coords_polygone_debris <- reactiveVal(data.frame(x = numeric(), y = numeric()))
  
  # Capturer les clics de souris sur le plot morphologique
  observeEvent(input$click_debris, {
    clic <- input$click_debris
    if (!is.null(clic)) {
      df_actuel <- coords_polygone_debris()
      nouvelle_ligne <- data.frame(x = clic$x, y = clic$y)
      coords_polygone_debris(rbind(df_actuel, nouvelle_ligne))
    }
  })
  
  # Bouton Reset pour effacer le tracé actuel
  observeEvent(input$reset_debris_gate, {
    coords_polygone_debris(data.frame(x = numeric(), y = numeric()))
    showNotification("Tracé du polygone effacé.", type = "message")
  })
  
  # Plot de morphologie (FSC vs SSC) avec superposition de la gate manuelle
  output$plot_debris <- renderPlot({
    # On récupère le flowFrame automatiquement aiguillé selon la source choisie (brutes/peacoqc/flowai)
    ff <- fcs_courant_debris()
    req(ff)
    
    # Extraction des données sous forme de dataframe pour ggplot
    donnees <- as.data.frame(flowCore::exprs(ff))
    
    # Détection automatique des canaux morphologiques FSC-A et SSC-A
    canal_fsc <- grep("FSC-A", colnames(donnees), ignore.case = TRUE, value = TRUE)[1]
    canal_ssc <- grep("SSC-A", colnames(donnees), ignore.case = TRUE, value = TRUE)[1]
    
    if(is.na(canal_fsc)) canal_fsc <- colnames(donnees)[1]
    if(is.na(canal_ssc)) canal_ssc <- colnames(donnees)[2]
    
    # Sous-échantillonnage pour conserver la fluidité
    if (nrow(donnees) > 10000) {
      set.seed(123)
      donnees <- donnees[sample(nrow(donnees), 10000), ]
    }
    
    # Création du graphique de base (densité morphologique)
    p <- ggplot2::ggplot(donnees, ggplot2::aes(x = .data[[canal_fsc]], y = .data[[canal_ssc]])) +
      ggpointdensity::geom_pointdensity(size = 0.5, alpha = 0.5) +
      ggplot2::scale_color_gradientn(colours = c("darkblue", "blue", "cyan", "greenyellow", "yellow", "darkorange", "red")) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste("Sélection Morphologique :", input$echantillon_visu_debris),
        subtitle = paste("Source :", input$source_donnees_debris, "| Événements affichés :", format(nrow(donnees), big.mark=" ")),
        x = canal_fsc,
        y = canal_ssc
      ) +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(face="bold"))
    
    # AJOUT DYNAMIQUE : Si l'utilisateur clique, on dessine le polygone rouge par-dessus
    pts <- coords_polygone_debris()
    if (nrow(pts) > 0) {
      p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y), color = "red", size = 2.5, inherit.aes = FALSE)
      if (nrow(pts) > 1) {
        p <- p + ggplot2::geom_path(data = pts, ggplot2::aes(x = x, y = y), color = "red", linewidth = 0.9, inherit.aes = FALSE)
        
        # CORRECTION DU GEOM_SEGMENT POUR EVITER LES AVERTISSEMENTS EN BOUCLE
        p <- p + ggplot2::geom_segment(
          data = data.frame(x = pts$x[nrow(pts)], y = pts$y[nrow(pts)], xend = pts$x[1], yend = pts$y[1]),
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          color = "red", linetype = "dashed", linewidth = 0.9, inherit.aes = FALSE
        )
      }
    }
    
    return(p)
  })
  
  # Bouton de validation globale de la gate Débris
  observeEvent(input$valider_debris, {
    req(r_obj_R6(), input$echantillon_visu_debris, input$source_donnees_debris)
    
    obj <- r_obj_R6()
    pts <- coords_polygone_debris()
    
    # Validation de sécurité sur l'UI avant l'envoi au R6
    if (is.null(pts) || nrow(pts) < 3) {
      showNotification("⚠️ Tracé invalide : Veuillez définir au moins 3 points sur le graphique pour former un polygone.", type = "error")
      return(NULL)
    }
    
    # Détection des canaux actuels pour l'affectation correcte des noms de colonnes dans la gate
    ff_temp <- fcs_courant_debris()
    noms_canaux <- colnames(flowCore::exprs(ff_temp))
    canal_fsc <- grep("FSC-A", noms_canaux, ignore.case = TRUE, value = TRUE)[1]
    canal_ssc <- grep("SSC-A", noms_canaux, ignore.case = TRUE, value = TRUE)[1]
    if(is.na(canal_fsc)) canal_fsc <- noms_canaux[1]
    if(is.na(canal_ssc)) canal_ssc <- noms_canaux[2]
    
    withProgress(message = "Application du filtre Débris (Batch global)...", value = 0, {
      tryCatch({
        incProgress(0.3, detail = "Transmission des coordonnées géométriques...")
        
        # APPEL REEL DE TA METHODE R6 COORDONNEE :
        # L'argument nom_echantillon est mis à NULL pour propager automatiquement cette gate sur toute la cohorte
        obj$retirer_les_debris(
          matrice_points   = pts,
          canal_x          = canal_fsc,
          canal_y          = canal_ssc,
          nom_echantillon  = NULL, 
          source_nettoyage = input$source_donnees_debris
        )
        
        # Sauvegarde et mise à jour de l'état réactif Shiny
        r_obj_R6(obj)
        showNotification("✅ Filtre morphologique appliqué avec succès sur tous les échantillons !", type = "message")
        
        # Optionnel : réinitialiser le calque des coordonnées après validation réussie
        coords_polygone_debris(data.frame(x = numeric(), y = numeric()))
        
      }, error = function(e) {
        showNotification(paste("Erreur lors du filtrage R6 :", e$message), type = "error", duration = NULL)
      })
    })
  })
  
  
  
  # 1. Génération dynamique du bouton "Bordures" (Toggle Button)
  output$ui_bouton_bordures <- renderUI({
    req(input$echantillon_visu_debris)
    nom_ech <- input$echantillon_visu_debris
    est_actif <- !is.null(status_bordures()[[nom_ech]]) && status_bordures()[[nom_ech]] == TRUE
    
    if (est_actif) {
      actionButton("toggle_bordures", "Bordures : ACTIVÉES", class = "btn-success", width = "100%", icon = shiny::icon("check"))
    } else {
      actionButton("toggle_bordures", "Bordures : DESACTIVÉES", class = "btn-default", width = "100%", icon = shiny::icon("times"))
    }
  })
  
  # Intercepter le clic sur le bouton Bordures pour inverser son statut (Toggle)
  observeEvent(input$toggle_bordures, {
    req(input$echantillon_visu_debris)
    nom_ech <- input$echantillon_visu_debris
    
    liste_etats <- status_bordures()
    # On inverse l'état actuel (si c'était TRUE ça devient FALSE, et inversement)
    actuel <- !is.null(liste_etats[[nom_ech]]) && liste_etats[[nom_ech]] == TRUE
    liste_etats[[nom_ech]] <- !actuel
    status_bordures(liste_etats)
  })
  
  # 2. LOGIQUE D'AIGUILLAGE EN CHAINE (L'ordre n'a plus d'importance)
  # Ce réactif contient les données prêtes à être AFFICHÉES sur le graphique de droite
  fcs_courant_debris <- reactive({
    req(r_obj_R6(), input$echantillon_visu_debris, input$source_donnees_debris)
    obj <- r_obj_R6()
    nom_ech <- input$echantillon_visu_debris
    
    # Étape A : On part de la source brute sélectionnée (Brute, PeacoQC ou flowAI)
    if (input$source_donnees_debris == "peacoqc" && !is.null(obj$post_PeacoQC[[nom_ech]])) {
      ff_base <- obj$post_PeacoQC[[nom_ech]]
    } else if (input$source_donnees_debris == "flowai" && !is.null(obj$post_flowAI[[nom_ech]])) {
      ff_base <- obj$post_flowAI[[nom_ech]]
    } else {
      ff_base <- if(length(obj$echantillons_compenses) > 0) obj$echantillons_compenses[[nom_ech]] else obj$fcs_compenses[[nom_ech]]
    }
    
    req(ff_base)
    
    # Étape B : Si l'interrupteur bordure est actif, on applique le retrait des bordures à la volée
    est_actif_bordures <- !is.null(status_bordures()[[nom_ech]]) && status_bordures()[[nom_ech]] == TRUE
    if (est_actif_bordures) {
      donnees_mat <- flowCore::exprs(ff_base)
      canal_fsc <- grep("FSC-A", colnames(donnees_mat), ignore.case = TRUE, value = TRUE)[1]
      canal_ssc <- grep("SSC-A", colnames(donnees_mat), ignore.case = TRUE, value = TRUE)[1]
      if(is.na(canal_fsc)) canal_fsc <- colnames(donnees_mat)[1]
      if(is.na(canal_ssc)) canal_ssc <- colnames(donnees_mat)[2]
      
      max_fsc <- max(donnees_mat[, canal_fsc]); min_fsc <- min(donnees_mat[, canal_fsc])
      max_ssc <- max(donnees_mat[, canal_ssc]); min_ssc <- min(donnees_mat[, canal_ssc])
      
      indices_valides <- which(
        donnees_mat[, canal_fsc] > min_fsc & donnees_mat[, canal_fsc] < (max_fsc - 100) &
          donnees_mat[, canal_ssc] > min_ssc & donnees_mat[, canal_ssc] < (max_ssc - 100)
      )
      ff_base <- ff_base[indices_valides, ]
    }
    
    return(ff_base)
  })
  
  # 3. INTERACTION RATON/SOURIS : CAPTURE DES CLICS
  observeEvent(input$click_debris, {
    clic <- input$click_debris
    req(clic)
    df_actuel <- coords_polygone_debris()
    coords_polygone_debris(rbind(df_actuel, data.frame(x = clic$x, y = clic$y)))
  })
  
  # 4. ACTION : EFFACER LE TRACÉ
  observeEvent(input$reset_debris_gate, {
    coords_polygone_debris(data.frame(x = numeric(), y = numeric()))
  })
  
  # 5. ACTION : COPIER LE PREMIER TRACÉ SUR TOUS LES ÉCHANTILLONS (GABARIT)
  gabarit_global_debris <- reactiveVal(NULL)
  
  # Modifie l'observeur correspondant :
  observeEvent(input$appliquer_tous_debris, {
    pts <- coords_polygone_debris()
    if (nrow(pts) < 3) {
      showNotification("⚠️ Veuillez d'abord dessiner un polygone complet (min. 3 points).", type = "warning")
      return(NULL)
    }
    gabarit_global_debris(pts) # Stockage dans le réactif Shiny plutôt que dans le R6 verrouillé
    showNotification("🎯 Ce polygone a été défini comme modèle pour toute la cohorte !", type = "default")
  })
  
  # Synchronisation automatique lors du changement d'échantillon inspecté
  observeEvent(input$echantillon_visu_debris, {
    obj <- r_obj_R6()
    req(obj)
    nom_ech <- input$echantillon_visu_debris
    
    # Si ce tube a déjà une gate enregistrée dans l'objet R6
    if (!is.null(obj$gate_debris[[nom_ech]])) {
      gate_obj <- obj$gate_debris[[nom_ech]]
      
      # CORRECTION ICI : Extraction propre et universelle des coordonnées
      if (.hasSlot(gate_obj, "boundaries")) {
        matrice_gate <- gate_obj@boundaries
      } else if (.hasSlot(gate_obj, ".gate")) {
        matrice_gate <- gate_obj@.gate
      } else {
        # Fallback si l'objet est structuré différemment
        matrice_gate <- flowCore::boundaries(gate_obj)
      }
      
      # Mise à jour des coordonnées pour le graphique
      coords_polygone_debris(data.frame(x = matrice_gate[,1], y = matrice_gate[,2]))
    } 
    # Sinon, si un gabarit global (dessiné sur un autre tube) existe, on le charge
    else if (!is.null(gabarit_global_debris())) { 
      coords_polygone_debris(gabarit_global_debris())
    } 
    # Sinon, on vide le graphique pour repartir à zéro
    else {
      coords_polygone_debris(data.frame(x = numeric(), y = numeric()))
    }
  })
  
  # 6. RENDU GRAPHIQUE ULTRA-DYNAMIQUE
  output$plot_debris <- renderPlot({
    ff <- fcs_courant_debris() # Il contient déjà le filtre bordure s'il est activé !
    req(ff)
    
    donnees <- as.data.frame(flowCore::exprs(ff))
    canal_fsc <- grep("FSC-A", colnames(donnees), ignore.case = TRUE, value = TRUE)[1]
    canal_ssc <- grep("SSC-A", colnames(donnees), ignore.case = TRUE, value = TRUE)[1]
    if(is.na(canal_fsc)) canal_fsc <- colnames(donnees)[1]
    if(is.na(canal_ssc)) canal_ssc <- colnames(donnees)[2]
    
    if (nrow(donnees) > 10000) {
      set.seed(123)
      donnees <- donnees[sample(nrow(donnees), 10000), ]
    }
    
    est_actif_bordures <- !is.null(status_bordures()[[input$echantillon_visu_debris]]) && status_bordures()[[input$echantillon_visu_debris]] == TRUE
    
    p <- ggplot2::ggplot(donnees, ggplot2::aes(x = .data[[canal_fsc]], y = .data[[canal_ssc]])) +
      ggpointdensity::geom_pointdensity(size = 0.5, alpha = 0.5) +
      ggplot2::scale_color_gradientn(colours = c("darkblue", "blue", "cyan", "greenyellow", "yellow", "darkorange", "red")) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste("FSC-A vs SSC-A :", input$echantillon_visu_debris),
        subtitle = paste("Filtre bordures :", if(est_actif_bordures) "ACTIF (Événements exclus)" else "INACTIF"),
        x = canal_fsc, y = canal_ssc
      ) + ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(face="bold"))
    
    # Dessin du polygone par-dessus
    pts <- coords_polygone_debris()
    if (nrow(pts) > 0) {
      p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y), color = "red", size = 2.5)
      if (nrow(pts) > 1) {
        p <- p + ggplot2::geom_path(data = pts, ggplot2::aes(x = x, y = y), color = "red", size = 0.9)
        p <- p + ggplot2::geom_segment(
          ggplot2::aes(x = pts$x[nrow(pts)], y = pts$y[nrow(pts)], xend = pts$x[1], yend = pts$y[1]),
          color = "red", linetype = "dashed", size = 0.9
        )
      }
    }
    return(p)
  })
  
  # 7. ACTION : VALIDER LA SÉLECTION POUR CE TUBE UNIQUE
  # Dans server.R / app.R
  observeEvent(input$valider_debris, {
    req(r_obj_R6(), input$echantillon_visu_debris)
    obj <- r_obj_R6()
    pts <- coords_polygone_debris()
    nom_ech <- input$echantillon_visu_debris
    
    if (is.null(pts) || nrow(pts) < 3) {
      showNotification("⚠️ Veuillez tracer un polygone valide (min. 3 points) avant de valider.", type = "error")
      return(NULL)
    }
    
    ff_a_filtrer <- fcs_courant_debris() 
    noms_canaux <- colnames(flowCore::exprs(ff_a_filtrer))
    canal_fsc <- grep("FSC-A", noms_canaux, ignore.case = TRUE, value = TRUE)[1]
    canal_ssc <- grep("SSC-A", noms_canaux, ignore.case = TRUE, value = TRUE)[1]
    if(is.na(canal_fsc)) canal_fsc <- noms_canaux[1]
    if(is.na(canal_ssc)) canal_ssc <- noms_canaux[2]
    
    withProgress(message = "Validation de l'échantillon...", {
      tryCatch({
        if (is.null(obj$post_bordures)) obj$post_bordures <- list()
        obj$post_bordures[[nom_ech]] <- ff_a_filtrer
        obj$retirer_les_debris(
          matrice_points      = pts,
          canal_x             = canal_fsc,
          canal_y             = canal_ssc,
          nom_echantillon     = nom_ech,
          flowframe_a_filtrer = ff_a_filtrer  
        )
        r_obj_R6(obj)
        showNotification(paste("✅ Nettoyage validé pour le tube :", nom_ech), type = "message")
      }, error = function(e) {
        showNotification(paste("Erreur lors du filtrage R6 :", e$message), type = "error", duration = NULL)
      })
    })
  })
  
  # 8. RENDU DU TABLEAU RÉCAPITULATIF (S'ADAPTE À TOUS LES ORDRES)
  output$table_stats_debris <- renderTable({
    req(r_obj_R6(), input$echantillon_visu_debris)
    obj <- r_obj_R6()
    nom_ech <- input$echantillon_visu_debris
    
    # A. Données initiales absolues arrivant à cet onglet
    if (input$source_donnees_debris == "peacoqc") ff_init <- obj$post_PeacoQC[[nom_ech]]
    else if (input$source_donnees_debris == "flowai") ff_init <- obj$post_flowAI[[nom_ech]]
    else ff_init <- if(length(obj$echantillons_compenses) > 0) obj$echantillons_compenses[[nom_ech]] else obj$fcs_compenses[[nom_ech]]
    
    req(ff_init)
    nb_initial <- nrow(flowCore::exprs(ff_init))
    
    # B. Événements après application du filtre Bordures (si actif)
    est_actif_bordures <- !is.null(status_bordures()[[nom_ech]]) && status_bordures()[[nom_ech]] == TRUE
    
    # On recalcule rapidement le volume théorique avec bordures retirées pour le tableau
    if (est_actif_bordures) {
      donnees_mat <- flowCore::exprs(ff_init)
      canal_fsc <- grep("FSC-A", colnames(donnees_mat), ignore.case = TRUE, value = TRUE)[1]
      canal_ssc <- grep("SSC-A", colnames(donnees_mat), ignore.case = TRUE, value = TRUE)[1]
      max_fsc <- max(donnees_mat[, canal_fsc]); min_fsc <- min(donnees_mat[, canal_fsc])
      max_ssc <- max(donnees_mat[, canal_ssc]); min_ssc <- min(donnees_mat[, canal_ssc])
      nb_bordures <- length(which(donnees_mat[, canal_fsc] > min_fsc & donnees_mat[, canal_fsc] < (max_fsc-100) & donnees_mat[, canal_ssc] > min_ssc & donnees_mat[, canal_ssc] < (max_ssc-100)))
    } else {
      nb_bordures <- nb_initial
    }
    
    # C. Événements après application définitive du polygone Débris
    nb_debris <- if(!is.null(obj$post_debris[[nom_ech]])) nrow(flowCore::exprs(obj$post_debris[[nom_ech]])) else NA
    
    pct_bordures <- round((nb_bordures / nb_initial) * 100, 1)
    pct_global_debris <- if(!is.na(nb_debris)) round((nb_debris / nb_initial) * 100, 1) else "-"
    
    data.frame(
      Étape = c("1. Données initiales reçues (Post-CQ)", "2. Si retrait des Bordures uniquement", "3. Résultat final après Polygone + Bordures"),
      `Événements Restants` = c(format(nb_initial, big.mark=" "), format(nb_bordures, big.mark=" "), if(!is.na(nb_debris)) format(nb_debris, big.mark=" ") else "En attente du clic 'Valider'..."),
      `% Conservation global` = c("100 %", paste(pct_bordures, "%"), paste(pct_global_debris, "%")),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, align = "c")
  
  
}



shinyApp(ui, server)