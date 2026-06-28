library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(shinyjs)
library(grid)
library(gridExtra)
library(flowCore)

source("pipeline_cytometrie.R")
options(shiny.maxRequestSize = 1024 * 1024 * 1024)

# ============================================================
# CONSTANTES GLOBALES
# ============================================================

CANAUX_CONNUS <- c(
  "", "FITC-A", "Alexa Fluor 488-A", "Alexa Fluor 700-A",
  "PE-A", "PE-Cy5-A", "PE-Cy7-A", "PerCP-Cy5-5-A", "GFP-A",
  "APC-A", "APC-R700-A", "APC-H7-A",
  "BV421-A", "AmCyan-A", "BV510-A", "BV605-A", "BV650-A",
  "Pacific Blue-A", "V450-A", "BV711-A", "BV786-A"
)

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


# ================================================================
# HELPERS
# ================================================================

get_labels_from_fcs <- function(fcs) {
  if (is.null(fcs)) return(list())
  pd <- tryCatch(flowCore::pData(flowCore::parameters(fcs)), error = function(e) NULL)
  if (is.null(pd)) return(list())
  setNames(lapply(seq_len(nrow(pd)), function(i) {
    d <- pd$desc[i]
    if (!is.null(d) && !is.na(d) && nchar(trimws(d)) > 0)
      paste0(pd$name[i], " | ", d)
    else pd$name[i]
  }), pd$name)
}

label_canal <- function(fcs, canal) {
  lbl <- get_labels_from_fcs(fcs)[[canal]]
  if (is.null(lbl)) canal else lbl
}

canaux_fluo_fcs <- function(fcs) {
  cx <- colnames(flowCore::exprs(fcs))
  cx[!grepl("^time$|^fsc|^ssc", cx, ignore.case = TRUE)]
}

choix_avec_labels <- function(fcs, canaux = NULL) {
  if (is.null(canaux)) canaux <- canaux_fluo_fcs(fcs)
  lbs <- get_labels_from_fcs(fcs)
  setNames(canaux, sapply(canaux, function(c) if (!is.null(lbs[[c]])) lbs[[c]] else c))
}

# ================================================================
# SERVER
# ================================================================

server <- function(input, output, session) {
  
  # ---------------------------------------------------------------
  # ÉTAT RÉACTIF CENTRAL
  # Tout est conservé dans rv — jamais effacé lors des changements
  # d'onglet.  Les triggers (entiers) forcent Shiny à ré-évaluer
  # les outputs qui dépendent de mutations internes d'un objet R6.
  # ---------------------------------------------------------------
  rv <- reactiveValues(
    # --- Données brutes ---
    carrot              = NULL,
    r_df_mono           = NULL,
    r_df_ech            = NULL,
    r_dictionnaire      = c(
      "", "FITC-A", "PE-A", "PE-Cy5-A", "PE-Cy7-A",
      "PerCP-Cy5-5-A", "APC-A", "APC-R700-A", "APC-H7-A",
      "BV421-A", "BV510-A", "BV605-A", "BV650-A",
      "BV711-A", "BV786-A", "AmCyan-A", "Pacific Blue-A",
      "V450-A", "GFP-A", "Alexa Fluor 488-A", "Alexa Fluor 700-A"
    ),
    groupes             = list(),
    
    # --- Historique undo/redo matrices (par échantillon) ---
    historique_matrices = list(),
    index_versions      = list(),
    
    # --- Triggers d'invalidation R6 → Shiny ---
    trig_fcs_charge     = 0L,
    trig_trans_applique = 0L,
    trig_gate_sauve     = 0L,
    trig_spillover      = 0L,
    trig_compensation   = 0L,
    trig_matrice_edit   = 0L,
    
    # --- Mode ---
    mode_cyto           = "Conventionnel",
    deja_compense       = FALSE   # TRUE si l'utilisateur indique que ses données sont déjà compensées
  )
  
  # --- Raccourcis réactifs ---
  carrot_ok      <- reactive({ rv$trig_fcs_charge;    !is.null(rv$carrot) })
  mono_ok        <- reactive({ rv$trig_fcs_charge;    carrot_ok() && length(rv$carrot$tubes_monomarques) > 0 })
  mono_trans_ok  <- reactive({ rv$trig_trans_applique; carrot_ok() && !is.null(rv$carrot$monomarques_trans) && length(rv$carrot$monomarques_trans) > 0 })
  spillover_ok   <- reactive({ rv$trig_spillover;      carrot_ok() && !is.null(rv$carrot$S_matrix) })
  comp_ok        <- reactive({ rv$trig_compensation;   carrot_ok() && length(rv$carrot$echantillons_traites) > 0 })
  
  # ================================================================
  #  ONGLET 1 — IMPORTATION
  # ================================================================
  
  # --- Chargement Monomarqués ---
  observeEvent(input$load_mono, {
    f             <- input$load_mono
    noms_proposes <- gsub("\\.fcs$", "", f$name, ignore.case = TRUE)
    dict          <- rv$r_dictionnaire
    
    canaux_choisis <- sapply(noms_proposes, function(np) {
      m <- dict[nchar(dict) > 0 & (
        sapply(dict, function(dc) grepl(dc, np, ignore.case = TRUE) || grepl(np, dc, ignore.case = TRUE))
      )]
      if (length(m) > 0) m[1] else ""
    })
    
    rv$r_df_mono <- data.frame(
      fichier = f$name,
      canal   = canaux_choisis,
      type    = ifelse(grepl("unstained", noms_proposes, ignore.case = TRUE), "Unstained", "Monomarque"),
      chemin  = f$datapath,
      stringsAsFactors = FALSE
    )
  })
  
  # --- Chargement Échantillons (lecture rapide des en-têtes) ---
  observeEvent(input$load_ech, {
    f <- input$load_ech
    withProgress(message = "Lecture des métadonnées FCS...", value = 0, {
      rows <- lapply(seq_len(nrow(f)), function(i) {
        incProgress(1 / nrow(f), detail = paste("Fichier", i, "/", nrow(f)))
        ch <- f$datapath[i]
        nb_events <- 0; tube_name <- "Inconnu"; exp_name <- "Inconnu"
        cytometre <- "Inconnu"; date_acq <- "Inconnu"
        tryCatch({
          hdr <- flowCore::read.FCS(ch, dataset = 1, which.lines = 1,
                                    transformation = FALSE, truncate_max_range = FALSE)
          kw        <- flowCore::keyword(hdr)
          nb_events <- ifelse(!is.null(kw[["$TOT"]]),           as.numeric(kw[["$TOT"]]), 0)
          tube_name <- ifelse(!is.null(kw[["TUBE NAME"]]),       kw[["TUBE NAME"]],        "Inconnu")
          exp_name  <- ifelse(!is.null(kw[["EXPERIMENT NAME"]]), kw[["EXPERIMENT NAME"]],  "Inconnu")
          cytometre <- ifelse(!is.null(kw[["$CYT"]]),            kw[["$CYT"]],             "Inconnu")
          date_acq  <- ifelse(!is.null(kw[["$DATE"]]),           kw[["$DATE"]],            "Inconnu")
        }, error = function(e) NULL)
        sz <- file.info(ch)$size
        poids <- if (is.na(sz)) "?" else if (sz >= 1e9) paste0(round(sz/1e9,2)," Go") else paste0(round(sz/1e6,1)," Mo")
        data.frame(tube_name=tube_name, fichier=f$name[i], chemin=ch,
                   nb_events=format(nb_events, big.mark="\u00a0"),
                   exp_name=exp_name, cytometre=cytometre, poids=poids, date=date_acq,
                   stringsAsFactors=FALSE)
      })
      rv$r_df_ech <- do.call(rbind, rows)
    })
  })
  
  # --- Table Monomarqués ---
  output$table_mono <- renderDT({
    req(rv$r_df_mono)
    df   <- rv$r_df_mono
    CANAUX <- rv$r_dictionnaire
    
    sel_canal <- sapply(seq_len(nrow(df)), function(i) {
      v    <- df$canal[i]
      opts <- paste0(sapply(CANAUX, function(o) {
        sel   <- if (identical(o, v)) " selected" else ""
        label <- if (o == "") "-- Choisir --" else o
        sprintf('<option value="%s"%s>%s</option>', o, sel, label)
      }), collapse = "")
      # Champ texte libre + datalist pour la recherche + select
      datalist_id <- paste0("dl_canal_", i)
      datalist_opts <- paste0(sapply(CANAUX[nchar(CANAUX)>0], function(o)
        sprintf('<option value="%s">', o)), collapse="")
      paste0(
        '<div style="display:flex;gap:4px;">',
        sprintf('<input type="text" class="dt-datalist-input" list="%s" id="txt_canal_%d"
          value="%s" placeholder="Rechercher ou saisir..."
          style="flex:1;"
          oninput=\'Shiny.setInputValue("change_table_canal",{row:%d,val:this.value},{priority:"event"})\'>',
                datalist_id, i, v, i),
        sprintf('<datalist id="%s">%s</datalist>', datalist_id, datalist_opts),
        '</div>'
      )
    })
    
    sel_type <- sapply(seq_len(nrow(df)), function(i) {
      opt_m <- if (df$type[i] == "Monomarque") " selected" else ""
      opt_u <- if (df$type[i] == "Unstained")  " selected" else ""
      sprintf('<select class="dt-select" onchange=\'Shiny.setInputValue("change_table_type",{row:%d,val:this.value},{priority:"event"})\'>
        <option value="Monomarque"%s>Monomarque</option>
        <option value="Unstained"%s>Unstained</option>
      </select>', i, opt_m, opt_u)
    })
    
    datatable(
      data.frame(Fichier=df$fichier, Canal=sel_canal, Type=sel_type, stringsAsFactors=FALSE),
      escape=FALSE, rownames=FALSE, selection="none",
      options=list(dom="t", paging=FALSE, ordering=FALSE,
                   preDrawCallback=JS("function(){Shiny.unbindAll(this.api().table().node());}"),
                   drawCallback=JS("function(){Shiny.bindAll(this.api().table().node());}"))
    )
  })
  
  observeEvent(input$change_table_canal, {
    info <- input$change_table_canal; req(rv$r_df_mono)
    rv$r_df_mono[info$row, "canal"] <- info$val
    if (nchar(info$val)>0 && !(info$val %in% rv$r_dictionnaire))
      rv$r_dictionnaire <- c(rv$r_dictionnaire, info$val)
  })
  observeEvent(input$change_table_type, {
    info <- input$change_table_type; req(rv$r_df_mono)
    rv$r_df_mono[info$row, "type"] <- info$val
  })
  
  # --- Table Échantillons ---
  output$table_ech <- renderDT({
    req(rv$r_df_ech)
    df  <- rv$r_df_ech
    cols <- c("tube_name","fichier","nb_events","exp_name","cytometre","poids","date")
    for (col in cols) if (!(col %in% names(df))) df[[col]] <- "—"
    datatable(
      df[, cols],
      colnames=c("Nom (modifiable)","Fichier","Évènements","Expérience","Cytomètre","Volume","Date"),
      editable=list(target="cell", disable=list(columns=1:6)),
      rownames=FALSE, selection="none",
      options=list(dom="t", paging=FALSE, scrollX=TRUE)
    )
  })
  observeEvent(input$table_ech_cell_edit, {
    info <- input$table_ech_cell_edit; req(rv$r_df_ech)
    rv$r_df_ech[info$row, "tube_name"] <- info$value
  })
  
  # --- Sélecteur cytomètre + option déjà compensé ---
  output$ui_cyto_select <- renderUI({
    tagList(
      selectInput("cyto_type_global", "Technologie Cytomètre :",
                  choices=c("Conventionnel","Spectral"), selected=rv$mode_cyto),
      hr(),
      radioButtons("fichiers_deja_compenses", "Vos fichiers sont-ils déjà compensés ?",
                   choices=c("Non — je veux faire la compensation"="non",
                             "Oui — données déjà compensées"="oui"),
                   selected=if(rv$deja_compense)"oui" else "non")
    )
  })
  
  observeEvent(input$cyto_type_global, {
    rv$mode_cyto <- input$cyto_type_global
    # Bloque/débloque les onglets compensation / unmixing
    if (input$cyto_type_global == "Conventionnel") {
      shinyjs::addClass(selector="li:has(a[data-value='unmixing_tab'])",  class="disabled-tab")
      shinyjs::removeClass(selector="li:has(a[data-value='compensation_tab'])", class="disabled-tab")
    } else {
      shinyjs::addClass(selector="li:has(a[data-value='compensation_tab'])", class="disabled-tab")
      shinyjs::removeClass(selector="li:has(a[data-value='unmixing_tab'])",  class="disabled-tab")
    }
  }, ignoreNULL=TRUE, ignoreInit=FALSE)
  
  observeEvent(input$fichiers_deja_compenses, {
    rv$deja_compense <- input$fichiers_deja_compenses == "oui"
  })
  
  # --- Matrice marqueurs (onglet 2) ---
  matrice_marqueurs_rv <- reactiveVal(NULL)
  observeEvent(rv$trig_fcs_charge, {
    req(carrot_ok(), length(rv$carrot$echantillons) > 0)
    noms_ech <- names(rv$carrot$echantillons)
    fcs_ref  <- rv$carrot$echantillons[[1]]
    cx_fluo  <- canaux_fluo_fcs(fcs_ref)
    lbs_ref  <- get_labels_from_fcs(fcs_ref)
    mat <- matrix("", nrow=length(noms_ech), ncol=length(cx_fluo),
                  dimnames=list(noms_ech, cx_fluo))
    for (nom in noms_ech) {
      lbs <- get_labels_from_fcs(rv$carrot$echantillons[[nom]])
      for (cx in cx_fluo) {
        lbl <- lbs[[cx]]
        mat[nom, cx] <- if (!is.null(lbl) && grepl(" \\| ", lbl))
          strsplit(lbl," \\| ")[[1]][2] else ""
      }
    }
    colnames(mat) <- sapply(cx_fluo, function(cx) { l <- lbs_ref[[cx]]; if (!is.null(l)) l else cx })
    matrice_marqueurs_rv(as.data.frame(mat, stringsAsFactors=FALSE))
  })
  output$table_matrice_marqueurs <- renderDT({
    req(matrice_marqueurs_rv())
    datatable(matrice_marqueurs_rv(), editable=list(target="cell"),
              rownames=TRUE, selection="none",
              options=list(scrollX=TRUE, pageLength=10, dom="t"))
  })
  observeEvent(input$table_matrice_marqueurs_cell_edit, {
    info <- input$table_matrice_marqueurs_cell_edit
    df <- matrice_marqueurs_rv(); df[info$row, info$col] <- info$value; matrice_marqueurs_rv(df)
  })
  
  # --- Groupes biologiques ---
  observeEvent(input$add_group, {
    req(input$group_name, nchar(trimws(input$group_name))>0)
    nom <- trimws(input$group_name)
    if (is.null(rv$groupes[[nom]])) {
      rv$groupes[[nom]] <- character(0)
      updateTextInput(session, "group_name", value="")
    } else showNotification("Ce groupe existe déjà.", type="warning")
  })
  output$ui_assignation_groupes <- renderUI({
    req(rv$r_df_ech, length(rv$groupes)>0)
    lapply(rv$r_df_ech$tube_name, function(tube) {
      current <- names(rv$groupes)[sapply(rv$groupes, function(g) tube %in% g)]
      selectInput(paste0("assign_",make.names(tube)), label=tube,
                  choices=c("(Non assigné)"="", names(rv$groupes)),
                  selected=if(length(current)) current[1] else "")
    })
  })
  output$table_resume_groupes <- renderDT({
    req(length(rv$groupes)>0, rv$r_df_ech)
    rows <- lapply(names(rv$groupes), function(g) {
      t <- rv$groupes[[g]]
      data.frame(Groupe=g,
                 Échantillons=if(length(t)) paste(t,collapse=", ") else "—",
                 N=length(t), stringsAsFactors=FALSE)
    })
    datatable(do.call(rbind,rows), rownames=FALSE, options=list(dom="t", pageLength=20))
  })
  observeEvent(input$save_cohort, {
    req(rv$r_df_ech, length(rv$groupes)>0)
    for (tube in rv$r_df_ech$tube_name) {
      grp <- input[[paste0("assign_",make.names(tube))]]
      if (!is.null(grp) && nchar(grp)>0) {
        for (g in names(rv$groupes)) rv$groupes[[g]] <- rv$groupes[[g]][rv$groupes[[g]]!=tube]
        rv$groupes[[grp]] <- c(rv$groupes[[grp]], tube)
      }
    }
    showNotification("Cohorte enregistrée.", type="message")
  })
  
  # --- Initialisation R6 (sans redirection forcée) ---
  observeEvent(input$init_r6, {
    req(rv$r_df_mono, rv$r_df_ech, input$cyto_type_global)
    withProgress(message="Chargement des fichiers FCS...", value=0.5, {
      tryCatch({
        obj <- CARROT$new(
          df_monomarques  = rv$r_df_mono,
          df_echantillons = rv$r_df_ech[, c("tube_name","chemin")],
          mode            = input$cyto_type_global
        )
        obj$charger_fcs()
        rv$r_df_ech        <- obj$chemins_echantillons
        rv$carrot          <- obj
        rv$mode_cyto       <- input$cyto_type_global
        rv$deja_compense   <- input$fichiers_deja_compenses == "oui"
        rv$trig_fcs_charge <- rv$trig_fcs_charge + 1L
        
        # Si déjà compensé : sauvegarder les FCS comme "traités" directement
        if (rv$deja_compense) {
          for (nom in names(obj$echantillons))
            obj$echantillons_traites[[nom]] <- obj$echantillons[[nom]]
          rv$trig_compensation <- rv$trig_compensation + 1L
          showNotification(
            paste0("✔ ", length(obj$echantillons),
                   " échantillon(s) chargés (déjà compensés, prêts pour la suite)."),
            type="message")
        } else {
          showNotification(
            paste0("✔ ", length(obj$echantillons), " échantillon(s) et ",
                   length(obj$tubes_monomarques), " contrôle(s) chargés."),
            type="message")
        }
        # L'utilisateur reste libre de naviguer où il veut
      }, error=function(e) {
        showNotification(paste("Erreur d'initialisation :", e$message),
                         type="error", duration=NULL)
      })
    })
  })
  
  
  # ================================================================
  #  ONGLET 2 — COMPENSATION
  # ================================================================
  
  # ---- Transformation ----------------------------------------
  
  output$ui_trans_file <- renderUI({
    rv$trig_fcs_charge
    req(mono_ok())
    selectInput("trans_file_sel", "Tube de prévisualisation :",
                choices=names(rv$carrot$tubes_monomarques))
  })
  
  canaux_trans_r <- reactive({
    rv$trig_fcs_charge; req(mono_ok())
    fcs <- rv$carrot$tubes_monomarques[[1]]
    cx  <- colnames(flowCore::exprs(fcs))
    cx[!grepl("^time$", cx, ignore.case=TRUE)]
  })
  
  output$ui_trans_cx <- renderUI({
    req(canaux_trans_r())
    fcs   <- rv$carrot$tubes_monomarques[[1]]
    choix <- choix_avec_labels(fcs, canaux_trans_r())
    selectInput("trans_cx", "Axe X :", choices=choix)
  })
  output$ui_trans_cy <- renderUI({
    req(canaux_trans_r())
    fcs   <- rv$carrot$tubes_monomarques[[1]]
    choix <- choix_avec_labels(fcs, canaux_trans_r())
    selectInput("trans_cy", "Axe Y :", choices=choix,
                selected=if(length(choix)>1) choix[2] else choix[1])
  })
  
  # Application arcsinh (bouton explicite)
  observeEvent(input$btn_apply_trans, {
    req(mono_ok(), input$cofacteur)
    cofac <- as.numeric(input$cofacteur)
    withProgress(message="Application de la transformation arcsinh...", value=0.5, {
      tryCatch({
        fcs_ref <- rv$carrot$tubes_monomarques[[1]]
        cx_fluo <- canaux_fluo_fcs(fcs_ref)
        tl <- setNames(lapply(cx_fluo, function(cx)
          flowCore::arcsinhTransform(paste0("arcsinh_",cx), a=1/cofac, b=1, c=0)), cx_fluo)
        tf <- flowCore::transformList(names(tl), tl)
        rv$carrot$trans_list <- tl
        rv$carrot$monomarques_trans <- lapply(rv$carrot$tubes_monomarques,
                                              function(fcs) flowCore::transform(fcs, tf))
        names(rv$carrot$monomarques_trans) <- names(rv$carrot$tubes_monomarques)
        for (cx in cx_fluo)
          rv$carrot$config_transformations[[cx]] <- list(type="arcsinh", cofacteur=cofac)
        rv$trig_trans_applique <- rv$trig_trans_applique + 1L
        showNotification(paste0("✔ Arcsinh appliqué (cofacteur = ", cofac, ")."),
                         type="message", duration=4)
      }, error=function(e) {
        showNotification(paste("Erreur transformation :", e$message), type="error")
      })
    })
  })
  
  # Prévisualisation live (se met à jour dès que cofacteur, cx, cy ou fichier changent)
  # sans relancer l'application sur les contrôles — juste un aperçu
  plot_trans_data <- reactive({
    rv$trig_fcs_charge
    req(mono_ok(), input$trans_cx, input$trans_cy, input$cofacteur, input$trans_file_sel)
    nom   <- input$trans_file_sel
    cx    <- input$trans_cx
    cy    <- input$trans_cy
    cofac <- as.numeric(input$cofacteur)
    req(nom %in% names(rv$carrot$tubes_monomarques))
    fcs_brut <- rv$carrot$tubes_monomarques[[nom]]
    canaux_ok <- colnames(flowCore::exprs(fcs_brut))
    req(cx %in% canaux_ok, cy %in% canaux_ok)
    
    tf_local <- flowCore::transformList(
      c(cx,cy),
      list(flowCore::arcsinhTransform("px", a=1/cofac, b=1, c=0),
           flowCore::arcsinhTransform("py", a=1/cofac, b=1, c=0))
    )
    fcs_t <- tryCatch(flowCore::transform(fcs_brut, tf_local), error=function(e) NULL)
    req(!is.null(fcs_t))
    mat <- as.data.frame(flowCore::exprs(fcs_t)[, c(cx,cy), drop=FALSE])
    n   <- min(nrow(mat), 5000)
    set.seed(rv$carrot$seed)
    list(
      mat   = mat[sample(seq_len(nrow(mat)), n), ],
      lbl_x = label_canal(fcs_brut, cx),
      lbl_y = label_canal(fcs_brut, cy),
      cofac = cofac,
      nom   = nom,
      cx    = cx, cy = cy
    )
  })
  
  output$plot_transformation <- renderPlot({
    d <- plot_trans_data(); req(!is.null(d))
    ggplot(d$mat, aes(x=.data[[d$cx]], y=.data[[d$cy]])) +
      ggpointdensity::geom_pointdensity(size=0.3, alpha=0.6) +
      scale_color_gradientn(colours=c("darkblue","blue","cyan","greenyellow","yellow","darkorange","red")) +
      theme_bw() +
      theme(legend.position="none", aspect.ratio=1,
            plot.title=element_text(face="bold")) +
      labs(title=paste("Prévisualisation arcsinh —", d$nom),
           subtitle=paste0("Cofacteur = ", d$cofac),
           x=d$lbl_x, y=d$lbl_y)
  })
  
  
  # ---- Gating interactif ------------------------------------
  
  output$ui_gate_canal <- renderUI({
    rv$trig_trans_applique; req(mono_trans_ok())
    canaux <- names(rv$carrot$monomarques_trans)
    canaux <- canaux[canaux != "TUBE_UNSTAINED"]
    req(length(canaux)>0)
    fcs_ref <- rv$carrot$tubes_monomarques[[canaux[1]]]
    lbs     <- get_labels_from_fcs(fcs_ref)
    choix   <- setNames(canaux, sapply(canaux, function(cx) { l <- lbs[[cx]]; if(!is.null(l)) l else cx }))
    selectInput("gate_canal", "Tube monomarqué :", choices=choix)
  })
  
  output$ui_gate_use_unstained_wrapper <- renderUI({
    rv$trig_trans_applique; req(mono_trans_ok())
    if ("TUBE_UNSTAINED" %in% names(rv$carrot$monomarques_trans))
      checkboxInput("gate_use_unstained","Référence négative : tube Unstained", value=TRUE)
    else helpText("Aucun Unstained détecté — population négative interne utilisée.")
  })
  
  # Bornes étendues pour les sliders : on laisse une marge pour permettre
  # à l'utilisateur de déplacer librement les curseurs au-delà des données
  bornes_canal <- reactive({
    rv$trig_trans_applique; req(mono_trans_ok(), input$gate_canal)
    canal <- input$gate_canal
    req(canal %in% names(rv$carrot$monomarques_trans))
    vals <- flowCore::exprs(rv$carrot$monomarques_trans[[canal]])[, canal]
    rng  <- range(vals, na.rm=TRUE)
    marge <- (rng[2] - rng[1]) * 0.3  # 30% de marge de chaque côté
    c(floor((rng[1]-marge)*10)/10, ceiling((rng[2]+marge)*10)/10)
  })
  
  output$ui_slider_neg <- renderUI({
    req(bornes_canal())
    b <- bornes_canal()
    span <- b[2] - b[1]
    sliderInput("gate_slider_neg", NULL,
                min=b[1], max=b[2],
                value=c(b[1], b[1]+span*0.25),
                step=0.01, width="100%")
  })
  output$ui_slider_pos <- renderUI({
    req(bornes_canal())
    b <- bornes_canal()
    span <- b[2] - b[1]
    sliderInput("gate_slider_pos", NULL,
                min=b[1], max=b[2],
                value=c(b[1]+span*0.65, b[2]),
                step=0.01, width="100%")
  })
  
  # Clic sur le graphique → déplace la borne la plus proche
  observeEvent(input$click_gating, {
    req(input$gate_slider_neg, input$gate_slider_pos)
    x_val <- input$click_gating$x
    req(!is.null(x_val), length(x_val)==1, is.finite(x_val))
    lim_n <- input$gate_slider_neg
    lim_p <- input$gate_slider_pos
    b     <- bornes_canal()
    bornes <- c(n1=lim_n[1], n2=lim_n[2], p1=lim_p[1], p2=lim_p[2])
    cible  <- names(which.min(abs(bornes - x_val)))
    if (length(cible)==0) return()
    # Étend la borne du slider si le clic est hors plage
    new_min <- min(b[1], x_val - 0.1)
    new_max <- max(b[2], x_val + 0.1)
    if (cible=="n1") updateSliderInput(session,"gate_slider_neg", min=new_min, max=new_max,
                                       value=c(x_val, max(x_val+0.05, lim_n[2])))
    else if (cible=="n2") updateSliderInput(session,"gate_slider_neg", min=new_min, max=new_max,
                                            value=c(min(lim_n[1], x_val-0.05), x_val))
    else if (cible=="p1") updateSliderInput(session,"gate_slider_pos", min=new_min, max=new_max,
                                            value=c(x_val, max(x_val+0.05, lim_p[2])))
    else if (cible=="p2") updateSliderInput(session,"gate_slider_pos", min=new_min, max=new_max,
                                            value=c(min(lim_p[1], x_val-0.05), x_val))
  })
  
  output$plot_gating <- renderPlot({
    rv$trig_trans_applique; req(mono_trans_ok(), input$gate_canal)
    canal   <- input$gate_canal
    lim_n   <- if(!is.null(input$gate_slider_neg)) input$gate_slider_neg else c(0,2)
    lim_p   <- if(!is.null(input$gate_slider_pos)) input$gate_slider_pos else c(4,7)
    use_uns <- if(!is.null(input$gate_use_unstained)) input$gate_use_unstained else TRUE
    tryCatch({
      rv$carrot$graphiques_gates(
        nom_canal=canal, shiny_neg=lim_n, shiny_pos=lim_p,
        afficher_unstained_neg=use_uns)
    }, error=function(e) {
      ggplot() +
        annotate("text",x=0.5,y=0.5,
                 label=paste0("Appliquez d'abord la transformation.\n\nDétail : ",e$message),
                 size=4,color="#605ca8",hjust=0.5) + theme_void()
    })
  })
  
  observeEvent(input$save_gate, {
    req(mono_trans_ok(), input$gate_canal, input$gate_slider_neg, input$gate_slider_pos)
    canal   <- input$gate_canal
    lim_n   <- input$gate_slider_neg
    lim_p   <- input$gate_slider_pos
    use_uns <- if(!is.null(input$gate_use_unstained)) input$gate_use_unstained else TRUE
    tryCatch({
      rv$carrot$definir_et_extraire(
        nom_canal=canal, intervalle_gate_negatif=lim_n,
        intervalle_gate_positif=lim_p, utiliser_unstained=use_uns)
      rv$carrot$plots_gates[[canal]] <- rv$carrot$graphiques_gates(
        nom_canal=canal, shiny_neg=lim_n, shiny_pos=lim_p,
        afficher_unstained_neg=use_uns)
      rv$trig_gate_sauve <- rv$trig_gate_sauve + 1L
      showNotification(
        paste0("✔ Gate enregistré : ", canal,
               "  Nég [",round(lim_n[1],2),"–",round(lim_n[2],2),"]",
               "  Pos [",round(lim_p[1],2),"–",round(lim_p[2],2),"]"),
        type="message", duration=5)
    }, error=function(e) {
      showNotification(paste("Erreur gate :", e$message), type="error")
    })
  })
  
  # Calcul automatique spillover quand tous les gates sont OK
  observeEvent(rv$trig_gate_sauve, {
    req(carrot_ok(), !is.null(rv$carrot$canaux), length(rv$carrot$canaux)>0)
    tous_ok <- all(sapply(rv$carrot$canaux, function(cx)
      !is.null(rv$carrot$gates_positifs[[cx]]) &&
        !is.null(rv$carrot$gates_negatifs[[cx]])))
    if (tous_ok && is.null(rv$carrot$S_matrix)) {
      tryCatch({
        rv$carrot$calculer_spillover()
        rv$trig_spillover <- rv$trig_spillover + 1L
        showNotification("✔ Matrice de spillover calculée automatiquement.",
                         type="message", duration=5)
      }, error=function(e) NULL)
    }
  })
  
  
  # ---- Matrice de spillover + Édition -----------------------
  
  output$ui_select_ech_matrice_origine <- renderUI({
    rv$trig_spillover; rv$trig_matrice_edit; req(spillover_ok())
    choix <- c("Matrice globale"="__global__")
    if (length(rv$carrot$S_matrices_par_echantillon)>0)
      choix <- c(choix, names(rv$carrot$S_matrices_par_echantillon))
    selectInput("sel_ech_matrice","Afficher la matrice de :", choices=choix)
  })
  
  matrice_courante_r <- reactive({
    rv$trig_spillover; rv$trig_matrice_edit; req(spillover_ok())
    sel <- input$sel_ech_matrice
    if (!is.null(sel) && sel!="__global__" &&
        !is.null(rv$carrot$S_matrices_par_echantillon[[sel]]))
      rv$carrot$S_matrices_par_echantillon[[sel]]
    else rv$carrot$S_matrix
  })
  
  output$matrice_origine_affichage <- renderTable({
    req(matrice_courante_r())
    as.data.frame(round(matrice_courante_r(),4))
  }, rownames=TRUE, striped=TRUE, hover=TRUE, bordered=TRUE)
  
  output$ui_edit_axes_c1 <- renderUI({
    rv$trig_spillover; req(spillover_ok())
    selectInput("edit_c1","Canal source (ligne) :", choices=rownames(rv$carrot$S_matrix))
  })
  output$ui_edit_axes_c2 <- renderUI({
    rv$trig_spillover; req(spillover_ok(), input$edit_c1)
    selectInput("edit_c2","Spillover vers (colonne) :",
                choices=setdiff(colnames(rv$carrot$S_matrix), input$edit_c1))
  })
  
  observeEvent(c(input$edit_c1, input$edit_c2, input$sel_ech_matrice), {
    req(spillover_ok(), input$edit_c1, input$edit_c2)
    mat <- matrice_courante_r()
    if (input$edit_c1 %in% rownames(mat) && input$edit_c2 %in% colnames(mat))
      updateNumericInput(session,"edit_axes_val",
                         value=round(mat[input$edit_c1,input$edit_c2],5))
  })
  
  observeEvent(input$valider_changement_axe, {
    req(spillover_ok(), input$edit_c1, input$edit_c2,
        input$edit_axes_val, input$sel_ech_matrice)
    nom <- input$sel_ech_matrice; req(nom!="__global__")
    tryCatch({
      if (is.null(rv$historique_matrices[[nom]])) {
        mat_base <- if(!is.null(rv$carrot$S_matrices_par_echantillon[[nom]]))
          rv$carrot$S_matrices_par_echantillon[[nom]] else rv$carrot$S_matrix
        rv$historique_matrices[[nom]] <- list(mat_base)
        rv$index_versions[[nom]]      <- 1L
      }
      rv$carrot$modifier_spillover(nom, input$edit_c1, input$edit_c2, input$edit_axes_val)
      idx  <- rv$index_versions[[nom]]
      hist <- rv$historique_matrices[[nom]][seq_len(idx)]
      hist <- c(hist, list(rv$carrot$S_matrices_par_echantillon[[nom]]))
      rv$historique_matrices[[nom]] <- hist
      rv$index_versions[[nom]]      <- length(hist)
      rv$trig_matrice_edit          <- rv$trig_matrice_edit + 1L
      showNotification(paste0("[",input$edit_c1," → ",input$edit_c2,"] = ",
                              round(input$edit_axes_val,5)),
                       type="message", duration=3)
    }, error=function(e) showNotification(paste("Erreur :",e$message), type="error"))
  })
  
  observeEvent(input$btn_undo_matrice, {
    nom <- input$sel_ech_matrice; req(!is.null(nom), nom!="__global__")
    idx <- rv$index_versions[[nom]]; req(!is.null(idx), idx>1)
    rv$index_versions[[nom]] <- idx-1L
    rv$carrot$S_matrices_par_echantillon[[nom]] <- rv$historique_matrices[[nom]][[idx-1L]]
    rv$trig_matrice_edit <- rv$trig_matrice_edit + 1L
    showNotification("Undo appliqué.", type="message", duration=2)
  })
  
  observeEvent(input$btn_redo_matrice, {
    nom <- input$sel_ech_matrice; req(!is.null(nom), nom!="__global__")
    idx <- rv$index_versions[[nom]]; tot <- length(rv$historique_matrices[[nom]])
    req(!is.null(idx), idx<tot)
    rv$index_versions[[nom]] <- idx+1L
    rv$carrot$S_matrices_par_echantillon[[nom]] <- rv$historique_matrices[[nom]][[idx+1L]]
    rv$trig_matrice_edit <- rv$trig_matrice_edit + 1L
    showNotification("Redo appliqué.", type="message", duration=2)
  })
  
  output$texte_version_matrice <- renderText({
    rv$trig_matrice_edit
    nom <- input$sel_ech_matrice
    if (is.null(nom) || nom=="__global__") return("Sélectionnez un échantillon pour éditer.")
    idx <- rv$index_versions[[nom]]; tot <- length(rv$historique_matrices[[nom]])
    if (is.null(idx)||tot==0) "Aucune modification enregistrée."
    else paste0("Version ", idx, " / ", tot)
  })
  
  observeEvent(input$sauvegarder_matrice_echantillon, {
    req(spillover_ok(), input$sel_ech_matrice)
    nom <- input$sel_ech_matrice; req(nom!="__global__")
    tryCatch({
      mat_perso <- rv$carrot$S_matrices_par_echantillon[[nom]]; req(!is.null(mat_perso))
      fcs_brut  <- rv$carrot$echantillons[[nom]]
      cx_com    <- intersect(colnames(mat_perso), flowCore::colnames(fcs_brut))
      rv$carrot$echantillons_traites[[nom]] <-
        flowCore::compensate(fcs_brut, mat_perso[cx_com,cx_com,drop=FALSE])
      rv$trig_compensation <- rv$trig_compensation + 1L
      showNotification(paste0("✔ Matrice personnalisée appliquée à : ",nom),
                       type="message", duration=4)
    }, error=function(e) showNotification(paste("Erreur :",e$message), type="error"))
  })
  
  
  # ---- Figures & Validations --------------------------------
  
  output$ui_fig_ech <- renderUI({
    rv$trig_fcs_charge; req(carrot_ok(), length(rv$carrot$echantillons)>0)
    selectInput("fig_ech_sel","Échantillon :", choices=names(rv$carrot$echantillons))
  })
  
  output$ui_conditionnel_axes <- renderUI({
    rv$trig_fcs_charge; rv$trig_trans_applique; req(carrot_ok(), input$fig_mode)
    fcs   <- rv$carrot$echantillons[[1]]
    choix <- choix_avec_labels(fcs)
    if (input$fig_mode=="pairs") {
      tagList(helpText("Une paire à la fois — sélectionnez ci-dessous."),
              uiOutput("ui_pairs_selector"))
    } else {
      tagList(
        selectInput("fig_cx","Axe X :", choices=choix),
        selectInput("fig_cy","Axe Y :", choices=choix,
                    selected=if(length(choix)>1) choix[2] else choix[1])
      )
    }
  })
  
  output$ui_pairs_selector <- renderUI({
    rv$trig_trans_applique; req(mono_trans_ok())
    cx_mono <- names(rv$carrot$monomarques_trans)
    cx_mono <- cx_mono[cx_mono!="TUBE_UNSTAINED"]
    if (length(cx_mono)<2) return(helpText("Pas assez de canaux."))
    paires <- sapply(seq_len(length(cx_mono)-1), function(i)
      paste0(cx_mono[i]," \u2192 ",cx_mono[i+1]))
    selectInput("fig_paire_sel","Paire de canaux :", choices=paires)
  })
  
  axes_selectionnes <- reactive({
    req(input$fig_mode)
    if (input$fig_mode=="free") {
      req(input$fig_cx, input$fig_cy)
      list(cx=input$fig_cx, cy=input$fig_cy)
    } else {
      req(input$fig_paire_sel, mono_trans_ok())
      cx_mono <- names(rv$carrot$monomarques_trans)
      cx_mono <- cx_mono[cx_mono!="TUBE_UNSTAINED"]
      parties <- strsplit(input$fig_paire_sel," \u2192 ")[[1]]
      list(cx=trimws(parties[1]), cy=trimws(parties[2]))
    }
  })
  
  # Graphique via la méthode R6 visualiser_compensation (carrés, axes partagés)
  output$plot_validation_comp <- renderPlot({
    rv$trig_compensation; rv$trig_matrice_edit; rv$trig_fcs_charge
    req(carrot_ok(), input$fig_ech_sel)
    
    ax <- tryCatch(axes_selectionnes(), error=function(e) NULL)
    req(!is.null(ax))
    
    ech      <- input$fig_ech_sel
    max_pts  <- if(!is.null(input$fig_max_points)) input$fig_max_points else 5000
    mode_aff <- if(is.null(input$fig_view_type)) "Both" else {
      v <- input$fig_view_type
      if (grepl("apr.s",v,ignore.case=TRUE) && grepl("avant",v,ignore.case=TRUE)) "Both"
      else if (grepl("apr.s",v,ignore.case=TRUE)) "After compensation only"
      else "Before compensation only"
    }
    
    fcs_brut <- rv$carrot$echantillons[[ech]]; req(!is.null(fcs_brut))
    cx_ok    <- colnames(flowCore::exprs(fcs_brut))
    req(ax$cx %in% cx_ok, ax$cy %in% cx_ok)
    
    # Utilise la méthode R6 qui garantit carrés + axes partagés
    tryCatch({
      result <- rv$carrot$visualiser_compensation(
        nom_echantillon = ech,
        canal_x         = ax$cx,
        canal_y         = ax$cy,
        max_points      = max_pts,
        affichage       = mode_aff
      )
      if (!is.null(result)) grid::grid.draw(result)
    }, error=function(e) {
      # Fallback ggplot si la méthode R6 échoue
      mat <- as.data.frame(flowCore::exprs(fcs_brut)[,c(ax$cx,ax$cy),drop=FALSE])
      if (!is.null(rv$carrot$config_transformations[[ax$cx]])) {
        cofac_x <- rv$carrot$config_transformations[[ax$cx]]$cofacteur
        cofac_y <- rv$carrot$config_transformations[[ax$cy]]$cofacteur
        if(!is.null(cofac_x)) mat[[ax$cx]] <- asinh(mat[[ax$cx]]/cofac_x)
        if(!is.null(cofac_y)) mat[[ax$cy]] <- asinh(mat[[ax$cy]]/cofac_y)
      }
      n   <- min(nrow(mat), max_pts)
      idx <- sample(seq_len(nrow(mat)), n)
      mat <- mat[idx,,drop=FALSE]
      p <- ggplot(mat, aes(x=.data[[ax$cx]], y=.data[[ax$cy]])) +
        ggpointdensity::geom_pointdensity(size=0.2, alpha=0.5) +
        scale_color_gradientn(colours=c("darkblue","blue","cyan","greenyellow","yellow","darkorange","red")) +
        theme_bw() + theme(legend.position="none", aspect.ratio=1) +
        labs(title=paste("Avant compensation —",ech),
             x=label_canal(fcs_brut,ax$cx), y=label_canal(fcs_brut,ax$cy))
      print(p)
    })
  }, height=450)
  
  # Lancer la compensation globale
  observeEvent(input$run_compensation, {
    req(carrot_ok())
    withProgress(message="Compensation en cours...", value=0.5, {
      tryCatch({
        if (is.null(rv$carrot$S_matrix)) {
          rv$carrot$calculer_spillover()
          rv$trig_spillover <- rv$trig_spillover + 1L
        }
        rv$carrot$compenser()
        rv$trig_compensation <- rv$trig_compensation + 1L
        showNotification(
          paste0("✔ Compensation appliquée sur ",
                 length(rv$carrot$echantillons_traites)," échantillon(s)."),
          type="message", duration=5)
      }, error=function(e) {
        showNotification(paste("Erreur :",e$message), type="error", duration=8)
      })
    })
  })
  
  
  # ---- Exports ----------------------------------------------
  
  output$ui_select_fcs_export <- renderUI({
    rv$trig_compensation; req(comp_ok())
    checkboxGroupInput("fcs_export_sel","Sélectionner les échantillons :",
                       choices=names(rv$carrot$echantillons_traites),
                       selected=names(rv$carrot$echantillons_traites))
  })
  
  observe({
    req(input$fcs_export_sel)
    n <- length(input$fcs_export_sel)
    if (n==1)     { shinyjs::show("wrapper_single_fcs"); shinyjs::hide("wrapper_zip_fcs") }
    else if (n>1) { shinyjs::hide("wrapper_single_fcs"); shinyjs::show("wrapper_zip_fcs") }
    else          { shinyjs::hide("wrapper_single_fcs"); shinyjs::hide("wrapper_zip_fcs") }
  })
  
  output$btn_fcs_unique <- downloadHandler(
    filename=function() paste0(input$fcs_export_sel[1],"_compense.fcs"),
    content=function(file) {
      req(comp_ok(), length(input$fcs_export_sel)==1)
      fcs_obj <- rv$carrot$echantillons_traites[[input$fcs_export_sel[1]]]
      if (!is.null(rv$carrot$S_matrix))
        flowCore::keyword(fcs_obj)[["$SPILLOVER"]] <- rv$carrot$S_matrix
      flowCore::write.FCS(fcs_obj, filename=file)
    }
  )
  
  output$btn_fcs_zip <- downloadHandler(
    filename=function() "FCS_Compenses.zip",
    content=function(file) {
      req(comp_ok(), length(input$fcs_export_sel)>1)
      tmp <- tempfile(); dir.create(tmp)
      for (nom in input$fcs_export_sel) {
        fcs_obj <- rv$carrot$echantillons_traites[[nom]]
        if (is.null(fcs_obj)) next
        if (!is.null(rv$carrot$S_matrix))
          flowCore::keyword(fcs_obj)[["$SPILLOVER"]] <- rv$carrot$S_matrix
        flowCore::write.FCS(fcs_obj,
                            file.path(tmp, paste0(gsub("[^a-zA-Z0-9_]","_",nom),"_compense.fcs")))
      }
      owd <- setwd(tmp); on.exit(setwd(owd))
      utils::zip(file, files=list.files(tmp))
    }
  )
  
  output$download_session_rds <- downloadHandler(
    filename=function() paste0("Session_Compensation_",Sys.Date(),".rds"),
    content=function(file) {
      req(carrot_ok())
      saveRDS(list(
        meta=list(date=Sys.time(), canaux=rv$carrot$canaux, mode=rv$carrot$mode),
        configuration_technique=list(
          trans_list=rv$carrot$trans_list,
          matrice_spillover=rv$carrot$S_matrix,
          matrices_par_ech=rv$carrot$S_matrices_par_echantillon
        ),
        gating=list(
          gates_positifs=rv$carrot$gates_positifs,
          gates_negatifs=rv$carrot$gates_negatifs,
          bornes_gates_pos=rv$carrot$bornes_gates_pos,
          bornes_gates_neg=rv$carrot$bornes_gates_neg
        ),
        visualisations=list(
          plots_gates=rv$carrot$plots_gates,
          plots_compensation=rv$carrot$plots_compensation
        ),
        echantillons_traites=rv$carrot$echantillons_traites
      ), file=file)
    }
  )
}

shinyApp(ui, server)