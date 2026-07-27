library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)

RES_PIXELS_GATING    <- 200 # Résolution de la grille de densité pour le gating interactif (plotly) — rapprochée de celle des figures statiques (400) pour un rendu visuel cohérent, tout en restant fluide en direct dans le navigateur.
TAILLE_PIXEL_GATING  <- 3 # Taille des marqueurs plotly représentant chaque pixel de densité (réduite en conséquence de la résolution plus fine, pour éviter tout chevauchement excessif)

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
    
    tabsetPanel(
      id = ns("pretrait_steps"),
      
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
                       tags$b("Double-clic"), " pour fermer le polygone. Le gate affiché est ",
                       "celui de l'échantillon sélectionné : chaque échantillon peut avoir sa propre forme."),
                   br(),
                   
                   fluidRow(
                     column(6, actionButton(ns("btn_undo_sommet"),
                                            tagList(icon("undo"), " Annuler sommet"),
                                            class = "btn-default btn-sm", style = "width:100%;")),
                     column(6, actionButton(ns("btn_reset_gate"),
                                            tagList(icon("eraser"), " Effacer"),
                                            class = "btn-default btn-sm", style = "width:100%;"))
                   )
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE CENTRALE — Graphiques
          # ──────────────────────────────────────────────────────────────────
          column(width = 6,
                 box(title = tagList(icon("draw-polygon"), " Dessin interactif du gate"),
                     width = NULL, status = "info", solidHeader = TRUE,
                     plotlyOutput(ns("plot_gate_dessin"), height = "480px"),
                     br(),
                     actionButton(ns("btn_save_gate"),
                                  tagList(icon("check-circle"), " Enregistrer ce gate"),
                                  class = "btn-success", style = "width:100%; font-weight:bold;")
                 ),
                 
                 box(title = uiOutput(ns("titre_resultat")), width = NULL,
                     status = "success", solidHeader = TRUE,
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
                     uiOutput(ns("ui_statut_doublet_auto")),
                     br(),
                     actionButton(ns("btn_appliquer_doublet_auto"),
                                  tagList(icon("play"), " Appliquer le retrait automatique"),
                                  class = "btn-success", style = "width:100%; font-weight:bold;")
                   ),
                   
                   # ── Méthode par gating manuel ────────────────────────────
                   conditionalPanel(
                     condition = paste0("input['", ns("methode_doublet"), "'] == 'gate'"),
                     div(class = "pretrait-instr",
                         icon("hand-pointer"),
                         " Cliquez sur le graphique pour ajouter des sommets. ",
                         tags$b("Double-clic"), " pour fermer le polygone. Le gate affiché est ",
                         "celui de l'échantillon sélectionné : chaque échantillon peut avoir sa propre forme."),
                     br(),
                     fluidRow(
                       column(6, actionButton(ns("btn_undo_sommet_doublet"),
                                              tagList(icon("undo"), " Annuler sommet"),
                                              class = "btn-default btn-sm", style = "width:100%;")),
                       column(6, actionButton(ns("btn_reset_gate_doublet"),
                                              tagList(icon("eraser"), " Effacer"),
                                              class = "btn-default btn-sm", style = "width:100%;"))
                     )
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
                       plotlyOutput(ns("plot_gate_doublet_dessin"), height = "480px"),
                       br(),
                       actionButton(ns("btn_save_gate_doublet"),
                                    tagList(icon("check-circle"), " Enregistrer ce gate"),
                                    class = "btn-success", style = "width:100%; font-weight:bold;")
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
                     div(class = "pretrait-instr",
                         icon("hand-pointer"),
                         " Cliquez sur le graphique pour ajouter des sommets autour des cellules ",
                         tags$b("vivantes"), ". ", tags$b("Double-clic"), " pour fermer le polygone. ",
                         "Le gate affiché est celui de l'échantillon sélectionné : chaque échantillon ",
                         "peut avoir sa propre forme."),
                     br(),
                     fluidRow(
                       column(6, actionButton(ns("btn_undo_sommet_viabilite"),
                                              tagList(icon("undo"), " Annuler sommet"),
                                              class = "btn-default btn-sm", style = "width:100%;")),
                       column(6, actionButton(ns("btn_reset_gate_viabilite"),
                                              tagList(icon("eraser"), " Effacer"),
                                              class = "btn-default btn-sm", style = "width:100%;"))
                     )
                   )
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE CENTRALE — Graphiques
          # ──────────────────────────────────────────────────────────────────
          column(width = 6,
                 box(title = tagList(icon("draw-polygon"), " Dessin interactif du gate"),
                     width = NULL, status = "info", solidHeader = TRUE,
                     plotlyOutput(ns("plot_gate_viabilite_dessin"), height = "480px"),
                     br(),
                     actionButton(ns("btn_save_gate_viabilite"),
                                  tagList(icon("check-circle"), " Enregistrer ce gate"),
                                  class = "btn-success", style = "width:100%; font-weight:bold;")
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
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET — EXPORT
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("file-export"), " Export"),
        fluidRow(
          # --- Colonne de Gauche : Sélection & options ---
          column(width = 4,
                 wellPanel(
                   h4("Export des résultats de prétraitement"),
                   
                   div(class = "pretrait-instr",
                       icon("info-circle"),
                       " Sélectionnez la ou les étapes de prétraitement à exporter, ainsi que les échantillons souhaités, ",
                       "puis téléchargez ", tags$b("les fichiers FCS"), " correspondants et/ou ",
                       tags$b("une session RDS"), " regroupant tous les paramètres utilisés (gates, seuils, canaux)."),
                   hr(),
                   
                   checkboxGroupInput(ns("export_etapes_pretrait"), "Étapes à inclure :",
                                      choices = c("Post débris"          = "debris",
                                                  "Post doublets"        = "doublets",
                                                  "Post cellules mortes" = "viabilite"),
                                      selected = c("debris", "doublets", "viabilite")),
                   
                   uiOutput(ns("ui_select_echantillons_export_pretrait")),
                   
                   hr(),
                   
                   downloadButton(ns("dl_export_fcs_pretrait"), tagList(icon("file-arrow-down"), " Télécharger les FCS (.zip)"),
                                  class = "btn-warning", style = "width:100%; font-weight:bold; margin-bottom:15px;"),
                   
                   hr(),
                   
                   textInput(ns("pretrait_rds_filename"), "Nom du fichier RDS :",
                             value = "Pretraitement_Session_Complete.rds",
                             placeholder = "nom_session.rds"),
                   downloadButton(ns("dl_export_rds_pretrait"), tagList(icon("file-arrow-down"), " Télécharger la session RDS"),
                                  class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          
          # --- Colonne de Droite : Récapitulatif de ce qui sera exporté ---
          column(width = 8,
                 box(title = "Aperçu des données disponibles pour l'export", width = NULL,
                     status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_recap_export_pretrait"))
                 )
          )
        )
      )
      
    ) # /tabsetPanel
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
    # Le gate "actif" (affiché en résultat / préchargé au changement d'échantillon)
    # est simplement celui dont le nom est actuellement saisi dans le champ
    # "Nom du gate" — plus besoin de suivre un état séparé.
    nom_gate_actuel <- reactive({
      nm <- trimws(input$nom_gate %||% "")
      if (nchar(nm) == 0) NULL else nm
    })
    
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
      # "brutes" (ou repli) : on ne remonte JAMAIS au-delà de la source d'entrée
      # du gating débris (bordures > PeacoQC/flowAI > compensées/démixées).
      # Utiliser p$obtenir_derniere_source() ici serait incorrect : une fois qu'un
      # gate de débris existe déjà pour au moins un échantillon, cette méthode
      # renvoie post_debris lui-même (résultat du filtrage), ce qui bouclerait
      # sur sa propre sortie et n'afficherait plus que les événements déjà
      # conservés par le gate, au lieu du nuage complet à gater.
      if (!is.null(p$post_retrait_bordures) && length(p$post_retrait_bordures) > 0) return(p$post_retrait_bordures)
      if (!is.null(p$post_PeacoQC) && length(p$post_PeacoQC) > 0) return(p$post_PeacoQC)
      if (!is.null(p$post_flowAI)  && length(p$post_flowAI)  > 0) return(p$post_flowAI)
      return(p$echantillons_traites)
    }
    
    # ════════════════════════════════════════════════════════════════════════
    # BORDURES (MARGINS)
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_canaux_bordures <- renderUI({
      p <- carrot_obj()
      src <- p$obtenir_derniere_source()
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
      # Préserve la sélection en cours si elle est toujours valide, plutôt que
      # de revenir systématiquement au premier échantillon à chaque
      # reconstruction de ce sélecteur (par ex. après l'enregistrement d'un
      # gate) — sans quoi l'utilisateur se retrouverait renvoyé sur un autre
      # échantillon sans s'en rendre compte, en plein milieu d'une édition.
      actuel <- input$sel_ech_debris
      sel <- if (!is.null(actuel) && actuel %in% names(src)) actuel else names(src)[1]
      selectInput(ns("sel_ech_debris"), "Échantillon :", choices = names(src), selected = sel)
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
    
    # Calcule, pour un échantillon et des canaux donnés, la forme du gate à
    # précharger (celle déjà enregistrée sous le nom de gate courant, si elle
    # existe et correspond aux mêmes canaux X/Y) — ou une liste vide sinon.
    # Utilisée à la fois pour synchroniser sommets_rv() (édition interactive)
    # et pour dessiner directement la forme dans le rendu initial du graphique
    # (afin que le tracé soit toujours correct dès le premier rendu, sans
    # dépendre d'une mise à jour ultérieure via plotlyProxy).
    # Construit, à partir d'une liste de sommets et du data.frame affiché, les
    # coordonnées de la trace fermée (polygone) et les formes de poignées
    # (cercles) correspondantes — utilisé à la fois pour le rendu initial du
    # graphique et pour la mise à jour via plotlyProxy lors de l'édition.
    # Quand l'utilisateur déplace un cercle-sommet (shape), Plotly renvoie via
    # plotly_relayout les 4 coins de sa boîte englobante (shapes[i].x0/x1/y0/y1)
    # simultanément dans le même événement — jamais le centre directement. Le
    # centre réel du cercle (donc la position du sommet) est leur moyenne ;
    # utiliser x0/y0 seuls décale systématiquement le sommet à chaque
    # déplacement (d'un montant égal au rayon du cercle), ce qui le fait
    # dériver de plus en plus loin du curseur au fil des mouvements.
    # Renvoie la liste de sommets mise à jour, ou NULL si rien n'a changé.
    appliquer_deplacement_shapes <- function(ev, soms) {
      indices <- regmatches(names(ev), regexpr("(?<=shapes\\[)[0-9]+(?=\\]\\.)", names(ev), perl = TRUE))
      indices <- unique(as.integer(indices[nchar(indices) > 0]))
      if (length(indices) == 0) return(NULL)
      
      modifie <- FALSE
      for (idx0 in indices) {
        i <- idx0 + 1L   # 0-indexé en JS → 1-indexé en R
        if (i < 1 || i > length(soms)) next
        prefix <- paste0("shapes[", idx0, "].")
        x0 <- ev[[paste0(prefix, "x0")]]; x1 <- ev[[paste0(prefix, "x1")]]
        y0 <- ev[[paste0(prefix, "y0")]]; y1 <- ev[[paste0(prefix, "y1")]]
        if (!is.null(x0) && !is.null(x1)) { soms[[i]]$x <- (x0 + x1) / 2; modifie <- TRUE }
        if (!is.null(y0) && !is.null(y1)) { soms[[i]]$y <- (y0 + y1) / 2; modifie <- TRUE }
      }
      if (modifie) soms else NULL
    }
    
    construire_trace_et_shapes_gate <- function(soms, df) {
      if (length(soms) == 0) {
        return(list(x = numeric(0), y = numeric(0), shapes = list()))
      }
      xs <- c(vapply(soms, `[[`, numeric(1), "x"), soms[[1]]$x)   # ferme visuellement
      ys <- c(vapply(soms, `[[`, numeric(1), "y"), soms[[1]]$y)
      
      shapes <- list()
      if (length(soms) >= 2) {
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
      }
      list(x = xs, y = ys, shapes = shapes)
    }
    
    calculer_soms_preload_debris <- function(p, nom, cx, cy) {
      ga <- isolate(nom_gate_actuel())
      if (is.null(ga) || is.null(nom)) return(list())
      infos <- p$gates_history[[ga]][[nom]]
      if (!is.null(infos) && !is.null(infos$polygone) &&
          identical(infos$canal_x, cx) && identical(infos$canal_y, cy)) {
        mat <- infos$polygone
        lapply(seq_len(nrow(mat)), function(i) list(x = mat[i, 1], y = mat[i, 2]))
      } else {
        list()
      }
    }
    
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
    
    # Recharge (ou réinitialise) le polygone à chaque changement pertinent :
    # échantillon, canaux, source ou nom du gate. Un seul observer gère tous
    # les cas pour éviter toute course entre deux observers concurrents :
    #  - si un gate de ce nom existe déjà pour CET échantillon ET sur ces
    #    MÊMES canaux X/Y, on recharge sa forme (modifiable) ;
    #  - sinon (canaux différents, nom inconnu, ou aucun gate pour cet
    #    échantillon), le canevas repart vide pour un nouveau tracé.
    observeEvent(list(input$sel_ech_debris, input$canal_x_debris,
                      input$canal_y_debris, input$source_debris, input$nom_gate), {
                        p   <- carrot_obj()
                        nom <- input$sel_ech_debris
                        cx  <- input$canal_x_debris
                        cy  <- input$canal_y_debris
                        req(nom)
                        
                        sommets_rv(calculer_soms_preload_debris(p, nom, cx, cy))
                      }, ignoreInit = FALSE)
    
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
      lbl_x <- p_obj$obtenir_label(fcs, input$canal_x_debris)
      lbl_y <- p_obj$obtenir_label(fcs, input$canal_y_debris)
      
      # Densité par binning raster (rapide, indépendant du nombre d'événements affichés).
      # Rendue en petits ronds colorés par densité (mode "scatter" classique, même famille
      # que le tracé du gate, pour éviter tout ralentissement dû au mélange SVG/WebGL) plutôt
      # qu'en heatmap plotly : on évite ainsi tout lissage/interpolation qui brouillerait les
      # frontières entre zones.
      lim_x <- range(df$X, na.rm = TRUE)
      lim_y <- range(df$Y, na.rm = TRUE)
      dens  <- calculer_densite_raster(df$X, df$Y, lim_x, lim_y, res = RES_PIXELS_GATING, lissage = FALSE)
      
      plt <- plot_ly(source = ns("plot_gate_dessin"))
      if (!is.null(dens)) {
        plt <- plt %>%
          add_trace(data = dens, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    marker = list(symbol = "circle", size = TAILLE_PIXEL_GATING,
                                  color = ~densite, colorscale = COLORSCALE_DENSITE_PLOTLY,
                                  line = list(width = 0)),
                    hoverinfo = "none")
      } else {
        plt <- plt %>%
          add_trace(data = df, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    name = "Cellules",
                    marker = list(size = 2, color = "#0077b6", opacity = 0.4),
                    hoverinfo = "none")
      }
      
      # Précharge, dès le rendu initial, la forme du gate déjà enregistrée pour
      # cet échantillon (si elle existe) — plutôt que de partir d'une trace
      # vide en attendant une mise à jour ultérieure via plotlyProxy, ce qui
      # évite toute course entre le redessin complet et l'injection du gate.
      soms_init <- calculer_soms_preload_debris(p_obj, input$sel_ech_debris, input$canal_x_debris, input$canal_y_debris)
      trace_init <- construire_trace_et_shapes_gate(soms_init, df)
      
      plt %>%
        add_trace(x = trace_init$x, y = trace_init$y, type = "scatter",
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
          shapes = trace_init$shapes
        ) %>%
        config(displayModeBar = TRUE, editable = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE,
               # Désactive le double-clic natif de Plotly (reset des axes / des
               # formes) : ce plot a son propre gestionnaire de double-clic
               # (fermeture du polygone) et le comportement natif interférait
               # avec lui en effaçant le tracé du gate en cours d'édition.
               doubleClick = FALSE)
    })
    
    # Pousse les sommets courants vers le graphique déjà rendu, sans redraw complet
    observeEvent(sommets_rv(), {
      soms  <- sommets_rv()
      proxy <- plotlyProxy(ns("plot_gate_dessin"), session)
      df    <- isolate(donnees_plot())
      trace <- construire_trace_et_shapes_gate(soms, df)
      
      plotlyProxyInvoke(proxy, "restyle", list(x = list(trace$x), y = list(trace$y)), list(1))
      plotlyProxyInvoke(proxy, "relayout", list(shapes = trace$shapes))
    }, ignoreNULL = FALSE)
    
    # Capture les déplacements de shapes (sommets éditables)
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_dessin")), {
      ev   <- event_data("plotly_relayout", source = ns("plot_gate_dessin"))
      req(ev)
      
      # Sécurité : si sommets_rv() n'a pas (encore) été synchronisé avec la
      # forme réellement affichée (ex. juste après un changement d'échantillon,
      # décalage de timing réactif), on la recharge avant d'appliquer le
      # déplacement — sinon le glissé d'un gate déjà enregistré n'aurait
      # silencieusement aucun effet tant qu'un enregistrement n'a pas forcé un
      # nouveau rendu complet.
      soms_actuels <- sommets_rv()
      if (length(soms_actuels) == 0) {
        p <- carrot_obj()
        soms_actuels <- calculer_soms_preload_debris(p, input$sel_ech_debris, input$canal_x_debris, input$canal_y_debris)
      }
      
      soms_maj <- appliquer_deplacement_shapes(ev, soms_actuels)
      if (!is.null(soms_maj)) sommets_rv(soms_maj)
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
      # (l'enveloppe convexe garantit un polygone propre même si les sommets
      # n'ont pas été cliqués dans un ordre parfaitement séquentiel).
      xs <- sapply(soms, `[[`, "x")
      ys <- sapply(soms, `[[`, "y")
      cx <- input$canal_x_debris
      cy <- input$canal_y_debris
      
      idx_hull <- grDevices::chull(xs, ys)
      mat_poly <- cbind(xs[idx_hull], ys[idx_hull])
      colnames(mat_poly) <- c(cx, cy)
      
      # Ce nom de gate a-t-il déjà été enregistré pour au moins un échantillon ?
      # Si non, l'enregistrement qui suit est le tout premier : il s'appliquera
      # par défaut à toute la cohorte (cf. appliquer_gate_nomme) — on le
      # signale explicitement à l'utilisateur pour lever toute ambiguïté.
      premier_enregistrement <- is.null(p$gates_history[[nom_g]]) || length(p$gates_history[[nom_g]]) == 0
      
      withProgress(message = paste0("Application du gate '", nom_g, "'..."), value = 0.4, {
        tryCatch({
          p$appliquer_gate_nomme(
            nom_gate       = nom_g,
            matrice_points = mat_poly,
            canal_x        = cx,
            canal_y        = cy,
            source_nettoyage = input$source_debris %||% "brutes",
            nom_echantillon  = input$sel_ech_debris
          )
          pipeline(p)
          gates_trigger(gates_trigger() + 1L)
          
          if (premier_enregistrement) {
            nb <- length(p$gates_history[[nom_g]])
            showNotification(
              paste0("Gate '", nom_g, "' créé et appliqué par défaut à ", nb, " échantillon(s). ",
                     "Changez d'échantillon pour ajuster sa forme individuellement, puis ré-enregistrez."),
              type = "message", duration = 8
            )
          } else {
            showNotification(paste0("Gate '", nom_g, "' mis à jour pour ", input$sel_ech_debris, " uniquement."),
                             type = "message")
          }
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # PANNEAU RÉSULTAT (DÉBRIS)
    # ════════════════════════════════════════════════════════════════════════
    
    output$titre_resultat <- renderUI({
      gates_trigger()
      ga <- nom_gate_actuel()
      if (is.null(ga)) "Résultat après filtration"
      else paste0("Résultat — gate : ", ga)
    })
    
    # Graphique résultat (ggplot via visualiser_debris)
    output$plot_gate_resultat <- renderPlot({
      gates_trigger()
      p  <- carrot_obj()
      ga <- nom_gate_actuel()
      req(input$sel_ech_debris)
      
      if (is.null(ga) || is.null(p$gates_history[[ga]])) {
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
        if (!is.null(p$post_retrait_bordures) && length(p$post_retrait_bordures) > 0) return(p$post_retrait_bordures)
        if (!is.null(p$post_PeacoQC) && length(p$post_PeacoQC) > 0) return(p$post_PeacoQC)
        if (!is.null(p$post_flowAI)  && length(p$post_flowAI)  > 0) return(p$post_flowAI)
        return(p$echantillons_traites)
      }
      if (!is.null(p$post_debris) && length(p$post_debris) > 0) return(p$post_debris)
      if (!is.null(p$post_retrait_bordures) && length(p$post_retrait_bordures) > 0) return(p$post_retrait_bordures)
      if (!is.null(p$post_PeacoQC) && length(p$post_PeacoQC) > 0) return(p$post_PeacoQC)
      if (!is.null(p$post_flowAI)  && length(p$post_flowAI)  > 0) return(p$post_flowAI)
      return(p$echantillons_traites)
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
      # Préserve la sélection en cours si elle est toujours valide (cf. même
      # correctif que pour les débris : éviter un retour silencieux au
      # premier échantillon après l'enregistrement d'un gate).
      actuel <- input$sel_ech_doublet
      sel <- if (!is.null(actuel) && actuel %in% names(src)) actuel else names(src)[1]
      selectInput(ns("sel_ech_doublet"), "Échantillon :",
                  choices = names(src), selected = sel)
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
      # Pas de sous-échantillonnage : la densité par binning raster reste rapide
      # indépendamment du nombre d'événements.
      df  <- as.data.frame(mat)
      colnames(df) <- c("X", "Y")
      attr(df, "canal_x") <- unname(canaux["x"])
      attr(df, "canal_y") <- unname(canaux["y"])
      df
    })
    
    # Calcule, pour l'échantillon/étape/axe courants, la forme du gate manuel
    # (polygone) déjà enregistrée — ou une liste vide sinon. Utilisée à la
    # fois pour synchroniser sommets_doublet_rv() et pour dessiner directement
    # la forme dans le rendu initial du graphique.
    calculer_soms_preload_doublet <- function(p, etape, axe, nom) {
      if (is.null(nom)) return(list())
      infos_gate <- if (etape == "FSC") p$gate_doublets_FSC[[nom]] else p$gate_doublets_SSC[[nom]]
      
      src <- obtenir_source_doublet(p, etape)
      canaux_actuels <- if (length(src) > 0) {
        resoudre_canaux_doublet(flowCore::colnames(src[[1]]), etape, axe)
      } else {
        c(x = NA, y = NA)
      }
      
      if (!is.null(infos_gate) && identical(infos_gate$type, "poly") && !is.null(infos_gate$gate) &&
          identical(unname(infos_gate$channels), unname(canaux_actuels))) {
        mat <- infos_gate$gate@boundaries
        lapply(seq_len(nrow(mat)), function(i) list(x = mat[i, 1], y = mat[i, 2]))
      } else {
        list()
      }
    }
    
    # Recharge (ou réinitialise) le polygone à chaque changement pertinent :
    # étape, axe, méthode ou échantillon. Un seul observer gère les deux cas
    # pour éviter toute course entre deux observers concurrents.
    observeEvent(list(input$etape_doublet, input$axe_doublet,
                      input$methode_doublet, input$sel_ech_doublet), {
                        p     <- carrot_obj()
                        etape <- input$etape_doublet %||% "FSC"
                        axe   <- input$axe_doublet %||% "H_A"
                        nom   <- input$sel_ech_doublet
                        req(nom)
                        
                        sommets_doublet_rv(calculer_soms_preload_doublet(p, etape, axe, nom))
                      }, ignoreInit = FALSE)
    
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
      lbl_x <- p_obj$obtenir_label(fcs, attr(df, "canal_x"))
      lbl_y <- p_obj$obtenir_label(fcs, attr(df, "canal_y"))
      
      # Densité par binning raster (rapide, indépendant du nombre d'événements affichés).
      # Rendue en petits ronds colorés par densité (mode "scatter" classique, même famille
      # que le tracé du gate, pour éviter tout ralentissement dû au mélange SVG/WebGL) plutôt
      # qu'en heatmap plotly : on évite ainsi tout lissage/interpolation qui brouillerait les
      # frontières entre zones.
      lim_x <- range(df$X, na.rm = TRUE)
      lim_y <- range(df$Y, na.rm = TRUE)
      dens  <- calculer_densite_raster(df$X, df$Y, lim_x, lim_y, res = RES_PIXELS_GATING, lissage = FALSE)
      
      plt <- plot_ly(source = ns("plot_gate_doublet_dessin"))
      if (!is.null(dens)) {
        plt <- plt %>%
          add_trace(data = dens, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    marker = list(symbol = "circle", size = TAILLE_PIXEL_GATING,
                                  color = ~densite, colorscale = COLORSCALE_DENSITE_PLOTLY,
                                  line = list(width = 0)),
                    hoverinfo = "none")
      } else {
        plt <- plt %>%
          add_trace(data = df, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    name = "Cellules",
                    marker = list(size = 2, color = "#0077b6", opacity = 0.4),
                    hoverinfo = "none")
      }
      
      # Précharge, dès le rendu initial, la forme du gate déjà enregistrée pour
      # cet échantillon (si elle existe), afin d'éviter toute course entre le
      # redessin complet et l'injection ultérieure via plotlyProxy.
      soms_init  <- calculer_soms_preload_doublet(p_obj, etape, input$axe_doublet %||% "H_A", input$sel_ech_doublet)
      trace_init <- construire_trace_et_shapes_gate(soms_init, df)
      
      plt %>%
        add_trace(x = trace_init$x, y = trace_init$y, type = "scatter",
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
          shapes = trace_init$shapes
        ) %>%
        config(displayModeBar = TRUE, editable = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE,
               # Désactive le double-clic natif de Plotly (reset des axes / des
               # formes) : ce plot a son propre gestionnaire de double-clic
               # (fermeture du polygone) et le comportement natif interférait
               # avec lui en effaçant le tracé du gate en cours d'édition.
               doubleClick = FALSE)
    })
    
    observeEvent(sommets_doublet_rv(), {
      soms  <- sommets_doublet_rv()
      proxy <- plotlyProxy(ns("plot_gate_doublet_dessin"), session)
      df    <- isolate(donnees_plot_doublet())
      trace <- construire_trace_et_shapes_gate(soms, df)
      
      plotlyProxyInvoke(proxy, "restyle", list(x = list(trace$x), y = list(trace$y)), list(1))
      plotlyProxyInvoke(proxy, "relayout", list(shapes = trace$shapes))
    }, ignoreNULL = FALSE)
    
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_doublet_dessin")), {
      ev   <- event_data("plotly_relayout", source = ns("plot_gate_doublet_dessin"))
      req(ev)
      
      # Sécurité : recharge la forme déjà enregistrée si sommets_doublet_rv()
      # n'a pas encore été synchronisée (voir même commentaire dans le
      # gestionnaire équivalent des débris).
      soms_actuels <- sommets_doublet_rv()
      if (length(soms_actuels) == 0) {
        p     <- carrot_obj()
        etape <- input$etape_doublet %||% "FSC"
        axe   <- input$axe_doublet %||% "H_A"
        soms_actuels <- calculer_soms_preload_doublet(p, etape, axe, input$sel_ech_doublet)
      }
      
      soms_maj <- appliquer_deplacement_shapes(ev, soms_actuels)
      if (!is.null(soms_maj)) sommets_doublet_rv(soms_maj)
    }, ignoreInit = TRUE)
    
    # Affiche le seuil (facteur/valeur) déjà enregistré pour l'échantillon
    # actuellement sélectionné, pour référence avant un nouvel enregistrement.
    output$ui_statut_doublet_auto <- renderUI({
      doublet_trigger()
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      nom   <- input$sel_ech_doublet
      req(nom)
      
      infos_gate <- if (etape == "FSC") p$gate_doublets_FSC[[nom]] else p$gate_doublets_SSC[[nom]]
      
      if (!is.null(infos_gate) && identical(infos_gate$type, "stat")) {
        div(class = "source-info-box",
            icon("check-circle"),
            " Seuil déjà enregistré pour cet échantillon : facteur = ", tags$b(infos_gate$facteur),
            " (seuil = ", round(infos_gate$seuil, 3), ").")
      } else if (!is.null(infos_gate)) {
        div(class = "source-info-box",
            icon("info-circle"),
            " Cet échantillon utilise actuellement un gate manuel (polygone) pour cette étape.")
      } else {
        NULL
      }
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # APPLICATION — MÉTHODE AUTOMATIQUE (STATISTIQUE)
    # ════════════════════════════════════════════════════════════════════════
    
    observeEvent(input$btn_appliquer_doublet_auto, {
      p     <- carrot_obj()
      etape <- input$etape_doublet %||% "FSC"
      fs    <- input$facteur_sensibilite_doublet %||% 4
      axe   <- input$axe_doublet %||% "H_A"
      req(input$sel_ech_doublet)
      
      gate_store <- if (etape == "FSC") p$gate_doublets_FSC else p$gate_doublets_SSC
      premier_enregistrement <- length(gate_store) == 0
      
      withProgress(message = paste0("Retrait automatique des doublets ", etape, "..."), value = 0.4, {
        tryCatch({
          if (etape == "FSC") {
            p$retirer_doublets_FSC(facteur_sensibilite = fs, axe_discrimination = axe,
                                   nom_echantillon = input$sel_ech_doublet)
          } else {
            p$retirer_doublets_SSC(facteur_sensibilite = fs, axe_discrimination = axe,
                                   nom_echantillon = input$sel_ech_doublet)
          }
          pipeline(p)
          doublet_trigger(doublet_trigger() + 1L)
          
          if (premier_enregistrement) {
            nb <- length(if (etape == "FSC") p$gate_doublets_FSC else p$gate_doublets_SSC)
            showNotification(
              paste0("Seuil ", etape, " (facteur = ", fs, ") appliqué par défaut à ", nb, " échantillon(s). ",
                     "Changez d'échantillon pour ajuster le seuil individuellement, puis ré-enregistrez."),
              type = "message", duration = 8
            )
          } else {
            showNotification(paste0("Seuil ", etape, " (facteur = ", fs, ") mis à jour pour ", input$sel_ech_doublet, " uniquement."),
                             type = "message")
          }
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
      req(input$sel_ech_doublet)
      
      p <- carrot_obj()
      
      xs <- sapply(soms, `[[`, "x")
      ys <- sapply(soms, `[[`, "y")
      idx_hull <- grDevices::chull(xs, ys)
      mat_poly <- cbind(xs[idx_hull], ys[idx_hull])
      
      gate_store <- if (etape == "FSC") p$gate_doublets_FSC else p$gate_doublets_SSC
      premier_enregistrement <- length(gate_store) == 0
      
      withProgress(message = paste0("Application du gate de doublets ", etape, "..."), value = 0.4, {
        tryCatch({
          if (etape == "FSC") {
            p$gate_les_doublets_FSC(points_utilisateur = mat_poly, axe_discrimination = axe,
                                    nom_echantillon = input$sel_ech_doublet)
          } else {
            p$gate_les_doublets_SSC(points_utilisateur = mat_poly, axe_discrimination = axe,
                                    nom_echantillon = input$sel_ech_doublet)
          }
          pipeline(p)
          doublet_trigger(doublet_trigger() + 1L)
          
          if (premier_enregistrement) {
            nb <- length(if (etape == "FSC") p$gate_doublets_FSC else p$gate_doublets_SSC)
            showNotification(
              paste0("Gate de doublets ", etape, " créé et appliqué par défaut à ", nb, " échantillon(s). ",
                     "Changez d'échantillon pour ajuster sa forme individuellement, puis ré-enregistrez."),
              type = "message", duration = 8
            )
          } else {
            showNotification(paste0("Gate de doublets ", etape, " mis à jour pour ", input$sel_ech_doublet, " uniquement."), type = "message")
          }
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
                              type_analyse     = etape),
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
    # (le polygone de gating est réinitialisé/rechargé séparément ci-dessous).
    observeEvent(input$canal_viabilite, {
      canal_transforme_rv(NULL)
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
      # Préserve la sélection en cours si elle est toujours valide (cf. même
      # correctif que pour les débris et les doublets).
      actuel <- input$sel_ech_viabilite
      noms   <- names(p$post_transformation)
      sel <- if (!is.null(actuel) && actuel %in% noms) actuel else noms[1]
      selectInput(ns("sel_ech_viabilite"), "Échantillon :",
                  choices = noms, selected = sel)
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
      # Pas de sous-échantillonnage : la densité par binning raster reste rapide
      # indépendamment du nombre d'événements.
      df  <- as.data.frame(mat)
      colnames(df) <- c("X", "Y")
      attr(df, "canal_x") <- cx
      attr(df, "canal_y") <- canal
      df
    })
    
    # Calcule, pour l'échantillon/canal courants, la forme du gate de
    # viabilité déjà enregistrée — ou une liste vide sinon. Utilisée à la fois
    # pour synchroniser sommets_viabilite_rv() et pour dessiner directement la
    # forme dans le rendu initial du graphique.
    calculer_soms_preload_viabilite <- function(p, nom, canal) {
      if (is.null(nom) || is.null(canal)) return(list())
      gate_sample <- p$gate_viabilite[[nom]]
      if (!is.null(gate_sample) && !is.null(gate_sample@boundaries) &&
          canal %in% colnames(gate_sample@boundaries)) {
        mat <- gate_sample@boundaries
        lapply(seq_len(nrow(mat)), function(i) list(x = mat[i, 1], y = mat[i, 2]))
      } else {
        list()
      }
    }
    
    # Recharge (ou réinitialise) le polygone à chaque changement pertinent :
    # canal ou échantillon. Un seul observer gère les deux cas pour éviter
    # toute course entre deux observers concurrents.
    observeEvent(list(input$canal_viabilite, input$sel_ech_viabilite), {
      p     <- carrot_obj()
      nom   <- input$sel_ech_viabilite
      canal <- input$canal_viabilite
      req(nom, canal)
      
      sommets_viabilite_rv(calculer_soms_preload_viabilite(p, nom, canal))
    }, ignoreInit = FALSE)
    
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
      lbl_x <- p_obj$obtenir_label(fcs, attr(df, "canal_x"))
      lbl_y <- p_obj$obtenir_label(fcs, attr(df, "canal_y"))
      
      # Densité par binning raster (rapide, indépendant du nombre d'événements affichés).
      # Rendue en petits ronds colorés par densité (mode "scatter" classique, même famille
      # que le tracé du gate, pour éviter tout ralentissement dû au mélange SVG/WebGL) plutôt
      # qu'en heatmap plotly : on évite ainsi tout lissage/interpolation qui brouillerait les
      # frontières entre zones.
      lim_x <- range(df$X, na.rm = TRUE)
      lim_y <- range(df$Y, na.rm = TRUE)
      dens  <- calculer_densite_raster(df$X, df$Y, lim_x, lim_y, res = RES_PIXELS_GATING, lissage = FALSE)
      
      plt <- plot_ly(source = ns("plot_gate_viabilite_dessin"))
      if (!is.null(dens)) {
        plt <- plt %>%
          add_trace(data = dens, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    marker = list(symbol = "circle", size = TAILLE_PIXEL_GATING,
                                  color = ~densite, colorscale = COLORSCALE_DENSITE_PLOTLY,
                                  line = list(width = 0)),
                    hoverinfo = "none")
      } else {
        plt <- plt %>%
          add_trace(data = df, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    name = "Cellules",
                    marker = list(size = 2, color = "#0077b6", opacity = 0.4),
                    hoverinfo = "none")
      }
      
      # Précharge, dès le rendu initial, la forme du gate déjà enregistrée pour
      # cet échantillon (si elle existe), afin d'éviter toute course entre le
      # redessin complet et l'injection ultérieure via plotlyProxy.
      soms_init  <- calculer_soms_preload_viabilite(p_obj, input$sel_ech_viabilite, input$canal_viabilite)
      trace_init <- construire_trace_et_shapes_gate(soms_init, df)
      
      plt %>%
        add_trace(x = trace_init$x, y = trace_init$y, type = "scatter",
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
          shapes = trace_init$shapes
        ) %>%
        config(displayModeBar = TRUE, editable = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE,
               # Désactive le double-clic natif de Plotly (reset des axes / des
               # formes) : ce plot a son propre gestionnaire de double-clic
               # (fermeture du polygone) et le comportement natif interférait
               # avec lui en effaçant le tracé du gate en cours d'édition.
               doubleClick = FALSE)
    })
    
    observeEvent(sommets_viabilite_rv(), {
      soms  <- sommets_viabilite_rv()
      proxy <- plotlyProxy(ns("plot_gate_viabilite_dessin"), session)
      df    <- isolate(donnees_plot_viabilite())
      trace <- construire_trace_et_shapes_gate(soms, df)
      
      plotlyProxyInvoke(proxy, "restyle", list(x = list(trace$x), y = list(trace$y)), list(1))
      plotlyProxyInvoke(proxy, "relayout", list(shapes = trace$shapes))
    }, ignoreNULL = FALSE)
    
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_viabilite_dessin")), {
      ev   <- event_data("plotly_relayout", source = ns("plot_gate_viabilite_dessin"))
      req(ev)
      
      # Sécurité : recharge la forme déjà enregistrée si sommets_viabilite_rv()
      # n'a pas encore été synchronisée (voir même commentaire dans le
      # gestionnaire équivalent des débris).
      soms_actuels <- sommets_viabilite_rv()
      if (length(soms_actuels) == 0) {
        p <- carrot_obj()
        soms_actuels <- calculer_soms_preload_viabilite(p, input$sel_ech_viabilite, input$canal_viabilite)
      }
      
      soms_maj <- appliquer_deplacement_shapes(ev, soms_actuels)
      if (!is.null(soms_maj)) sommets_viabilite_rv(soms_maj)
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
      req(input$sel_ech_viabilite)
      
      p   <- carrot_obj()
      fcs <- p$post_transformation[[input$sel_ech_viabilite]]
      req(!is.null(fcs))
      canal_fsc <- resoudre_canal_fsc_viabilite(flowCore::colnames(fcs))
      req(!is.na(canal_fsc))
      
      xs <- sapply(soms, `[[`, "x")
      ys <- sapply(soms, `[[`, "y")
      idx_hull <- grDevices::chull(xs, ys)
      mat_poly <- cbind(xs[idx_hull], ys[idx_hull])
      
      premier_enregistrement <- length(p$gate_viabilite) == 0
      
      withProgress(message = "Application du gate de viabilité...", value = 0.4, {
        tryCatch({
          p$retirer_les_cellules_mortes(
            canal_fsc          = canal_fsc,
            marqueur_viabilite  = canal,
            points_utilisateur  = mat_poly,
            nom_echantillon     = input$sel_ech_viabilite
          )
          pipeline(p)
          viabilite_trigger(viabilite_trigger() + 1L)
          
          if (premier_enregistrement) {
            nb <- length(p$gate_viabilite)
            showNotification(
              paste0("Gate de cellules vivantes créé et appliqué par défaut à ", nb, " échantillon(s). ",
                     "Changez d'échantillon pour ajuster sa forme individuellement, puis ré-enregistrez."),
              type = "message", duration = 8
            )
          } else {
            showNotification(paste0("Gate de cellules vivantes mis à jour pour ", input$sel_ech_viabilite, " uniquement."), type = "message")
          }
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
        p$visualiser_viabilite(nom_echantillon = input$sel_ech_viabilite),
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
    
    # Reconstruit, pour un échantillon donné, la chaîne complète des étapes de
    # prétraitement qu'il a traversées, dans l'ORDRE réel du pipeline
    # (compensées -> PeacoQC/flowAI -> bordures -> gates de débris nommés,
    # empilés dans l'ordre de création -> doublets FSC -> doublets SSC ->
    # viabilité), en ne retenant que les étapes réellement effectuées pour CET
    # échantillon (une étape absente est simplement omise, pas affichée à 0).
    # Chaque étape mémorise son propre effectif (n) ET l'effectif de référence
    # (ref) par rapport auquel calculer son pourcentage de conservation — cette
    # référence est toujours l'étape immédiatement précédente dans la chaîne
    # réellement traversée, pas la population totale d'origine (sinon un
    # pourcentage élevé à une étape tardive masquerait des pertes déjà
    # subies aux étapes précédentes).
    construire_recap <- function(p, nom) {
      etapes <- list()
      
      if (!is.null(p$echantillons_traites[[nom]])) {
        n <- nrow(flowCore::exprs(p$echantillons_traites[[nom]]))
        etapes[["Compensées"]] <- list(n = n, ref = n) # Point de départ : sa propre référence est lui-même (100%), puisqu'aucune étape antérieure n'existe
      }
      ref <- etapes[["Compensées"]]$n %||% NA
      
      if (!is.null(p$post_PeacoQC[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_PeacoQC[[nom]]))
        etapes[["PeacoQC"]] <- list(n = n, ref = ref)
      }
      if (!is.null(p$post_flowAI[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_flowAI[[nom]]))
        etapes[["flowAI"]] <- list(n = n, ref = ref) # PeacoQC et flowAI utilisent la MÊME référence (Compensées) : ce sont deux méthodes de QC alternatives, potentiellement appliquées toutes les deux indépendamment, pas chaînées l'une à l'autre
      }
      if (!is.null(p$post_retrait_bordures[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_retrait_bordures[[nom]]))
        etapes[["Bordures"]] <- list(n = n, ref = ref)
      }
      
      # Gates de débris nommés (empilés)
      ref_gate <- if (length(etapes) > 0) etapes[[length(etapes)]]$n else ref # Repart de la DERNIÈRE étape déjà enregistrée (quelle qu'elle soit) : à partir d'ici, chaque étape est chaînée séquentiellement à la précédente, pas toutes à la même référence d'origine
      for (nm_gate in names(p$gates_history)) { # Parcourt les gates de débris DANS L'ORDRE DE CRÉATION (l'ordre des noms dans la liste), pour respecter l'ordre réel d'application successive
        infos <- p$gates_history[[nm_gate]][[nom]]
        if (!is.null(infos)) {
          etapes[[nm_gate]] <- list(n = infos$n_apres, ref = infos$n_avant) # n_avant/n_apres déjà mémorisés lors de la création du gate (voir appliquer_gate_nomme() dans pipeline_cytometrie.R) : évite de devoir relire le flowFrame juste pour compter les lignes
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
        etapes[["Doublets SSC"]] <- list(n = n, ref = ref_gate) # Chaîné après FSC (pas après les débris directement) : SSC s'applique sur le résultat déjà filtré par FSC quand les deux sont utilisés
        ref_gate <- n
      }
      
      # Viabilité, chaînée à la suite des doublets/débris
      if (!is.null(p$post_viabilite[[nom]])) {
        n <- nrow(flowCore::exprs(p$post_viabilite[[nom]]))
        etapes[["Viabilité"]] <- list(n = n, ref = ref_gate) # ref_gate pointe ici vers la DERNIÈRE étape de doublets/débris réellement effectuée, quelle qu'elle soit (FSC seul, SSC seul, les deux, ou aucun)
        ref_gate <- n
      }
      
      etapes
    }
    
    # Transforme la liste d'étapes (construire_recap()) en affichage HTML : une
    # ligne par étape, avec effectif et pourcentage de conservation par rapport
    # à l'étape précédente, coloré selon le même code couleur d'alerte que le
    # reste de l'application (vert ≥90%, orange ≥70%, rouge <70%).
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
    
    # ════════════════════════════════════════════════════════════════════════
    # SECTION EXPORT (FCS post-débris / post-doublets / post-cellules mortes + session RDS)
    # ════════════════════════════════════════════════════════════════════════
    
    # Sélecteur multiple des échantillons disponibles pour l'export (union débris + doublets + viabilité, filtrée par les étapes cochées)
    output$ui_select_echantillons_export_pretrait <- renderUI({
      bordures_trigger(); gates_trigger(); doublet_trigger(); viabilite_trigger()
      p <- carrot_obj()
      
      etapes <- input$export_etapes_pretrait %||% c("debris", "doublets", "viabilite")
      noms_disponibles <- unique(c(
        if ("debris"    %in% etapes) names(p$post_debris)        else character(0),
        if ("doublets"  %in% etapes) names(p$post_doublets_final) else character(0),
        if ("viabilite" %in% etapes) names(p$post_viabilite)      else character(0)
      ))
      
      if (length(noms_disponibles) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"),
                   " Aucun échantillon disponible pour les étapes sélectionnées."))
      }
      
      selectizeInput(ns("export_echantillons_pretrait"), "Échantillons à exporter :",
                     choices = noms_disponibles, selected = noms_disponibles,
                     multiple = TRUE, options = list(plugins = list("remove_button")))
    })
    
    # Tableau récapitulatif de ce qui sera inclus dans l'export
    output$ui_recap_export_pretrait <- renderUI({
      bordures_trigger(); gates_trigger(); doublet_trigger(); viabilite_trigger()
      p <- carrot_obj()
      
      noms <- unique(c(names(p$post_debris), names(p$post_doublets_final), names(p$post_viabilite)))
      if (length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"),
                   " Aucun résultat de prétraitement disponible. Exécutez d'abord le retrait des débris, des doublets et/ou des cellules mortes."))
      }
      
      icone_ok   <- icon("circle-check", style = "color:#2e7d32;")
      icone_vide <- icon("circle-xmark", style = "color:#c62828;")
      
      lignes <- lapply(noms, function(nom) {
        a_debris    <- !is.null(p$post_debris[[nom]])
        a_doublets  <- !is.null(p$post_doublets_final[[nom]])
        a_viabilite <- !is.null(p$post_viabilite[[nom]])
        tags$tr(
          tags$td(strong(nom)),
          tags$td(style = "text-align:center;", if (a_debris) icone_ok else icone_vide),
          tags$td(style = "text-align:center;", if (a_doublets) icone_ok else icone_vide),
          tags$td(style = "text-align:center;", if (a_viabilite) icone_ok else icone_vide)
        )
      })
      
      tags$table(class = "table table-condensed", style = "font-size:12px;",
                 tags$thead(tags$tr(tags$th("Échantillon"), tags$th("Post débris"), tags$th("Post doublets"), tags$th("Post cellules mortes"))),
                 tags$tbody(lignes)
      )
    })
    
    # Export ZIP des fichiers FCS de prétraitement (étape(s) au choix de l'utilisateur)
    output$dl_export_fcs_pretrait <- downloadHandler(
      filename = function() paste0("Pretraitement_FCS_", format(Sys.Date(), "%Y%m%d"), ".zip"),
      content = function(file) {
        p <- carrot_obj()
        req(length(input$export_etapes_pretrait) > 0, length(input$export_echantillons_pretrait) > 0)
        
        dossier_temp <- file.path(tempdir(), paste0("export_fcs_pretrait_", as.integer(Sys.time())))
        dir.create(dossier_temp, recursive = TRUE)
        
        tryCatch({
          p$exporter_fcs_pretraitement(
            noms_echantillons = input$export_echantillons_pretrait,
            etapes            = input$export_etapes_pretrait,
            dossier_export    = dossier_temp
          )
        }, error = function(e) {
          showNotification(paste("Erreur export FCS :", conditionMessage(e)), type = "error")
        })
        
        fichiers <- list.files(dossier_temp, full.names = FALSE)
        req(length(fichiers) > 0)
        
        ancien_wd <- setwd(dossier_temp)
        on.exit(setwd(ancien_wd), add = TRUE)
        utils::zip(zipfile = file, files = fichiers)
      },
      contentType = "application/zip"
    )
    
    # Export RDS de la session complète de prétraitement (gates, seuils, canaux, figures)
    output$dl_export_rds_pretrait <- downloadHandler(
      filename = function() {
        nm <- trimws(input$pretrait_rds_filename)
        if (nchar(nm) == 0) nm <- "Pretraitement_Session_Complete.rds"
        if (!grepl("\\.rds$", nm, ignore.case = TRUE)) nm <- paste0(nm, ".rds")
        nm
      },
      content = function(file) {
        p <- carrot_obj()
        withProgress(message = "Sérialisation de la session...", value = 0.5, {
          tryCatch({
            p$sauvegarder_session_pretraitement_rds(nom_fichier = file)
          }, error = function(e) {
            showNotification(paste("Erreur export RDS :", conditionMessage(e)), type = "error")
          })
        })
      },
      contentType = "application/octet-stream"
    )
    
  })
}

# Opérateur null-coalescing utilitaire
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}