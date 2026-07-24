library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(DT)

RES_PIXELS_GATING_ANALYSE   <- 200 # Résolution de la grille de densité pour le dessin interactif des gates (même logique que dans le module Prétraitement)
TAILLE_PIXEL_GATING_ANALYSE <- 3   # Taille des marqueurs plotly représentant chaque pixel de densité

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

analyse_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    
    tags$style(HTML("
      .analyse-instr {
        background:#eef6fb; border-left:3px solid #0077b6;
        padding:8px 12px; border-radius:4px;
        font-size:12px; color:#444; margin-top:8px;
      }
      .analyse-instr b { color:#0077b6; }
    ")),

    tabsetPanel(
      id = ns("tabs_analyse"),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 1 — GATING
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("draw-polygon"), " Gating"),
        fluidRow(
          # ──────────────────────────────────────────────────────────────────
          # COLONNE GAUCHE — Paramètres du nouveau gate
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 wellPanel(
                   h4("Nouveau gate"),
                   textInput(ns("nom_gate"), "Nom du gate :", placeholder = "ex: Lymphocytes, CD4+..."),
                   
                   selectInput(ns("type_gate"), "Type de gate :",
                               choices = c("Polygone (tracé libre)" = "polygon", "Rectangle (bornes numériques)" = "rectangle")),
                   
                   uiOutput(ns("ui_gate_parent")),
                   div(class = "analyse-instr", icon("sitemap"),
                       " Choisir un parent permet de créer un gate sur une ",
                       tags$b("sous-population"), " (gating hiérarchique) au lieu de repartir de toutes les cellules."),
                   
                   hr(),
                   uiOutput(ns("ui_echantillon_ref_gating")),
                   fluidRow(
                     column(6, uiOutput(ns("ui_canal_x_gating"))),
                     column(6, uiOutput(ns("ui_canal_y_gating")))
                   ),
                   numericInput(ns("cofacteur_gating"), "Cofacteur (transformation Arcsinh) :",
                                value = 150, min = 1, step = 10),
                   div(class = "analyse-instr", icon("wave-square"),
                       " Le gating se fait sur des données transformées (Arcsinh) : ajustez le cofacteur pour bien séparer les populations sur les canaux choisis."),
                   
                   hr(),
                   conditionalPanel(
                     condition = paste0("input['", ns("type_gate"), "'] == 'polygon'"),
                     div(class = "analyse-instr", icon("hand-pointer"),
                         " Cliquez sur le graphique pour ajouter des sommets. ",
                         tags$b("Double-clic"), " pour fermer le polygone."),
                     br(),
                     fluidRow(
                       column(6, actionButton(ns("btn_undo_sommet_gate"),
                                              tagList(icon("undo"), " Annuler sommet"),
                                              class = "btn-default btn-sm", style = "width:100%;")),
                       column(6, actionButton(ns("btn_reset_gate_dessin"),
                                              tagList(icon("eraser"), " Effacer"),
                                              class = "btn-default btn-sm", style = "width:100%;"))
                     )
                   ),
                   conditionalPanel(
                     condition = paste0("input['", ns("type_gate"), "'] == 'rectangle'"),
                     div(class = "analyse-instr", icon("info-circle"),
                         " Le rectangle prévisualisé sur le graphique se met à jour en direct selon les bornes ci-dessous."),
                     fluidRow(
                       column(6, numericInput(ns("rect_xmin"), "X min :", value = 0)),
                       column(6, numericInput(ns("rect_xmax"), "X max :", value = 1))
                     ),
                     fluidRow(
                       column(6, numericInput(ns("rect_ymin"), "Y min :", value = 0)),
                       column(6, numericInput(ns("rect_ymax"), "Y max :", value = 1))
                     )
                   ),
                   
                   hr(),
                   actionButton(ns("btn_creer_gate"), tagList(icon("check-circle"), " Créer le gate"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE CENTRALE — Dessin interactif
          # ──────────────────────────────────────────────────────────────────
          column(width = 6,
                 box(title = tagList(icon("draw-polygon"), " Dessin interactif"), width = NULL,
                     status = "info", solidHeader = TRUE,
                     plotlyOutput(ns("plot_gate_dessin_analyse"), height = "520px")
                 )
          ),
          
          # ──────────────────────────────────────────────────────────────────
          # COLONNE DROITE — Gates créés
          # ──────────────────────────────────────────────────────────────────
          column(width = 3,
                 box(title = "Gates créés", width = NULL, status = "success", solidHeader = TRUE,
                     DTOutput(ns("table_gates")),
                     br(),
                     actionButton(ns("btn_supprimer_gate"), tagList(icon("trash"), " Supprimer le gate sélectionné"),
                                  class = "btn-danger btn-sm", style = "width:100%;"),
                     hr(),
                     h5("Résumé de la sélection"),
                     DTOutput(ns("table_resume_gate"))
                 )
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 2 — PCA
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("chart-scatter"), " PCA"),
        fluidRow(
          column(width = 3,
                 wellPanel(
                   uiOutput(ns("ui_gate_select_pca")),
                   uiOutput(ns("ui_canaux_select_pca")),
                   numericInput(ns("pca_n_components"), "Nombre de composantes :", value = 2, min = 2, max = 10, step = 1),
                   checkboxInput(ns("pca_centrer"), "Centrer les données", TRUE),
                   checkboxInput(ns("pca_reduire"), "Réduire (variance unitaire)", TRUE),
                   numericInput(ns("pca_max_cells"), "Sous-échantillonnage max (vide = aucun) :", value = NA, min = 1000, step = 1000),
                   numericInput(ns("pca_cofacteur"), "Cofacteur (transformation Arcsinh) :", value = 150, min = 1, step = 10),
                   actionButton(ns("btn_lancer_pca"), tagList(icon("play"), " Lancer la PCA"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          column(width = 9,
                 box(title = "Projection PCA", width = NULL, status = "info", solidHeader = TRUE,
                     plotlyOutput(ns("plot_pca"), height = "500px")),
                 box(title = "Variance expliquée", width = NULL, status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_variance_pca")))
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 3 — UMAP
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("braille"), " UMAP"),
        fluidRow(
          column(width = 3,
                 wellPanel(
                   uiOutput(ns("ui_gate_select_umap")),
                   uiOutput(ns("ui_canaux_select_umap")),
                   numericInput(ns("umap_n_neighbors"), "n_neighbors :", value = 15, min = 2, step = 1),
                   numericInput(ns("umap_min_dist"), "min_dist :", value = 0.1, min = 0, max = 1, step = 0.05),
                   selectInput(ns("umap_metric"), "Métrique :",
                               choices = c("euclidean", "cosine", "manhattan", "correlation", "hamming")),
                   numericInput(ns("umap_max_cells"), "Sous-échantillonnage max :", value = 50000, min = 1000, step = 1000),
                   numericInput(ns("umap_cofacteur"), "Cofacteur (transformation Arcsinh) :", value = 150, min = 1, step = 10),
                   actionButton(ns("btn_lancer_umap"), tagList(icon("play"), " Lancer l'UMAP"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          column(width = 9,
                 box(title = "Projection UMAP", width = NULL, status = "info", solidHeader = TRUE,
                     uiOutput(ns("ui_couleur_umap")),
                     plotlyOutput(ns("plot_umap"), height = "550px"))
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 4 — t-SNE
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("braille"), " t-SNE"),
        fluidRow(
          column(width = 3,
                 wellPanel(
                   uiOutput(ns("ui_gate_select_tsne")),
                   uiOutput(ns("ui_canaux_select_tsne")),
                   numericInput(ns("tsne_perplexity"), "Perplexité :", value = 30, min = 5, step = 1),
                   numericInput(ns("tsne_theta"), "Theta (Barnes-Hut) :", value = 0.5, min = 0, max = 1, step = 0.05),
                   numericInput(ns("tsne_max_iter"), "Itérations max :", value = 1000, min = 100, step = 100),
                   numericInput(ns("tsne_max_cells"), "Sous-échantillonnage max :", value = 20000, min = 1000, step = 1000),
                   numericInput(ns("tsne_cofacteur"), "Cofacteur (transformation Arcsinh) :", value = 150, min = 1, step = 10),
                   actionButton(ns("btn_lancer_tsne"), tagList(icon("play"), " Lancer le t-SNE"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          column(width = 9,
                 box(title = "Projection t-SNE", width = NULL, status = "info", solidHeader = TRUE,
                     uiOutput(ns("ui_couleur_tsne")),
                     plotlyOutput(ns("plot_tsne"), height = "550px"))
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 5 — CLUSTERING (FlowSOM)
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("diagram-project"), " Clustering (FlowSOM)"),
        fluidRow(
          column(width = 3,
                 wellPanel(
                   uiOutput(ns("ui_gate_select_cluster")),
                   uiOutput(ns("ui_reutiliser_donnees_cluster")),
                   conditionalPanel(
                     condition = paste0("input['", ns("reutiliser_donnees_cluster"), "'] == 'AUCUN'"),
                     uiOutput(ns("ui_canaux_select_cluster")),
                     numericInput(ns("fsom_max_cells"), "Sous-échantillonnage max (vide = aucun) :", value = NA, min = 1000, step = 1000),
                     numericInput(ns("fsom_cofacteur"), "Cofacteur (transformation Arcsinh) :", value = 150, min = 1, step = 10)
                   ),
                   conditionalPanel(
                     condition = paste0("input['", ns("reutiliser_donnees_cluster"), "'] != 'AUCUN'"),
                     div(class = "analyse-instr", icon("check-circle"),
                         " Le clustering portera exactement sur les mêmes cellules (mêmes canaux, même sous-échantillonnage) que la projection choisie, ",
                         "garantissant un alignement parfait pour la superposition ci-dessous.")
                   ),
                   fluidRow(
                     column(6, numericInput(ns("fsom_xdim"), "Largeur grille SOM :", value = 10, min = 2, step = 1)),
                     column(6, numericInput(ns("fsom_ydim"), "Hauteur grille SOM :", value = 10, min = 2, step = 1))
                   ),
                   numericInput(ns("fsom_n_metaclusters"), "Nombre de métaclusters :", value = 10, min = 2, step = 1),
                   actionButton(ns("btn_lancer_cluster"), tagList(icon("play"), " Lancer le clustering"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          column(width = 9,
                 box(title = "Répartition des métaclusters", width = NULL, status = "info", solidHeader = TRUE,
                     plotlyOutput(ns("plot_cluster_repartition"), height = "350px")),
                 box(title = "Heatmap d'expression médiane par métacluster", width = NULL, status = "primary", solidHeader = TRUE,
                     div(class = "analyse-instr", icon("info-circle"),
                         " Expression médiane (données transformées) de chaque canal, par métacluster — z-scorée par colonne pour comparer les canaux entre eux. ",
                         "Sert à identifier biologiquement chaque métacluster avant de le renommer ci-dessous."),
                     plotlyOutput(ns("plot_heatmap_clusters"), height = "400px")),
                 box(title = "Annoter les métaclusters", width = NULL, status = "success", solidHeader = TRUE,
                     div(class = "analyse-instr", icon("tag"),
                         " Donnez un nom biologique à chaque métacluster une fois son identité déterminée (ex: \u00ab Lymphocytes T CD4+ \u00bb)."),
                     uiOutput(ns("ui_annotation_metaclusters")),
                     actionButton(ns("btn_enregistrer_annotations"), tagList(icon("save"), " Enregistrer les noms"),
                                  class = "btn-success btn-sm")),
                 box(title = "Projection colorée par métacluster", width = NULL, status = "warning", solidHeader = TRUE,
                     uiOutput(ns("ui_choix_embedding_cluster")),
                     plotlyOutput(ns("plot_cluster_embedding"), height = "450px"))
          )
        )
      ),
      
      # ══════════════════════════════════════════════════════════════════════
      # ONGLET 6 — COMPARAISON GROUPES/CONDITIONS
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = tagList(icon("scale-balanced"), " Comparaison Groupes"),
        fluidRow(
          column(width = 3,
                 wellPanel(
                   h4("Groupes/conditions"),
                   div(class = "analyse-instr", icon("info-circle"),
                       " Les groupes sont assignés dans l'onglet ", tags$b("Import > 3. Groupes d'Échantillons"),
                       " et réutilisés directement ici — inutile de les recréer."),
                   uiOutput(ns("ui_resume_groupes_existants"))
                 )
          ),
          column(width = 3,
                 wellPanel(
                   h4("Comparaison"),
                   uiOutput(ns("ui_gate_select_comparaison")),
                   selectInput(ns("comparaison_variable"), "Variable à comparer :",
                               choices = c("% du parent (population du gate)" = "pourcentage",
                                           "MFI d'un canal (au sein du gate)" = "MFI")),
                   conditionalPanel(
                     condition = paste0("input['", ns("comparaison_variable"), "'] == 'MFI'"),
                     uiOutput(ns("ui_canal_comparaison"))
                   ),
                   actionButton(ns("btn_lancer_comparaison"), tagList(icon("play"), " Comparer les groupes"),
                                class = "btn-success", style = "width:100%; font-weight:bold;")
                 )
          ),
          column(width = 6,
                 box(title = "Résultat de la comparaison", width = NULL, status = "info", solidHeader = TRUE,
                     uiOutput(ns("ui_resultat_comparaison")),
                     plotlyOutput(ns("plot_comparaison_groupes"), height = "450px"))
          )
        )
      )
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════

analyse_server <- function(id, pipeline, pipeline_version) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    carrot_obj <- reactive({
      pipeline_version()
      pipeline()
    })
    
    # ── Signaux locaux ──────────────────────────────────────────────────────
    gates_trigger   <- reactiveVal(0L)
    pca_trigger     <- reactiveVal(0L)
    umap_trigger    <- reactiveVal(0L)
    tsne_trigger    <- reactiveVal(0L)
    cluster_trigger <- reactiveVal(0L)
    
    sommets_gate_rv <- reactiveVal(list())
    
    # ════════════════════════════════════════════════════════════════════════
    # ONGLET GATING
    # ════════════════════════════════════════════════════════════════════════
    
    # Résout la source de données à afficher/gater selon le parent choisi
    # (population totale si aucun parent, sinon la sous-population déjà
    # sélectionnée par le gate parent — gating hiérarchique).
    obtenir_source_gating <- function(p, parent) {
      if (is.null(parent) || identical(parent, "AUCUN") || nchar(parent) == 0) {
        p$get_derniere_source()
      } else {
        tryCatch(p$resoudre_population_gate(parent), error = function(e) list())
      }
    }
    
    output$ui_gate_parent <- renderUI({
      gates_trigger()
      p <- carrot_obj()
      noms <- names(p$gates_personnalisees)
      selectInput(ns("gate_parent_select"), "Population parente :",
                  choices = c("Population totale (aucun parent)" = "AUCUN", noms))
    })
    
    output$ui_echantillon_ref_gating <- renderUI({
      gates_trigger()
      p <- carrot_obj()
      src <- obtenir_source_gating(p, input$gate_parent_select)
      if (length(src) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Aucune cellule disponible pour cette population parente."))
      }
      selectInput(ns("echantillon_ref_gating"), "Échantillon (référence pour le tracé) :", choices = names(src))
    })
    
    output$ui_canal_x_gating <- renderUI({
      p   <- carrot_obj()
      src <- obtenir_source_gating(p, input$gate_parent_select)
      req(length(src) > 0, input$echantillon_ref_gating %in% names(src))
      fcs_ref <- src[[input$echantillon_ref_gating]]
      cols    <- flowCore::colnames(fcs_ref)
      choix   <- stats::setNames(cols, vapply(cols, function(c) p$get_label(fcs_ref, c), character(1))) # Affiche "canal | marqueur" dans le menu, valeur = nom technique du canal
      selectInput(ns("canal_x_gating"), "Canal X :", choices = choix, selected = cols[1])
    })
    
    output$ui_canal_y_gating <- renderUI({
      p   <- carrot_obj()
      src <- obtenir_source_gating(p, input$gate_parent_select)
      req(length(src) > 0, input$echantillon_ref_gating %in% names(src))
      fcs_ref <- src[[input$echantillon_ref_gating]]
      cols    <- flowCore::colnames(fcs_ref)
      choix   <- stats::setNames(cols, vapply(cols, function(c) p$get_label(fcs_ref, c), character(1)))
      selectInput(ns("canal_y_gating"), "Canal Y :", choices = choix, selected = cols[min(2, length(cols))])
    })
    
    # Le gate "actif" (préchargé au changement d'échantillon, ciblé à
    # l'enregistrement) est celui dont le nom est actuellement saisi dans le
    # champ "Nom du gate" — même principe que dans module_pretraitement.R.
    nom_gate_actuel <- reactive({
      nm <- trimws(input$nom_gate %||% "")
      if (nchar(nm) == 0) NULL else nm
    })
    
    donnees_plot_gating <- reactive({
      gates_trigger()
      p <- carrot_obj()
      req(input$echantillon_ref_gating, input$canal_x_gating, input$canal_y_gating, input$cofacteur_gating)
      src <- obtenir_source_gating(p, input$gate_parent_select)
      req(length(src) > 0, input$echantillon_ref_gating %in% names(src))
      
      fcs <- src[[input$echantillon_ref_gating]]
      cx  <- input$canal_x_gating
      cy  <- input$canal_y_gating
      req(cx %in% flowCore::colnames(fcs), cy %in% flowCore::colnames(fcs))
      
      # Le gating se fait sur des données transformées (Arcsinh) : on applique ici
      # la même transformation (et le même cofacteur) que celle utilisée par
      # resoudre_population_gate() lors de la création du gate, pour que ce qui
      # est dessiné corresponde exactement à ce qui sera filtré.
      trans_axes <- flowCore::transformList(c(cx, cy), flowCore::arcsinhTransform(a = 0, b = 1 / input$cofacteur_gating, c = 0))
      fcs_trans  <- flowCore::transform(fcs, trans_axes)
      
      mat <- flowCore::exprs(fcs_trans)[, c(cx, cy), drop = FALSE]
      df  <- as.data.frame(mat)
      colnames(df) <- c("X", "Y")
      df
    })
    
    # Construit, à partir d'une liste de sommets et du data.frame affiché, les
    # coordonnées de la trace fermée (polygone) et les formes de poignées
    # (cercles) correspondantes — utilisé à la fois pour le rendu initial du
    # graphique et pour la mise à jour via plotlyProxy lors de l'édition.
    construire_trace_et_shapes_gate_analyse <- function(soms, df) {
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
    
    # Calcule, pour l'échantillon et le gate actif courants, la forme déjà
    # enregistrée à précharger (polygone uniquement — le rectangle a son
    # propre mécanisme de bornes numériques) — ou une liste vide sinon.
    calculer_soms_preload_gate_analyse <- function(p, nom) {
      ga <- isolate(nom_gate_actuel())
      if (is.null(ga) || is.null(nom) || is.null(p$gates_personnalisees[[ga]])) return(list())
      if (!identical(p$gates_personnalisees[[ga]]$type, "polygon")) return(list())
      
      mat_forme <- p$obtenir_forme_gate(ga, nom)
      if (is.null(mat_forme)) return(list())
      lapply(seq_len(nrow(mat_forme)), function(i) list(x = mat_forme[i, 1], y = mat_forme[i, 2]))
    }
    
    # Recharge (ou réinitialise) le polygone à chaque changement pertinent :
    # échantillon, canaux, source parente, nom de gate ou cofacteur. Un seul
    # observer gère les deux cas pour éviter toute course entre observers
    # concurrents (même principe que module_pretraitement.R) :
    #  - si un gate polygonal de ce nom existe déjà pour CET échantillon, on
    #    recharge sa forme (modifiable) ;
    #  - sinon, le canevas repart vide pour un nouveau tracé.
    observeEvent(list(input$echantillon_ref_gating, input$canal_x_gating, input$canal_y_gating,
                      input$gate_parent_select, input$type_gate, input$cofacteur_gating, input$nom_gate), {
                        p   <- carrot_obj()
                        nom <- input$echantillon_ref_gating
                        req(nom)
                        sommets_gate_rv(calculer_soms_preload_gate_analyse(p, nom))
                      }, ignoreInit = FALSE)
    
    output$plot_gate_dessin_analyse <- renderPlotly({
      df <- donnees_plot_gating()
      req(nrow(df) > 0)
      
      p   <- carrot_obj()
      src <- obtenir_source_gating(p, input$gate_parent_select)
      fcs_ref <- src[[input$echantillon_ref_gating]]
      lbl_x <- p$get_label(fcs_ref, input$canal_x_gating) # Combine canal technique et marqueur biologique (ex: "PE-A | CD3"), cohérent avec le reste de l'application
      lbl_y <- p$get_label(fcs_ref, input$canal_y_gating)
      
      lim_x <- range(df$X, na.rm = TRUE)
      lim_y <- range(df$Y, na.rm = TRUE)
      dens  <- calculer_densite_raster(df$X, df$Y, lim_x, lim_y, res = RES_PIXELS_GATING_ANALYSE, lissage = FALSE)
      
      plt <- plot_ly(source = ns("plot_gate_dessin_analyse"))
      if (!is.null(dens)) {
        plt <- plt %>%
          add_trace(data = dens, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    marker = list(symbol = "circle", size = TAILLE_PIXEL_GATING_ANALYSE,
                                  color = ~densite, colorscale = COLORSCALE_DENSITE_PLOTLY,
                                  line = list(width = 0)),
                    hoverinfo = "none")
      } else {
        plt <- plt %>%
          add_trace(data = df, x = ~X, y = ~Y, type = "scatter", mode = "markers",
                    marker = list(size = 2, color = "#0077b6", opacity = 0.4), hoverinfo = "none")
      }
      
      # Précharge, dès le rendu initial, la forme déjà enregistrée pour cet
      # échantillon (le cas échéant) — évite toute course entre le redessin
      # complet et l'injection ultérieure via plotlyProxy.
      soms_init  <- calculer_soms_preload_gate_analyse(p, input$echantillon_ref_gating)
      trace_init <- construire_trace_et_shapes_gate_analyse(soms_init, df)
      
      # Aperçu du rectangle en direct si ce type est sélectionné
      shapes_init <- trace_init$shapes
      if (identical(input$type_gate, "rectangle") &&
          !is.null(input$rect_xmin) && !is.null(input$rect_xmax) &&
          !is.null(input$rect_ymin) && !is.null(input$rect_ymax)) {
        shapes_init <- list(list(type = "rect",
                                 x0 = input$rect_xmin, x1 = input$rect_xmax,
                                 y0 = input$rect_ymin, y1 = input$rect_ymax,
                                 line = list(color = "#e65100", width = 2),
                                 fillcolor = "rgba(230,101,0,0.15)"))
      }
      
      plt %>%
        add_trace(x = trace_init$x, y = trace_init$y, type = "scatter", mode = "lines+markers",
                  name = "Gate en cours",
                  line   = list(color = "#e65100", width = 2, dash = "dot"),
                  marker = list(size = 5, color = "#e65100", symbol = "circle"),
                  hoverinfo = "x+y") %>%
        event_register("plotly_click") %>%
        event_register("plotly_doubleclick") %>%
        event_register("plotly_relayout") %>%
        layout(
          dragmode = "zoom",
          xaxis  = list(title = lbl_x),
          yaxis  = list(title = lbl_y),
          legend = list(orientation = "h", y = -0.15),
          margin = list(b = 60),
          shapes = shapes_init
        ) %>%
        config(displayModeBar = TRUE, editable = TRUE,
               modeBarButtonsToRemove = list("lasso2d", "select2d"),
               displaylogo = FALSE, doubleClick = FALSE)
    })
    
    # Ajout de sommet (clic), uniquement en mode polygone
    observeEvent(event_data("plotly_click", source = ns("plot_gate_dessin_analyse")), {
      req(identical(input$type_gate, "polygon"))
      ev <- event_data("plotly_click", source = ns("plot_gate_dessin_analyse"))
      req(ev)
      soms <- sommets_gate_rv()
      soms[[length(soms) + 1]] <- list(x = ev$x, y = ev$y)
      sommets_gate_rv(soms)
    })
    
    observeEvent(input$btn_undo_sommet_gate, {
      soms <- sommets_gate_rv()
      if (length(soms) > 0) sommets_gate_rv(soms[-length(soms)])
    })
    observeEvent(input$btn_reset_gate_dessin, { sommets_gate_rv(list()) })
    
    # Pousse les sommets courants vers le graphique déjà rendu, sans redraw complet
    observeEvent(sommets_gate_rv(), {
      soms  <- sommets_gate_rv()
      proxy <- plotlyProxy(ns("plot_gate_dessin_analyse"), session)
      df    <- isolate(donnees_plot_gating())
      trace <- construire_trace_et_shapes_gate_analyse(soms, df)
      
      plotlyProxyInvoke(proxy, "restyle", list(x = list(trace$x), y = list(trace$y)), list(1))
      plotlyProxyInvoke(proxy, "relayout", list(shapes = trace$shapes))
    }, ignoreNULL = FALSE)
    
    # Déplacement d'un sommet (drag) : le centre réel du cercle est la moyenne
    # de ses deux coins opposés (x0+x1)/2, (y0+y1)/2 — jamais x0/y0 seuls.
    observeEvent(event_data("plotly_relayout", source = ns("plot_gate_dessin_analyse")), {
      ev <- event_data("plotly_relayout", source = ns("plot_gate_dessin_analyse"))
      req(ev)
      soms <- sommets_gate_rv()
      
      indices <- regmatches(names(ev), regexpr("(?<=shapes\\[)[0-9]+(?=\\]\\.)", names(ev), perl = TRUE))
      indices <- unique(as.integer(indices[nchar(indices) > 0]))
      req(length(indices) > 0)
      
      modifie <- FALSE
      for (idx0 in indices) {
        i <- idx0 + 1L
        if (i < 1 || i > length(soms)) next
        prefix <- paste0("shapes[", idx0, "].")
        x0 <- ev[[paste0(prefix, "x0")]]; x1 <- ev[[paste0(prefix, "x1")]]
        y0 <- ev[[paste0(prefix, "y0")]]; y1 <- ev[[paste0(prefix, "y1")]]
        if (!is.null(x0) && !is.null(x1)) { soms[[i]]$x <- (x0 + x1) / 2; modifie <- TRUE }
        if (!is.null(y0) && !is.null(y1)) { soms[[i]]$y <- (y0 + y1) / 2; modifie <- TRUE }
      }
      if (modifie) sommets_gate_rv(soms)
    }, ignoreInit = TRUE)
    
    # Rafraîchit l'aperçu du rectangle en direct quand les bornes numériques changent
    observeEvent(list(input$rect_xmin, input$rect_xmax, input$rect_ymin, input$rect_ymax), {
      req(identical(input$type_gate, "rectangle"))
      proxy <- plotlyProxy(ns("plot_gate_dessin_analyse"), session)
      plotlyProxyInvoke(proxy, "relayout", list(shapes = list(list(
        type = "rect",
        x0 = input$rect_xmin, x1 = input$rect_xmax,
        y0 = input$rect_ymin, y1 = input$rect_ymax,
        line = list(color = "#e65100", width = 2),
        fillcolor = "rgba(230,101,0,0.15)"
      ))))
    }, ignoreInit = TRUE)
    
    observeEvent(input$btn_creer_gate, {
      nom <- trimws(input$nom_gate %||% "")
      if (nchar(nom) == 0) {
        showNotification("Donnez un nom au gate avant de le créer.", type = "error")
        return(invisible(NULL))
      }
      
      p      <- carrot_obj()
      parent <- if (is.null(input$gate_parent_select) || identical(input$gate_parent_select, "AUCUN")) NULL else input$gate_parent_select
      cx     <- input$canal_x_gating
      cy     <- input$canal_y_gating
      req(cx, cy, input$cofacteur_gating, input$echantillon_ref_gating)
      
      # Un gate déjà créé (au moins une forme enregistrée) ne verra sa forme mise
      # à jour QUE pour l'échantillon actuellement affiché (ajustement par
      # sous-population/échantillon) ; la toute première création s'applique par
      # défaut à tous les échantillons — logique gérée directement dans creer_gate().
      premiere_creation <- is.null(p$gates_personnalisees[[nom]]) || length(p$gates_personnalisees[[nom]]$formes) == 0
      
      tryCatch({
        if (identical(input$type_gate, "polygon")) {
          soms <- sommets_gate_rv()
          if (length(soms) < 3) stop("Tracez au moins 3 sommets pour définir un polygone.")
          xs <- sapply(soms, `[[`, "x")
          ys <- sapply(soms, `[[`, "y")
          idx_hull <- grDevices::chull(xs, ys)
          mat_poly <- cbind(xs[idx_hull], ys[idx_hull])
          p$creer_gate(nom_gate = nom, type = "polygon", axes = c(cx, cy), points = mat_poly,
                       gate_parent = parent, cofacteur = input$cofacteur_gating,
                       nom_echantillon = input$echantillon_ref_gating)
        } else {
          req(input$rect_xmin, input$rect_xmax, input$rect_ymin, input$rect_ymax)
          if (input$rect_xmin >= input$rect_xmax || input$rect_ymin >= input$rect_ymax) {
            stop("Les bornes minimales doivent être strictement inférieures aux bornes maximales.")
          }
          p$creer_gate(nom_gate = nom, type = "rectangle", axes = c(cx, cy),
                       points = c(input$rect_xmin, input$rect_xmax, input$rect_ymin, input$rect_ymax),
                       gate_parent = parent, cofacteur = input$cofacteur_gating,
                       nom_echantillon = input$echantillon_ref_gating)
        }
        pipeline(p)
        gates_trigger(gates_trigger() + 1L)
        
        if (premiere_creation) {
          nb <- length(p$gates_personnalisees[[nom]]$formes)
          showNotification(paste0("Gate '", nom, "' créé et appliqué par défaut à ", nb, " échantillon(s). ",
                                  "Changez d'échantillon pour ajuster sa forme individuellement, puis ré-enregistrez."),
                           type = "message", duration = 8)
        } else {
          showNotification(paste0("Gate '", nom, "' mis à jour pour ", input$echantillon_ref_gating, " uniquement."), type = "message")
        }
      }, error = function(e) showNotification(conditionMessage(e), type = "error"))
    })
    
    output$table_gates <- DT::renderDT({
      gates_trigger()
      p     <- carrot_obj()
      infos <- p$gates_personnalisees
      
      if (length(infos) == 0) {
        return(DT::datatable(data.frame(Message = "Aucun gate créé pour l'instant."),
                             rownames = FALSE, options = list(dom = "t")))
      }
      
      df <- do.call(rbind, lapply(names(infos), function(nm) {
        g <- infos[[nm]]
        data.frame(
          Nom    = nm,
          Axes   = paste(g$axes, collapse = " / "),
          Type   = if (identical(g$type, "polygon")) "Polygone" else "Rectangle",
          Parent = if (is.null(g$gate_parent)) "\u2014" else g$gate_parent,
          stringsAsFactors = FALSE
        )
      }))
      
      DT::datatable(df, rownames = FALSE, selection = "single",
                    options = list(dom = "t", paging = FALSE, scrollX = TRUE))
    })
    
    observeEvent(input$btn_supprimer_gate, {
      sel <- input$table_gates_rows_selected
      req(sel)
      p    <- carrot_obj()
      noms <- names(p$gates_personnalisees)
      req(sel <= length(noms))
      
      nom_a_supprimer <- noms[sel]
      tryCatch({
        p$supprimer_gate(nom_a_supprimer, cascade = TRUE)
        pipeline(p)
        gates_trigger(gates_trigger() + 1L)
        showNotification(paste0("Gate '", nom_a_supprimer, "' supprimé."), type = "message")
      }, error = function(e) showNotification(conditionMessage(e), type = "error"))
    })
    
    output$table_resume_gate <- DT::renderDT({
      gates_trigger()
      p    <- carrot_obj()
      sel  <- input$table_gates_rows_selected
      noms <- names(p$gates_personnalisees)
      
      if (is.null(sel) || sel > length(noms)) {
        return(DT::datatable(data.frame(Message = "Sélectionnez un gate dans le tableau ci-dessus."),
                             rownames = FALSE, options = list(dom = "t")))
      }
      
      resume <- tryCatch(p$resumer_gate(noms[sel]), error = function(e) NULL)
      req(!is.null(resume))
      DT::datatable(resume, rownames = FALSE,
                    colnames = c("Échantillon", "N. événements", "% du parent"),
                    options = list(dom = "t", paging = FALSE, scrollX = TRUE))
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # ONGLET PCA
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_gate_select_pca <- renderUI({
      gates_trigger()
      p    <- carrot_obj()
      noms <- names(p$gates_personnalisees)
      if (length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Créez d'abord un gate dans l'onglet Gating."))
      }
      selectInput(ns("gate_select_pca"), "Gate à analyser :", choices = noms)
    })
    
    output$ui_canaux_select_pca <- renderUI({
      p <- carrot_obj()
      req(input$gate_select_pca)
      pop <- tryCatch(p$resoudre_population_gate(input$gate_select_pca), error = function(e) list())
      req(length(pop) > 0)
      cols <- flowCore::colnames(pop[[1]])
      fluo <- cols[!grepl("FSC|SSC|Time", cols, ignore.case = TRUE)]
      selectizeInput(ns("canaux_select_pca"), "Canaux à utiliser (vide = tous) :",
                     choices = fluo, multiple = TRUE)
    })
    
    observeEvent(input$btn_lancer_pca, {
      p <- carrot_obj()
      req(input$gate_select_pca)
      canaux    <- if (length(input$canaux_select_pca) == 0) NULL else input$canaux_select_pca
      max_cells <- if (is.null(input$pca_max_cells) || is.na(input$pca_max_cells)) NULL else input$pca_max_cells
      
      withProgress(message = "Calcul de la PCA...", value = 0.4, {
        tryCatch({
          p$projection_PCA(nom_gate = input$gate_select_pca, canaux = canaux,
                           n_components = input$pca_n_components, centrer = input$pca_centrer,
                           reduire = input$pca_reduire, sous_echantillonnage_max = max_cells,
                           cofacteur = input$pca_cofacteur)
          pipeline(p)
          pca_trigger(pca_trigger() + 1L)
          showNotification("PCA calculée.", type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    output$plot_pca <- renderPlotly({
      pca_trigger()
      p <- carrot_obj()
      req(input$gate_select_pca)
      res <- p$analyses_pca[[input$gate_select_pca]]
      validate(need(!is.null(res), "Lancez la PCA pour afficher la projection."))
      
      df <- as.data.frame(res$embedding)
      df$Echantillon <- res$echantillon_origine
      noms_axes <- colnames(res$embedding)
      
      plot_ly(df, x = df[[1]], y = df[[2]], color = ~Echantillon,
              type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
        layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
    })
    
    output$ui_variance_pca <- renderUI({
      pca_trigger()
      p <- carrot_obj()
      req(input$gate_select_pca)
      res <- p$analyses_pca[[input$gate_select_pca]]
      req(!is.null(res))
      
      noms_axes <- colnames(res$embedding)
      tagList(
        tags$ul(lapply(seq_along(res$variance_expliquee), function(i) {
          tags$li(paste0(noms_axes[i], " : ", round(res$variance_expliquee[i] * 100, 1), "% de variance expliquée"))
        })),
        tags$p(tags$b(paste0("Total sur les composantes conservées : ", round(sum(res$variance_expliquee) * 100, 1), "%")))
      )
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # ONGLET UMAP
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_gate_select_umap <- renderUI({
      gates_trigger()
      p    <- carrot_obj()
      noms <- names(p$gates_personnalisees)
      if (length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Créez d'abord un gate dans l'onglet Gating."))
      }
      selectInput(ns("gate_select_umap"), "Gate à analyser :", choices = noms)
    })
    
    output$ui_canaux_select_umap <- renderUI({
      p <- carrot_obj()
      req(input$gate_select_umap)
      pop <- tryCatch(p$resoudre_population_gate(input$gate_select_umap), error = function(e) list())
      req(length(pop) > 0)
      cols <- flowCore::colnames(pop[[1]])
      fluo <- cols[!grepl("FSC|SSC|Time", cols, ignore.case = TRUE)]
      selectizeInput(ns("canaux_select_umap"), "Canaux à utiliser (vide = tous) :",
                     choices = fluo, multiple = TRUE)
    })
    
    observeEvent(input$btn_lancer_umap, {
      p <- carrot_obj()
      req(input$gate_select_umap)
      canaux <- if (length(input$canaux_select_umap) == 0) NULL else input$canaux_select_umap
      
      withProgress(message = "Calcul de l'UMAP en cours (peut prendre quelques minutes)...", value = 0.3, {
        tryCatch({
          p$projection_UMAP(nom_gate = input$gate_select_umap, canaux = canaux,
                            n_neighbors = input$umap_n_neighbors, min_dist = input$umap_min_dist,
                            metric = input$umap_metric, sous_echantillonnage_max = input$umap_max_cells,
                            cofacteur = input$umap_cofacteur)
          pipeline(p)
          umap_trigger(umap_trigger() + 1L)
          showNotification("UMAP calculée.", type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    output$ui_couleur_umap <- renderUI({
      umap_trigger(); cluster_trigger()
      p <- carrot_obj()
      req(input$gate_select_umap)
      res <- p$analyses_umap[[input$gate_select_umap]]
      req(!is.null(res))
      
      choix <- c("Échantillon" = "ECHANTILLON", stats::setNames(res$canaux, res$canaux))
      if (!is.null(p$clusters_flowsom[[input$gate_select_umap]])) {
        choix <- c(choix, "Métacluster (FlowSOM)" = "METACLUSTER")
      }
      selectInput(ns("couleur_umap"), "Colorer par :", choices = choix)
    })
    
    output$plot_umap <- renderPlotly({
      umap_trigger()
      p <- carrot_obj()
      req(input$gate_select_umap)
      res <- p$analyses_umap[[input$gate_select_umap]]
      validate(need(!is.null(res), "Lancez l'UMAP pour afficher la projection."))
      
      df <- as.data.frame(res$embedding)
      noms_axes <- colnames(res$embedding)
      choix_couleur <- input$couleur_umap %||% "ECHANTILLON"
      
      if (identical(choix_couleur, "ECHANTILLON")) {
        df$Couleur <- res$echantillon_origine
        plot_ly(df, x = df[[1]], y = df[[2]], color = ~Couleur,
                type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
          layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
      } else if (identical(choix_couleur, "METACLUSTER")) {
        clust <- p$clusters_flowsom[[input$gate_select_umap]]
        validate(need(!is.null(clust) && length(clust$metaclusters) == nrow(df),
                      "Le clustering doit être calculé avec le même sous-échantillonnage que cette UMAP pour être superposé."))
        df$Metacluster <- as.factor(clust$metaclusters)
        plot_ly(df, x = df[[1]], y = df[[2]], color = ~Metacluster,
                type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
          layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
      } else {
        req(choix_couleur %in% colnames(res$expression))
        df$Marqueur <- res$expression[, choix_couleur]
        plot_ly(df, x = df[[1]], y = df[[2]], color = ~Marqueur, colors = PALETTE_DENSITE,
                type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
          layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
      }
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # ONGLET t-SNE
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_gate_select_tsne <- renderUI({
      gates_trigger()
      p    <- carrot_obj()
      noms <- names(p$gates_personnalisees)
      if (length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Créez d'abord un gate dans l'onglet Gating."))
      }
      selectInput(ns("gate_select_tsne"), "Gate à analyser :", choices = noms)
    })
    
    output$ui_canaux_select_tsne <- renderUI({
      p <- carrot_obj()
      req(input$gate_select_tsne)
      pop <- tryCatch(p$resoudre_population_gate(input$gate_select_tsne), error = function(e) list())
      req(length(pop) > 0)
      cols <- flowCore::colnames(pop[[1]])
      fluo <- cols[!grepl("FSC|SSC|Time", cols, ignore.case = TRUE)]
      selectizeInput(ns("canaux_select_tsne"), "Canaux à utiliser (vide = tous) :",
                     choices = fluo, multiple = TRUE)
    })
    
    observeEvent(input$btn_lancer_tsne, {
      p <- carrot_obj()
      req(input$gate_select_tsne)
      canaux <- if (length(input$canaux_select_tsne) == 0) NULL else input$canaux_select_tsne
      
      withProgress(message = "Calcul du t-SNE en cours (peut prendre quelques minutes)...", value = 0.3, {
        tryCatch({
          p$projection_tSNE(nom_gate = input$gate_select_tsne, canaux = canaux,
                            perplexity = input$tsne_perplexity, theta = input$tsne_theta,
                            max_iter = input$tsne_max_iter, sous_echantillonnage_max = input$tsne_max_cells,
                            cofacteur = input$tsne_cofacteur)
          pipeline(p)
          tsne_trigger(tsne_trigger() + 1L)
          showNotification("t-SNE calculé.", type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    output$ui_couleur_tsne <- renderUI({
      tsne_trigger(); cluster_trigger()
      p <- carrot_obj()
      req(input$gate_select_tsne)
      res <- p$analyses_tsne[[input$gate_select_tsne]]
      req(!is.null(res))
      
      choix <- c("Échantillon" = "ECHANTILLON", stats::setNames(res$canaux, res$canaux))
      if (!is.null(p$clusters_flowsom[[input$gate_select_tsne]])) {
        choix <- c(choix, "Métacluster (FlowSOM)" = "METACLUSTER")
      }
      selectInput(ns("couleur_tsne"), "Colorer par :", choices = choix)
    })
    
    output$plot_tsne <- renderPlotly({
      tsne_trigger()
      p <- carrot_obj()
      req(input$gate_select_tsne)
      res <- p$analyses_tsne[[input$gate_select_tsne]]
      validate(need(!is.null(res), "Lancez le t-SNE pour afficher la projection."))
      
      df <- as.data.frame(res$embedding)
      noms_axes <- colnames(res$embedding)
      choix_couleur <- input$couleur_tsne %||% "ECHANTILLON"
      
      if (identical(choix_couleur, "ECHANTILLON")) {
        df$Couleur <- res$echantillon_origine
        plot_ly(df, x = df[[1]], y = df[[2]], color = ~Couleur,
                type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
          layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
      } else if (identical(choix_couleur, "METACLUSTER")) {
        clust <- p$clusters_flowsom[[input$gate_select_tsne]]
        validate(need(!is.null(clust) && length(clust$metaclusters) == nrow(df),
                      "Le clustering doit être calculé avec le même sous-échantillonnage que ce t-SNE pour être superposé."))
        df$Metacluster <- as.factor(clust$metaclusters)
        plot_ly(df, x = df[[1]], y = df[[2]], color = ~Metacluster,
                type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
          layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
      } else {
        req(choix_couleur %in% colnames(res$expression))
        df$Marqueur <- res$expression[, choix_couleur]
        plot_ly(df, x = df[[1]], y = df[[2]], color = ~Marqueur, colors = PALETTE_DENSITE,
                type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
          layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
      }
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # ONGLET CLUSTERING (FlowSOM)
    # ════════════════════════════════════════════════════════════════════════
    
    output$ui_gate_select_cluster <- renderUI({
      gates_trigger()
      p    <- carrot_obj()
      noms <- names(p$gates_personnalisees)
      if (length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Créez d'abord un gate dans l'onglet Gating."))
      }
      selectInput(ns("gate_select_cluster"), "Gate à analyser :", choices = noms)
    })
    
    # Permet de reprendre exactement les mêmes cellules qu'une UMAP/t-SNE déjà
    # calculée pour ce gate, afin de garantir un alignement parfait entre le
    # clustering et la projection (voir creer_clusters(reutiliser_donnees_de=)).
    output$ui_reutiliser_donnees_cluster <- renderUI({
      umap_trigger(); tsne_trigger()
      p <- carrot_obj()
      req(input$gate_select_cluster)
      
      choix <- c("Extraction indépendante (nouveau sous-échantillonnage)" = "AUCUN")
      if (!is.null(p$analyses_umap[[input$gate_select_cluster]])) choix <- c(choix, "Mêmes cellules que l'UMAP existante" = "umap")
      if (!is.null(p$analyses_tsne[[input$gate_select_cluster]])) choix <- c(choix, "Mêmes cellules que le t-SNE existant" = "tsne")
      
      selectInput(ns("reutiliser_donnees_cluster"), "Cellules à utiliser :", choices = choix)
    })
    
    output$ui_canaux_select_cluster <- renderUI({
      p <- carrot_obj()
      req(input$gate_select_cluster)
      pop <- tryCatch(p$resoudre_population_gate(input$gate_select_cluster), error = function(e) list())
      req(length(pop) > 0)
      cols <- flowCore::colnames(pop[[1]])
      fluo <- cols[!grepl("FSC|SSC|Time", cols, ignore.case = TRUE)]
      selectizeInput(ns("canaux_select_cluster"), "Canaux à utiliser (vide = tous) :",
                     choices = fluo, multiple = TRUE)
    })
    
    observeEvent(input$btn_lancer_cluster, {
      p <- carrot_obj()
      req(input$gate_select_cluster)
      reutiliser <- if (is.null(input$reutiliser_donnees_cluster) || identical(input$reutiliser_donnees_cluster, "AUCUN")) NULL else input$reutiliser_donnees_cluster
      canaux     <- if (length(input$canaux_select_cluster) == 0) NULL else input$canaux_select_cluster
      max_cells  <- if (is.null(input$fsom_max_cells) || is.na(input$fsom_max_cells)) NULL else input$fsom_max_cells
      
      withProgress(message = "Clustering FlowSOM en cours...", value = 0.4, {
        tryCatch({
          p$creer_clusters(nom_gate = input$gate_select_cluster, canaux = canaux,
                           xdim = input$fsom_xdim, ydim = input$fsom_ydim,
                           n_metaclusters = input$fsom_n_metaclusters, sous_echantillonnage_max = max_cells,
                           reutiliser_donnees_de = reutiliser, cofacteur = input$fsom_cofacteur)
          pipeline(p)
          cluster_trigger(cluster_trigger() + 1L)
          showNotification("Clustering calculé.", type = "message")
        }, error = function(e) showNotification(conditionMessage(e), type = "error"))
      })
    })
    
    # Applique les noms personnalisés (renommer_metacluster()) aux identifiants
    # de métacluster, pour un affichage plus lisible partout où ils apparaissent.
    etiqueter_metaclusters <- function(res, valeurs) {
      vapply(as.character(valeurs), function(id) {
        nom_perso <- if (!is.null(res$noms_metaclusters)) res$noms_metaclusters[[id]] else NULL
        if (!is.null(nom_perso) && nchar(trimws(nom_perso)) > 0) paste0(id, " - ", nom_perso) else id
      }, character(1), USE.NAMES = FALSE)
    }
    
    output$plot_cluster_repartition <- renderPlotly({
      cluster_trigger()
      p <- carrot_obj()
      req(input$gate_select_cluster)
      res <- p$clusters_flowsom[[input$gate_select_cluster]]
      validate(need(!is.null(res), "Lancez le clustering pour afficher la répartition des métaclusters."))
      
      df <- as.data.frame(table(Metacluster = etiqueter_metaclusters(res, res$metaclusters)))
      plot_ly(df, x = ~Metacluster, y = ~Freq, type = "bar", marker = list(color = "#0077b6")) %>%
        layout(xaxis = list(title = "Métacluster"), yaxis = list(title = "Nombre de cellules"))
    })
    
    # ── Heatmap d'expression médiane par métacluster ──────────────────────────
    output$plot_heatmap_clusters <- renderPlotly({
      cluster_trigger()
      p <- carrot_obj()
      req(input$gate_select_cluster)
      
      medianes <- tryCatch(p$resumer_expression_clusters(input$gate_select_cluster, niveau = "metacluster"), error = function(e) NULL)
      validate(need(!is.null(medianes), "Lancez le clustering pour afficher la heatmap d'expression."))
      
      # Z-score par colonne (canal) : compare la position relative de chaque
      # métacluster sur un canal donné, plutôt que les valeurs brutes dont
      # l'échelle diffère d'un marqueur à l'autre.
      medianes_z <- scale(medianes)
      medianes_z[is.na(medianes_z)] <- 0
      
      plot_ly(x = colnames(medianes_z), y = rownames(medianes_z), z = medianes_z,
              type = "heatmap", colorscale = "RdBu", reversescale = TRUE,
              colorbar = list(title = "Z-score")) %>%
        layout(xaxis = list(title = "", tickangle = -45), yaxis = list(title = ""))
    })
    
    # ── Annotation manuelle des métaclusters ───────────────────────────────────
    output$ui_annotation_metaclusters <- renderUI({
      cluster_trigger()
      p <- carrot_obj()
      req(input$gate_select_cluster)
      res <- p$clusters_flowsom[[input$gate_select_cluster]]
      validate(need(!is.null(res), "Lancez le clustering pour pouvoir annoter les métaclusters."))
      
      ids <- sort(unique(res$metaclusters))
      tagList(lapply(ids, function(id) {
        nom_existant <- if (!is.null(res$noms_metaclusters)) res$noms_metaclusters[[as.character(id)]] else ""
        textInput(ns(paste0("nom_metacluster_", id)), paste0("Métacluster ", id, " :"),
                  value = nom_existant %||% "", placeholder = "ex: Lymphocytes T CD4+")
      }))
    })
    
    observeEvent(input$btn_enregistrer_annotations, {
      p <- carrot_obj()
      req(input$gate_select_cluster)
      res <- p$clusters_flowsom[[input$gate_select_cluster]]
      req(!is.null(res))
      
      ids <- sort(unique(res$metaclusters))
      tryCatch({
        for (id in ids) {
          valeur <- input[[paste0("nom_metacluster_", id)]]
          if (!is.null(valeur)) p$renommer_metacluster(input$gate_select_cluster, id, valeur)
        }
        pipeline(p)
        cluster_trigger(cluster_trigger() + 1L)
        showNotification("Noms des métaclusters enregistrés.", type = "message")
      }, error = function(e) showNotification(conditionMessage(e), type = "error"))
    })
    
    output$ui_choix_embedding_cluster <- renderUI({
      cluster_trigger(); umap_trigger(); tsne_trigger()
      p <- carrot_obj()
      req(input$gate_select_cluster)
      
      choix <- c()
      if (!is.null(p$analyses_umap[[input$gate_select_cluster]])) choix <- c(choix, "UMAP" = "UMAP")
      if (!is.null(p$analyses_tsne[[input$gate_select_cluster]])) choix <- c(choix, "t-SNE" = "TSNE")
      
      if (length(choix) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("info-circle"),
                   " Calculez d'abord une UMAP ou un t-SNE pour ce même gate afin de visualiser les métaclusters dessus."))
      }
      selectInput(ns("choix_embedding_cluster"), "Projection à colorer :", choices = choix)
    })
    
    output$plot_cluster_embedding <- renderPlotly({
      cluster_trigger()
      p <- carrot_obj()
      req(input$gate_select_cluster, input$choix_embedding_cluster)
      
      type_embed <- if (identical(input$choix_embedding_cluster, "UMAP")) "umap" else "tsne"
      
      res_cluster <- p$clusters_flowsom[[input$gate_select_cluster]]
      validate(need(!is.null(res_cluster), "Lancez le clustering pour afficher cette projection."))
      
      res_embed <- if (identical(type_embed, "umap")) p$analyses_umap[[input$gate_select_cluster]] else p$analyses_tsne[[input$gate_select_cluster]]
      validate(need(!is.null(res_embed), "Cette projection n'est pas disponible pour ce gate."))
      
      # Alignement garanti si le clustering a été explicitement calculé sur les
      # mêmes cellules que cette projection (option "Mêmes cellules que..." dans
      # l'onglet Clustering) ; sinon, on se rabat sur une simple vérification du
      # nombre de lignes (moins fiable : un même nombre ne garantit pas les mêmes
      # cellules dans le même ordre si les deux extractions étaient indépendantes).
      aligne_confirme <- identical(res_cluster$aligne_sur, type_embed)
      if (!aligne_confirme) {
        validate(need(nrow(res_embed$embedding) == length(res_cluster$metaclusters),
                      paste0("Le clustering et la projection ", input$choix_embedding_cluster,
                             " ont été calculés indépendamment et n'ont pas le même nombre de cellules. ",
                             "Dans l'onglet Clustering, choisissez \u00ab Mêmes cellules que l'", input$choix_embedding_cluster,
                             " \u00bb puis relancez le clustering pour garantir un alignement correct.")))
      }
      
      df <- as.data.frame(res_embed$embedding)
      noms_axes <- colnames(res_embed$embedding)
      df$Metacluster <- as.factor(etiqueter_metaclusters(res_cluster, res_cluster$metaclusters))
      
      plot_ly(df, x = df[[1]], y = df[[2]], color = ~Metacluster,
              type = "scatter", mode = "markers", marker = list(size = 3, opacity = 0.6)) %>%
        layout(xaxis = list(title = noms_axes[1]), yaxis = list(title = noms_axes[2]))
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # ONGLET COMPARAISON GROUPES/CONDITIONS
    # ════════════════════════════════════════════════════════════════════════
    
    # Les groupes sont assignés dans module_import_R.R (onglet "3. Groupes
    # d'Échantillons") puis répercutés dans l'objet pipeline partagé
    # (p$groupes_echantillons, échantillon -> groupe) : ce module se contente
    # de les lire — aucune création ni assignation en double ici.
    output$ui_resume_groupes_existants <- renderUI({
      p <- carrot_obj() # Dépend déjà de pipeline_version(), incrémenté par module_import_R.R lors de l'enregistrement des groupes
      
      if (is.null(p$groupes_echantillons) || length(p$groupes_echantillons) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"),
                   " Aucun groupe assigné pour l'instant. Rendez-vous dans l'onglet Import > ",
                   tags$b("3. Groupes d'Échantillons"), " pour créer des groupes et y assigner vos échantillons."))
      }
      
      lignes <- lapply(names(p$groupes_echantillons), function(nom) {
        tags$tr(tags$td(nom), tags$td(tags$b(p$groupes_echantillons[[nom]])))
      })
      tags$table(class = "table table-condensed", style = "font-size:12px;",
                 tags$thead(tags$tr(tags$th("Échantillon"), tags$th("Groupe"))),
                 tags$tbody(lignes))
    })
    
    
    output$ui_gate_select_comparaison <- renderUI({
      gates_trigger()
      p    <- carrot_obj()
      noms <- names(p$gates_personnalisees)
      if (length(noms) == 0) {
        return(div(class = "alert alert-warning", style = "font-size:12px; padding:8px;",
                   icon("exclamation-triangle"), " Créez d'abord un gate dans l'onglet Gating."))
      }
      selectInput(ns("gate_select_comparaison"), "Gate à comparer :", choices = noms)
    })
    
    output$ui_canal_comparaison <- renderUI({
      p <- carrot_obj()
      req(input$gate_select_comparaison)
      pop <- tryCatch(p$resoudre_population_gate(input$gate_select_comparaison), error = function(e) list())
      req(length(pop) > 0)
      cols <- flowCore::colnames(pop[[1]])
      fluo <- cols[!grepl("FSC|SSC|Time", cols, ignore.case = TRUE)]
      fcs_ref <- pop[[1]]
      choix <- stats::setNames(fluo, vapply(fluo, function(c) p$get_label(fcs_ref, c), character(1)))
      selectInput(ns("canal_comparaison"), "Canal :", choices = choix)
    })
    
    resultat_comparaison_rv <- reactiveVal(NULL)
    
    observeEvent(input$btn_lancer_comparaison, {
      p <- carrot_obj()
      req(input$gate_select_comparaison, input$comparaison_variable)
      canal <- if (identical(input$comparaison_variable, "MFI")) input$canal_comparaison else NULL
      
      tryCatch({
        res <- p$comparer_groupes(nom_gate = input$gate_select_comparaison, variable = input$comparaison_variable, canal = canal)
        resultat_comparaison_rv(res)
        showNotification("Comparaison effectuée.", type = "message")
      }, error = function(e) {
        resultat_comparaison_rv(NULL)
        showNotification(conditionMessage(e), type = "error")
      })
    })
    
    output$ui_resultat_comparaison <- renderUI({
      res <- resultat_comparaison_rv()
      req(!is.null(res))
      
      tagList(
        tags$p(tags$b("Test utilisé : "), res$methode),
        tags$p(tags$b("p-value : "), signif(res$p_value, 3),
               if (!is.na(res$p_value) && res$p_value < 0.05) tags$span(" (significatif à 5%)", style = "color:#2e7d32; font-weight:bold;") else NULL),
        tags$p(tags$b("Groupes comparés : "), paste(res$groupes, collapse = " vs "))
      )
    })
    
    output$plot_comparaison_groupes <- renderPlotly({
      res <- resultat_comparaison_rv()
      req(!is.null(res))
      
      df <- res$donnees
      lbl_y <- if (identical(res$variable, "MFI")) paste0("MFI (", res$canal, ")") else "% du parent"
      
      plot_ly(df, x = ~groupe, y = ~valeur, type = "box", boxpoints = "all", jitter = 0.4,
              pointpos = 0, marker = list(color = "#0077b6"), text = ~echantillon, hoverinfo = "text+y") %>%
        layout(xaxis = list(title = "Groupe"), yaxis = list(title = lbl_y))
    })
    
  })
}

# Opérateur null-coalescing utilitaire
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}