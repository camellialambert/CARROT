library(shiny)
library(DT)
library(flowCore)
library(tools)

# ══════════════════════════════════════════════════════════════════════════════
# UI — Module d'import : point d'entrée de toute l'application. 3 onglets :
#   1. Importation : type de cytomètre, fichiers contrôles/échantillons
#   2. Configuration Marqueurs : annotation canal -> marqueur biologique
#   3. Groupes d'Échantillons : assignation échantillon -> groupe/condition
#      (lue ensuite par le module Analyses pour les comparaisons statistiques)
# ══════════════════════════════════════════════════════════════════════════════

import_data_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    tabPanel("1. Importation", br(),
             fluidRow(
               column(4,
                      box(width = 12, status = "primary", solidHeader = TRUE,
                          title = "Console d'Importation",
                          h4("1. Choix du Cytomètre"),
                          selectInput(ns("cyto_type"), "Technologie :", choices = c("Conventionnel", "Spectral")), # Conditionne tout l'aiguillage ultérieur du pipeline (compensation classique vs démixage spectral AutoSpectral)
                          radioButtons(ns("deja_compense"), "Données déjà compensées/démixées ?",
                                       choices = c("Non" = "non", "Oui" = "oui"), inline = TRUE), # Si "Oui" : p$charger_fcs() copiera directement les échantillons dans echantillons_traites, sans repasser par la compensation/l'unmixing (voir pipeline_cytometrie.R)
                          hr(),
                          h4("2. Chargement des fichiers"),
                          fileInput(ns("files_controls"), "Tubes Monomarqués / Unstained (optionnel)", multiple = TRUE, accept = ".fcs"), # Optionnel : absent si "déjà compensé/démixé" ou en spectral (les contrôles se gèrent alors dans l'onglet Démixage)
                          fileInput(ns("files_samples"),  "Échantillons Biologiques", multiple = TRUE, accept = ".fcs"), # Seul champ réellement obligatoire pour pouvoir initialiser le pipeline
                          uiOutput(ns("aide_controles_ui")), # Message d'aide contextuel qui change selon cyto_type/deja_compense (voir server ci-dessous)
                          hr(),
                          actionButton(ns("init_r6"), "Initialiser l'objet",
                                       class = "btn-success",
                                       style = "width:100%; height:30px; font-weight:bold; font-size:16px;") # Crée/réinitialise l'objet pipeline R6 partagé (CARROT$new()) et charge tous les FCS en mémoire
                      )
               ),
               column(8,
                      box(width = 12, status = "warning", solidHeader = TRUE,
                          title = "Annotation & Visualisation",
                          tabsetPanel(
                            tabPanel("Contrôles",    br(), DTOutput(ns("table_controls"))), # Tableau éditable : assignation canal (conventionnel) et type Monomarque/Unstained par tube contrôle
                            tabPanel("Échantillons", br(), DTOutput(ns("table_samples")))   # Tableau éditable : nom de chaque échantillon (seule colonne modifiable), métadonnées en lecture seule
                          )
                      )
               )
             )
    ),
    tabPanel("2. Configuration Marqueurs", br(),
             fluidRow(
               column(12,
                      h3("Matrice de Mapping Marqueurs / Canaux"),
                      p(style = "color:gray;", "Double-cliquez sur une case pour corriger manuellement un marqueur."),
                      div(style = "overflow-x:auto; background:white; padding:10px;",
                          DTOutput(ns("table_matrice_marqueurs"))), # Une ligne par échantillon, une colonne par canal fluorescent : permet de corriger/compléter le marqueur biologique associé à chaque canal, échantillon par échantillon si besoin
                      br(),
                      actionButton(ns("save_marker_config"), "Enregistrer la configuration",
                                   class = "btn-info", style = "width:200px;") # Pousse la table dans p$definir_config_marqueurs() : reconstruit le dictionnaire canal -> marqueur utilisé PARTOUT ensuite pour les libellés d'axes (p$obtenir_label())
               )
             )
    ),
    tabPanel("3. Groupes d'Échantillons", br(),
             fluidRow(
               column(width = 4,
                      box(title = "Gestion des Groupes", width = 12, status = "info", solidHeader = TRUE,
                          textInput(ns("group_name"), "Créer un nouveau groupe", placeholder = "ex: Contrôle..."),
                          actionButton(ns("add_group"), "Ajouter le groupe", class = "btn-info", style = "width:100%;"),
                          hr(),
                          h5("Groupes existants :"),
                          uiOutput(ns("list_groups_ui")) # Simple liste à puces des noms de groupes créés (pas encore l'assignation elle-même, voir colonne de droite)
                      )
               ),
               column(width = 8,
                      box(title = "Assignation des Tubes", width = 12, status = "primary", solidHeader = TRUE,
                          uiOutput(ns("ui_assignation_groupes")), # Un menu déroulant par échantillon (groupe assigné ou "Non assigné")
                          actionButton(ns("save_cohort"), "Enregistrer l'assignation", class = "btn-success") # Répercute l'assignation dans le pipeline partagé (p$groupes_echantillons) : c'est CE bouton qui rend les groupes visibles par le module Analyses
                      ),
                      box(title = "Résumé de la Cohorte", width = 12, status = "warning", solidHeader = TRUE,
                          DTOutput(ns("table_resume_groupes")) # Récapitulatif : un groupe par ligne, avec la liste des tubes qui lui sont assignés et leur effectif
                      )
               )
             )
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════

import_data_server <- function(id, pipeline, pipeline_version,canaux) {
  moduleServer(id, function(input, output, session) {
    
    # État purement local à ce module (tables encore en cours d'édition, avant
    # "Initialiser l'objet"/"Enregistrer..." qui les poussent dans le pipeline
    # R6 partagé). rv$groupes est indexé PAR GROUPE (groupe -> vecteur de noms
    # de tubes), à l'inverse de p$groupes_echantillons qui est indexé PAR
    # ÉCHANTILLON (échantillon -> nom de groupe) — les deux représentations
    # sont converties l'une en l'autre dans save_cohort ci-dessous.
    rv <- reactiveValues(
      r_df_mono            = NULL, # data.frame des tubes contrôles (fichier, canal, type, chemin), NULL si aucun contrôle chargé
      r_df_ech             = NULL, # data.frame des échantillons (tube_name, fichier, chemin, métadonnées FCS)
      r_dictionnaire      = canaux, # Liste de référence des noms de canaux connus (CANAUX_CONNUS, injectée depuis shiny.R), utilisée pour la saisie assistée du canal de chaque tube contrôle
      r_matrice_marqueurs  = NULL, # Non utilisé directement (voir matrice_marqueurs_rv, un reactiveVal séparé, plus bas) — conservé pour compatibilité éventuelle
      groupes              = list() # Groupes créés localement, avant enregistrement définitif dans le pipeline (voir save_cohort)
    )
    
    # ── Texte d'aide contextuel (dépend du type de cytomètre et de la réponse
    # "déjà compensées/démixées ?"), affiché sous les fileInput ────────────
    output$aide_controles_ui <- renderUI({
      if (identical(input$deja_compense, "oui")) { # Cas 1 : données déjà traitées en amont, quel que soit le type de cytomètre
        div(style = "color:gray; font-size:12px; margin-top:-8px;",
            icon("info-circle"),
            " Vos échantillons sont déclarés déjà compensés (conventionnel) ou déjà démixés (spectral) : ",
            "vous pouvez ne charger que les échantillons ci-dessus et passer directement aux étapes de prétraitement.")
      } else if (identical(input$cyto_type, "Spectral")) { # Cas 2 : cytomètre spectral, données brutes à démixer
        div(style = "color:gray; font-size:12px; margin-top:-8px;",
            icon("info-circle"),
            " En spectral, il n'est pas nécessaire d'indiquer de canal pour les tubes contrôles : ",
            "le fichier de contrôle AutoSpectral (fluorophores/marqueurs) sera configuré dans l'onglet Démixage.")
      } else { # Cas 3 (par défaut) : cytomètre conventionnel, données brutes à compenser
        div(style = "color:gray; font-size:12px; margin-top:-8px;",
            icon("info-circle"),
            " Si vos échantillons sont déjà compensés, vous pouvez ne charger que les échantillons ",
            "et passer directement aux étapes suivantes.")
      }
    })
    
    # Dès que des fichiers contrôles sont déposés, construit le tableau
    # d'annotation initial (canal encore vide, à saisir manuellement ensuite
    # dans table_controls) — devine le type Unstained/Monomarque à partir du
    # nom de fichier, à titre de valeur par défaut modifiable.
    observeEvent(input$files_controls, {
      f <- input$files_controls
      noms_proposes <- gsub("\\.fcs$", "", f$name, ignore.case = TRUE) # Nom de fichier sans l'extension, utilisé uniquement pour la détection heuristique du type ci-dessous
      
      rv$r_df_mono <- data.frame(
        fichier = f$name,
        canal   = "", # Vide par défaut : à renseigner par l'utilisateur (conventionnel) ou laissé vide (spectral, voir table_controls)
        type    = ifelse(grepl("unstained", noms_proposes, ignore.case = TRUE), "Unstained", "Monomarque"), # Devine le type depuis le nom de fichier ; reste modifiable ensuite via le menu déroulant de la table
        chemin  = f$datapath, # Chemin temporaire local où Shiny a stocké le fichier uploadé
        stringsAsFactors = FALSE
      )
    })
    
    # Tableau éditable des tubes contrôles : la colonne Type est toujours un
    # menu déroulant HTML natif (via injection de <select> avec un onchange
    # JS ciblant Shiny.setInputValue) ; la colonne Canal n'apparaît qu'en mode
    # Conventionnel (voir aide_controles_ui ci-dessus pour la justification).
    output$table_controls <- renderDT({
      req(rv$r_df_mono)
      df <- rv$r_df_mono
      CANAUX <- rv$r_dictionnaire
      
      sel_type <- sapply(seq_len(nrow(df)), function(i) { # Construit, pour chaque ligne, le HTML brut d'un menu déroulant Monomarque/Unstained (DT n'a pas de widget natif d'édition par liste déroulante, d'où cette injection manuelle)
        opt_m <- if (df$type[i] == "Monomarque") " selected" else ""
        opt_u <- if (df$type[i] == "Unstained")  " selected" else ""
        sprintf('<select class="dt-select" onchange=\'Shiny.setInputValue("%s",{row:%d,val:this.value},{priority:"event"})\'>
            <option value="Monomarque"%s>Monomarque</option>
            <option value="Unstained"%s>Unstained</option></select>',
                session$ns("change_table_type"), i, opt_m, opt_u) # Le onchange pousse {row, val} vers input$change_table_type (voir observeEvent plus bas), avec priority="event" pour garantir le déclenchement même si la valeur "semble" identique
      })
      
      # ── Cytomètre Conventionnel : le canal d'acquisition doit être renseigné
      # manuellement pour chaque tube monomarqué (utilisé ensuite pour la
      # compensation par gates positifs/négatifs, cf. pipeline_cytometrie.R).
      # ── Cytomètre Spectral : le canal n'est pas requis à cette étape. Le
      # rattachement fluorophore/marqueur est configuré séparément dans
      # l'onglet Démixage (fichier de contrôle AutoSpectral), donc on
      # n'affiche que le nom du fichier et le type (Monomarque/Unstained).
      if (identical(input$cyto_type, "Spectral")) {
        datatable(data.frame(Fichier = df$fichier, Type = sel_type, stringsAsFactors = FALSE),
                  escape = FALSE, rownames = FALSE, selection = "none", # escape=FALSE : indispensable pour que le HTML des <select> soit interprété plutôt qu'affiché comme texte brut
                  options = list(dom = "t", paging = FALSE, ordering = FALSE,
                                 preDrawCallback = JS("function(){Shiny.unbindAll(this.api().table().node());}"), # Détache les anciens écouteurs Shiny avant chaque redessin DT (évite les doublons d'événements après un tri/filtre)
                                 drawCallback    = JS("function(){Shiny.bindAll(this.api().table().node());}")))  # Ré-attache les écouteurs Shiny sur les nouveaux éléments HTML injectés après le redessin
      } else {
        sel_canal <- sapply(seq_len(nrow(df)), function(i) { # Construit, pour chaque ligne, un champ texte avec autocomplétion (HTML <datalist>) proposant les canaux connus, tout en acceptant une saisie libre
          v <- df$canal[i]
          datalist_id <- session$ns(paste0("dl_canal_", i))
          input_id    <- session$ns(paste0("txt_canal_", i))
          datalist_opts <- paste0(sapply(CANAUX[nchar(CANAUX) > 0], function(o) sprintf('<option value="%s">', o)), collapse = "") # Ignore les entrées vides du dictionnaire (ex: le "" utilisé comme choix placeholder ailleurs)
          
          # Modification : on remplace oninput par onchange
          paste0('<div style="display:flex;gap:4px;">',
                 sprintf('<input type="text" class="dt-datalist-input" list="%s" id="%s" value="%s" placeholder="Choisir le canal" style="flex:1;" onchange=\'Shiny.setInputValue("%s",{row:%d,val:this.value},{priority:"event"})\'>',
                         datalist_id, input_id, v, session$ns("change_table_canal"), i), # onchange (et non oninput) : ne déclenche la mise à jour Shiny qu'une fois la saisie/sélection terminée, pas à chaque frappe de touche
                 sprintf('<datalist id="%s">%s</datalist>', datalist_id, datalist_opts), '</div>')
        })
        
        datatable(data.frame(Fichier = df$fichier, Canal = sel_canal, Type = sel_type, stringsAsFactors = FALSE),
                  escape = FALSE, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE, ordering = FALSE,
                                 preDrawCallback = JS("function(){Shiny.unbindAll(this.api().table().node());}"),
                                 drawCallback    = JS("function(){Shiny.bindAll(this.api().table().node());}")))
      }
    })
    
    # Répercute dans rv$r_df_mono la modification d'un canal saisie dans le
    # champ texte injecté ci-dessus (voir onchange -> Shiny.setInputValue).
    observeEvent(input$change_table_canal, {
      req(rv$r_df_mono)
      info <- input$change_table_canal # {row, val} envoyé par le JS de la table
      rv$r_df_mono$canal[info$row] <- trimws(info$val)
    })
    
    # Répercute dans rv$r_df_mono la modification du type Monomarque/Unstained
    # sélectionnée dans le menu déroulant injecté ci-dessus.
    observeEvent(input$change_table_type, {
      req(rv$r_df_mono)
      info <- input$change_table_type
      rv$r_df_mono$type[info$row] <- info$val
    })
    
    # Bouton "Initialiser l'objet" : crée l'objet pipeline R6 (ou le
    # réinitialise entièrement si déjà existant) et charge tous les FCS en
    # mémoire. C'est le SEUL endroit de toute l'application où p$initialize()
    # et p$charger_fcs() sont appelés.
    observeEvent(input$init_r6, {
      # Seuls les échantillons biologiques sont obligatoires : les tubes contrôles
      # (monomarqués / unstained) sont optionnels, notamment quand les échantillons
      # sont déjà compensés (conventionnel) ou déjà unmixés (spectral).
      req(rv$r_df_ech)
      
      p <- pipeline()
      p$initialize(
        df_monomarques  = rv$r_df_mono,       # NULL si aucun tube contrôle n'a été chargé
        df_echantillons = rv$r_df_ech,
        mode            = input$cyto_type,
        deja_traite     = identical(input$deja_compense, "oui")
      )
      
      withProgress(message = "Chargement...", value = 0.2, {
        tryCatch({
          p$charger_fcs() # Lit effectivement tous les fichiers FCS depuis le disque (contrôles + échantillons) et peuple self$echantillons/self$tubes_monomarques
          pipeline(p) # Republie l'objet R6 modifié dans le reactiveVal partagé, pour que tous les autres modules le voient
          pipeline_version(pipeline_version() + 1L) # Incrémente le compteur global : force la réinvalidation de tous les éléments réactifs des autres modules qui en dépendent
          
          if (is.null(rv$r_df_mono) && isTRUE(p$deja_traite) && p$mode == "Conventionnel") { # Cas particulier : conventionnel, pas de contrôles fournis, données annoncées déjà compensées -> tente une récupération automatique de la matrice de spillover depuis les métadonnées FCS (voir p$importer_spillover_fcs(), appelé en interne par charger_fcs())
            if (!is.null(p$S_matrix)) {
              showNotification("Initialisation effectuée : matrice de compensation détectée dans les métadonnées FCS (onglet Compensation).", type = "message", duration = 8)
            } else {
              showNotification("Initialisation effectuée. Aucune matrice de compensation n'a été trouvée dans les métadonnées FCS : vous pourrez la définir manuellement dans l'onglet Compensation.", type = "warning", duration = 8)
            }
          } else {
            showNotification("Initialisation effectuée.", type = "message")
          }
        }, error = function(e) {
          showNotification(paste("Erreur :", conditionMessage(e)), type = "error")
        })
      })
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # UPLOAD ÉCHANTILLONS
    # Lecture des métadonnées FCS (TUBE NAME, $TOT, $CYT, $DATE...)
    # Fallback sur le nom du fichier si TUBE NAME absent ou vide
    # ════════════════════════════════════════════════════════════════════════
    observeEvent(input$files_samples, {
      f <- input$files_samples
      withProgress(message = "Lecture des métadonnées FCS...", value = 0, {
        rows <- lapply(seq_len(nrow(f)), function(i) {
          incProgress(1 / nrow(f), detail = paste("Fichier", i, "/", nrow(f))) # Barre de progression : utile car la lecture d'en-tête de nombreux gros fichiers FCS peut prendre plusieurs secondes
          ch        <- f$datapath[i]
          nb_events <- 0; tube_name <- "Inconnu"; exp_name  <- "Inconnu" # Valeurs de repli si la lecture des métadonnées échoue (fichier corrompu, format non standard...)
          cytometre <- "Inconnu"; date_acq  <- "Inconnu"
          
          tryCatch({
            hdr       <- flowCore::read.FCS(ch, dataset = 1, which.lines = 1,
                                            transformation = FALSE, truncate_max_range = FALSE) # which.lines=1 : ne lit qu'UN seul événement (juste pour accéder aux métadonnées d'en-tête), évite de charger tout le fichier en mémoire à ce stade
            kw        <- flowCore::keyword(hdr)
            nb_events <- if (!is.null(kw[["$TOT"]])) as.numeric(kw[["$TOT"]]) else 0 # $TOT = nombre total d'événements réellement présents dans le fichier (mot-clé FCS standard)
            tube_name <- if (!is.null(kw[["TUBE NAME"]]) && nchar(kw[["TUBE NAME"]]) > 0)
              kw[["TUBE NAME"]] # Priorité au nom de tube tel que saisi sur le cytomètre à l'acquisition
            else tools::file_path_sans_ext(f$name[i]) # Repli : nom de fichier sans extension, si TUBE NAME est absent ou vide
            exp_name  <- if (!is.null(kw[["EXPERIMENT NAME"]])) kw[["EXPERIMENT NAME"]] else "Inconnu"
            cytometre <- if (!is.null(kw[["$CYT"]])) kw[["$CYT"]] else "Inconnu" # $CYT = modèle de cytomètre déclaré dans le FCS (mot-clé standard)
            date_acq  <- if (!is.null(kw[["$DATE"]])) kw[["$DATE"]] else "Inconnu"
          }, error = function(e) NULL) # Échec silencieux : garde les valeurs de repli "Inconnu" plutôt que de bloquer tout l'import pour un seul fichier problématique
          
          sz    <- file.info(ch)$size
          poids <- if (is.na(sz)) "?"
          else if (sz >= 1e9) paste0(round(sz / 1e9, 2), " Go") # Affiche en Go au-delà de 1 milliard d'octets, sinon en Mo (lisibilité pour l'utilisateur)
          else paste0(round(sz / 1e6, 1), " Mo")
          
          data.frame(tube_name = tube_name, fichier = f$name[i], chemin = ch,
                     nb_events = format(nb_events, big.mark = "\u00a0"), # Espace insécable comme séparateur de milliers (évite un retour à la ligne indésirable dans la table)
                     exp_name  = exp_name, cytometre = cytometre,
                     poids = poids, date = date_acq,
                     stringsAsFactors = FALSE)
        })
        rv$r_df_ech <- do.call(rbind, rows) # Empile toutes les lignes (une par fichier) en un seul data.frame
      })
    })
    
    # ── Table des échantillons ─────────────────────────────────────────────
    # Seule la colonne "Nom" (tube_name, colonne d'index 0) est éditable ;
    # toutes les métadonnées lues automatiquement du FCS sont verrouillées
    # (disable = columns 1 à 6) pour éviter toute incohérence avec le fichier
    # réel sur le disque.
    output$table_samples <- renderDT({
      req(rv$r_df_ech)
      df   <- rv$r_df_ech
      cols <- c("tube_name", "fichier", "nb_events", "exp_name", "cytometre", "poids", "date")
      datatable(
        df[, cols],
        colnames = c("Nom", "Fichier", "Évènements",
                     "Expérience", "Cytomètre", "Volume", "Date"),
        editable = list(target = "cell", disable = list(columns = c(1, 2, 3, 4, 5, 6))),
        rownames = FALSE, selection = "none",
        options  = list(dom = "t", paging = FALSE, scrollX = TRUE)
      )
    })
    
    # Gestionnaire du renommage d'échantillon (édition de cellule dans la
    # colonne "Nom" ci-dessus). Sans ce gestionnaire, l'édition resterait
    # purement visuelle côté client (DT) et ne serait jamais répercutée dans
    # rv$r_df_ech : le nouveau nom se perdrait silencieusement au moindre
    # nouveau rendu de la table.
    observeEvent(input$table_samples_cell_edit, {
      info <- input$table_samples_cell_edit
      req(rv$r_df_ech, info$col == 0) # Sécurité : n'agit que sur la colonne tube_name, la seule éditable
      
      nouveau_nom <- trimws(info$value)
      ancien_nom  <- rv$r_df_ech$tube_name[info$row]
      
      if (nchar(nouveau_nom) == 0) {
        showNotification("Le nom de l'échantillon ne peut pas être vide.", type = "error")
        return(invisible(NULL))
      }
      if (nouveau_nom != ancien_nom && nouveau_nom %in% rv$r_df_ech$tube_name) { # Empêche deux échantillons de porter le même nom (les noms servent de clés partout dans le pipeline : une collision créerait un écrasement silencieux)
        showNotification(paste0("Le nom '", nouveau_nom, "' est déjà utilisé par un autre échantillon."), type = "error")
        return(invisible(NULL))
      }
      
      df <- rv$r_df_ech
      df$tube_name[info$row] <- nouveau_nom
      rv$r_df_ech <- df
      
      # Répercute le renommage dans l'objet pipeline partagé, s'il a déjà été
      # initialisé ("Initialiser l'objet" déjà cliqué) : sans cela, tous les
      # autres modules (compensation, QC, prétraitement, analyses) continuent
      # d'afficher l'ancien nom, puisqu'ils lisent le pipeline, pas cette table.
      p <- pipeline()
      if (inherits(p, "R6") && !is.null(p$echantillons) && length(p$echantillons) > 0) {
        p$renommer_echantillon(ancien_nom, nouveau_nom) # Renomme la clé de l'échantillon dans TOUTES les structures du pipeline qui la référencent (echantillons, post_*, gates_history, groupes_echantillons...) — voir pipeline_cytometrie.R pour la liste exhaustive
        pipeline(p)
        pipeline_version(pipeline_version() + 1L)
      }
      
      showNotification(paste0("Échantillon renommé : '", ancien_nom, "' → '", nouveau_nom, "'."), type = "message")
    })
    
    # ── Matrice marqueurs ──────────────────────────────────────────────────
    # IMPORTANT : les colonnes de matrice_marqueurs_rv restent indexées par le nom de
    # CANAL brut (ex: "PE-A"), pas par le libellé affiché, afin de pouvoir reconstruire
    # un dictionnaire canal -> marqueur exploitable par obtenir_label() une fois enregistré.
    matrice_marqueurs_rv     <- reactiveVal(NULL)
    canal_labels_affiches_rv <- reactiveVal(NULL) # libellés "jolis" (ex: "PE-A | CD4") utilisés uniquement pour l'affichage DT
    
    # Dès que la liste des échantillons est disponible, pré-remplit la matrice
    # marqueurs à partir des descriptions déjà présentes dans les métadonnées
    # FCS de chaque fichier (le cas échéant), pour éviter à l'utilisateur de
    # tout ressaisir manuellement s'il a déjà correctement annoté ses fichiers
    # à l'acquisition.
    observeEvent(rv$r_df_ech, {
      req(rv$r_df_ech)
      tryCatch({
        fcs_ref <- flowCore::read.FCS(rv$r_df_ech$chemin[1],
                                      transformation = FALSE, truncate_max_range = FALSE) # Sert uniquement de référence pour connaître la liste des canaux fluorescents (cx_fluo) : suppose que tous les échantillons partagent le même panel
        cx_fluo <- canaux_fluo_fcs(fcs_ref) # Utilitaire de utils.R : renvoie tous les noms de colonnes (canaux) du flowFrame
        mat     <- matrix("", nrow = nrow(rv$r_df_ech), ncol = length(cx_fluo)) # Grille vide (échantillons x canaux), remplie ci-dessous
        colnames(mat) <- cx_fluo
        rownames(mat) <- rv$r_df_ech$tube_name
        
        for (i in seq_len(nrow(rv$r_df_ech))) { # Relit CHAQUE fichier individuellement (pas seulement la référence) : le marqueur annoté sur un canal peut varier d'un échantillon à l'autre dans certains protocoles
          fcs_obj <- flowCore::read.FCS(rv$r_df_ech$chemin[i],
                                        transformation = FALSE, truncate_max_range = FALSE)
          lbs <- get_labels_from_fcs(fcs_obj) # utils.R : "canal | marqueur" si une description existe dans les métadonnées, sinon juste le canal
          for (cx in cx_fluo) {
            lbl              <- lbs[[cx]]
            marqueur_extrait <- if (!is.null(lbl) && grepl(" \\| ", lbl))
              strsplit(lbl, " \\| ")[[1]][2] else "" # N'extrait que la partie "marqueur" du libellé combiné (après le " | "), pas le canal technique
            mat[i, cx] <- if (nchar(trimws(marqueur_extrait)) > 0) marqueur_extrait else cx # Si aucun marqueur n'a été trouvé dans les métadonnées, pré-remplit avec le nom du canal lui-même (valeur de repli neutre, à corriger manuellement par l'utilisateur si besoin)
          }
        }
        lbs_ref <- get_labels_from_fcs(fcs_ref)
        # Les colonnes de la matrice restent les canaux bruts (cx_fluo) ; on ne conserve
        # les libellés "jolis" que pour l'affichage (colnames = ... dans renderDT).
        canal_labels_affiches_rv(sapply(cx_fluo, function(cx) { l <- lbs_ref[[cx]]; if (!is.null(l)) l else cx }))
        matrice_marqueurs_rv(as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)) # check.names=FALSE : préserve les noms de canaux exacts (ex: "FSC-A") sans que R ne les modifie (il remplacerait sinon le tiret par un point)
      }, error = function(e) NULL) # Échec silencieux (ex: fichier illisible) : la table reste vide, l'utilisateur pourra quand même la remplir manuellement
    })
    
    output$table_matrice_marqueurs <- renderDT({
      df <- matrice_marqueurs_rv()
      req(df)
      labels_affiches <- canal_labels_affiches_rv()
      # IMPORTANT : DT interprète un vecteur de colnames NOMMÉ comme un mapping
      # "Nouveau nom" = "Nom existant dans les données" (voir doc DT::datatable).
      # Comme labels_affiches est un vecteur nommé (noms = canaux), il faut le
      # "dénommer" (unname) avant de le passer, sinon DT tente de faire correspondre
      # ces noms aux colonnes réelles du data.frame et échoue avec l'erreur
      # "Some column names in the 'escape' argument not found in data".
      entetes <- if (!is.null(labels_affiches)) unname(c("Tube", labels_affiches)) else colnames(df) # "Tube" en premier : correspond à rownames(df) affiché comme première colonne visuelle (rownames=TRUE plus bas)
      datatable(df, editable = list(target = "cell"),
                rownames = TRUE, selection = "none",
                colnames = entetes,
                options  = list(scrollX = TRUE, pageLength = 20, dom = "ft")) # dom="ft" : affiche uniquement le champ de recherche (f) et le tableau (t), pas de pagination ni d'info superflue
    })
    
    # Répercute dans matrice_marqueurs_rv toute correction manuelle d'une
    # cellule (marqueur mal détecté ou absent des métadonnées d'origine).
    observeEvent(input$table_matrice_marqueurs_cell_edit, {
      info <- input$table_matrice_marqueurs_cell_edit
      df   <- matrice_marqueurs_rv()
      df[info$row, info$col] <- info$value
      matrice_marqueurs_rv(df)
    })
    
    # Bouton "Enregistrer la configuration" : pousse la matrice marqueurs
    # définitive dans le pipeline partagé.
    observeEvent(input$save_marker_config, {
      data_to_save <- matrice_marqueurs_rv()
      req(data_to_save)
      p <- pipeline()
      if (inherits(p, "R6")) {
        # p$definir_config_marqueurs() enregistre la table ET reconstruit le dictionnaire
        # canal -> marqueur utilisé ensuite par obtenir_label() pour tous les libellés d'axes
        # (biplots de compensation, gates, QC, prétraitement...).
        p$definir_config_marqueurs(data_to_save)
        pipeline(p)
        pipeline_version(pipeline_version() + 1L)
        showNotification("Configuration enregistrée : les libellés des graphiques utiliseront désormais ces annotations.", type = "message")
      } else {
        showNotification("Erreur : Pipeline non valide.", type = "error")
      }
    })
    
    # ── Groupes ───────────────────────────────────────────────────────────
    # Création d'un nouveau groupe (nom libre, ex: "Contrôle", "Traité", "J7"...).
    observeEvent(input$add_group, {
      req(input$group_name, nchar(trimws(input$group_name)) > 0)
      nom <- trimws(input$group_name)
      if (!(nom %in% names(rv$groupes))) { # Empêche la création d'un doublon de nom de groupe
        temp <- rv$groupes; temp[[nom]] <- character(0); rv$groupes <- temp # Nouveau groupe initialement vide (aucun tube assigné) ; passe par une variable temporaire pour ne déclencher qu'UNE seule invalidation réactive de rv$groupes, pas deux
        updateTextInput(session, "group_name", value = "") # Vide le champ de saisie après création, prêt pour le groupe suivant
      } else showNotification("Ce groupe existe déjà.", type = "warning")
    })
    
    output$list_groups_ui <- renderUI({
      if (length(rv$groupes) == 0) return(p("Aucun groupe créé.", style = "color:gray;"))
      tags$ul(lapply(names(rv$groupes), function(g) tags$li(tags$b(g))))
    })
    
    # Un menu déroulant par échantillon (groupe assigné ou "Non assigné").
    output$ui_assignation_groupes <- renderUI({
      req(rv$r_df_ech, length(rv$groupes) > 0)
      ns_func <- session$ns
      fluidRow(lapply(rv$r_df_ech$tube_name, function(tube) {
        input_id <- paste0("assign_", make.names(tube)) # make.names() : garantit un ID Shiny valide même si le nom de tube contient des espaces/caractères spéciaux
        
        # Préserve le choix déjà fait dans le menu déroulant (même si pas encore
        # enregistré via "save_cohort") s'il reste valide parmi les groupes
        # actuels. Sans cela, la création d'un nouveau groupe reconstruit tous
        # les sélecteurs et recalcule "selected" uniquement à partir de
        # rv$groupes (qui ne connaît pas encore ce choix en attente) : chaque
        # assignation non sauvegardée se retrouvait donc silencieusement
        # effacée dès qu'on ajoutait un groupe supplémentaire.
        choix_en_attente <- input[[input_id]]
        
        valeur_selectionnee <- if (!is.null(choix_en_attente) && choix_en_attente %in% c("", names(rv$groupes))) {
          choix_en_attente # Le choix en attente reste valide (le groupe existe toujours) : on le garde tel quel
        } else {
          grp <- names(rv$groupes)[sapply(rv$groupes, function(g) tube %in% g)] # Repli : cherche si ce tube est déjà formellement assigné à un groupe (assignation déjà enregistrée lors d'un save_cohort précédent)
          if (length(grp)) grp[1] else "" # Un tube ne devrait appartenir qu'à un seul groupe ; grp[1] par sécurité si jamais plusieurs groupes le référençaient
        }
        
        column(6,
               selectInput(ns_func(input_id), label = tube,
                           choices  = c("(Non assigné)" = "", names(rv$groupes)),
                           selected = valeur_selectionnee)
        )
      }))
    })
    
    output$table_resume_groupes <- renderDT({
      req(length(rv$groupes) > 0)
      datatable(
        do.call(rbind, lapply(names(rv$groupes), function(g) {
          data.frame(Groupe = g, Fichiers = paste(rv$groupes[[g]], collapse = ", "),
                     N = length(rv$groupes[[g]]), stringsAsFactors = FALSE)
        })),
        options = list(dom = "t"), rownames = FALSE
      )
    })
    
    # Bouton "Enregistrer l'assignation" : lit la valeur actuelle de CHAQUE
    # menu déroulant (input$assign_*), reconstruit rv$groupes (groupe -> tubes)
    # à partir de ça, ET répercute la même information dans le pipeline partagé
    # au format inverse (échantillon -> groupe), qui est le format lu par le
    # module Analyses pour ses comparaisons statistiques.
    observeEvent(input$save_cohort, {
      req(rv$r_df_ech)
      new_groupes <- lapply(rv$groupes, function(x) character(0)) # Repart d'une liste de groupes tous vides (mêmes noms que rv$groupes), pour re-remplir depuis zéro à partir des sélections actuelles des menus déroulants
      for (tube in rv$r_df_ech$tube_name) {
        val <- input[[paste0("assign_", make.names(tube))]]
        if (!is.null(val) && val != "") new_groupes[[val]] <- c(new_groupes[[val]], tube) # Ajoute ce tube à la liste des tubes du groupe choisi (val != "" exclut les tubes "Non assigné")
      }
      rv$groupes <- new_groupes
      
      # Répercute l'assignation dans l'objet pipeline partagé (p$groupes_echantillons,
      # au format échantillon -> groupe) : c'est ce que lit directement le module
      # Analyses (onglet Comparaison Groupes) pour ses comparaisons statistiques,
      # évitant ainsi de devoir recréer les groupes une seconde fois là-bas.
      p <- pipeline()
      if (inherits(p, "R6")) {
        for (tube in rv$r_df_ech$tube_name) {
          val <- input[[paste0("assign_", make.names(tube))]]
          p$definir_groupe_echantillon(tube, if (is.null(val)) "" else val) # Une valeur vide retire l'échantillon de tout groupe (voir p$definir_groupe_echantillon() côté pipeline)
        }
        pipeline(p)
        pipeline_version(pipeline_version() + 1L)
      }
      
      showNotification("Cohorte enregistrée !", type = "message")
    })
    
  })
}