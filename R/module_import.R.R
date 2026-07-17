library(shiny)
library(DT)
library(flowCore)
library(tools)

import_data_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    tabPanel("1. Importation", br(),
             fluidRow(
               column(4,
                      box(width = 12, status = "primary", solidHeader = TRUE,
                          title = "Console d'Importation",
                          h4("1. Chargement des fichiers"),
                          fileInput(ns("files_controls"), "Tubes Monomarqués / Unstained (optionnel)", multiple = TRUE, accept = ".fcs"),
                          fileInput(ns("files_samples"),  "Échantillons Biologiques", multiple = TRUE, accept = ".fcs"),
                          div(style = "color:gray; font-size:12px; margin-top:-8px;",
                              icon("info-circle"),
                              " Si vos échantillons sont déjà compensés (conventionnel) ou déjà unmixés (spectral), ",
                              "vous pouvez ne charger que les échantillons et passer directement aux étapes suivantes."),
                          hr(),
                          h4("2. Choix du Cytomètre"),
                          selectInput(ns("cyto_type"), "Technologie :", choices = c("Conventionnel", "Spectral")),
                          radioButtons(ns("deja_compense"), "Données déjà compensées ?",
                                       choices = c("Non" = "non", "Oui" = "oui"), inline = TRUE),
                          hr(),
                          actionButton(ns("init_r6"), "Initialiser l'objet",
                                       class = "btn-success",
                                       style = "width:100%; height:30px; font-weight:bold; font-size:16px;")
                      )
               ),
               column(8,
                      box(width = 12, status = "warning", solidHeader = TRUE,
                          title = "Annotation & Visualisation",
                          tabsetPanel(
                            tabPanel("Contrôles",    br(), DTOutput(ns("table_controls"))),
                            tabPanel("Échantillons", br(), DTOutput(ns("table_samples")))
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
                          DTOutput(ns("table_matrice_marqueurs"))),
                      br(),
                      actionButton(ns("save_marker_config"), "Enregistrer la configuration",
                                   class = "btn-info", style = "width:200px;")
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
                          uiOutput(ns("list_groups_ui"))
                      )
               ),
               column(width = 8,
                      box(title = "Assignation des Tubes", width = 12, status = "primary", solidHeader = TRUE,
                          uiOutput(ns("ui_assignation_groupes")),
                          actionButton(ns("save_cohort"), "Enregistrer l'assignation", class = "btn-success")
                      ),
                      box(title = "Résumé de la Cohorte", width = 12, status = "warning", solidHeader = TRUE,
                          DTOutput(ns("table_resume_groupes"))
                      )
               )
             )
    )
  )
}

import_data_server <- function(id, pipeline, pipeline_version,canaux) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(
      r_df_mono            = NULL,
      r_df_ech             = NULL,
      r_dictionnaire      = canaux,
      r_matrice_marqueurs  = NULL,
      groupes              = list()
    )
    
    observeEvent(input$files_controls, {
      f <- input$files_controls
      noms_proposes <- gsub("\\.fcs$", "", f$name, ignore.case = TRUE)
      
      rv$r_df_mono <- data.frame(
        fichier = f$name,
        canal   = "",
        type    = ifelse(grepl("unstained", noms_proposes, ignore.case = TRUE), "Unstained", "Monomarque"),
        chemin  = f$datapath,
        stringsAsFactors = FALSE
      )
    })
    
    output$table_controls <- renderDT({
      req(rv$r_df_mono)
      df <- rv$r_df_mono
      CANAUX <- rv$r_dictionnaire
      
      sel_canal <- sapply(seq_len(nrow(df)), function(i) {
        v <- df$canal[i]
        datalist_id <- session$ns(paste0("dl_canal_", i))
        input_id    <- session$ns(paste0("txt_canal_", i))
        datalist_opts <- paste0(sapply(CANAUX[nchar(CANAUX) > 0], function(o) sprintf('<option value="%s">', o)), collapse = "")
        
        # Modification : on remplace oninput par onchange
        paste0('<div style="display:flex;gap:4px;">',
               sprintf('<input type="text" class="dt-datalist-input" list="%s" id="%s" value="%s" placeholder="Choisir le canal" style="flex:1;" onchange=\'Shiny.setInputValue("%s",{row:%d,val:this.value},{priority:"event"})\'>',
                       datalist_id, input_id, v, session$ns("change_table_canal"), i),
               sprintf('<datalist id="%s">%s</datalist>', datalist_id, datalist_opts), '</div>')
      })
      
      sel_type <- sapply(seq_len(nrow(df)), function(i) {
        opt_m <- if (df$type[i] == "Monomarque") " selected" else ""
        opt_u <- if (df$type[i] == "Unstained")  " selected" else ""
        sprintf('<select class="dt-select" onchange=\'Shiny.setInputValue("%s",{row:%d,val:this.value},{priority:"event"})\'>
            <option value="Monomarque"%s>Monomarque</option>
            <option value="Unstained"%s>Unstained</option></select>',
                session$ns("change_table_type"), i, opt_m, opt_u)
      })
      
      datatable(data.frame(Fichier = df$fichier, Canal = sel_canal, Type = sel_type, stringsAsFactors = FALSE),
                escape = FALSE, rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, ordering = FALSE,
                               preDrawCallback = JS("function(){Shiny.unbindAll(this.api().table().node());}"),
                               drawCallback    = JS("function(){Shiny.bindAll(this.api().table().node());}")))
    })
    
    observeEvent(input$change_table_canal, {
      req(rv$r_df_mono)
      info <- input$change_table_canal
      rv$r_df_mono$canal[info$row] <- trimws(info$val)
    })
    
    observeEvent(input$change_table_type, {
      req(rv$r_df_mono)
      info <- input$change_table_type
      rv$r_df_mono$type[info$row] <- info$val
    })
    
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
          p$charger_fcs()
          pipeline(p)
          pipeline_version(pipeline_version() + 1L)
          
          if (is.null(rv$r_df_mono) && isTRUE(p$deja_traite) && p$mode == "Conventionnel") {
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
          incProgress(1 / nrow(f), detail = paste("Fichier", i, "/", nrow(f)))
          ch        <- f$datapath[i]
          nb_events <- 0; tube_name <- "Inconnu"; exp_name  <- "Inconnu"
          cytometre <- "Inconnu"; date_acq  <- "Inconnu"
          
          tryCatch({
            hdr       <- flowCore::read.FCS(ch, dataset = 1, which.lines = 1,
                                            transformation = FALSE, truncate_max_range = FALSE)
            kw        <- flowCore::keyword(hdr)
            nb_events <- if (!is.null(kw[["$TOT"]])) as.numeric(kw[["$TOT"]]) else 0
            tube_name <- if (!is.null(kw[["TUBE NAME"]]) && nchar(kw[["TUBE NAME"]]) > 0)
              kw[["TUBE NAME"]]
            else tools::file_path_sans_ext(f$name[i])
            exp_name  <- if (!is.null(kw[["EXPERIMENT NAME"]])) kw[["EXPERIMENT NAME"]] else "Inconnu"
            cytometre <- if (!is.null(kw[["$CYT"]])) kw[["$CYT"]] else "Inconnu"
            date_acq  <- if (!is.null(kw[["$DATE"]])) kw[["$DATE"]] else "Inconnu"
          }, error = function(e) NULL)
          
          sz    <- file.info(ch)$size
          poids <- if (is.na(sz)) "?"
          else if (sz >= 1e9) paste0(round(sz / 1e9, 2), " Go")
          else paste0(round(sz / 1e6, 1), " Mo")
          
          data.frame(tube_name = tube_name, fichier = f$name[i], chemin = ch,
                     nb_events = format(nb_events, big.mark = "\u00a0"),
                     exp_name  = exp_name, cytometre = cytometre,
                     poids = poids, date = date_acq,
                     stringsAsFactors = FALSE)
        })
        rv$r_df_ech <- do.call(rbind, rows)
      })
    })
    
    # ── Table des échantillons ─────────────────────────────────────────────
    output$table_samples <- renderDT({
      req(rv$r_df_ech)
      df   <- rv$r_df_ech
      cols <- c("tube_name", "fichier", "nb_events", "exp_name", "cytometre", "poids", "date")
      datatable(
        df[, cols],
        colnames = c("Nom (modifiable)", "Fichier", "Évènements",
                     "Expérience", "Cytomètre", "Volume", "Date"),
        editable = list(target = "cell", disable = list(columns = c(1, 2, 3, 4, 5, 6))),
        rownames = FALSE, selection = "none",
        options  = list(dom = "t", paging = FALSE, scrollX = TRUE)
      )
    })
    
    # ── Matrice marqueurs ──────────────────────────────────────────────────
    # IMPORTANT : les colonnes de matrice_marqueurs_rv restent indexées par le nom de
    # CANAL brut (ex: "PE-A"), pas par le libellé affiché, afin de pouvoir reconstruire
    # un dictionnaire canal -> marqueur exploitable par get_label() une fois enregistré.
    matrice_marqueurs_rv     <- reactiveVal(NULL)
    canal_labels_affiches_rv <- reactiveVal(NULL) # libellés "jolis" (ex: "PE-A | CD4") utilisés uniquement pour l'affichage DT
    
    observeEvent(rv$r_df_ech, {
      req(rv$r_df_ech)
      tryCatch({
        fcs_ref <- flowCore::read.FCS(rv$r_df_ech$chemin[1],
                                      transformation = FALSE, truncate_max_range = FALSE)
        cx_fluo <- canaux_fluo_fcs(fcs_ref)
        mat     <- matrix("", nrow = nrow(rv$r_df_ech), ncol = length(cx_fluo))
        colnames(mat) <- cx_fluo
        rownames(mat) <- rv$r_df_ech$tube_name
        
        for (i in seq_len(nrow(rv$r_df_ech))) {
          fcs_obj <- flowCore::read.FCS(rv$r_df_ech$chemin[i],
                                        transformation = FALSE, truncate_max_range = FALSE)
          lbs <- get_labels_from_fcs(fcs_obj)
          for (cx in cx_fluo) {
            lbl              <- lbs[[cx]]
            marqueur_extrait <- if (!is.null(lbl) && grepl(" \\| ", lbl))
              strsplit(lbl, " \\| ")[[1]][2] else ""
            mat[i, cx] <- if (nchar(trimws(marqueur_extrait)) > 0) marqueur_extrait else cx
          }
        }
        lbs_ref <- get_labels_from_fcs(fcs_ref)
        # Les colonnes de la matrice restent les canaux bruts (cx_fluo) ; on ne conserve
        # les libellés "jolis" que pour l'affichage (colnames = ... dans renderDT).
        canal_labels_affiches_rv(sapply(cx_fluo, function(cx) { l <- lbs_ref[[cx]]; if (!is.null(l)) l else cx }))
        matrice_marqueurs_rv(as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE))
      }, error = function(e) NULL)
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
      entetes <- if (!is.null(labels_affiches)) unname(c("Tube", labels_affiches)) else colnames(df)
      datatable(df, editable = list(target = "cell"),
                rownames = TRUE, selection = "none",
                colnames = entetes,
                options  = list(scrollX = TRUE, pageLength = 20, dom = "ft"))
    })
    
    observeEvent(input$table_matrice_marqueurs_cell_edit, {
      info <- input$table_matrice_marqueurs_cell_edit
      df   <- matrice_marqueurs_rv()
      df[info$row, info$col] <- info$value
      matrice_marqueurs_rv(df)
    })
    
    observeEvent(input$save_marker_config, {
      data_to_save <- matrice_marqueurs_rv()
      req(data_to_save)
      p <- pipeline()
      if (inherits(p, "R6")) {
        # p$definir_config_marqueurs() enregistre la table ET reconstruit le dictionnaire
        # canal -> marqueur utilisé ensuite par get_label() pour tous les libellés d'axes
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
    observeEvent(input$add_group, {
      req(input$group_name, nchar(trimws(input$group_name)) > 0)
      nom <- trimws(input$group_name)
      if (!(nom %in% names(rv$groupes))) {
        temp <- rv$groupes; temp[[nom]] <- character(0); rv$groupes <- temp
        updateTextInput(session, "group_name", value = "")
      } else showNotification("Ce groupe existe déjà.", type = "warning")
    })
    
    output$list_groups_ui <- renderUI({
      if (length(rv$groupes) == 0) return(p("Aucun groupe créé.", style = "color:gray;"))
      tags$ul(lapply(names(rv$groupes), function(g) tags$li(tags$b(g))))
    })
    
    output$ui_assignation_groupes <- renderUI({
      req(rv$r_df_ech, length(rv$groupes) > 0)
      ns_func <- session$ns
      fluidRow(lapply(rv$r_df_ech$tube_name, function(tube) {
        column(6,
               selectInput(ns_func(paste0("assign_", make.names(tube))), label = tube,
                           choices  = c("(Non assigné)" = "", names(rv$groupes)),
                           selected = { grp <- names(rv$groupes)[sapply(rv$groupes, function(g) tube %in% g)]
                           if (length(grp)) grp[1] else "" })
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
    
    observeEvent(input$save_cohort, {
      req(rv$r_df_ech)
      new_groupes <- lapply(rv$groupes, function(x) character(0))
      for (tube in rv$r_df_ech$tube_name) {
        val <- input[[paste0("assign_", make.names(tube))]]
        if (!is.null(val) && val != "") new_groupes[[val]] <- c(new_groupes[[val]], tube)
      }
      rv$groupes <- new_groupes
      showNotification("Cohorte enregistrée !", type = "message")
    })
    
  })
}