library(R6)
library(ggplot2)
library(ggpointdensity)
library(flowCore)

EchantillonCompense <- R6Class(
  classname = "EchantillonCompense",
  public = list(
    mode = "Conventionnel", # l'utilisateur précise si ses données sont issues d'un cytomètre conventionnel ou spectral
    canaux = NULL, # liste des canaux (PE-A, FITC-A, etc) présents dans le fichier
    chemins_monomarques = NULL, # chemins utilisés pour accéder aux fichiers controles (monomarqués et unstained)
    tubes_monomarques = list(), # liste contenant le nom du fichier, le canal, le type (Monomarqué ou Unstained), et le chemin
    chemins_echantillons = NULL, # chemins utilisés pour accéder aux fichiers échantillons
    echantillons = list(), # liste des noms des fichiers échantillons
    dossier_racine = NULL, # nom du dossier dans lequel l'utilisateur souhaite accueillir les sorties d'AutoSpectral
    asp_control_file = "fcs_control_file.csv", # nom du fichier d'entrée d'AutoSpectral
    asp_config = NULL, # paramètres de configuration d'AutoSpectral
    flow.control = NULL, # liste des fichiers controles après l'étape de nettoyage
    gates = list(), # liste des gates réealisés
    spectra = NULL, # résultat de la fonction get.spectra d'AutoSpectral
    af_spectra = NULL, # résultat de la fonction get.af_spectra d'AutoSpectral
    variants = list(), # liste des variants 
    plots_unmixing = list(), # regroupe les dots plots réalisés
    echantillons_traites = list(), # regroupe tous les fichiers échantillons après la compensation ou l'unmixing
    
    initialize = function(df_monomarques = NULL, df_echantillons = NULL, chemin_racine = NULL, mode = "Conventionnel") { # initialiser la classe R6 
      self$mode <- mode # Enregistre le mode utilisé (Conventionnel ou Spectral)
      self$chemins_monomarques  <- df_monomarques # Enregistre les chemins et métadonnées des fichiers contrôles (monomarqués et unstained)
      self$chemins_echantillons <- df_echantillons # Enregistre les chemins et métadonnées des fichiers échantillons
      
      if (self$mode == "Conventionnel") { # Si le mode utilisé est Conventionnel
        self$dossier_racine   <- NULL # Alors on n'a pas besoin de dossier_racine puisqu'AutoSpectral ne sera pas utilisé
        self$asp_control_file <- NULL # Alors on n'a pas besoin de configurer le fichier CSV dont a besoin AutoSpectral
        
      } else if (self$mode == "Spectral") { # Si le mode utilisé est Spectral
        if (!is.null(chemin_racine)) { # Si un chemin racine est fourni (pour lire/déposer les résultats d'AutoSpectral)
          if (!dir.exists(chemin_racine)) { # Mais que ce chemin n'existe pas physiquement sur le disque dur
            stop("Le chemin fourni n'existe pas ou il y a une erreur de frappe dans : ", chemin_racine) # Message alertant l'utilisateur de l'erreur
          }
          self$dossier_racine <- chemin_racine # Enregistre le chemin racine validé dans l'objet R6
        } else {
          self$dossier_racine <- NULL # Si aucun chemin n'est transmis on l'initialise à NULL
        }
        self$asp_control_file <- "fcs_control_file.csv" # Configure le nom par défaut du fichier de métadonnées requis par AutoSpectral
      }
      
      if (!is.null(self$chemins_monomarques) && any(colnames(self$chemins_monomarques) == "type")) { # Vérifie si le tableau des monomarqués existe et si la colonne "type" est présente parmi les noms de colonnes
        lignes_monomarques <- self$chemins_monomarques[self$chemins_monomarques$type == "Monomarque", ]  # Filtre et extrait uniquement les lignes correspondant aux tubes de contrôles monomarqués
        
        if (nrow(lignes_monomarques) > 0 && any(colnames(lignes_monomarques) == "canal")) { # Vérifie qu'il y a au moins une ligne de monomarque ET que la colonne "canal" est bien présente
          self$canaux <- unique(lignes_monomarques$canal) # Récupère la liste unique de tous les canaux/détecteurs configurés
          self$canaux <- self$canaux[!is.na(self$canaux) & self$canaux != ""] # Nettoie la liste en retirant les valeurs manquantes (NA) ou vides ("")
        }
      } else { # Si le tableau n'est pas encore chargé
        self$canaux <- NULL # la liste des canaux reste initialisée à NULL
      }
    },
  
    
    charger_fcs = function() { # Méthode qui permet de charger les fichiers FCS en mémoire
      if (is.null(self$chemins_monomarques) || is.null(self$chemins_echantillons)) { # Si aucun fichier (contrôles ou échanillons) n'est configuré dans l'objet R6
        return(invisible(NULL)) # Alors on arrête immédiatement la fonction 
      }
      
      self$tubes_monomarques <- lapply(seq_len(nrow(self$chemins_monomarques)), function(i) { # Boucle sur chaque ligne du tableau des fichiers contrôles
        row <- self$chemins_monomarques[i, ] # Extrait la ligne actuelle (les métadonnées du tube courant)
        flowCore::read.FCS(row$chemin, transformation = FALSE, truncate_max_range = FALSE) # Lecture du fichier FCS associé 
      }) 
      
      noms_tubes <- sapply(seq_len(nrow(self$chemins_monomarques)), function(i) { # Boucle pour générer un vecteur d'identifiants textuels pour les contrôles
        row <- self$chemins_monomarques[i, ] # Extrait la ligne actuelle pour analyser ses propriétés
        if (!is.na(row$type) && row$type == "Unstained") return("TUBE_UNSTAINED") # Pour le cas du témoin négatif, on lui donne l'étiquette fixe TUBE_UNSTAINED
        return(row$canal) # Pour les autres, on utilise le nom de leur canal/détecteur (ex: FITC-A, PE-A) comme nom unique
      }) 
      names(self$tubes_monomarques) <- noms_tubes # Assigne les étiquettes générées aux éléments de la liste des objets flowFrame des controles
      
      if (self$mode == "Conventionnel") { # Si le mode de l'expérience est configuré sur Conventionnel
        self$echantillons <- lapply(self$chemins_echantillons$chemin, function(f) { # Boucle sur chaque chemin d'accès des échantillons
          flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE) # Lecture physique et brute de chaque fichier FCS échantillon
        }) # Fin du chargement de la cohorte conventionnelle
        names(self$echantillons) <- self$chemins_echantillons$tube_name # Assigne les vrais noms des tubes (ex: Patient_1) aux éléments chargés
        
      } else if (self$mode == "Spectral") { # Sinon, si le mode de l'expérience est configuré sur Spectral
        if (is.null(self$dossier_racine) || !dir.exists(self$dossier_racine)) { # Vérifie si le dossier requis par AutoSpectral est absent ou introuvable
          stop("Impossible de charger les fichiers : le dossier racine AutoSpectral est introuvable ou non configuré.") # Alerte et bloque le script si le dossier n'existe pas
        } # Fin de la vérification du dossier racine
        self$echantillons <- lapply(self$chemins_echantillons$chemin, function(f) { # Boucle sur chaque chemin d'accès de la cohorte spectrale
          flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE) # Lecture physique et brute de chaque fichier FCS échantillon
        }) # Fin du chargement de la cohorte spectrale
        names(self$echantillons) <- self$chemins_echantillons$tube_name # Assigne les vrais noms des tubes échantillons aux éléments spectralement chargés
      } 
    },
    
    get_label = function(fcs, canal) { # Permet d'extraire à partir des métadonnées du fichier FCS le nom des marqueurs biologiques
      if (is.null(fcs) || is.na(canal) || canal == "") return(canal) # Si on n'a ni fichier FCS ni canal donné, on retourne le nom du canal par défaut
      param_data <- flowCore::pData(flowCore::parameters(fcs)) # On accède aux métadonnées du fichier FCS, spécifiquement dans la table des paramètres
      if (is.null(param_data) || !any(colnames(param_data) == "name") || !any(colnames(param_data) == "desc")) { # Vérifie que la table des paramètres n'est pas vide et que les colonnes indispensables "name" et "desc" existent bien
        return(canal) # Si la structure des métadonnées est incomplète, on retourne le nom du canal brut par sécurité
      }
      
      desc <- param_data$desc # Extrait la liste des descriptions des marqueurs (ex: CD4, CD8, CD3) => correspond spécifiquement aux "$PXS"
      nom  <- param_data$name # Extrait la liste des noms techniques des détecteurs/canaux (ex: FITC-A, V3-A) => correspond spécifiquement aux "$PXN"
      idx  <- which(nom == canal) # Cherche la position de la ligne correspondant au canal demandé dans le tableau
      
      if (length(idx) == 0) return(canal) # Si le canal demandé n'est pas trouvé du tout dans le fichier FCS, on retourne le canal brut
      marqueur <- desc[idx[1]] # Récupère la description (le marqueur biologique) située sur la première ligne correspondante trouvée
      if (is.null(marqueur) || is.na(marqueur) || marqueur == "" || is.nan(marqueur)) { # Vérifie si le nom du marqueur extrait est vide, manquant (NA) ou invalide
        return(canal) # Si le marqueur n'est pas renseigné dans le fichier FCS, on retourne uniquement le nom du canal brut
      }
      return(paste0(canal, " | ", marqueur)) # Concatène et retourne le résultat sous un format propre (ex: "V3-A | CD4")
    },
    
    update_pipeline = function(etape, nom_echantillon = NULL) { # méthode qui permet de savoir à quelle étape on se situe
      if (is.null(self$pipeline_status)) self$pipeline_status <- list() # si on a une étape du pipeline réalisée alors on l'ajoute a la liste de statut du pipeline
      horodatage <- format(Sys.time(), "%H:%M:%S") # donne l'horodatage
      if (is.null(nom_echantillon)) { # si on a pas de 
        self$pipeline_status[[etape]] <- list(
          statut = "Terminé", # assigne le statut à Terminé
          date = horodatage, 
          portee = "Tous les échantillons"
        )
      } else {
        if (is.null(self$pipeline_status[[etape]])) self$pipeline_status[[etape]] <- list()
        self$pipeline_status[[etape]][[nom_echantillon]] <- list(
          statut = "Terminé", 
          date = horodatage
        )
      }
    },
    
    # DEBUT DES ETAPES SPECIFIQUES A AUTOSPECTRAL 
    
    lancer_asp = function(type_cytometre = "aurora") {
      
      if (self$mode != "Spectral") {
        stop("Cette méthode nécessite le mode 'Spectral'.")
      }
      
      self$asp_config <- AutoSpectral::get.autospectral.param(cytometer = type_cytometre)
      control_dir <- path.expand(dirname(self$chemins_monomarques$chemin[1]))
      fichier_csv <- file.path(path.expand(self$dossier_racine), "fcs_control_file")
      
      if (file.exists(paste0(fichier_csv, ".csv"))) {
        file.remove(paste0(fichier_csv, ".csv"))
      }
      
      AutoSpectral::create.control.file(
        control.dir = control_dir,
        asp = self$asp_config,
        filename = fichier_csv
      )
    },
    
    verifier_asp = function(warning = 5000, error = 1000) {
      
      if (is.null(self$asp_config)) {
        stop("Erreur : asp_config est NULL. Lancez d'abord lancer_asp().")
      }
      chemin_csv_complet <- file.path(path.expand(self$dossier_racine), "fcs_control_file.csv")
      dossier_fcs <- path.expand(dirname(self$chemins_monomarques$chemin[1]))
      
      if (!file.exists(chemin_csv_complet)) {
        stop("Fichier de contrôle introuvable à : ", chemin_csv_complet)
      }
      
      verification <- AutoSpectral::check.control.file(
        control.dir = dossier_fcs, 
        control.def.file = chemin_csv_complet, 
        asp = self$asp_config,
        min.event.warning = warning, 
        min.event.error = error
      )
      
      if (is.null(verification)) {
      } else {
        print(verification)
      }
      return(verification)
    },
    
    definir_gates_landmarks = function(control_name, n.cells = 2000, percentile = 70, 
                                       grid.n = 100, bandwidth.factor = 1, 
                                       fsc.channel = NULL, ssc.channel = NULL) {
      
      if (is.null(self$asp_config)) {
        stop("La configuration ASP n'est pas initialisée. Appelez d'abord lancer_asp().")
      }
      
      old_wd <- getwd()
      setwd(path.expand(self$dossier_racine))
      on.exit(setwd(old_wd))
      
      output_dir <- file.path(self$dossier_racine, "figure_gate")
      if (!dir.exists(output_dir)) dir.create(output_dir)
      
      gate_result <- AutoSpectral::define.gate.landmarks(
        control.file = "fcs_control_file.csv", 
        control.dir = path.expand(dirname(self$chemins_monomarques$chemin[1])),
        asp = self$asp_config,
        gate.name = control_name,
        n.cells = n.cells,
        percentile = percentile,
        grid.n = grid.n,
        bandwidth.factor = bandwidth.factor,
        fsc.channel = fsc.channel,
        ssc.channel = ssc.channel,
        output.dir = output_dir
      )
      
      if (is.null(self$gates)) {
        self$gates <- list()
      }
      self$gates[[control_name]] <- gate_result
      return(gate_result)
    },
    
    definir_gates_density = function(control_name, n.cells = 2000, grid.n = 100, 
                                     bandwidth.factor = 1, fsc.channel = NULL, 
                                     ssc.channel = NULL) {
      
      if (is.null(self$asp_config)) {
        stop("La configuration ASP n'est pas initialisée. Appelez d'abord lancer_asp().")
      }
      
      old_wd <- getwd()
      setwd(path.expand(self$dossier_racine))
      on.exit(setwd(old_wd))
      
      output_dir <- file.path(self$dossier_racine, "figure_gate")
      if (!dir.exists(output_dir)) dir.create(output_dir)
      
      gate_result <- AutoSpectral::define.gate.density(
        control.file = "fcs_control_file.csv",
        control.dir = path.expand(dirname(self$chemins_monomarques$chemin[1])),
        asp = self$asp_config,
        gate.name = control_name,
        n.cells = n.cells,
        grid.n = grid.n,
        bandwidth.factor = bandwidth.factor,
        fsc.channel = fsc.channel,
        ssc.channel = ssc.channel,
        output.dir = output_dir
      )
      
      if (is.null(self$gates)) {
        self$gates <- list()
      }
      
      self$gates[[control_name]] <- gate_result
      return(gate_result)
    },
    
    definir_tune_gates = function(gate.name, n_cells = 2000, percentile = 70, bandwidth = 1) {
      
      csv_file <- file.path(path.expand(self$dossier_racine), "fcs_control_file.csv")
      
      if (!file.exists(csv_file)) {
        stop("Fichier CSV introuvable à : ", csv_file, ". Avez-vous bien lancé lancer_asp() ?")
      }
      
      dossier_fcs <- path.expand(dirname(self$chemins_monomarques$chemin[1]))
      output_dir <- file.path(self$dossier_racine, "figure_gate_tuning")
      if (!dir.exists(output_dir)) dir.create(output_dir)
      
      gate_tuned <- AutoSpectral::tune.gate(
        control.file = csv_file, 
        control.dir = dossier_fcs, 
        asp = self$asp_config, 
        gate.name = gate.name, 
        n.cells = n_cells,
        percentiles = percentile,
        bandwidth.factor = bandwidth,
        output.dir = output_dir,
        filename = paste0("tuned_", gate.name)
      )
      
      self$gates[[gate.name]] <- gate_tuned
      return(invisible(gate_tuned))
    },
    
    charger_et_nettoyer = function() {
      # 1. Vérification optionnelle : on ne bloque plus si gates est vide, 
      # mais on informe l'utilisateur.
      if (is.null(self$gates) || length(self$gates) == 0) {
        message("⚠️ Aucune gate définie. Les fichiers seront traités sans filtrage spatial.")
      }
      
      dossier_fcs <- path.expand(dirname(self$chemins_monomarques$chemin[1]))
      dossier_figures <- path.expand(file.path(self$dossier_racine, "figure_clean_controls"))
      
      if (!dir.exists(dossier_figures)) {
        dir.create(dossier_figures, recursive = TRUE)
      }
      
      old_wd <- getwd()
      setwd(dossier_figures)
      on.exit(setwd(old_wd))
      
      # 2. Préparation du contrôle : 
      # AutoSpectral::define.flow.control gérera les fichiers avec gate.define = FALSE 
      # automatiquement si votre CSV est bien configuré.
      flow_res <- AutoSpectral::define.flow.control(
        control.dir = dossier_fcs, 
        control.def.file = file.path(self$dossier_racine, self$asp_control_file),
        asp = self$asp_config,
        gate.list = if(length(self$gates) > 0) self$gates else NULL
      )
      
      # 3. Nettoyage :
      # Si aucune gate n'est définie, clean.controls traitera les fichiers "bruts" 
      # (autofluorescence uniquement).
      flow_cleaned <- AutoSpectral::clean.controls(
        flow.control = flow_res,
        asp = self$asp_config,
        main.figures = TRUE 
      )
      
      self$flow.control <- flow_cleaned
      message("✅ Chargement et nettoyage terminés avec succès.")
      return(invisible(self$flow.control))
    },
    
    extraire_fluorophore_spectre = function() {
      if (is.null(self$flow.control)) {
        stop("Erreur : flow.control n'est pas chargé.")
      }
      
      old_wd <- getwd()
      setwd(path.expand(self$dossier_racine))
      
      spectra_result <- AutoSpectral::get.fluorophore.spectra(
        flow.control = self$flow.control,
        asp = self$asp_config
      )
      
      setwd(old_wd)
      self$spectra <- spectra_result
      return(invisible(self$spectra))
    },
    
    # =========================
    # étapes optionnelles
    # =========================
    extraire_spectre_af = function(unstained_fcs_path, tissue_name, refine = TRUE) {
      
      if (is.null(self$spectra)) {
        stop("Erreur : Les spectres fluorophores n'ont pas été extraits. Lancez extraire_fluorophore_spectre() d'abord.")
      }
      
      dossier_figures <- path.expand(file.path(self$dossier_racine, "figure_autofluorescence"))
      dossier_tables <- path.expand(file.path(self$dossier_racine, "table_autofluorescence"))
      if (!dir.exists(dossier_figures)) dir.create(dossier_figures, recursive = TRUE)
      if (!dir.exists(dossier_tables)) dir.create(dossier_tables, recursive = TRUE)
      
      af_result <- AutoSpectral::get.af.spectra(
        unstained.sample = path.expand(unstained_fcs_path),
        asp = self$asp_config,
        spectra = self$spectra,
        refine = refine,
        figures = TRUE,
        plot.dir = dossier_figures,  # Forcé ici
        table.dir = dossier_tables,  # Forcé ici
        title = paste("Autofluorescence -", tissue_name)
      )
      
      if (is.null(self$af_spectra)) {
        self$af_spectra <- list()
      }
      
      self$af_spectra[[tissue_name]] <- af_result
      return(invisible(af_result))
    },
    
    preparer_variants_spectraux = function(tissue_af_name = NULL, refine = TRUE) {
      
      nom_tissu <- if (!is.null(tissue_af_name)) {
        tissue_af_name
      } else if (length(self$af_spectra) > 0) {
        names(self$af_spectra)[1]
      } else {
        stop("Erreur : Aucune AF trouvée dans self$af_spectra. Lancez extraire_spectre_af() d'abord.")
      }
      
      if (is.null(self$spectra) || is.null(self$af_spectra[[nom_tissu]])) {
        stop("Erreur : Spectres ou AF non trouvés pour le tissu : ", nom_tissu)
      }
      
      dossier_variants <- path.expand(file.path(self$dossier_racine, "figure_spectral_variants"))
      if (!dir.exists(dossier_variants)) dir.create(dossier_variants, recursive = TRUE)
      
      chemin_dossier_controles <- path.expand(dirname(self$chemins_monomarques$chemin[1]))
      chemin_fichier_csv <- path.expand(file.path(self$dossier_racine, self$asp_control_file))
      
      variants_result <- AutoSpectral::get.spectral.variants(
        control.dir = chemin_dossier_controles,
        control.def.file = chemin_fichier_csv,
        asp = self$asp_config,
        spectra = self$spectra,
        af.spectra = self$af_spectra[[nom_tissu]],
        refine = refine,
        figures = TRUE, 
        output.dir = dossier_variants 
      )
      
      if (is.null(self$variants)) self$variants <- list()
      self$variants[[nom_tissu]] <- variants_result
      return(invisible(variants_result))
    },
    
    # ===========================   
    # unmixing 
    # ===========================
    
    unmix_fcs = function(fcs_file_path, tissue_name = NULL, method = "AutoSpectral", speed = "slow") {
      
      dossier_sortie <- file.path(self$dossier_racine, "AutoSpectral_unmixed")
      if (!dir.exists(dossier_sortie)) dir.create(dossier_sortie)
      
      n_detectors <- ncol(self$spectra[[1]]) 
      af_s <- if (!is.null(tissue_name) && !is.null(self$af_spectra[[tissue_name]])) {
        self$af_spectra[[tissue_name]]
      } else if (length(self$af_spectra) > 0) {
        self$af_spectra[[1]]
      } else {
        matrix(nrow = 0, ncol = n_detectors) 
      }
      
      var_s <- if (!is.null(tissue_name) && !is.null(self$variants[[tissue_name]])) {
        self$variants[[tissue_name]]
      } else {
        NULL
      }
      
      AutoSpectral::unmix.fcs(
        fcs.file = fcs_file_path,
        spectra = self$spectra,
        asp = self$asp_config,
        flow.control = self$flow.control,
        method = method,
        af.spectra = af_s,
        spectra.variants = var_s,
        speed = speed,
        output.dir = dossier_sortie,  
        parallel = TRUE
      )
      
    },
    
    unmix_folder = function(folder_path, tissue_name = NULL, method = "AutoSpectral", speed = "slow") {
      
      dossier_sortie <- file.path(self$dossier_racine, "AutoSpectral_unmixed")
      if (!dir.exists(dossier_sortie)) dir.create(dossier_sortie, recursive = TRUE)
      n_detectors <- ncol(self$spectra[[1]]) 
      
      af_s <- if (!is.null(tissue_name) && !is.null(self$af_spectra[[tissue_name]])) {
        self$af_spectra[[tissue_name]]
      } else if (length(self$af_spectra) > 0) {
        self$af_spectra[[1]]
      } else {
        matrix(nrow = 0, ncol = n_detectors) 
      }
      
      var_s <- if (!is.null(tissue_name) && !is.null(self$variants[[tissue_name]])) {
        self$variants[[tissue_name]]
      } else {
        NULL
      }
      
      AutoSpectral::unmix.folder(
        fcs.dir = folder_path,
        spectra = self$spectra,
        asp = self$asp_config,
        flow.control = self$flow.control,
        method = method,
        af.spectra = af_s,
        spectra.variants = var_s,
        speed = speed,
        output.dir = dossier_sortie, 
        parallel = TRUE
      )
    },
    
    verifier_qualite_unmix = function(fluorophore, single_stained_fcs, unstained_fcs, cytometer = "aurora", gate = TRUE) {
      dossier_gates <- file.path(obj$dossier_racine, "figure_gate")
      if (!dir.exists(dossier_gates)) dir.create(dossier_gates, recursive = TRUE)
      if (is.null(self$spectra)) {
        stop("Erreur : Les spectres n'ont pas été extraits. Lancez extraire_fluorophore_spectre() d'abord.")
      }
      
      if (!fluorophore %in% rownames(self$spectra)) {
        stop("Erreur : Fluorophore introuvable dans les spectres.")
      }
      
      dossier_figures <- path.expand(file.path(self$dossier_racine, "figure_compare_unmix"))
      if (!dir.exists(dossier_figures)) dir.create(dossier_figures, recursive = TRUE)
      spectre_cible <- self$spectra[fluorophore, ]
      
      qc_result <- AutoSpectral::compare.unmix(
        single.stained.fcs = path.expand(single_stained_fcs),
        unstained.fcs      = path.expand(unstained_fcs),
        fluorophore        = fluorophore,
        spectra            = self$spectra,
        ref.spectrum       = spectre_cible,
        test.spectrum      = spectre_cible,
        cytometer          = cytometer,
        gate               = gate, 
        plot.dir           = dossier_figures
      )
      
      return(invisible(qc_result))
    },
    
    # exportation
    
    charger_fcs_unmixes = function(dossier = "AutoSpectral_unmixed") {
      chemin_complet <- file.path(self$dossier_racine, dossier)
      if(!dir.exists(chemin_complet)) stop("Dossier introuvable : ", chemin_complet)
      
      fichiers <- list.files(chemin_complet, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
      
      for (f in fichiers) {
        nom_cle <- basename(f)
        self$echantillons_traites[[nom_cle]] <- flowCore::read.FCS(f, truncate_max_range = FALSE)
      }
      message("Chargement terminé : ", length(fichiers), " échantillons importés depuis ", dossier)
    },
  
    
    get_chemins_figures = function(control_name) {
      dossier <- file.path(self$dossier_racine, "figure_gate")
      if (!dir.exists(dossier)) return(NULL)
      list.files(dossier, pattern=paste0(control_name, ".*\\.png$"), full.names=TRUE)
    },
    
    visualiser_unmixing = function(nom_fichier_fcs, canal_x, canal_y, cofacteur = 150, max_points = 10000) {
      fcs_unmixed <- self$echantillons_traites[[nom_fichier_fcs]]
      if (is.null(fcs_unmixed)) stop("Fichier introuvable en mémoire.")
      trans_list <- flowCore::transformList(c(canal_x, canal_y), flowCore::arcsinhTransform(a = 0, b = 1/cofacteur, c = 0))
      mat <- flowCore::exprs(flowCore::transform(fcs_unmixed, trans_list))[, c(canal_x, canal_y)]
      indices <- sample(seq_len(nrow(mat)), min(nrow(mat), max_points))
      df <- as.data.frame(mat[indices, ])
      colnames(df) <- c("Axe_X", "Axe_Y")
      
      ggplot(df, aes(x = Axe_X, y = Axe_Y)) +
        ggpointdensity::geom_pointdensity(size = 0.2, alpha = 0.5) +
        scale_color_gradientn(
          colours = c("darkblue", "blue", "cyan", "greenyellow", "yellow", "darkorange", "red"),
          name = "Densité"
        ) +
        theme_bw() + 
        labs(
          title = paste("Résultat après Unmixing :", nom_fichier_fcs), 
          x = self$get_label(fcs_unmixed, canal_x), 
          y = self$get_label(fcs_unmixed, canal_y)
        )
    },
    
    visualiser_figures = function(dossier_nom) {
      chemin_dossier <- file.path(self$dossier_racine, dossier_nom)
      fichiers <- list.files(chemin_dossier, pattern="\\.(jpg|jpeg|png)$", full.names=TRUE, ignore.case=TRUE)
      if (length(fichiers) == 0) return(message("Aucune image."))
      
      html_elements <- sapply(fichiers, function(f) {
        mime <- ifelse(grepl("\\.(jpg|jpeg)", f, ignore.case=TRUE), "image/jpeg", "image/png")
        paste0("<div><h3>", basename(f), "</h3><img src='", base64enc::dataURI(file=f, mime=mime), "' style='max-width:100%'></div>")
      })
      
      temp_html <- tempfile(fileext=".html")
      writeLines(c("<html><body>", html_elements, "</body></html>"), temp_html)
      rstudioapi::viewer(temp_html)
    }
    
  ),
  private = list(df_control_file = NULL)
)
