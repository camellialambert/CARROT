library(R6)
library(ggplot2)
library(ggpointdensity)
library(raster)
library(flowCore)
library(base64enc)

calculer_densite_raster <- function(x, y, xlim, ylim, res = 800, lissage = FALSE) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NULL)
  
  if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5)
  if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
  
  r_vide <- raster::raster(xmn = xlim[1], xmx = xlim[2],
                           ymn = ylim[1], ymx = ylim[2],
                           nrows = res, ncols = res)
  r_densite <- raster::rasterize(cbind(x, y), r_vide, fun = "count", background = 0)
  
  if (isTRUE(lissage)) {
    r_densite <- raster::focal(r_densite, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE, pad = TRUE)
  }
  
  df <- as.data.frame(r_densite, xy = TRUE)
  names(df) <- c("X", "Y", "densite")
  df <- df[!is.na(df$densite) & df$densite > 0, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df
}

PALETTE_DENSITE <- c("darkblue", "blue", "cyan", "greenyellow", "yellow", "darkorange", "red")

COLORSCALE_DENSITE_PLOTLY <- local({
  n <- length(PALETTE_DENSITE)
  lapply(seq_len(n), function(i) list((i - 1) / (n - 1), PALETTE_DENSITE[i]))
})


generer_image_densite_base64 <- function(x, y, xlim, ylim, res = 800, largeur_px = 600, hauteur_px = 600) {
  df <- calculer_densite_raster(x, y, xlim, ylim, res = res, lissage = FALSE)
  if (is.null(df)) return(NULL)
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x = X, y = Y, fill = densite)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE) +
    ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none",
                   plot.margin = ggplot2::margin(0, 0, 0, 0))
  
  fichier_tmp <- tempfile(fileext = ".png")
  on.exit(unlink(fichier_tmp), add = TRUE)
  ggplot2::ggsave(fichier_tmp, plot = g, width = largeur_px / 96, height = hauteur_px / 96,
                  dpi = 96, bg = "transparent")
  
  base64enc::dataURI(file = fichier_tmp, mime = "image/png")
}

CARROT <- R6Class(
  classname = "CARROT",
  public = list(
    # Importation des fichiers
    mode = "Conventionnel", # l'utilisateur précise si ses données sont issues d'un cytomètre conventionnel ou spectral
    mapping_canal_fichier = NULL,
    canaux = NULL, # liste des canaux (PE-A, FITC-A, etc) présents dans le fichier
    chemins_monomarques = NULL, # chemins utilisés pour accéder aux fichiers controles (monomarqués et unstained)
    tubes_monomarques = list(), # liste contenant le nom du fichier, le canal, le type (Monomarqué ou Unstained), et le chemin
    chemins_echantillons = NULL, # chemins utilisés pour accéder aux fichiers échantillons
    echantillons = list(), # liste des noms des fichiers échantillons
    seed = 123,
    # Variables de Compensation
    config_transformations = NULL, # contient les cofacteurs de transformation utilisés pour chaque échantillon
    trans_list = NULL, # Contient la liste des fonctions de transformation (ex: Arcsinh) ou des paramètres appliqués aux canaux
    monomarques_trans = NULL, # Contient la liste des fichiers FCS des contrôles (monomarqués et unstained) après application de la transformation
    bornes_gates_pos = list(), # contient les bornes des gates positifs pour chacun des tubes monomarqués
    bornes_gates_neg = list(), # contient les bornes des gates négatifs pour chacun des tubes monomarqués
    source_neg_utilisee = list(), # Enregistre la source du signal négatif pour chaque canal (soit le tube Unstained global, soit la population négative interne du tube)
    gates_positifs = list(), # contient tous les évènements faisant partie du gate positif, pour chacun des tubes monomarqués
    gates_negatifs = list(), # contient tous les évènements faisant partie du gate négatif, pour chacun des tubes monomarqués
    plots_gates = list(), # contient toutes les figures des gates (courbes de densité)
    S_matrix = NULL, # contient la matrice de spillover
    S_matrices_par_echantillon = list(), # contient les matrices de spillover par échantillon
    comparaison_medianes = list(), # contient, par échantillon (ou "GLOBAL"), les médianes brutes positive/négative (et leur écart) mesurées avant/après compensation
    plots_compensation = list(), # contient toutes les figures biplots, avant et/ou après compensation
    
    # Variable de l'unmixing
    dossier_racine = NULL, # nom du dossier dans lequel l'utilisateur souhaite accueillir les sorties d'AutoSpectral
    dossier_monomarques = NULL,
    dossier_echantillons = NULL,
    dossier_controles_asp = NULL, # dossier de travail isolé ne contenant QUE les fichiers monomarqués/unstained réellement sélectionnés à l'import (voir preparer_dossier_controles_asp) : empêche AutoSpectral de scanner tout un dossier source et d'y confondre échantillons, contrôles non sélectionnés ou fichiers d'une expérience précédente
    asp_control_file = "fcs_control_file.csv", # nom du fichier d'entrée d'AutoSpectral
    asp_config = NULL, # paramètres de configuration d'AutoSpectral
    asp_correspondance_incomplete = FALSE, # TRUE si AutoSpectral a signalé des fluorophores/marqueurs "No match" lors de create.control.file() (non bloquant : à corriger manuellement dans le tableau puis à enregistrer)
    flow.control = NULL, # liste des fichiers controles après l'étape de nettoyage
    gates = list(), # liste des gates réealisés
    spectra = NULL, # résultat de la fonction get.spectra d'AutoSpectral
    af_spectra = NULL, # résultat de la fonction get.af_spectra d'AutoSpectral
    variants = list(), # liste des variants 
    plots_unmixing = list(), # regroupe les dots plots réalisés
    
    #pour compensation OU unmixing
    echantillons_traites = list(), # regroupe tous les fichiers échantillons après la compensation ou l'unmixing 
    
    # Variables de Prétraitement 
    pipeline_status = list(), # Enregistre l'état d'avancement du nettoyage (ex: TRUE/FALSE ou étapes validées) pour chaque échantillon
    canaux_bordures = NULL, # Liste technique des canaux ou détecteurs ciblés pour l'analyse et le retrait des événements saturés aux bordures
    post_PeacoQC = list(), # Contient les données des fichiers FCS après le contrôle qualité de PeacoQC (retrait des anomalies de débit et d'instabilité du signal)
    post_flowAI = list(), # Contient les données des fichiers FCS après le contrôle qualité alternatif flowAI (vérification du débit, de la lueur et0 de la stabilité)
    rapports_flowai = list(), # Contient, par échantillon, le chemin du dossier temporaire regroupant le rapport natif flowAI (HTML, TXT, PNG) généré par flow_auto_qc
    parametres_peacoqc_utilises = NULL, # Mémorise les derniers réglages PeacoQC appliqués à la cohorte (pour affichage dans le résumé PDF)
    parametres_flowai_utilises = NULL, # Mémorise les derniers réglages flowAI appliqués à la cohorte (pour affichage dans le résumé PDF)
    post_retrait_bordures = list(), # Stocke la matrice d'expression des échantillons nettoyée des signaux saturés (valeurs maximales ou minimales des détecteurs)
    gate_debris = list(), # Contient les coordonnées et les structures géométriques des fenêtres (gates) de sélection des cellules (retrait des débris en FSC vs SSC)
    post_debris = list(), # Stocke les données des échantillons filtrées où seuls les événements correspondants aux cellules (hors débris) ont été conservés
    gates_history = list(), # Historique ordonné des gates nommés appliqués : list(nom_gate -> list(nom_echantillon -> list(polygone, cx, cy, post_data, n_avant, n_apres)))
    gate_doublets_FSC = list(), # Contient les coordonnées du gate de discrimination des doublets basé sur les paramètres du Forward Scatter (ex: FSC-A vs FSC-H)
    gate_doublets_SSC = list(), # Contient les coordonnées du gate de discrimination des doublets basé sur les paramètres du Side Scatter (ex: SSC-A vs SSC-H)
    post_doublets_FSC = list(), # Stocke les données des échantillons après l'élimination des doublets par le filtre FSC
    post_doublets_SSC = list(), # Stocke les données des échantillons après l'élimination des doublets par le filtre SSC
    post_doublets_final = list(), # Stocke les données des échantillons totalement débarrassées de tous les agrégats cellulaires (doublets FSC et SSC combinés)
    gate_viabilite = list(), # Contient les coordonnées de la fenêtre de sélection des cellules vivantes (basée sur le canal du marqueur de viabilité / "live-dead")
    post_viabilite = list(), # Stocke le produit final du prétraitement : les données des cellules vivantes, singulets et de haute qualité, prêtes pour l'analyse immunologique
    
    # Variables de Visualisation
    plots_peacoqc = list(), # contient toutes les figures après application de PeacoQC
    plots_peacoqc_natif = list(), # contient les chemins des PNG diagnostiques natifs générés par PeacoQC::PeacoQC (bins retirés par canal)
    plots_flowai = list(), # contient toutes les figures après application de flowAI (archivées pour réutilisation dans le résumé PDF d'export)
    plots_debris = list(), # contient toutes les figures après application du gate de débris
    plots_doublets = list(), # contient toutes les figures après application du gate de doublets
    plots_viabilite = list(), # contient toutes les figures après application du gate de viabilité
    
    # variables d'analyses
    post_transformation = list(), # Stocke les données des échantillons après application de la transformation (ex: arcsinh)
    cofactor_transformation = NULL,
    canaux_transformes = list(),
    
    # variable 
    pipeline_callback = NULL,
    
    #pour siny
    config_marqueurs = NULL, # data.frame brut (lignes = tubes/échantillons, colonnes = canaux) des annotations marqueurs saisies par l'utilisateur
    dictionnaire_marqueurs = NULL, # vecteur nommé canal -> marqueur biologique, dérivé de config_marqueurs, utilisé par get_label() pour tous les libellés d'axes
    
    #pour l'import sans fichiers contrôles (échantillons déjà compensés / unmixés)
    sans_controles = FALSE, # TRUE si l'utilisateur n'a fourni aucun tube monomarqué/unstained
    deja_traite = FALSE, # TRUE si l'utilisateur indique que ses échantillons sont déjà compensés (conventionnel) ou unmixés (spectral)
    cofacteur_defaut = 150, # cofacteur Arcsinh utilisé par défaut pour construire trans_list automatiquement en l'absence de contrôles (ajustable ensuite via l'onglet Transformation)
    
    
    #initialisation de la classe
    initialize = function(df_monomarques = NULL, df_echantillons = NULL, chemin_racine = NULL, mode = "Conventionnel", deja_traite = FALSE) { # initialiser la classe R6 
      self$mode <- mode # Enregistre le mode utilisé (Conventionnel ou Spectral)
      # Un data.frame de 0 ligne est traité comme une absence de contrôles
      if (!is.null(df_monomarques) && nrow(df_monomarques) == 0) df_monomarques <- NULL
      self$chemins_monomarques  <- df_monomarques # Enregistre les chemins et métadonnées des fichiers contrôles (monomarqués et unstained), NULL si aucun contrôle fourni
      self$chemins_echantillons <- df_echantillons # Enregistre les chemins et métadonnées des fichiers échantillons
      self$sans_controles <- is.null(self$chemins_monomarques) # Mémorise si l'utilisateur travaille sans tubes monomarqués/unstained
      self$deja_traite    <- isTRUE(deja_traite) # Mémorise si les échantillons sont annoncés comme déjà compensés/unmixés
      
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
    
    
    charger_fcs = function() {
      # ── Tubes contrôles (monomarqués / unstained) : optionnels ──────────────
      # Si l'utilisateur n'a fourni aucun fichier contrôle (échantillons déjà
      # compensés ou déjà unmixés), on saute simplement cette étape.
      if (!is.null(self$chemins_monomarques) && nrow(self$chemins_monomarques) > 0) {
        self$tubes_monomarques <- lapply(seq_len(nrow(self$chemins_monomarques)), function(i) {
          row <- self$chemins_monomarques[i, ]
          flowCore::read.FCS(row$chemin, transformation = FALSE, truncate_max_range = FALSE)
        })
        noms_tubes <- sapply(seq_len(nrow(self$chemins_monomarques)), function(i) {
          row <- self$chemins_monomarques[i, ]
          if (!is.na(row$type) && row$type == "Unstained") return("TUBE_UNSTAINED")
          return(row$canal)
        })
        names(self$tubes_monomarques) <- noms_tubes
        self$sans_controles <- FALSE
      } else {
        self$tubes_monomarques <- list()
        self$sans_controles <- TRUE
      }
      
      # ── Échantillons biologiques : toujours requis ───────────────────────────
      self$echantillons = lapply(self$chemins_echantillons$chemin, function(f) {
        flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
      }) 
      names(self$echantillons) <- self$chemins_echantillons$tube_name
      
      if (self$sans_controles && self$mode == "Conventionnel" && self$deja_traite) {
        invisible(tryCatch(self$importer_spillover_fcs(), error = function(e) NULL))
      }
      
      if (self$sans_controles) {
        invisible(tryCatch(self$transformer_fcs(cofacteur = self$cofacteur_defaut),
                           error = function(e) NULL))
      }
      
      # ── Transition directe vers echantillons_traites ────────────────────────
      # Si l'utilisateur a répondu "Oui" à "Données déjà compensées/démixées ?",
      # les fichiers importés sont, par définition, déjà dans leur état final
      # (compensation conventionnelle ou unmixing spectral déjà appliqués en
      # amont, hors de l'application). On ne doit donc PAS repasser par
      # self$compenser() (qui réappliquerait une matrice de compensation sur des
      # données déjà compensées) ni par le pipeline AutoSpectral : on copie
      # simplement self$echantillons tel quel dans self$echantillons_traites,
      # qu'il y ait ou non des tubes contrôles associés. C'est ce champ que lisent
      # ensuite les modules QC et Prétraitement, quel que soit le type de
      # cytomètre (Conventionnel ou Spectral).
      if (isTRUE(self$deja_traite)) {
        self$echantillons_traites <- self$echantillons
      }
    },
    
    importer_spillover_fcs = function(nom_echantillon = NULL) {
      if (is.null(self$echantillons) || length(self$echantillons) == 0) {
        stop("Aucun échantillon chargé.")
      }
      nom <- if (!is.null(nom_echantillon)) nom_echantillon else names(self$echantillons)[1]
      fcs <- self$echantillons[[nom]]
      if (is.null(fcs)) stop("Échantillon introuvable : ", nom)
      
      mat <- extraire_spillover_depuis_fcs(fcs)
      
      if (is.null(mat)) {
        cles_trouvees <- lister_cles_spillover_fcs(fcs)
        if (length(cles_trouvees) == 0) {
          return(list(succes = FALSE,
                      message = "Aucun mot-clé de compensation (SPILL/SPILLOVER/COMP) n'a été trouvé dans les métadonnées de ce fichier FCS."))
        } else {
          return(list(succes = FALSE,
                      message = paste0("Mot-clé(s) détecté(s) mais non exploitable(s) : ",
                                       paste(cles_trouvees, collapse = ", "),
                                       ". Vérifiez le format (attendu : \"n,chan1,...,chanN,v11,...,vNN\").")))
        }
      }
      
      # On ne conserve que les canaux réellement présents dans l'échantillon
      canaux_valides <- intersect(rownames(mat), flowCore::colnames(fcs))
      if (length(canaux_valides) == 0) {
        return(list(succes = FALSE,
                    message = paste0("Une matrice a été trouvée (canaux : ",
                                     paste(rownames(mat), collapse = ", "),
                                     ") mais aucun de ces canaux ne correspond aux paramètres de l'échantillon (",
                                     paste(flowCore::colnames(fcs), collapse = ", "), ").")))
      }
      
      mat <- mat[canaux_valides, canaux_valides, drop = FALSE]
      
      self$S_matrix <- mat
      self$canaux   <- canaux_valides
      return(list(succes = TRUE,
                  message = paste0("Matrice importée pour ", length(canaux_valides), " canaux.")))
    },
    
    get_label = function(fcs, canal) { # Permet d'extraire le nom des marqueurs biologiques à afficher sur les axes des graphiques
      if (is.na(canal) || canal == "") return(canal) # Si aucun canal n'est donné, on retourne tel quel
      
      # Priorité à l'annotation manuelle de l'utilisateur (onglet "Configuration Marqueurs"),
      # car elle prime sur la description brute embarquée dans le fichier FCS.
      if (!is.null(self$dictionnaire_marqueurs) && canal %in% names(self$dictionnaire_marqueurs)) {
        marqueur_manuel <- self$dictionnaire_marqueurs[[canal]]
        if (!is.null(marqueur_manuel) && !is.na(marqueur_manuel) && nchar(trimws(marqueur_manuel)) > 0) {
          return(paste0(canal, " | ", marqueur_manuel))
        }
      }
      
      if (is.null(fcs)) return(canal) # Sans fichier FCS de référence, impossible d'aller chercher plus loin
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
    
    
    definir_config_marqueurs = function(df_config) {
      self$config_marqueurs <- df_config
      if (is.null(df_config) || nrow(df_config) == 0 || ncol(df_config) == 0) {
        self$dictionnaire_marqueurs <- NULL
        return(invisible(NULL))
      }
      dict <- sapply(colnames(df_config), function(cx) {
        valeurs <- trimws(as.character(df_config[[cx]]))
        valeurs <- valeurs[!is.na(valeurs) & nchar(valeurs) > 0]
        if (length(valeurs) > 0) valeurs[1] else NA_character_
      })
      dict <- dict[!is.na(dict)]
      self$dictionnaire_marqueurs <- as.list(dict)
      invisible(self$dictionnaire_marqueurs)
    },
    
    update_pipeline = function(etape, nom_echantillon = NULL) { # méthode qui permet de savoir à quelle étape on se situe
      if (is.null(self$pipeline_status)) self$pipeline_status <- list() 
      horodatage <- format(Sys.time(), "%H:%M:%S") 
      if (is.null(nom_echantillon)) { 
        self$pipeline_status[[etape]] <- list(
          statut = "Terminé", 
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
      
      if (is.function(self$pipeline_callback)) {
        self$pipeline_callback(etape, nom_echantillon)
      }
    },
    
    # ===========================
    # COMPENSATION (Conventionnel)
    # ===========================
    
    transformer_fcs = function(cofacteur) {
      
      avec_controles <- !is.null(self$tubes_monomarques) && length(self$tubes_monomarques) > 0
      
      if (avec_controles) {
        # Recalcule les canaux fluo depuis le premier tube chargé
        # (plus fiable que self$canaux qui peut être NULL si init sans arguments)
        tous_canaux <- flowCore::colnames(self$tubes_monomarques[[1]])
      } else if (!is.null(self$echantillons) && length(self$echantillons) > 0) {
        # Pas de contrôles : on dérive les canaux directement des échantillons
        # (cas d'un import sans tubes monomarqués/unstained, ex: données déjà compensées)
        tous_canaux <- flowCore::colnames(self$echantillons[[1]])
      } else {
        stop("Aucun tube monomarqué ni échantillon chargé.")
      }
      canaux_fluo <- tous_canaux[!grepl("fsc|ssc|time", tous_canaux, ignore.case = TRUE)]
      
      arsinh_fun <- flowCore::arcsinhTransform(a = 0, b = 1/cofacteur, c = 0)
      
      if (avec_controles) {
        self$monomarques_trans <- lapply(self$tubes_monomarques, function(fcs) {
          canaux_presents <- intersect(canaux_fluo, flowCore::colnames(fcs))
          if (length(canaux_presents) > 0) {
            funs_locales <- lapply(seq_along(canaux_presents), function(x) arsinh_fun)
            local_trans  <- flowCore::transformList(
              from = canaux_presents,
              tfun = funs_locales,
              to   = canaux_presents
            )
            return(flowCore::transform(fcs, local_trans))
          } else {
            return(fcs)
          }
        })
        names(self$monomarques_trans) <- names(self$tubes_monomarques)
      } else {
        self$monomarques_trans <- NULL
      }
      
      # Met aussi à jour trans_list pour que controler_monomarques / visualiser_compensation puisse l'utiliser
      funs_globales    <- lapply(seq_along(canaux_fluo), function(x) arsinh_fun)
      self$trans_list  <- flowCore::transformList(
        from = canaux_fluo,
        tfun = funs_globales,
        to   = canaux_fluo
      )
    }, 
    
    definir_et_extraire = function(nom_canal, intervalle_gate_negatif, intervalle_gate_positif, utiliser_unstained = TRUE) {
      if (!nom_canal %in% self$canaux) stop("Le canal spécifié n'existe pas.")
      
      existe_unstained <- "TUBE_UNSTAINED" %in% names(self$monomarques_trans)
      nom_tube_neg <- if(utiliser_unstained && existe_unstained) "TUBE_UNSTAINED" else nom_canal
      self$source_neg_utilisee[[nom_canal]] <- nom_tube_neg
      
      source_trans_neg <- self$monomarques_trans[[nom_tube_neg]]
      source_brute_neg <- self$tubes_monomarques[[nom_tube_neg]]
      
      limites_negatif <- setNames(list(intervalle_gate_negatif), nom_canal)
      gate_negatif    <- flowCore::rectangleGate(filterId = paste0("Gate_Negatif_", nom_canal), .gate = limites_negatif)
      self$bornes_gates_neg[[nom_canal]] <- gate_negatif
      garde_evts_du_gate_negatif <- flowCore::filter(source_trans_neg, gate_negatif)
      self$gates_negatifs[[nom_canal]] <- source_brute_neg[garde_evts_du_gate_negatif@subSet, ]
      
      if (!is.null(intervalle_gate_positif)) {
        limites_positif <- setNames(list(intervalle_gate_positif), nom_canal)
        gate_positif    <- flowCore::rectangleGate(filterId = paste0("Gate_Positif_", nom_canal), .gate = limites_positif)
        self$bornes_gates_pos[[nom_canal]] <- gate_positif
        garde_evts_du_gate_positif <- flowCore::filter(self$monomarques_trans[[nom_canal]], gate_positif)
        self$gates_positifs[[nom_canal]]   <- self$tubes_monomarques[[nom_canal]][garde_evts_du_gate_positif@subSet, ]
      }
    },
    
    graphiques_gates = function(nom_canal = NULL, shiny_neg = NULL, shiny_pos = NULL, afficher_unstained_neg = TRUE) {
      canaux_a_generer <- if (is.null(nom_canal)) self$canaux else nom_canal
      
      plots <- lapply(canaux_a_generer, function(canal) {
        existe_unstained <- "TUBE_UNSTAINED" %in% names(self$monomarques_trans)
        tube_neg_a_tracer <- if(afficher_unstained_neg && existe_unstained) "TUBE_UNSTAINED" else canal
        
        fcs_brut  <- self$tubes_monomarques[[canal]]
        nom_bio   <- self$get_label(fcs_brut, canal) 
        
        df_trans_pos  <- as.data.frame(flowCore::exprs(self$monomarques_trans[[canal]]))
        df_trans_neg  <- as.data.frame(flowCore::exprs(self$monomarques_trans[[tube_neg_a_tracer]]))
        
        lim_n <- if (!is.null(shiny_neg)) shiny_neg else c(0, 2)
        lim_p <- if (!is.null(shiny_pos)) shiny_pos else c(4, 7)
        
        pct_n <- round(sum(df_trans_neg[[canal]] >= lim_n[1] & df_trans_neg[[canal]] <= lim_n[2], na.rm=TRUE) / nrow(df_trans_neg) * 100, 1)
        pct_p <- round(sum(df_trans_pos[[canal]] >= lim_p[1] & df_trans_pos[[canal]] <= lim_p[2], na.rm=TRUE) / nrow(df_trans_pos) * 100, 1)
        
        p <- ggplot() +
          geom_density(data = df_trans_pos, aes(x = .data[[canal]], y = after_stat(count)), fill = "#d90429", alpha = 0.4) +
          geom_vline(xintercept = lim_p, color = "#d90429", linetype = "solid", linewidth = 0.9) +
          annotate("text", x = mean(lim_p), y = Inf, label = paste0("Pos: ", pct_p, "%"), vjust = 1.5, color = "#d90429", fontface = "bold") +
          theme_bw() +
          labs(title = paste("Ajustement des Gates -", canal), x = nom_bio, y = "Nombre d'événements")
        
        if(tube_neg_a_tracer == "TUBE_UNSTAINED") {
          p <- p + geom_density(data = df_trans_neg, aes(x = .data[[canal]], y = after_stat(count)), fill = "#0077b6", alpha = 0.4) +
            geom_vline(xintercept = lim_n, color = "#0077b6", linetype = "dashed", linewidth = 0.9) +
            annotate("text", x = mean(lim_n), y = Inf, label = paste0("Unstained: ", pct_n, "%"), vjust = 3, color = "#0077b6", fontface = "bold") +
            labs(subtitle = "Bleu : Tube Unstained | Rouge : Tube Monomarqué")
        } else {
          p <- p + geom_vline(xintercept = lim_n, color = "#0077b6", linetype = "dashed", linewidth = 0.9) +
            annotate("text", x = mean(lim_n), y = Inf, label = paste0("Neg: ", pct_n, "%"), vjust = 3, color = "#0077b6", fontface = "bold") +
            labs(subtitle = "Distribution du tube monomarqué (Négatif interne)")
        }
        return(p)
      })
      
      names(plots) <- canaux_a_generer
      if (!is.null(nom_canal)) return(plots[[nom_canal]])
      return(plots)
    }, 
    
    calculer_spillover = function() { # Méthode calculant la matrice de spillover (compensation) à partir des médianes de fluorescence
      nombre_canaux <- length(self$canaux) # Évalue le nombre total de canaux d'acquisition actifs dans l'expérience
      S <- matrix(0, nrow = nombre_canaux, ncol = nombre_canaux, dimnames = list(self$canaux, self$canaux)) # Initialise une matrice carrée vide nommée avec les identifiants des canaux
      
      for (canal_principal in self$canaux) { # Boucle itérative sur chaque canal (chaque ligne de la future matrice)
        if (is.null(self$gates_positifs[[canal_principal]]) || is.null(self$gates_negatifs[[canal_principal]])) { # Vérifie que les populations cellulaires positives et négatives ont bien été isolées pour ce canal
          stop(paste("Gates manquants pour :", canal_principal)) # Bloque l'exécution et alerte l'utilisateur si une fenêtre de sélection est manquante
        } 
        exprs_pos <- flowCore::exprs(self$gates_positifs[[canal_principal]])[, self$canaux, drop = FALSE] # Extrait la table des intensités de fluorescence brute pour la population cellulaire positive
        exprs_neg <- flowCore::exprs(self$gates_negatifs[[canal_principal]])[, self$canaux, drop = FALSE] # Extrait la table des intensités de fluorescence brute pour la population cellulaire négative
        med_pos <- apply(exprs_pos, 2, median, na.rm = TRUE) # Calcule la valeur médiane de fluorescence de la population positive pour l'ensemble des canaux
        med_neg <- apply(exprs_neg, 2, median, na.rm = TRUE) # Calcule la valeur médiane de fluorescence de la population négative pour l'ensemble des canaux
        delta_signal <- med_pos - med_neg # Soustrait la médiane négative de la médiane positive pour isoler le signal spécifique (fluorescence nette)
        delta_signal <- pmax(delta_signal, 0)
        signal_principal <- delta_signal[canal_principal] # Isole l'intensité du signal net émis dans le canal primaire (le fluorophore correspondant au tube)
        if (is.na(signal_principal) || signal_principal <= 0) signal_principal <- 1e-5 # Sécurise le calcul en remplaçant un signal nul ou négatif par une valeur infime pour éviter une division par zéro
        S[canal_principal, ] <- delta_signal / signal_principal # Calcule le ratio de chevauchement spectral (spillover) pour tous les canaux par rapport au canal primaire
        S[canal_principal, canal_principal] <- 1 # Impose une valeur stricte de 1 (100%) sur la diagonale pour l'autofluorescence du canal primaire
      } 
      self$S_matrix <- S # Enregistre la matrice de spillover globale finalisée dans l'objet R6
      return(self$S_matrix) # Renvoie la matrice calculée pour pouvoir l'exploiter ou l'afficher dans Shiny
    }, 
    
    modifier_spillover = function(nom_echantillon, canal1, canal2, valeur) { # Méthode permettant d'ajuster manuellement un coefficient de compensation spécifique pour un échantillon donné
      if (canal1 == canal2) { # Si l'utilisateur tente de modifier le croisement d'un canal avec lui-même
        stop("Interdit de modifier la diagonale") # Bloque l'action car la diagonale théorique doit contractuellement rester égale à 1 (100%)
      } 
      if (is.null(self$S_matrices_par_echantillon)) { # Si la structure de stockage des matrices individualisées n'existe pas encore
        self$S_matrices_par_echantillon <- list() # Initialise une liste vide dédiée pour mémoriser les matrices par échantillon
      } 
      if (is.null(self$S_matrices_par_echantillon[[nom_echantillon]])) { # Si cet échantillon de patient précis ne possède pas encore sa propre matrice personnalisée
        self$S_matrices_par_echantillon[[nom_echantillon]] <- self$S_matrix # Duplique la matrice de spillover générale calculée pour lui servir de base de travail initiale
      } 
      
      self$S_matrices_par_echantillon[[nom_echantillon]][canal1, canal2] <- valeur # Remplace le coefficient ciblé (intersection canal donneur / canal receveur) par la nouvelle valeur saisie
      return(self$S_matrices_par_echantillon[[nom_echantillon]]) # Renvoie la matrice de compensation individualisée et mise à jour pour cet échantillon
    }, 
    
    comparer_medianes_spillover = function(nom_echantillon = NULL) { # Méthode calculant, canal par canal, les médianes brutes des populations positive et négative AVANT et APRÈS application de la matrice de compensation
      matrice_utilisee <- if (!is.null(nom_echantillon) && !is.null(self$S_matrices_par_echantillon[[nom_echantillon]])) { # Recherche si une matrice individualisée existe pour l'échantillon demandé
        self$S_matrices_par_echantillon[[nom_echantillon]] # Utilise la matrice personnalisée de cet échantillon si elle existe
      } else {
        self$S_matrix # Sinon utilise la matrice de spillover générale
      }
      if (is.null(matrice_utilisee)) stop("Veuillez d'abord calculer la matrice de spillover.") # Bloque le calcul si aucune matrice n'est disponible
      
      nombre_canaux <- length(self$canaux) # Évalue le nombre total de canaux d'acquisition actifs
      mat_pos_avant   <- matrix(NA, nrow = nombre_canaux, ncol = nombre_canaux, dimnames = list(self$canaux, self$canaux)) # Médiane de la population positive AVANT compensation
      mat_neg_avant   <- matrix(NA, nrow = nombre_canaux, ncol = nombre_canaux, dimnames = list(self$canaux, self$canaux)) # Médiane de la population négative AVANT compensation
      mat_pos_apres   <- matrix(NA, nrow = nombre_canaux, ncol = nombre_canaux, dimnames = list(self$canaux, self$canaux)) # Médiane de la population positive APRÈS compensation
      mat_neg_apres   <- matrix(NA, nrow = nombre_canaux, ncol = nombre_canaux, dimnames = list(self$canaux, self$canaux)) # Médiane de la population négative APRÈS compensation
      mat_ecart_avant <- matrix(NA, nrow = nombre_canaux, ncol = nombre_canaux, dimnames = list(self$canaux, self$canaux)) # Écart brut (médiane positive - médiane négative) AVANT compensation, pour vérifier l'égalité des médianes
      mat_ecart_apres <- matrix(NA, nrow = nombre_canaux, ncol = nombre_canaux, dimnames = list(self$canaux, self$canaux)) # Écart brut (médiane positive - médiane négative) APRÈS compensation, pour vérifier l'égalité des médianes
      
      for (canal_principal in self$canaux) { # Boucle sur chaque canal monomarqué disponible (une ligne de la matrice)
        if (is.null(self$gates_positifs[[canal_principal]]) || is.null(self$gates_negatifs[[canal_principal]])) { # Vérifie que les gates ont bien été validés pour ce canal
          stop(paste("Gates manquants pour :", canal_principal)) # Alerte si un gate est manquant
        }
        
        fcs_pos_brut <- self$gates_positifs[[canal_principal]] # Population positive brute (non compensée) pour ce canal
        fcs_neg_brut <- self$gates_negatifs[[canal_principal]] # Population négative brute (non compensée) pour ce canal
        fcs_pos_comp <- flowCore::compensate(fcs_pos_brut, matrice_utilisee) # Applique la matrice de compensation à la population positive
        fcs_neg_comp <- flowCore::compensate(fcs_neg_brut, matrice_utilisee) # Applique la matrice de compensation à la population négative
        
        calculer_medianes <- function(fcs_cible, fcs_ref_neg) { # Sous-fonction interne calculant les médianes brutes positive/négative, canal par canal
          mat_exprs_cible <- flowCore::exprs(fcs_cible)[, self$canaux, drop = FALSE] # Extrait la matrice des intensités de fluorescence de la population positive
          mat_exprs_neg   <- flowCore::exprs(fcs_ref_neg)[, self$canaux, drop = FALSE] # Extrait la matrice des intensités de fluorescence de la population négative
          med_pos <- apply(mat_exprs_cible, 2, median, na.rm = TRUE) # Médiane de la population positive, pour l'ensemble des canaux
          med_neg <- apply(mat_exprs_neg, 2, median, na.rm = TRUE) # Médiane de la population négative, pour l'ensemble des canaux
          list(pos = round(med_pos, 2), neg = round(med_neg, 2), ecart = round(med_pos - med_neg, 3)) # Renvoie les médianes brutes (non bornées) ainsi que leur écart, utile pour le test d'égalité
        }
        
        res_brut <- calculer_medianes(fcs_pos_brut, fcs_neg_brut) # Calcule les médianes brutes AVANT compensation
        res_comp <- calculer_medianes(fcs_pos_comp, fcs_neg_comp) # Calcule les médianes brutes APRÈS compensation
        
        mat_pos_avant[canal_principal, ] <- res_brut$pos # Remplit la ligne "médiane positive avant compensation" pour ce canal
        mat_neg_avant[canal_principal, ] <- res_brut$neg # Remplit la ligne "médiane négative avant compensation" pour ce canal
        mat_pos_apres[canal_principal, ] <- res_comp$pos # Remplit la ligne "médiane positive après compensation" pour ce canal
        mat_neg_apres[canal_principal, ] <- res_comp$neg # Remplit la ligne "médiane négative après compensation" pour ce canal
        
        mat_ecart_avant[canal_principal, ] <- res_brut$ecart # Écart brut médiane positive - médiane négative AVANT compensation (toutes unités confondues, positif ou négatif)
        mat_ecart_apres[canal_principal, ] <- res_comp$ecart # Écart brut médiane positive - médiane négative APRÈS compensation : une compensation correcte donne un écart proche de 0
        
        mat_ecart_avant[canal_principal, canal_principal] <- NA # La diagonale correspond au signal réel du marqueur (et non à un spillover résiduel) : elle n'est pas concernée par le test d'égalité des médianes
        mat_ecart_apres[canal_principal, canal_principal] <- NA # Idem après compensation : on neutralise la diagonale pour ne pas fausser la lecture du critère
      }
      
      resultat <- list(
        avant = list(pos = mat_pos_avant, neg = mat_neg_avant, ecart = mat_ecart_avant), # Médianes brutes positive/négative (+ écart) AVANT compensation
        apres = list(pos = mat_pos_apres, neg = mat_neg_apres, ecart = mat_ecart_apres)  # Médianes brutes positive/négative (+ écart) APRÈS compensation
      ) # Compile l'ensemble des résultats comparatifs
      
      cle <- if (!is.null(nom_echantillon)) nom_echantillon else "GLOBAL" # Clé de stockage : nom de l'échantillon ou "GLOBAL" si aucun échantillon n'est précisé
      self$comparaison_medianes[[cle]] <- resultat # Archive le résultat pour un accès ultérieur (ex: affichage Shiny)
      
      return(resultat) # Renvoie la liste des résultats comparatifs
    },
    
    controler_monomarques = function(fichier_monomarque, canal_x, canal_y, max_points = 200000) { # Méthode générant un graphique biplot comparatif avant/après compensation pour contrôler le spillover
      if (is.null(canal_x) || canal_x == "" || is.null(canal_y) || canal_y == "") return(NULL) # Intercepte et arrête la fonction si l'un des deux canaux d'acquisition n'est pas renseigné
      fcs_brut_original  <- self$tubes_monomarques[[fichier_monomarque]] # Récupère le fichier d'acquisition FCS brut d'origine correspondant au tube monomarqué
      fcs_compense_brut  <- flowCore::compensate(fcs_brut_original, self$S_matrix) # Applique la matrice de spillover calculée pour soustraire mathématiquement les fluorescences croisées
      nom_tube_neg <- self$source_neg_utilisee[[fichier_monomarque]] # Identifie le tube de référence négative qui a servi à calibrer ce canal
      fcs_trans_pour_pos <- self$monomarques_trans[[fichier_monomarque]] # Récupère les données transformées du tube positif pour réutiliser ses fenêtres de sélection
      fcs_trans_pour_neg <- self$monomarques_trans[[nom_tube_neg]] # Récupère les données transformées du tube négatif pour réutiliser ses fenêtres de sélection
      indices_neg <- flowCore::filter(fcs_trans_pour_neg, self$bornes_gates_neg[[fichier_monomarque]])@subSet # Extrait la liste des événements cellulaires inclus dans le gate négatif
      indices_pos <- flowCore::filter(fcs_trans_pour_pos, self$bornes_gates_pos[[fichier_monomarque]])@subSet # Extrait la liste des événements cellulaires inclus dans le gate positif
      
      calculer_spillover_pourcentage = function(fcs_cible, fcs_ref_neg = NULL) { # Sous-fonction interne calculant le pourcentage réel de débordement d'un canal sur un autre
        mat_exprs_cible <- flowCore::exprs(fcs_cible) # Extrait la matrice des intensités de fluorescence de l'échantillon cible analysé
        mat_exprs_neg   <- if(is.null(fcs_ref_neg)) mat_exprs_cible else flowCore::exprs(fcs_ref_neg) # Utilise la matrice de l'autocontrôle négatif approprié (interne ou externe)
        delta_observe <- median(mat_exprs_cible[indices_pos, canal_y], na.rm = TRUE) - median(mat_exprs_neg[indices_neg, canal_y], na.rm = TRUE) # Calcule le signal net reçu par le canal secondaire importun (canal Y)
        delta_source  <- median(mat_exprs_cible[indices_pos, canal_x], na.rm = TRUE) - median(mat_exprs_neg[indices_neg, canal_x], na.rm = TRUE) # Calcule le signal net reçu par le canal primaire légitime (canal X)
        delta_observe <- max(delta_observe, 0, na.rm = TRUE) # Empêche toute valeur négative : un spillover ne peut pas être inférieur à 0
        delta_source  <- max(delta_source, 0, na.rm = TRUE) # Idem pour le signal source, borné à 0
        if (is.na(delta_source) || delta_source == 0) return(0.00) # Évite une division par zéro si aucun signal n'est détecté dans le canal émetteur
        return(round((delta_observe / delta_source) * 100, 2)) # Renvoie le ratio de spillover exprimé en pourcentage, arrondi à deux décimales, toujours ≥ 0
      }
      
      tube_unstained_brut <- if(nom_tube_neg == "TUBE_UNSTAINED") self$tubes_monomarques[["TUBE_UNSTAINED"]] else NULL # Récupère le tube Unstained brut s'il s'agit de la référence négative globale
      tube_unstained_comp <- if(nom_tube_neg == "TUBE_UNSTAINED") flowCore::compensate(tube_unstained_brut, self$S_matrix) else NULL # Applique la compensation au tube Unstained brut pour le calcul résiduel
      valeur_spill_initial  <- calculer_spillover_pourcentage(fcs_brut_original, tube_unstained_brut) # Évalue le pourcentage de chevauchement spectral initial (avant compensation)
      valeur_spill_residuel <- calculer_spillover_pourcentage(fcs_compense_brut, tube_unstained_comp) # Évalue le pourcentage de chevauchement spectral résiduel (après application de la matrice)
      fcs_avant_trans    <- self$monomarques_trans[[fichier_monomarque]] # Récupère l'objet flowFrame des contrôles originaux déjà transformés
      fcs_apres_trans    <- flowCore::transform(fcs_compense_brut, self$trans_list) # Applique la transformation mathématique Arcsinh sur le fichier nouvellement compensé
      mat_avant <- flowCore::exprs(fcs_avant_trans)[, c(canal_x, canal_y), drop = FALSE] # Isoles l'intensité de fluorescence des deux canaux d'intérêt avant compensation
      mat_apres <- flowCore::exprs(fcs_apres_trans)[, c(canal_x, canal_y), drop = FALSE] # Isoles l'intensité de fluorescence des deux canaux d'intérêt après compensation
      nb_evenements <- nrow(mat_avant) # Compte le nombre total de cellules (événements) présentes dans le fichier FCS
      taille_echantillon <- min(nb_evenements, max_points) # Restreint le nombre d'événements à afficher pour optimiser la vitesse du rendu graphique
      indices_sub <- sample(seq_len(nb_evenements), taille_echantillon) # Tire au sort de manière aléatoire les indices des cellules à afficher
      df_avant <- as.data.frame(mat_avant[indices_sub, , drop = FALSE]) # Crée un tableau R contenant le sous-échantillon d'événements non compensés
      df_apres <- as.data.frame(mat_apres[indices_sub, , drop = FALSE]) # Crée un tableau R contenant le sous-échantillon d'événements compensés
      label_x_explicite <- self$get_label(fcs_brut_original, canal_x) # Génère le libellé biologique complet pour l'axe X (ex: "V3-A | CD4")
      label_y_explicite <- self$get_label(fcs_brut_original, canal_y) # Génère le libellé biologique complet pour l'axe Y (ex: "B1-A | CD8")
      limite_x <- range(c(df_avant[[canal_x]], df_apres[[canal_x]]), na.rm = TRUE) + c(-0.5, 0.5) # Calcule des limites d'affichage identiques en abscisse pour les deux graphiques
      limite_y <- range(c(df_avant[[canal_y]], df_apres[[canal_y]]), na.rm = TRUE) + c(-0.5, 0.5) # Calcule des limites d'affichage identiques en ordonnée pour les deux graphiques
      
      creer_affichage_comparatif = function(df_points, titre_plot, valeur_spill, lab_x, lab_y, cx, cy) { # Sous-fonction standardisant le style visuel des figures ggplot2
        # Densité calculée par binning raster (rapide, indépendant du nombre d'événements)
        # plutôt que par ggpointdensity (densité par point, coûteuse avec beaucoup d'événements).
        df_densite <- calculer_densite_raster(df_points[[cx]], df_points[[cy]], limite_x, limite_y)
        if (is.null(df_densite)) {
          return(ggplot() + theme_bw() +
                   labs(title = titre_plot, subtitle = "Pas assez d'événements pour tracer la densité", x = lab_x, y = lab_y))
        }
        ggplot(df_densite, aes(x = X, y = Y, fill = densite)) + # Initialise la figure biplot à partir de la grille de densité
          geom_raster(interpolate = TRUE) + # Affiche la densité sous forme d'image raster (rapide, même avec beaucoup d'événements)
          scale_fill_gradientn(colours = PALETTE_DENSITE) + # Applique une palette de couleurs pseudo-spectrale allant du bleu (faible densité) au rouge (forte densité)
          coord_cartesian(xlim = limite_x, ylim = limite_y) + # Verrouille les fenêtres d'affichage pour aligner visuellement les deux graphiques côte à côte
          theme_bw() + theme(legend.position = "none") + # Applique un arrière-plan blanc et masque la légende de l'échelle colorimétrique
          labs(title = titre_plot, subtitle = paste0("Spillover : ", valeur_spill, " %"), x = lab_x, y = lab_y) # Assigne le titre, la valeur calculée du spillover et les axes
      }
      
      plot_avant <- creer_affichage_comparatif(df_avant, "Avant Compensation", valeur_spill_initial, label_x_explicite, label_y_explicite, canal_x, canal_y) # Génère la figure représentative des données brutes croisées
      plot_apres <- creer_affichage_comparatif(df_apres, "Après Compensation", valeur_spill_residuel, label_x_explicite, label_y_explicite, canal_x, canal_y) # Génère la figure représentative des données nettoyées et corrigées
      return(gridExtra::grid.arrange(plot_avant, plot_apres, ncol = 2)) # Assemble et retourne les deux graphiques côte à côte au sein d'une seule figure combinée
    }, 
    
    compenser = function() {
      self$echantillons_traites <- lapply(names(self$echantillons), function(nom) {
        fcs <- self$echantillons[[nom]]
        
        # Sélectionner la matrice : spécifique si elle existe, sinon globale
        matrice_cible <- if (!is.null(self$S_matrices_par_echantillon[[nom]])) {
          self$S_matrices_par_echantillon[[nom]]
        } else {
          self$S_matrix
        }
        
        canaux_communs <- intersect(colnames(matrice_cible), flowCore::colnames(fcs))
        sub_matrix <- matrice_cible[canaux_communs, canaux_communs, drop = FALSE]
        
        return(flowCore::compensate(fcs, sub_matrix))
      })
      names(self$echantillons_traites) <- names(self$echantillons)
    },
    
    visualiser_compensation = function(nom_echantillon, canal_x, canal_y, affichage = "Both") {
      local_trans <- self$trans_list
      fcs_source  <- self$echantillons[[nom_echantillon]]
      
      # ───────────────────────────────────────────────
      # Cas "ALL" : toutes les combinaisons fluo
      # ───────────────────────────────────────────────
      if (canal_x == "ALL" || canal_y == "ALL") {
        # "Both" n'est pas autorisé en vue d'ensemble : on ne peut afficher ici que
        # l'ensemble des plots "avant" OU l'ensemble des plots "après", jamais les
        # deux combinés par paire (l'interface Shiny empêche déjà cette sélection,
        # ce garde-fou couvre tout appel direct à la méthode).
        if (identical(affichage, "Both")) {
          stop("Le mode 'Vue d'ensemble' ne permet pas l'affichage 'Both' : choisissez 'Before compensation only' ou 'After compensation only'.")
        }
        
        tous_canaux  <- flowCore::colnames(fcs_source)
        canaux_fluo  <- tous_canaux[!grepl("FSC|SSC|Time", tous_canaux, ignore.case = TRUE)]
        combinaisons <- expand.grid(x = canaux_fluo, y = canaux_fluo, stringsAsFactors = FALSE)
        # x < y : ne garde qu'une seule orientation par paire de canaux (évite le
        # doublon "A en X / B en Y" ET "B en X / A en Y", qui montrent la même
        # information biologique retournée)
        combinaisons <- combinaisons[combinaisons$x < combinaisons$y, ]
        
        plots_list <- lapply(seq_len(nrow(combinaisons)), function(i) {
          self$visualiser_compensation(
            nom_echantillon,
            combinaisons$x[i],
            combinaisons$y[i],
            affichage = affichage
          )
        })
        
        plots_list <- Filter(Negate(is.null), plots_list)
        if (length(plots_list) == 0) return(NULL)
        
        return(structure(plots_list, class = c("carrot_plots_list", "list")))
      }
      
      # ───────────────────────────────────────────────
      # Cas simple : une paire de canaux
      # ───────────────────────────────────────────────
      label_x <- self$get_label(fcs_source, canal_x)
      label_y <- self$get_label(fcs_source, canal_y)
      
      extraire_matrice <- function(fcs, cx, cy) {
        if (is.null(fcs)) return(NULL)
        if (!(cx %in% flowCore::colnames(fcs)) || !(cy %in% flowCore::colnames(fcs))) return(NULL)
        
        fcs_trans <- flowCore::transform(fcs, local_trans)
        mat       <- flowCore::exprs(fcs_trans)[, c(cx, cy), drop = FALSE]
        
        # → ICI : on garde TOUS les évènements
        as.data.frame(mat)
      }
      
      df_avant <- if (affichage %in% c("Both", "Before compensation only"))
        extraire_matrice(self$echantillons[[nom_echantillon]], canal_x, canal_y) else NULL
      
      df_apres <- if (affichage %in% c("Both", "After compensation only"))
        extraire_matrice(self$echantillons_traites[[nom_echantillon]], canal_x, canal_y) else NULL
      
      # Limites communes
      all_data <- rbind(df_avant, df_apres)
      if (!is.null(all_data) && nrow(all_data) > 0) {
        lim_x <- range(all_data[, 1], na.rm = TRUE)
        lim_y <- range(all_data[, 2], na.rm = TRUE)
      } else {
        lim_x <- c(0, 1)
        lim_y <- c(0, 1)
      }
      
      # ───────────────────────────────────────────────
      # Fonction interne : densité → raster → ggplot
      # ───────────────────────────────────────────────
      creer_plot_raster <- function(df, titre) {
        if (is.null(df) || nrow(df) == 0) return(NULL)
        colnames(df) <- c("X", "Y")
        
        # Résolution de la grille (pixels plus petits si ↑)
        resolution <- 400
        
        x_breaks <- seq(lim_x[1], lim_x[2], length.out = resolution + 1)
        y_breaks <- seq(lim_y[1], lim_y[2], length.out = resolution + 1)
        
        df_binned <- df |>
          dplyr::mutate(
            x_bin = cut(X, breaks = x_breaks, include.lowest = TRUE),
            y_bin = cut(Y, breaks = y_breaks, include.lowest = TRUE)
          ) |>
          dplyr::count(x_bin, y_bin, name = "densite") |>
          tidyr::drop_na()
        
        # centres des bins
        x_centers <- (head(x_breaks, -1) + tail(x_breaks, -1)) / 2
        y_centers <- (head(y_breaks, -1) + tail(y_breaks, -1)) / 2
        
        df_binned <- df_binned |>
          dplyr::mutate(
            X = x_centers[as.integer(x_bin)],
            Y = y_centers[as.integer(y_bin)]
          )
        
        ggplot2::ggplot(df_binned, ggplot2::aes(x = X, y = Y, fill = densite)) +
          ggplot2::geom_raster(interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(
            colours = PALETTE_DENSITE
          ) +
          ggplot2::coord_cartesian(xlim = lim_x, ylim = lim_y) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
          ggplot2::labs(title = titre, x = label_x, y = label_y)
      }
      
      # Création des deux plots
      plot_avant <- creer_plot_raster(df_avant, paste(nom_echantillon, "- Avant"))
      plot_apres <- creer_plot_raster(df_apres, paste(nom_echantillon, "- Après"))
      
      # Rendu final
      if (affichage == "Both") {
        plots_valides <- Filter(Negate(is.null), list(plot_avant, plot_apres))
        if (length(plots_valides) == 0) return(NULL)
        if (length(plots_valides) == 1) return(plots_valides[[1]])
        return(gridExtra::arrangeGrob(grobs = plots_valides, ncol = 2))
      }
      
      if (affichage == "Before compensation only") return(plot_avant)
      if (affichage == "After compensation only")  return(plot_apres)
      
      return(NULL)
    },
    
    exporter_fcs_compenses = function(noms_echantillons = "all", dossier_export = ".") { # Méthode permettant d'écrire et de sauvegarder les fichiers FCS compensés sur le disque dur
      if (length(self$echantillons_traites) == 0) { # Vérifie si la liste contenant les échantillons corrigés est vide
        stop("Aucun échantillon compensé n'est disponible dans l'objet R6.") # Bloque l'exécution si aucun traitement de compensation n'a été préalablement réalisé
      } 
      if (length(noms_echantillons) == 1 && noms_echantillons == "all") { # Si l'utilisateur demande explicitement d'exporter la totalité de la cohorte
        tubes_a_exporter <- names(self$echantillons_traites) # Sélectionne l'intégralité des identifiants des patients stockés en mémoire
      } else { # Sinon, si une liste restreinte de noms de patients spécifiques a été fournie
        tubes_a_exporter <- intersect(noms_echantillons, names(self$echantillons_traites)) # Identifie par intersection les échantillons demandés qui existent réellement en mémoire
        if (length(tubes_a_exporter) == 0) { # Si aucun des noms de la liste fournie par l'utilisateur ne correspond aux fichiers disponibles
          stop("Aucun des échantillons spécifiés n'a été trouvé dans les données compensées.") # Interrompt le processus et renvoie un message d'erreur explicite
        } 
      } 
      if (!dir.exists(dossier_export)) { # Vérifie si le dossier de destination spécifié n'existe pas encore physiquement sur le disque
        dir.create(dossier_export, recursive = TRUE) # Crée automatiquement l'arborescence des dossiers manquants pour éviter une erreur d'écriture
      } 
      
      for (nom in tubes_a_exporter) { # Boucle de traitement itérative pour exporter chaque fichier sélectionné
        fcs_obj <- self$echantillons_traites[[nom]] # Extrait l'objet flowFrame de l'échantillon compensé courant
        if (is.null(fcs_obj)) next # Sécurité : passe immédiatement au fichier suivant si l'objet en mémoire est anormalement vide
        
        if (!is.null(self$S_matrix)) { # Si une matrice de compensation valide est configurée dans l'objet R6
          flowCore::keyword(fcs_obj)[["$SPILLOVER"]] <- self$S_matrix # Intègre de manière standardisée la matrice de spillover au sein des métadonnées (mots-clés) du fichier FCS
        } 
        nom_fichier_propre <- paste0(gsub("[^a-zA-Z0-9_]", "_", nom), "_compense.fcs") # Nettoie le nom de l'échantillon en remplaçant les caractères spéciaux par des tirets bas pour la compatibilité système
        nom_fcs <- file.path(dossier_export, nom_fichier_propre) # Concatène le chemin du dossier et le nom du fichier pour obtenir l'adresse d'écriture finale
        flowCore::write.FCS(fcs_obj, filename = nom_fcs) # Enregistre physiquement le nouvel objet flowFrame au format binaire FCS 3.0 ou 3.1 standard sur le disque dur
      } 
      message(paste(length(tubes_a_exporter), "fichier(s) FCS compensé(s) exporté(s) dans :", dossier_export)) # Affiche un message de confirmation récapitulatif dans la console de commande
    },
    
    exporter_fcs_qc = function(noms_echantillons = "all", sources = c("peacoqc", "flowai"), dossier_export = ".") { # Méthode permettant d'écrire sur le disque les fichiers FCS ayant subi le contrôle qualité PeacoQC et/ou flowAI
      sources <- intersect(sources, c("peacoqc", "flowai")) # Filtre les sources demandées pour ne garder que les valeurs reconnues ("peacoqc" et/ou "flowai")
      if (length(sources) == 0) {
        stop("Aucune source de contrôle qualité valide sélectionnée (attendu : 'peacoqc' et/ou 'flowai').")
      }
      
      liste_sources <- list(peacoqc = self$post_PeacoQC, flowai = self$post_flowAI) # Regroupe les deux listes de résultats QC disponibles pour un accès homogène
      suffixes      <- list(peacoqc = "_PeacoQC", flowai = "_flowAI") # Associe à chaque source le suffixe de nommage utilisé pour distinguer les fichiers exportés
      
      noms_disponibles <- unique(unlist(lapply(sources, function(s) names(liste_sources[[s]])))) # Recense l'ensemble des échantillons disponibles, toutes sources sélectionnées confondues
      if (length(noms_disponibles) == 0) {
        stop("Aucun résultat PeacoQC ou flowAI disponible : exécutez d'abord le contrôle qualité correspondant.")
      }
      
      tubes_a_exporter <- if (length(noms_echantillons) == 1 && noms_echantillons == "all") { # Si l'utilisateur souhaite exporter la totalité des échantillons traités
        noms_disponibles # Sélectionne l'intégralité des échantillons disponibles pour les sources demandées
      } else { # Sinon, si une liste restreinte de noms a été fournie par l'utilisateur
        intersect(noms_echantillons, noms_disponibles) # Identifie par intersection les échantillons demandés qui existent réellement en mémoire
      }
      if (length(tubes_a_exporter) == 0) {
        stop("Aucun des échantillons spécifiés n'a été trouvé dans les résultats PeacoQC/flowAI.")
      }
      
      if (!dir.exists(dossier_export)) { # Vérifie si le dossier de destination spécifié n'existe pas encore physiquement sur le disque
        dir.create(dossier_export, recursive = TRUE) # Crée automatiquement l'arborescence des dossiers manquants pour éviter une erreur d'écriture
      }
      
      fichiers_ecrits <- character(0) # Accumule au fil de la boucle les chemins des fichiers FCS effectivement écrits sur le disque
      
      for (nom in tubes_a_exporter) { # Boucle de traitement itérative pour chaque échantillon sélectionné
        for (source in sources) { # Boucle interne pour chaque source de contrôle qualité demandée (PeacoQC et/ou flowAI)
          fcs_obj <- liste_sources[[source]][[nom]] # Extrait l'objet flowFrame nettoyé correspondant à cet échantillon et cette source
          if (is.null(fcs_obj)) next # Sécurité : passe au suivant si cet échantillon n'a pas de résultat pour cette source précise
          
          nom_fichier_propre <- paste0(gsub("[^a-zA-Z0-9_]", "_", nom), suffixes[[source]], ".fcs") # Construit un nom de fichier sûr, distinguant la source d'origine du nettoyage
          chemin_fcs <- file.path(dossier_export, nom_fichier_propre) # Concatène le chemin du dossier et le nom du fichier pour obtenir l'adresse d'écriture finale
          flowCore::write.FCS(fcs_obj, filename = chemin_fcs) # Enregistre physiquement l'objet flowFrame nettoyé au format binaire FCS standard sur le disque dur
          fichiers_ecrits <- c(fichiers_ecrits, chemin_fcs) # Ajoute le chemin du fichier fraîchement écrit à la liste de suivi
        }
      }
      
      message(paste(length(fichiers_ecrits), "fichier(s) FCS post-QC exporté(s) dans :", dossier_export)) # Affiche un message de confirmation récapitulatif dans la console de commande
      invisible(fichiers_ecrits) # Renvoie de manière invisible la liste des chemins écrits (utile pour zipper ensuite depuis Shiny)
    },
    
    generer_pdf_resume_qc = function(chemin_pdf, noms_echantillons = NULL, sources = c("peacoqc", "flowai")) { # Méthode générant un PDF multi-pages regroupant toutes les figures de contrôle qualité et les paramètres utilisés
      sources <- intersect(sources, c("peacoqc", "flowai")) # Filtre les sources demandées pour ne garder que les valeurs reconnues
      if (length(sources) == 0) {
        stop("Aucune source de contrôle qualité valide sélectionnée (attendu : 'peacoqc' et/ou 'flowai').")
      }
      
      if (is.null(noms_echantillons) || (length(noms_echantillons) == 1 && noms_echantillons == "all")) { # Si aucune sélection précise n'est fournie
        noms_echantillons <- unique(c(names(self$post_PeacoQC), names(self$post_flowAI))) # Reprend l'ensemble des échantillons ayant un résultat PeacoQC et/ou flowAI
      }
      if (length(noms_echantillons) == 0) {
        stop("Aucun échantillon PeacoQC ou flowAI disponible pour générer le résumé.")
      }
      
      formater_parametres <- function(titre, params) { # Fonction utilitaire transformant une liste de réglages nommés en un bloc de texte lisible pour la page de garde des paramètres
        if (is.null(params) || length(params) == 0) return(paste0(titre, "\n(paramètres non disponibles)")) # Sécurité si aucun réglage n'a encore été archivé
        lignes <- vapply(names(params), function(n) { # Parcourt chaque réglage pour produire une ligne "nom : valeur"
          paste0("  - ", n, " : ", paste(params[[n]], collapse = ", ")) # Concatène le nom du paramètre et sa valeur (gère aussi les vecteurs comme ChExcludeFS)
        }, character(1))
        paste(c(titre, lignes), collapse = "\n") # Assemble le titre et l'ensemble des lignes de paramètres en un seul bloc de texte
      }
      
      grDevices::pdf(chemin_pdf, width = 11, height = 8.5) # Ouvre un nouveau périphérique graphique PDF multi-pages à l'emplacement demandé
      on.exit(grDevices::dev.off(), add = TRUE) # Garantit la fermeture propre du périphérique PDF même en cas d'erreur en cours de génération
      
      # Page de garde récapitulant les paramètres utilisés pour chaque méthode sélectionnée
      grid::grid.newpage() # Initialise une page blanche pour la synthèse textuelle des réglages
      blocs_params <- character(0) # Accumule les blocs de texte des paramètres pour chaque source demandée
      if ("peacoqc" %in% sources) {
        blocs_params <- c(blocs_params, formater_parametres("Paramètres PeacoQC utilisés :", self$parametres_peacoqc_utilises))
      }
      if ("flowai" %in% sources) {
        blocs_params <- c(blocs_params, formater_parametres("Paramètres flowAI utilisés :", self$parametres_flowai_utilises))
      }
      texte_page_garde <- paste(c(
        paste0("Résumé du contrôle qualité — ", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
        paste0("Échantillons inclus (", length(noms_echantillons), ") : ", paste(noms_echantillons, collapse = ", ")), "",
        blocs_params
      ), collapse = "\n") # Assemble titre, liste des échantillons et blocs de paramètres en un seul texte multi-lignes
      grid::grid.text(texte_page_garde, x = 0.05, y = 0.95, just = c("left", "top"), # Affiche le texte de synthèse en haut à gauche de la page PDF
                      gp = grid::gpar(fontsize = 10, fontfamily = "mono")) # Utilise une police à chasse fixe pour un alignement propre des paramètres
      
      # Une page par échantillon et par source sélectionnée, contenant la figure de diagnostic correspondante
      for (nom in noms_echantillons) { # Boucle sur chaque échantillon à documenter dans le résumé
        if ("peacoqc" %in% sources && !is.null(self$post_PeacoQC[[nom]])) { # Si un résultat PeacoQC existe pour cet échantillon
          graphique <- self$plots_peacoqc[[nom]] # Récupère la figure déjà archivée si l'utilisateur l'a préalablement visualisée dans l'onglet PeacoQC
          if (is.null(graphique)) { # Si la figure n'a encore jamais été générée pour cet échantillon
            graphique <- tryCatch(self$visualiser_peacoqc(nom_echantillon = nom), error = function(e) NULL) # Génère la figure à la volée en cas d'absence
          }
          if (!is.null(graphique)) print(graphique) # Insère la figure PeacoQC dans une nouvelle page du PDF
          
          chemin_png <- self$plots_peacoqc_natif[[nom]] # Récupère le chemin du rapport diagnostique natif PNG généré par PeacoQC::PeacoQC
          if (!is.null(chemin_png) && file.exists(chemin_png)) { # Si ce PNG natif est disponible sur le disque
            image_native <- png::readPNG(chemin_png) # Charge le PNG natif en mémoire sous forme de matrice de pixels
            grid::grid.newpage() # Ouvre une nouvelle page dédiée au rapport natif
            grid::grid.raster(image_native) # Dessine l'image du rapport natif PeacoQC sur toute la page
          }
        }
        
        if ("flowai" %in% sources && !is.null(self$post_flowAI[[nom]])) { # Si un résultat flowAI existe pour cet échantillon
          graphique <- self$plots_flowai[[nom]] # Récupère la figure déjà archivée si l'utilisateur l'a préalablement visualisée dans l'onglet flowAI
          if (is.null(graphique)) { # Si la figure n'a encore jamais été générée pour cet échantillon
            graphique <- tryCatch(self$visualiser_flowai(nom_echantillon = nom), error = function(e) NULL) # Génère la figure à la volée en cas d'absence
          }
          if (!is.null(graphique)) print(graphique) # Insère la figure flowAI dans une nouvelle page du PDF
        }
      }
      
      invisible(chemin_pdf) # Renvoie de manière invisible le chemin du PDF généré, prêt à être proposé au téléchargement depuis Shiny
    },
    
    sauvegarder_session_rds = function(nom_fichier = "Compensation_Session_Complete.rds") { # Méthode permettant de sauvegarder l'intégralité de l'état et des résultats de la session d'analyse dans un fichier binaire R (.rds)
      sauvegarde <- list( # Initialise une structure de liste imbriquée pour regrouper de manière organisée tous les éléments à archiver
        meta = list( # Sous-liste dédiée aux informations générales et de traçabilité de l'expérience
          date_export   = Sys.time(), # Enregistre l'horodatage exact (date et heure) de la création de la sauvegarde
          canaux_utilises = self$canaux # Mémorise la liste des canaux d'acquisition actifs durant cette session d'analyse
        ), 
        configuration_technique = list( # Sous-liste dédiée aux paramètres et opérateurs mathématiques appliqués aux données
          trans_list        = self$trans_list, # Archive l'opérateur global contenant les fonctions de transformation Arcsinh des canaux
          matrice_spillover = self$S_matrix # Archive la matrice de spillover (coefficients de compensation) calculée ou ajustée
        ), 
        gating = list( # Sous-liste dédiée aux populations cellulaires filtrées et extraites lors de l'analyse
          gates_positifs = self$gates_positifs, # Archive les structures flowFrame contenant les événements des populations cellulaires positives isolées
          gates_negatifs = self$gates_negatifs # Archive les structures flowFrame contenant les événements des populations cellulaires négatives isolées
        ), 
        visualisations = list( # Sous-liste dédiée à l'archivage des rendus graphiques produits pour l'assurance qualité
          plots_gates        = self$plots_gates, # Archive l'historique des figures de densité affichant le positionnement des fenêtres de sélection
          plots_compensation = self$plots_compensation # Archive l'historique des figures biplots comparatifs avant/après compensation générés pour les patients
        ) 
      ) 
      
      saveRDS(sauvegarde, file = nom_fichier) # Sérialise et enregistre l'objet liste complet sous forme de fichier binaire compressé (.rds) sur le stockage local
      message(paste("Session et paramètres sauvegardés avec succès dans :", nom_fichier)) # Génère un message de confirmation explicite au sein de la console R de commande
    }, 
    
    # =======================
    #   SECTION UNMIXING 
    # =======================
    
    # type cytometre = a8, s8, a5se, aurora, id7000, mosaic, opteon, xenith, minimal, auroraNL
    
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
      dossier_gates <- file.path(self$dossier_racine, "figure_gate")
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
    }, 
    
    
    # ============================================================
    #       ️ SECTION PRÉ-TRAITEMENT
    # ============================================================
    
    get_derniere_source = function() {
      # 1. Étapes de transformation et gating final
      if (!is.null(self$post_transformation) && length(self$post_transformation) > 0) return(self$post_transformation)
      if (!is.null(self$post_viabilite) && length(self$post_viabilite) > 0) return(self$post_viabilite)
      
      # 2. Étapes de nettoyage des doublets et débris
      if (!is.null(self$post_doublets_final) && length(self$post_doublets_final) > 0) return(self$post_doublets_final)
      if (!is.null(self$post_debris) && length(self$post_debris) > 0) return(self$post_debris)
      
      # 3. Étapes de contrôle qualité (QC)
      if (!is.null(self$post_retrait_bordures) && length(self$post_retrait_bordures) > 0) return(self$post_retrait_bordures)
      if (!is.null(self$post_PeacoQC) && length(self$post_PeacoQC) > 0) return(self$post_PeacoQC)
      if (!is.null(self$post_flowAI) && length(self$post_flowAI) > 0) return(self$post_flowAI)
      
      # 4. Sources initiales (fallback)
      if (!is.null(self$echantillons_traites) && length(self$echantillons_traites) > 0) return(self$echantillons_traites)
      
      return(self$echantillons) # Données brutes chargées
    },
    
    appliquer_peacoqc = function(dossier_rapports = NULL, reglages_specifiques = list()) {
      
      # Paramètres par défaut
      parametres_par_defaut <- list(
        determine_good_cells = "all",
        min_cells = 150,
        max_bins = 100,
        step = 500,
        MAD = 6,
        IT_limit = 0.6,
        consecutive_bins = 5,
        remove_zeros = FALSE,
        force_IT = 150,
        peak_removal = 1/3,
        min_nr_bins_peakdetection = 10
      )
      
      # Fusion des paramètres
      config_qc <- utils::modifyList(parametres_par_defaut, reglages_specifiques)
      
      # Sorties
      self$post_PeacoQC <- list()
      self$plots_peacoqc <- list()
      self$plots_peacoqc_natif <- list()
      
      # Boucle échantillons
      for (nom in names(self$echantillons_traites)) {
        
        message("Exécution PeacoQC sur : ", nom)
        
        ff_actuel <- self$echantillons_traites[[nom]]
        
        # Détection du canal temporel
        vrai_canal_temps <- grep("time", colnames(ff_actuel), value = TRUE, ignore.case = TRUE)[1]
        if (is.na(vrai_canal_temps)) vrai_canal_temps <- "Time"
        
        # Détermination des canaux à analyser (doc officielle)
        tous_canaux <- colnames(ff_actuel)
        canaux_fluo <- tous_canaux[!grepl("fsc|ssc|time", tous_canaux, ignore.case = TRUE)]
        
        # Dossier temporaire pour les PNG
        dossier_temp_plot <- file.path(
          tempdir(),
          paste0("PeacoQC_", gsub("[^a-zA-Z0-9_]", "_", nom), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        )
        dir.create(dossier_temp_plot, recursive = TRUE, showWarnings = FALSE)
        
        # Exécution PeacoQC
        res <- PeacoQC::PeacoQC(
          ff                     = ff_actuel,
          channels               = canaux_fluo,
          determine_good_cells   = config_qc$determine_good_cells,
          transform              = FALSE,
          time_channel_parameter = vrai_canal_temps,
          save_fcs               = FALSE,
          min_cells              = config_qc$min_cells,
          max_bins               = config_qc$max_bins,
          step                   = config_qc$step,
          MAD                    = config_qc$MAD,
          IT_limit               = config_qc$IT_limit,
          consecutive_bins       = config_qc$consecutive_bins,
          remove_zeros           = config_qc$remove_zeros,
          force_IT               = config_qc$force_IT,
          peak_removal           = config_qc$peak_removal,
          min_nr_bins_peakdetection = config_qc$min_nr_bins_peakdetection,
          plot                   = TRUE,
          output_directory       = dossier_temp_plot,
          report                 = FALSE
        )
        
        # Stockage du flowFrame filtré
        self$post_PeacoQC[[nom]] <- res$FinalFF
        
        # Stockage du plot ggplot
        self$plots_peacoqc[[nom]] <- res$Plot
        
        # Stockage du PNG natif
        # ✓ CORRIGÉ : PeacoQC::PeacoQC() n'écrit PAS le PNG directement dans "output_directory".
        # En interne, il crée un sous-dossier "output_directory/<name_directory>/" (name_directory
        # vaut "PeacoQC_results" par défaut, non redéfini ici) et c'est LÀ qu'atterrit le PNG natif.
        # Un list.files() non récursif sur dossier_temp_plot ne le trouvait donc jamais, laissant
        # plots_peacoqc_natif toujours à NULL et empêchant tout affichage côté interface.
        pngs <- list.files(dossier_temp_plot, pattern = "\\.png$", full.names = TRUE, recursive = TRUE)
        self$plots_peacoqc_natif[[nom]] <- if (length(pngs) > 0) pngs[1] else NULL # Un seul PNG "montage" est généré par échantillon : on n'en garde qu'un chemin scalaire
      }
      
      # Archive des paramètres
      self$parametres_peacoqc_utilises <- config_qc
      
      # Mise à jour pipeline
      if (!is.null(self$update_pipeline)) self$update_pipeline("PeacoQC")
    },
    
    appliquer_flowai = function(reglages_specifiques = list()) {
      
      if (length(self$echantillons_traites) == 0) {
        stop("Aucun échantillon traité disponible pour flowAI.")
      }
      
      # Paramètres sécurisés
      remove_from       <- reglages_specifiques$remove_from       %||% "all"
      timeCh            <- reglages_specifiques$timeCh            %||% NULL
      second_fractionFR <- reglages_specifiques$second_fractionFR %||% 0.1
      alphaFR           <- reglages_specifiques$alphaFR           %||% 0.01
      decompFR          <- if (isTRUE(reglages_specifiques$decompFR)) "cffilter" else "loess"
      
      ChExcludeFS <- reglages_specifiques$ChExcludeFS
      if (is.null(ChExcludeFS) || length(ChExcludeFS) == 0) ChExcludeFS <- NULL
      
      outlier_binsFS <- reglages_specifiques$outlier_binsFS %||% FALSE
      pen_valueFS    <- reglages_specifiques$pen_valueFS    %||% 500
      max_cptFS      <- reglages_specifiques$max_cptFS      %||% 3
      
      ChExcludeFM <- reglages_specifiques$ChExcludeFM
      if (is.null(ChExcludeFM) || length(ChExcludeFM) == 0) ChExcludeFM <- NULL
      
      sideFM       <- reglages_specifiques$sideFM       %||% "both"
      neg_valuesFM <- reglages_specifiques$neg_valuesFM %||% 1
      
      # Sorties
      self$post_flowAI <- list()
      self$rapports_flowai <- list()
      
      # Boucle échantillons
      for (nom in names(self$echantillons_traites)) {
        
        fcs_obj <- self$echantillons_traites[[nom]]
        
        # Dossier unique par échantillon
        dossier_tmp <- file.path(tempdir(), paste0("flowAI_", nom, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
        dir.create(dossier_tmp, recursive = TRUE, showWarnings = FALSE)
        
        res <- tryCatch({
          
          flowAI::flow_auto_qc(
            fcsfiles            = fcs_obj,
            remove_from         = remove_from,
            output              = 1,
            timeCh              = timeCh,
            second_fractionFR   = second_fractionFR,
            alphaFR             = alphaFR,
            decompFR            = decompFR,
            ChExcludeFS         = ChExcludeFS,
            outlier_binsFS      = outlier_binsFS,
            pen_valueFS         = pen_valueFS,
            max_cptFS           = max_cptFS,
            ChExcludeFM         = ChExcludeFM,
            sideFM              = sideFM,
            neg_valuesFM        = neg_valuesFM,
            html_report         = "_QC",
            mini_report         = "_QCmini",
            fcs_QC              = "_QC",
            fcs_highQ           = FALSE,
            fcs_lowQ            = FALSE,
            folder_results      = dossier_tmp
          )
          
        }, error = function(e) {
          stop(paste("Erreur flowAI pour l'échantillon", nom, ":", conditionMessage(e)))
        })
        
        # res est un flowSet → on extrait le flowFrame
        if (inherits(res, "flowSet")) {
          res <- res[[1]]
        }
        
        self$post_flowAI[[nom]] <- res
        self$rapports_flowai[[nom]] <- dossier_tmp
      }
      
      self$update_pipeline("flowAI")
      invisible(TRUE)
    },
    
    retirer_les_bordures = function(canal1, canal2, nom_echantillon = NULL) {
      self$canaux_bordures <- c(canal1, canal2)
      liste_source <- self$get_derniere_source()
      noms_a_traiter <- if (is.null(nom_echantillon)) names(liste_source) else nom_echantillon
      
      if (is.null(noms_a_traiter) || length(noms_a_traiter) == 0) {
        warning("Aucun échantillon trouvé à traiter pour le retrait des bordures.")
        return(NULL)
      }
      
      for (nom in noms_a_traiter) {
        if (!is.null(liste_source[[nom]])) {
          message("Retrait des bordures (Margins) sur : ", nom)
          self$post_retrait_bordures[[nom]] <- PeacoQC::RemoveMargins(
            ff       = liste_source[[nom]], 
            channels = self$canaux_bordures
          )
        }
      }
      if (!is.null(self$update_pipeline)) self$update_pipeline("bordures", nom_echantillon)
    },
    
    retirer_les_debris = function(matrice_points, canal_x, canal_y, nom_echantillon = NULL, source_nettoyage = "brutes") {
      if (is.null(matrice_points)) stop("Aucune coordonnée de gate fournie.")
      if (!is.matrix(matrice_points) && !is.data.frame(matrice_points)) {
        stop("La structure de la gate doit être une matrice ou un data.frame.")
      }
      if (nrow(matrice_points) < 3) stop("Un polygone de filtrage nécessite au moins 3 points.")
      
      liste_source <- NULL
      if (source_nettoyage == "peacoqc") {
        liste_source <- self$post_PeacoQC
      } else if (source_nettoyage == "flowai") {
        liste_source <- self$post_flowAI
      }
      
      if (is.null(liste_source) || length(liste_source) == 0) {
        if (source_nettoyage %in% c("peacoqc", "flowai")) {
          message("⚠️ [Warning] Source '", source_nettoyage, "' introuvable ou vide. Repli sur les données compensées brutes.")
        }
        liste_source <- self$echantillons_traites
      }
      
      if (is.null(liste_source) || length(liste_source) == 0) {
        stop("Aucune donnée disponible (brute, PeacoQC ou flowAI) pour appliquer le filtre débris.")
      }
      
      matrice_points <- as.matrix(matrice_points[, 1:2])
      colnames(matrice_points) <- c(canal_x, canal_y)
      polygone_debris <- flowCore::polygonGate(.gate = matrice_points, filterId = "Gate_Debris")
      if (is.null(self$gate_debris)) self$gate_debris <- list()
      if (is.null(self$post_debris)) self$post_debris <- list()
      
      appliquer_le_filtrage = function(nom) {
        flowframe_entree <- liste_source[[nom]]
        if (is.null(flowframe_entree)) return(NULL)
        
        message("Application du filtre Débris (PolygonGate) via source '", source_nettoyage, "' sur : ", nom)
        resultat_filtre <- flowCore::filter(flowframe_entree, polygone_debris)
        self$gate_debris[[nom]] <- polygone_debris
        self$post_debris[[nom]] <- flowframe_entree[resultat_filtre@subSet, ]
      }
      
      noms_a_traiter <- if (is.null(nom_echantillon)) names(liste_source) else nom_echantillon
      for (nom in noms_a_traiter) { 
        appliquer_le_filtrage(nom) 
      }
      
      if (!is.null(self$update_pipeline)) {
        self$update_pipeline("debris", nom_echantillon)
      }
      
      return(invisible(self))
    },
    
    # Applique un gate polygonal nommé séquentiellement sur tous les échantillons.
    # Chaque gate est empilé sur le précédent (le résultat du gate N-1 alimente le gate N).
    # L'entrée du 1er gate est déterminée par source_nettoyage.
    # gates_history stocke, pour chaque nom de gate, les résultats par échantillon.
    appliquer_gate_nomme = function(nom_gate, matrice_points, canal_x, canal_y,
                                    source_nettoyage = "brutes") {
      if (is.null(nom_gate) || nchar(trimws(nom_gate)) == 0) stop("Le gate doit avoir un nom.")
      if (is.null(matrice_points) || nrow(matrice_points) < 3) stop("Polygone invalide (< 3 points).")
      
      # Source d'entrée : dernier gate validé si disponible, sinon source QC choisie
      if (length(self$gates_history) > 0) {
        dernier_gate <- self$gates_history[[length(self$gates_history)]]
        liste_source <- lapply(names(dernier_gate), function(n) dernier_gate[[n]]$post_data)
      } else {
        if (source_nettoyage == "peacoqc" && length(self$post_PeacoQC) > 0) {
          liste_source <- self$post_PeacoQC
        } else if (source_nettoyage == "flowai" && length(self$post_flowAI) > 0) {
          liste_source <- self$post_flowAI
        } else {
          liste_source <- self$get_derniere_source()
        }
      }
      
      if (is.null(liste_source) || length(liste_source) == 0) {
        stop("Aucune donnée disponible pour appliquer le gate.")
      }
      
      mat_gate <- as.matrix(matrice_points[, 1:2])
      colnames(mat_gate) <- c(canal_x, canal_y)
      polygone <- flowCore::polygonGate(.gate = mat_gate, filterId = nom_gate)
      
      resultats_gate <- list()
      for (nom in names(liste_source)) {
        ff_entree <- liste_source[[nom]]
        if (is.null(ff_entree)) next
        n_avant <- nrow(flowCore::exprs(ff_entree))
        message("Application gate '", nom_gate, "' sur : ", nom)
        res_filtre <- flowCore::filter(ff_entree, polygone)
        ff_apres   <- ff_entree[res_filtre@subSet, ]
        n_apres    <- nrow(flowCore::exprs(ff_apres))
        resultats_gate[[nom]] <- list(
          polygone  = mat_gate,
          canal_x   = canal_x,
          canal_y   = canal_y,
          post_data = ff_apres,
          n_avant   = n_avant,
          n_apres   = n_apres
        )
      }
      
      # Enregistre dans l'historique (écrase si même nom)
      self$gates_history[[nom_gate]] <- resultats_gate
      
      # Met aussi à jour post_debris pour compatibilité avec visualiser_debris et les étapes aval
      self$gate_debris <- lapply(resultats_gate, function(r) {
        flowCore::polygonGate(.gate = r$polygone, filterId = nom_gate)
      })
      self$post_debris <- lapply(resultats_gate, function(r) r$post_data)
      
      return(invisible(self))
    },
    
    retirer_doublets_FSC = function(facteur_sensibilite = 4, axe_discrimination = "H_A", nom_echantillon = NULL) { # Méthode statistique éliminant les agrégats cellulaires (doublets) sur l'axe Forward Scatter (FSC) en évaluant la linéarité géométrique des signaux
      liste_fcs_source <- if (!is.null(self$post_debris) && length(self$post_debris) > 0) self$post_debris else self$get_derniere_source() # Sélectionne par défaut les données issues du filtre débris ou applique le mécanisme de repli pyramidal sur le dernier état valide
      label_source <- if (!is.null(self$post_debris) && length(self$post_debris) > 0) "post_debris" else "source_parente" # Génère une étiquette textuelle décrivant le niveau d'origine des données pour la traçabilité du pipeline
      
      if (is.null(self$gate_doublets_FSC)) self$gate_doublets_FSC <- list()
      if (is.null(self$post_doublets_FSC)) self$post_doublets_FSC <- list()
      if (is.null(self$post_doublets_final)) self$post_doublets_final <- list()
      
      calculer_doublets = function(nom, lbl_src) { # Sous-fonction interne calculant de manière isolée les limites statistiques d'exclusion des agrégats pour un échantillon donné
        flowframe_entree <- liste_fcs_source[[nom]] # Extrait l'objet flowFrame d'entrée correspondant à l'indexation de la cohorte
        if (is.null(flowframe_entree)) return(NULL) # Quitte proprement la sous-routine si l'échantillon spécifié est introuvable en mémoire
        matrice_exprs <- flowCore::exprs(flowframe_entree) # Extrait la matrice bidimensionnelle des expressions d'intensités brutes de toutes les cellules de l'échantillon
        
        if (axe_discrimination == "H_A") {
          canal_x <- grep("FSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("FSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "W_A") {
          canal_x <- grep("FSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("FSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "H_W") {
          canal_x <- grep("FSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("FSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else {
          stop("axe_discrimination invalide. Choisir parmi 'H_A', 'W_A' ou 'H_W'")
        }
        
        if (is.na(canal_x) || is.na(canal_y)) { # Si l'un des deux signaux morphologiques fondamentaux indispensables à l'évaluation géométrique est manquant
          warning("Canaux FSC requis pour l'axe ", axe_discrimination, " introuvables pour ", nom, ". Étape ignorée.") # Émet un avertissement pour alerter sur l'incompatibilité de la structure du fichier FCS
          self$post_doublets_FSC[[nom]] <- flowframe_entree # Duplique les données d'entrée sans altération pour ne pas rompre la chaîne de traitement en aval
          return(NULL) # Interrompt la sous-routine de cet échantillon pour passer au tube suivant
        } 
        
        ratio_Y_X <- matrice_exprs[, canal_y] / (matrice_exprs[, canal_x] + 1e-6) # Calcule le ratio de linéarité cellulaire en ajoutant un epsilon régulateur pour interdire les divisions critiques par zéro
        val_mad <- stats::mad(ratio_Y_X, na.rm = TRUE) # Calcule l'écart absolu à la médiane (MAD), indicateur robuste de la dispersion de la population de cellules uniques (singlets)
        if (val_mad == 0) val_mad <- mean(ratio_Y_X, na.rm = TRUE) * 0.05 # Securité statistique : si la dispersion est nulle (artefact), calcule une variance de substitution basée sur la moyenne
        seuil_statistique <- stats::median(ratio_Y_X, na.rm = TRUE) + (facteur_sensibilite * val_mad) # Fixe la frontière critique d'exclusion au-delà de laquelle la déformation géométrique trahit un doublet
        
        self$gate_doublets_FSC[[nom]] <- list( # Sauvegarde la structure des paramètres de la coupure statistique au sein de l'environnement R6 pour les exports de métadonnées
          type = "stat", seuil = seuil_statistique, facteur = facteur_sensibilite, # Consigne la nature mathématique du filtre, le seuil calculé et le coefficient multiplicateur appliqué
          source = lbl_src, channels = c(canal_x, canal_y) # Enregistre la provenance des données d'entrée ainsi que le couple d'axes physiques utilisés
        ) 
        
        flowframe_filtre <- flowframe_entree[ratio_Y_X < seuil_statistique & is.finite(ratio_Y_X), ] # Filtre la matrice en excluant les événements au-dessus du seuil (doublets) ou présentant des instabilités mathématiques
        self$post_doublets_FSC[[nom]] <- flowframe_filtre # Enregistre l'échantillon purgé dans la mémoire de stockage intermédiaire dédiée aux filtres FSC
        self$post_doublets_final[[nom]] <- flowframe_filtre # Met à jour la mémoire de stockage terminale des cellules uniques validées (singlets)
      } 
      
      noms_a_traiter <- if (is.null(nom_echantillon)) names(liste_fcs_source) else nom_echantillon # Détermine la liste d'échantillons à traiter en lot (tous les fichiers ou l'identifiant exclusif fourni)
      for (nom in noms_a_traiter) { calculer_doublets(nom, label_source) } # Parcourt et traite séquentiellement via la boucle chaque échantillon configuré dans la liste cible
      self$update_pipeline("doublets_FSC", nom_echantillon) # Active la mise à jour des graphes d'état ou rafraîchit l'interface Shiny pour cette étape d'isolement
    },
    
    retirer_doublets_SSC = function(facteur_sensibilite = 4, axe_discrimination = "H_A", nom_echantillon = NULL) { # Méthode statistique éliminant les agrégats cellulaires (doublets) sur l'axe Side Scatter (SSC) en évaluant la linéarité géométrique des signaux de granularité
      liste_fcs_source <- if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) { # Évalue si la liste issue du filtrage des doublets FSC est préalablement disponible
        self$post_doublets_FSC # Privilégie le chaînage direct après l'isolement des singlets sur l'axe FSC
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) { # À défaut, évalue si la liste filtrée pour les débris est accessible
        self$post_debris # Se connecte en aval de l'étape de filtration des débris cellulaires
      } else { # Si aucune des structures de tri précédentes n'est peuplée ou initialisée
        self$get_derniere_source() # Active le mécanisme de repli automatique sur le dernier niveau de traitement valide du pipeline
      } 
      
      label_source <- if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) { # Détermine l'origine exacte des données pour documenter l'historique d'analyse
        "post_doublets_FSC" # Attribue le label indiquant un chaînage standard post-filtre FSC
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) { # Si l'étape FSC a été contournée ou omise volontairement
        "post_debris" # Attribue le label identifiant la provenance directe du filtre débris
      } else { # Si le script s'est replié sur une structure amont
        "source_parente" # Attribue le label générique de traçabilité par défaut
      } 
      
      if (is.null(self$gate_doublets_SSC)) self$gate_doublets_SSC <- list()
      if (is.null(self$post_doublets_SSC)) self$post_doublets_SSC <- list()
      if (is.null(self$post_doublets_final)) self$post_doublets_final <- list()
      
      calculer_doublets = function(nom, lbl_src) { # Sous-fonction interne isolant les limites de coupure statistique et triant les événements pour un échantillon individuel
        flowframe_entree <- liste_fcs_source[[nom]] # Extrait l'objet flowFrame d'entrée indexé pour la clé de l'itération active
        if (is.null(flowframe_entree)) return(NULL) # Quitte la sous-routine si l'échantillon de la cohorte est inexistant en mémoire vive
        matrice_exprs <- flowCore::exprs(flowframe_entree) # Extrait sous forme de matrice bivariée les intensités de numérisation de chaque cellule
        
        # === ROUTAGE DYNAMIQUE DES CANAUX SSC ===
        if (axe_discrimination == "H_A") {
          canal_x <- grep("SSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("SSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "W_A") {
          canal_x <- grep("SSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("SSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "H_W") {
          canal_x <- grep("SSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("SSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else {
          stop("axe_discrimination invalide. Choisir parmi 'H_A', 'W_A' ou 'H_W'")
        }
        
        if (is.na(canal_x) || is.na(canal_y)) { # Si l'un des deux signaux SSC fondamentaux pour discriminer la géométrie cellulaire est introuvable
          warning("Canaux SSC requis pour l'axe ", axe_discrimination, " introuvables pour ", nom, ". Étape ignorée.") # Alerte l'utilisateur via un avertissement concernant la non-conformité structurelle du fichier
          self$post_doublets_SSC[[nom]] <- flowframe_entree # Duplique les données brutes dans la structure cible pour préserver le flux algorithmique
          return(NULL) # Interrompt proprement la sous-routine de cet échantillon pour basculer sur le tube suivant
        } 
        
        ratio_Y_X <- matrice_exprs[, canal_y] / (matrice_exprs[, canal_x] + 1e-6) # Calcule le ratio de granularité en ajoutant un epsilon stabilisateur pour empêcher les divisions interdites par zéro
        val_mad <- stats::mad(ratio_Y_X, na.rm = TRUE) # Calcule l'écart absolu à la médiane (MAD), estimateur robuste de la dispersion de la population de cellules uniques (singlets)
        if (val_mad == 0) val_mad <- mean(ratio_Y_X, na.rm = TRUE) * 0.05 # Sécurité statistique : calcule une variance artificielle basée sur la moyenne si la dispersion réelle est nulle
        seuil_statistique <- stats::median(ratio_Y_X, na.rm = TRUE) + (facteur_sensibilite * val_mad) # Calcule le seuil discriminant au-delà duquel la déformation temporelle du signal valide la présence d'un agrégat
        
        self$gate_doublets_SSC[[nom]] <- list( # Consigne la totalité des variables de tri au sein de l'environnement R6 pour les besoins de traçabilité
          type = "stat", seuil = seuil_statistique, facteur = facteur_sensibilite, # Sauvegarde la nature du filtre, la valeur de la barrière calculée et la sensibilité appliquée
          source = lbl_src, channels = c(canal_x, canal_y) # Enregistre l'origine des données ainsi que le couple de paramètres SSC sollicités
        )
        
        flowframe_filtre <- flowframe_entree[ratio_Y_X < seuil_statistique & is.finite(ratio_Y_X), ] # Filtre la matrice en excluant les doublets (au-dessus du seuil) et les valeurs infinies aberrantes
        self$post_doublets_SSC[[nom]] <- flowframe_filtre # Sauvegarde l'échantillon épuré dans le compartiment intermédiaire dédié aux filtres SSC
        self$post_doublets_final[[nom]] <- flowframe_filtre # Met à jour la mémoire finale centralisée de l'objet contenant les cellules uniques (singlets) qualifiées
      } 
      
      noms_a_traiter <- if (is.null(nom_echantillon)) names(liste_fcs_source) else nom_echantillon # Cible la cohorte entière si l'identifiant est omis, sinon restreint l'exécution au fichier unique spécifié
      for (nom in noms_a_traiter) { calculer_doublets(nom, label_source) } # Parcourt et traite séquentiellement l'ensemble de la liste via une boucle d'exécution unitaire
      self$update_pipeline("doublets_SSC", nom_echantillon) # Déclenche directement la mise à jour des graphes de suivi ou actualise l'interface graphique UI Shiny
    },
    
    gate_les_doublets_FSC = function(points_utilisateur, axe_discrimination = "H_A", nom_echantillon = NULL) { # Méthode appliquant un fenêtrage polygonal manuel (PolygonGate) fourni par l'utilisateur pour discriminer et exclure les doublets sur les axes Forward Scatter
      liste_fcs_source <- if (!is.null(self$post_debris) && length(self$post_debris) > 0) self$post_debris else self$get_derniere_source() # Sélectionne en priorité les données issues du filtre débris ou active le mécanisme de repli hiérarchique pyramidal
      if (length(liste_fcs_source) == 0) stop("Aucune donnée source disponible pour le gating des doublets FSC.") # Sécurité : bloque l'exécution si aucune matrice cellulaire n'est localisée en mémoire vive
      
      premier_ff  <- liste_fcs_source[[1]] # Extrait le premier objet flowFrame disponible de la liste pour analyser sa structure technique
      tous_canaux <- flowCore::colnames(premier_ff) # Récupère la liste complète des étiquettes de colonnes (canaux physiques de numérisation) de l'appareil
      
      if (axe_discrimination == "H_A") {
        canal_x <- grep("FSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("FSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "W_A") {
        canal_x <- grep("FSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("FSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "H_W") {
        canal_x <- grep("FSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("FSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      }
      
      matrice_coords <- as.matrix(points_utilisateur[, 1:2]) # Force la conversion des deux premières colonnes de coordonnées utilisateur en une matrice R standard pour l'interface flowCore
      colnames(matrice_coords) <- c(canal_x, canal_y) # Aligne et synchronise obligatoirement les noms des colonnes de la matrice sur les canaux physiques cibles détectés
      poly_fsc <- flowCore::polygonGate(.gate = matrice_coords, filterId = "Gate_Doublets_FSC") # Instancie l'objet géométrique formel polygonGate définissant la barrière d'inclusion des cellules uniques
      
      appliquer_fsc = function(nom) { # Sous-fonction encapsulant le tri topologique de point dans un polygone pour un échantillon de la cohorte
        ff_entree  <- liste_fcs_source[[nom]] # Extrait l'objet flowFrame d'entrée correspondant à l'identifiant de l'itération active
        if (is.null(ff_entree)) return(NULL) # Quitte proprement la sous-routine si l'échantillon spécifié est introuvable ou mal chargé
        res_filtre <- flowCore::filter(ff_entree, poly_fsc) # Exécute le filtrage géométrique bidimensionnel pour évaluer l'appartenance de chaque événement cellulaire au polygone
        self$gate_doublets_FSC[[nom]] <- list(type = "poly", gate = poly_fsc, channels = c(canal_x, canal_y)) # Consigne le polygone et les métadonnées de tri dans l'environnement de l'objet R6 pour traçabilité
        ff_propre <- ff_entree[res_filtre@subSet, ] # Sous-échantillonne la matrice d'expression pour ne conserver que les cellules validées par l'indice logique TRUE (singlets)
        self$post_doublets_FSC[[nom]]    <- ff_propre # Enregistre l'échantillon nettoyé dans la mémoire intermédiaire dédiée aux structures FSC
        self$post_doublets_final[[nom]]  <- ff_propre # Met à jour la structure finale de stockage centralisant les cellules uniques qualifiées de l'expérience
      } 
      
      noms <- if (is.null(nom_echantillon)) names(liste_fcs_source) else nom_echantillon # Cible la cohorte complète si l'identifiant est omis, ou restreint l'exécution au fichier unique spécifié
      for (n in noms) { appliquer_fsc(n) } # Parcourt et traite séquentiellement l'ensemble de la liste via une boucle d'exécution unitaire
      if (!is.null(self$update_pipeline)) self$update_pipeline("doublets_FSC", nom_echantillon) # Déclenche la mise à jour des graphes de suivi ou actualise l'interface graphique Shiny
    },
    
    gate_les_doublets_SSC = function(points_utilisateur, axe_discrimination = "H_A", nom_echantillon = NULL) { # Méthode appliquant un fenêtrage polygonal manuel (PolygonGate) fourni par l'utilisateur pour discriminer et exclure les doublets sur les axes Side Scatter
      liste_fcs_source <- if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) { # Évalue si la liste issue du filtrage ou du gating des doublets FSC est disponible en mémoire
        self$post_doublets_FSC # Privilégie le chaînage direct des filtres en se connectant en aval des cellules uniques isolées sur l'axe FSC
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) { # À défaut, évalue si la liste filtrée pour les débris cellulaires est accessible
        self$post_debris # Connecte le flux en aval immédiat de l'étape de filtration des débris
      } else { # Si aucune des structures de tri précédentes n'est peuplée ou initialisée
        self$get_derniere_source() # Active le mécanisme de repli automatique sur le dernier niveau de traitement valide détecté dans le pipeline
      } 
      if (length(liste_fcs_source) == 0) stop("Aucune donnée source disponible pour le gating des doublets SSC.") # Sécurité : bloque l'exécution si aucune matrice cellulaire n'est localisée en mémoire vive
      
      premier_ff  <- liste_fcs_source[[1]] # Extrait le premier objet flowFrame disponible de la liste pour analyser sa structure technique
      tous_canaux <- flowCore::colnames(premier_ff) # Récupère la liste complète des étiquettes de colonnes (canaux physiques de numérisation) de l'appareil
      
      if (axe_discrimination == "H_A") {
        canal_x <- grep("SSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("SSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "W_A") {
        canal_x <- grep("SSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("SSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "H_W") {
        canal_x <- grep("SSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("SSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      }
      
      matrice_coords <- as.matrix(points_utilisateur[, 1:2]) # Force la conversion des deux premières colonnes de coordonnées utilisateur en une matrice R standard pour l'interface flowCore
      colnames(matrice_coords) <- c(canal_x, canal_y) # Aligne et synchronise obligatoirement les noms des colonnes de la matrice sur les canaux physiques cibles détectés
      poly_ssc <- flowCore::polygonGate(.gate = matrice_coords, filterId = "Gate_Doublets_SSC") # Instancie l'objet géométrique formel polygonGate définissant la barrière d'inclusion des cellules uniques
      
      appliquer_ssc = function(nom) { # Sous-fonction encapsulant le tri topologique de point dans un polygone pour un échantillon de la cohorte
        ff_entree <- liste_fcs_source[[nom]] # Extrait l'objet flowFrame d'entrée correspondant à l'identifiant de l'itération active
        if (is.null(ff_entree)) return(NULL) # Quitte proprement la sous-routine si l'échantillon spécifié est introuvable ou mal chargé
        res_filtre <- flowCore::filter(ff_entree, poly_ssc) # Exécute le filtrage géométrique bidimensionnel pour évaluer l'appartenance de chaque événement cellulaire au polygone
        self$gate_doublets_SSC[[nom]] <- list(type = "poly", gate = poly_ssc, channels = c(canal_x, canal_y)) # Consigne le polygone et les métadonnées de tri dans l'environnement de l'objet R6 pour traçabilité
        ff_propre <- ff_entree[res_filtre@subSet, ] # Sous-échantillonne la matrice d'expression pour ne conserver que les cellules validées par l'indice logique TRUE (singlets)
        self$post_doublets_SSC[[nom]]   <- ff_propre # Enregistre l'échantillon nettoyé dans la mémoire intermédiaire dédiée aux structures SSC
        self$post_doublets_final[[nom]] <- ff_propre # Met à jour la structure finale de stockage centralisant les cellules uniques qualifiées de l'expérience
      } 
      
      noms <- if (is.null(nom_echantillon)) names(liste_fcs_source) else nom_echantillon # Cible la cohorte complète si l'identifiant est omis, ou restreint l'exécution au fichier unique spécifié
      for (n in noms) { appliquer_ssc(n) } # Parcourt et traite séquentiellement l'ensemble de la liste via une boucle d'exécution unitaire
      if (!is.null(self$update_pipeline)) self$update_pipeline("doublets_SSC", nom_echantillon) # Déclenche la mise à jour des graphes de suivi ou actualise l'interface graphique Shiny
    },
    
    visualiser_peacoqc = function(nom_echantillon, max_points = 10000) { # Méthode générant un graphique cinétique comparatif pour évaluer visuellement l'efficacité du filtrage de bruit de flux opéré par PeacoQC
      if (is.null(self$post_PeacoQC[[nom_echantillon]])) { # Évalue si la structure ou le fichier cible nettoyé par PeacoQC est absent de la mémoire vive
        message("Pas de données PeacoQC pour ", nom_echantillon) # Notification d'avertissement en console si l'étape amont n'a pas été exécutée
        return(NULL) # Interrompt proprement la fonction et renvoie NULL pour ne pas provoquer de plantage de l'interface
      } 
      
      flowframe_initial <- self$echantillons_traites[[nom_echantillon]] # Récupère l'objet flowFrame initial d'origine (avant le processus PeacoQC)
      flowframe_nettoye <- self$post_PeacoQC[[nom_echantillon]] # Extrait l'objet flowFrame épuré contenant uniquement la population cellulaire conforme
      donnees_initiales <- as.data.frame(flowCore::exprs(flowframe_initial)) # Convertit la matrice des intensités d'origine en tableau de données exploitable par ggplot2
      donnees_nettoyees <- as.data.frame(flowCore::exprs(flowframe_nettoye)) # Convertit la matrice des intensités épurées en tableau de données exploitable par ggplot2
      canal_temps <- grep("time", colnames(donnees_initiales), value = TRUE, ignore.case = TRUE)[1] # Détecte par expression régulière le libellé de la colonne associée à la variable temporelle
      canal_fsc   <- grep("FSC", colnames(donnees_initiales), value = TRUE, ignore.case = TRUE)[1] # Détecte par expression régulière le libellé de la colonne associée au paramètre morphologique FSC
      if (is.na(canal_temps)) canal_temps <- colnames(donnees_initiales)[1] # Sécurité : force l'affectation de la première colonne si aucun paramètre temporel n'est trouvé
      if (is.na(canal_fsc))   canal_fsc   <- colnames(donnees_initiales)[2] # Sécurité : force l'affectation de la deuxième colonne si aucun paramètre morphologique n'est trouvé
      limites_temps <- range(donnees_initiales[[canal_temps]], na.rm = TRUE) # Calcule l'amplitude minimale et maximale du temps pour fixer les frontières absolues de l'axe X
      limites_fsc   <- range(donnees_initiales[[canal_fsc]], na.rm = TRUE) # Calcule l'amplitude minimale et maximale du signal FSC pour fixer les frontières absolues de l'axe Y
      total_evenements  <- nrow(donnees_initiales) # Dénombre le nombre total d'événements cellulaires bruts enregistrés à l'acquisition
      evenements_gardes <- nrow(donnees_nettoyees) # Dénombre le nombre de cellules conservées post-contrôle qualité statistique
      pourcentage_conservation <- if(total_evenements > 0) round((evenements_gardes / total_evenements) * 100, 1) else 0 # Calcule le rendement d'acquisition après filtration, arrondi au dixième
      
      if (!is.null(max_points) && total_evenements > max_points) { # Si la taille de l'échantillon outrepasse la limite maximale définie pour la fluidité du rendu graphique
        
        if (!is.null(self$seed)) set.seed(self$seed)
        donnees_init_visu  <- donnees_initiales[sample(seq_len(total_evenements), max_points), ] # Sous-échantillonne de manière aléatoire la matrice initiale pour alléger la charge graphique
        
        if (evenements_gardes > 0) { # S'il reste des cellules viables après le passage de l'algorithme PeacoQC
          nb_nettoye_visu   <- min(evenements_gardes, max_points) # Détermine la borne supérieure optimale de points à conserver pour le rendu des données épurées
          
          # === SÉCURITÉ SEED : Fixe la graine avant le sous-échantillonnage nettoyé ===
          if (!is.null(self$seed)) set.seed(self$seed)
          donnees_nett_visu  <- donnees_nettoyees[sample(seq_len(evenements_gardes), nb_nettoye_visu), ] # Sous-échantillonne aléatoirement la matrice nettoyée au même prorata visuel
        } else { # Si l'échantillon a été intégralement rejeté par le contrôle qualité
          donnees_nett_visu  <- donnees_nettoyees # Assigne la structure vide directement sans calcul d'échantillonnage
        } 
      } else { # Si le nombre total de cellules est inférieur au seuil max_points
        donnees_init_visu  <- donnees_initiales # Conserve l'intégralité de la matrice de départ pour la représentation graphique
        donnees_nett_visu  <- donnees_nettoyees # Conserve l'intégralité de la matrice épurée pour la représentation graphique
      } 
      
      lbl_x <- if (!is.null(self$get_label)) self$get_label(flowframe_initial, canal_temps) else canal_temps # Extrait le libellé biologique de l'axe X via get_label, ou utilise le nom technique brut
      lbl_y <- if (!is.null(self$get_label)) self$get_label(flowframe_initial, canal_fsc) else canal_fsc # Extrait le libellé biologique de l'axe Y via get_label, ou utilise le nom technique brut
      graphique_qc <- ggplot2::ggplot() + # Initialise l'objet graphique ggplot2 vide multicouche
        ggplot2::geom_point(data = donnees_init_visu, ggplot2::aes(x = .data[[canal_temps]], y = .data[[canal_fsc]]), # Ajoute la couche de fond représentant l'ensemble des points d'origine
                            size = 0.2, alpha = 0.2, color = "grey70") # Paramètre les points en gris clair et transparent pour matérialiser les événements exclus ou masqués
      if (nrow(donnees_nett_visu) > 0) { # Si la sous-matrice nettoyée contient des événements cellulaires à tracer
        graphique_qc <- graphique_qc + # Superpose une seconde couche d'événements par-dessus le bruit de fond gris
          ggplot2::geom_point(data = donnees_nett_visu, ggplot2::aes(x = .data[[canal_temps]], y = .data[[canal_fsc]]), # Spécifie les coordonnées cinétiques des cellules conformes validées
                              size = 0.2, alpha = 0.4, color = "darkblue") # Paramètre les cellules saines en bleu foncé contrasté pour mettre en évidence les zones d'instabilité supprimées
      } 
      graphique_qc <- graphique_qc + # Finalise la mise en forme structurelle et textuelle de la figure de diagnostic
        ggplot2::coord_cartesian(xlim = limites_temps, ylim = limites_fsc) + # Force un cadrage strict sur les limites initiales du fichier pour éviter tout effet de zoom déformant
        ggplot2::theme_bw() + # Applique un habillage blanc structuré et épuré facilitant la lecture des densités de points
        ggplot2::theme( # Ajuste les propriétés typographiques de la figure de contrôle
          plot.title = ggplot2::element_text(face = "bold"), # Renforce la visibilité du titre principal en l'affichant en caractères gras
          plot.subtitle = ggplot2::element_text(color = "darkblue", size = 11) # Distingue le sous-texte statistique en lui appliquant une coloration bleue
        ) + # Fin des ajustements de thème
        ggplot2::labs( # Définit l'ensemble des titres et des descriptions scientifiques entourant les axes
          title = paste("Contrôle qualité PeacoQC :", nom_echantillon), # Affiche dynamiquement le titre de la méthode couplé au nom du fichier FCS audité
          subtitle = paste0("Événements conservés : ", format(evenements_gardes, big.mark=" "), # Documente les métriques de tri incluant le décompte des cellules saines restantes
                            " | ", pourcentage_conservation, "% (Affichage max : ", max_points, " pts)"), # Affiche le rendement en pourcentage et précise le niveau de sous-échantillonnage graphique appliqué
          x = lbl_x, y = lbl_y # Attribue les libellés biologiques ou techniques finaux aux axes X et Y
        ) 
      
      if (is.null(self$plots_peacoqc)) self$plots_peacoqc <- list() # Initialise la structure de liste dédiée au stockage des graphiques PeacoQC si inexistante en mémoire
      self$plots_peacoqc[[nom_echantillon]] <- graphique_qc # Archive l'objet graphique au sein de l'environnement R6 pour permettre des exports en lot ultérieurs
      return(graphique_qc) # Renvoie l'objet graphique complet, prêt pour affichage à l'écran ou intégration dans une interface Shiny
    },
    
    visualiser_flowai = function(nom_echantillon, max_points = 10000) { # Méthode générant un graphique cinétique de contrôle qualité pour visualiser l'impact du nettoyage flowAI sur un échantillon
      if (is.null(self$post_flowAI) || is.null(self$post_flowAI[[nom_echantillon]])) { # Vérifie si la structure ou l'échantillon ciblé par le nettoyage flowAI est manquant en mémoire
        stop("Aucun résultat flowAI trouvé pour cet échantillon.") # Interrompt le script et exige l'exécution préalable de la méthode appliquer_flowai
      } 
      
      ff_nettoye <- self$post_flowAI[[nom_echantillon]] # Extrait l'objet flowFrame épuré (contenant uniquement les cellules validées par le QC)
      fcs_initial <- self$echantillons_traites[[nom_echantillon]] # Récupère l'échantillon compensé d'origine (avant QC)
      exprs_initiales <- as.data.frame(flowCore::exprs(fcs_initial)) # Extrait sous forme de tableau de données R la matrice d'expression de tous les événements cellulaires initiaux
      canal_temps <- grep("time", colnames(exprs_initiales), ignore.case = TRUE, value = TRUE)[1] # Détecte dynamiquement par expression régulière le nom du canal dédié au suivi du temps d'acquisition
      canal_taille <- grep("FSC", colnames(exprs_initiales), ignore.case = TRUE, value = TRUE)[1] # Détecte dynamiquement par expression régulière le nom du canal de taille cellulaire relative (Forward Scatter)
      if (is.na(canal_temps) || is.na(canal_taille)) { # Sécurité : si la détection automatique des canaux échoue ou si les libellés sont atypiques
        canal_temps <- colnames(exprs_initiales)[ncol(exprs_initiales)] # Assigne par défaut la toute dernière colonne de la matrice comme axe temporel
        canal_taille <- colnames(exprs_initiales)[1] # Assigne par défaut la toute première colonne de la matrice comme paramètre morphologique FSC
      } 
      exprs_nettoyees <- as.data.frame(flowCore::exprs(ff_nettoye)) # Extrait sous forme de tableau de données R la matrice d'expression des événements conservés post-QC
      exprs_initiales$Status <- "Éliminé (flowAI)" # Initialise par défaut l'état de chaque événement de la matrice d'origine comme étant rejeté par le filtre
      indices_conserves <- which(exprs_initiales[[canal_temps]] %in% exprs_nettoyees[[canal_temps]]) # Identifie par intersection les indices temporels des cellules ayant survécu au nettoyage de flowAI
      exprs_initiales$Status[indices_conserves] <- "Conservé" # Assigne le statut de conformité aux événements cellulaires validés par l'algorithme
      total_pts <- nrow(exprs_initiales) # Mémorise le nombre total d'événements cellulaires initialement présents dans le fichier d'acquisition
      
      if (total_pts > max_points) { # Si la taille de la matrice dépasse le seuil maximal de points fixé pour le tracé graphique
        if (!is.null(self$seed)) set.seed(self$seed)
        exprs_initiales <- exprs_initiales[sample(total_pts, max_points), ] # Échantillonne aléatoirement un nombre restreint de lignes pour optimiser le rendu graphique sans saturer la mémoire
      } 
      total_conserves <- length(indices_conserves) # Calcule le décompte absolu d'événements cellulaires conformes conservés post-QC
      pourcentage_conservation <- round((total_conserves / total_pts) * 100, 1) # Déduit le rendement d'acquisition exprimé en pourcentage de cellules saines conservées
      
      graphique_flowai <- ggplot2::ggplot(exprs_initiales, ggplot2::aes(x = .data[[canal_temps]], y = .data[[canal_taille]], color = Status)) + # Initialise le graphique biplot Temps vs FSC
        ggplot2::geom_point(size = 0.4, alpha = 0.6) + # Dessine le nuage de points cytométriques avec une taille fine et une légère transparence pour révéler l'empilement
        ggplot2::scale_color_manual(values = c("Conservé" = "#1f77b4", "Éliminé (flowAI)" = "#d62728")) + # Applique un code couleur binaire contrasté (bleu pour le signal sain, rouge pour les anomalies)
        ggplot2::theme_bw() + # Applique un habillage blanc structuré et épuré facilitant l'évaluation visuelle des coupures cinétiques
        ggplot2::labs( # Configure les textes et légendes scientifiques entourant la figure de diagnostic
          title = paste("Contrôle Qualité flowAI :", nom_echantillon), # Affiche le titre de l'analyse associé au nom du fichier FCS audité
          subtitle = paste0("Événements conservés : ", format(total_conserves, big.mark=" "), " / ", 
                            format(total_pts, big.mark=" "), " (", pourcentage_conservation, "%)"), # Affiche les métriques clés de rendement et d'élimination de la routine de QC
          x = paste("Axe du Temps :", canal_temps), # Documente le nom exact du canal temporel de l'axe des abscisses
          y = paste("Axe Morphologique :", canal_taille) # Documente le nom exact du canal de taille de l'axe des ordonnées
        ) + 
        ggplot2::theme( # Ajuste finement la typographie et la disposition des éléments de la figure
          legend.position = "bottom", # Positionne la légende des statuts sous le graphique pour maximiser la largeur de la zone de tracé
          plot.title    = ggplot2::element_text(face = "bold", size = 14), # Renforce l'accentuation visuelle du titre principal du diagnostic
          legend.title  = ggplot2::element_blank() # Masque l'intitulé de la légende devenu superflu grâce à l'explicitation des étiquettes
        ) 
      
      if (is.null(self$plots_flowai)) self$plots_flowai <- list() # Initialise la structure de liste dédiée au stockage des graphiques flowAI si inexistante en mémoire
      self$plots_flowai[[nom_echantillon]] <- graphique_flowai # Archive l'objet graphique au sein de l'environnement R6 pour permettre sa réutilisation dans le résumé PDF d'export
      return(graphique_flowai) # Renvoie l'objet graphique ggplot2 complet prêt pour affichage ou intégration UI Shiny
    },
    
    visualiser_debris = function(nom_echantillon) {
      
      # ───────────────────────────────────────────────
      # 1. Vérification des données disponibles
      # ───────────────────────────────────────────────
      if (is.null(self$post_debris[[nom_echantillon]])) {
        message("Pas de données Débris pour ", nom_echantillon)
        return(NULL)
      }
      
      # Sélection de la source "avant débris"
      flowframe_avant <- if (!is.null(self$post_retrait_bordures) && length(self$post_retrait_bordures) > 0) {
        self$post_retrait_bordures[[nom_echantillon]]
      } else if (!is.null(self$post_PeacoQC) && length(self$post_PeacoQC) > 0) {
        self$post_PeacoQC[[nom_echantillon]]
      } else if (!is.null(self$post_flowAI) && length(self$post_flowAI) > 0) {
        self$post_flowAI[[nom_echantillon]]
      } else {
        self$echantillons_traites[[nom_echantillon]]
      }
      
      if (is.null(flowframe_avant)) return(NULL)
      
      flowframe_apres <- self$post_debris[[nom_echantillon]]
      donnees_globales <- as.data.frame(flowCore::exprs(flowframe_avant))
      
      # ───────────────────────────────────────────────
      # 2. Extraction du polygone de gating
      # ───────────────────────────────────────────────
      gate_polygone    <- self$gate_debris[[nom_echantillon]]
      coordonnees_gate <- as.data.frame(gate_polygone@boundaries)
      colnames(coordonnees_gate) <- c("x", "y")
      
      params   <- flowCore::parameters(gate_polygone)
      canal_x  <- params[1]
      canal_y  <- params[2]
      
      total_evenements_avant <- nrow(donnees_globales)
      total_evenements_apres <- nrow(flowCore::exprs(flowframe_apres))
      pourcentage_conservation <- if (total_evenements_avant > 0)
        round((total_evenements_apres / total_evenements_avant) * 100, 1) else 0
      
      lbl_x <- self$get_label(flowframe_avant, canal_x)
      lbl_y <- self$get_label(flowframe_avant, canal_y)
      
      lim_x_globale <- c(0, max(donnees_globales[[canal_x]], na.rm = TRUE))
      lim_y_globale <- c(0, max(donnees_globales[[canal_y]], na.rm = TRUE))
      
      # ───────────────────────────────────────────────
      # 3. Densité (binning raster rapide, pixels petits, sans lissage pour
      #    bien distinguer les zones) — même helper que pour les débris/doublets/
      #    viabilité, afin de garder un rendu uniforme dans toute l'application.
      # ───────────────────────────────────────────────
      df_densite <- calculer_densite_raster(donnees_globales[[canal_x]], donnees_globales[[canal_y]], lim_x_globale, lim_y_globale)
      
      # ───────────────────────────────────────────────
      # 4. Construction du graphique (geom_raster)
      # ───────────────────────────────────────────────
      graphique_debris <- if (is.null(df_densite)) {
        ggplot2::ggplot() + ggplot2::theme_bw() +
          ggplot2::labs(title = paste("Nettoyage des débris :", nom_echantillon), subtitle = "Pas assez d'événements pour tracer la densité", x = lbl_x, y = lbl_y)
      } else {
        ggplot2::ggplot(df_densite, ggplot2::aes(x = X, y = Y, fill = densite)) +
          ggplot2::geom_raster() + # pixels réguliers et petits (grille haute résolution), rendu rapide même avec beaucoup d'événements
          ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE) +
          ggplot2::geom_polygon(
            data = coordonnees_gate,
            ggplot2::aes(x = x, y = y),
            fill = NA, color = "black", linewidth = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::coord_cartesian(xlim = lim_x_globale, ylim = lim_y_globale) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            legend.position = "none",
            plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(color = "darkblue", size = 11)
          ) +
          ggplot2::labs(
            title = paste("Nettoyage des débris :", nom_echantillon),
            subtitle = paste0(
              "Événements conservés : ",
              format(total_evenements_apres, big.mark = " "),
              " | ", pourcentage_conservation, "%"
            ),
            x = lbl_x, y = lbl_y
          )
      }
      
      # ───────────────────────────────────────────────
      # 5. Stockage dans l'objet CARROT
      # ───────────────────────────────────────────────
      if (is.null(self$plots_debris)) self$plots_debris <- list()
      self$plots_debris[[nom_echantillon]] <- graphique_debris
      
      return(graphique_debris)
    },
    
    visualiser_doublets = function(nom_echantillon, type_analyse = "FSC") { # Méthode générant un graphique bidimensionnel de densité (biplot) pour cartographier et valider l'exclusion des agrégats cellulaires (doublets) selon l'axe choisi (FSC ou SSC)
      infos_gate <- if (type_analyse == "FSC") self$gate_doublets_FSC[[nom_echantillon]] else self$gate_doublets_SSC[[nom_echantillon]] # Extrait les paramètres et métadonnées de la barrière de tri (statistique ou polygonale) correspondant au paramètre spécifié
      if (is.null(infos_gate)) return(NULL) # interrompt proprement la fonction si aucune information de gating n'est localisée pour cet échantillon
      
      ff_avant <- if (type_analyse == "FSC") { # Routage adaptatif de la source amont : si l'analyse porte sur les doublets Forward Scatter
        if (!is.null(self$post_debris) && length(self$post_debris) > 0) self$post_debris[[nom_echantillon]] else self$get_derniere_source()[[nom_echantillon]] # Récupère les données post-débris ou applique le mécanisme de repli hiérarchique sur la dernière source valide
      } else { # Sinon, si l'analyse porte sur les doublets Side Scatter (granularité)
        if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) { # Évalue si les cellules ont déjà subi l'exclusion des doublets sur l'axe FSC
          self$post_doublets_FSC[[nom_echantillon]] # Charge la matrice cellulaire résultant du premier niveau d'exclusion des doublets FSC
        } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) { # À défaut, vérifie si la structure filtrée pour les débris est accessible
          self$post_debris[[nom_echantillon]] # Connecte le flux en aval immédiat de l'étape de filtration des débris
        } else { # Si aucun filtre intermédiaire n'est présent en mémoire de l'objet
          self$get_derniere_source()[[nom_echantillon]] # Se replie automatiquement sur le dernier état de traitement valide disponible
        } 
      } 
      if (is.null(ff_avant)) return(NULL) # Sécurité critique : quitte proprement la routine si aucun flowFrame antécédent valide n'est identifié
      
      ff_apres <- if (type_analyse == "FSC") self$post_doublets_FSC[[nom_echantillon]] else self$post_doublets_SSC[[nom_echantillon]] # Extrait l'objet flowFrame épuré (contenant uniquement les cellules uniques validées) associé à l'axe d'analyse courant
      donnees_source <- as.data.frame(flowCore::exprs(ff_avant)) # Convertit la matrice des intensités de la source d'entrée en tableau de données exploitable par ggplot2
      
      nb_avant <- nrow(donnees_source) # Dénombre la population totale d'événements cellulaires entrant dans l'étape de discrimination des doublets
      # Densité calculée par binning raster : on garde systématiquement l'intégralité des
      # événements (pas de sous-échantillonnage), le coût du binning restant indépendant de N.
      donnees_visu <- donnees_source
      
      canal_x <- infos_gate$channels[1] # Extrait le nom technique du canal affecté à l'axe des abscisses (généralement le paramètre de Hauteur ou Largeur)
      canal_y <- infos_gate$channels[2] # Extrait le nom technique du canal affecté à l'axe des ordonnées (généralement le paramètre d'Aire)
      nb_apres    <- if (!is.null(ff_apres)) nrow(flowCore::exprs(ff_apres)) else 0 # Dénombre précisément la population de cellules uniques (singlets) conservée après application de la coupure
      pourcentage <- if (nb_avant > 0) round((nb_apres / nb_avant) * 100, 1) else 0 # Déduit le rendement de la filtration en pourcentage, arrondi au dixième
      limite_max  <- max(donnees_source[, c(canal_x, canal_y)], na.rm = TRUE) # Détermine la valeur maximale absolue observée sur les deux canaux pour harmoniser et proportionner les axes de la figure
      lbl_x <- if (!is.null(self$get_label)) self$get_label(ff_avant, canal_x) else canal_x # Extrait le libellé biologique ou l'étiquette de l'axe X via get_label, ou garde le nom technique brut
      lbl_y <- if (!is.null(self$get_label)) self$get_label(ff_avant, canal_y) else canal_y # Extrait le libellé biologique ou l'étiquette de l'axe Y via get_label, ou garde le nom technique brut
      
      # --- Détermination de l'appartenance conservé / retiré pour chaque évènement échantillonné ---
      if (!is.null(ff_apres) && nb_apres > 0) { # Seulement si un résultat de filtration existe réellement pour cet échantillon
        donnees_apres <- as.data.frame(flowCore::exprs(ff_apres)) # Convertit la matrice des cellules retenues (singlets) en tableau de données comparable
        
        colonnes_cle <- intersect(colnames(donnees_visu), colnames(donnees_apres)) # Sécurise la comparaison en ne retenant que les canaux communs aux deux tableaux
        cle_visu  <- do.call(paste, c(lapply(colonnes_cle, function(col) donnees_visu[[col]]),  sep = "\r")) # Construit une clé composite unique par évènement à partir de l'ensemble des canaux (réduit fortement le risque de collision)
        cle_apres <- do.call(paste, c(lapply(colonnes_cle, function(col) donnees_apres[[col]]), sep = "\r")) # Construit la même clé composite pour la population conservée en aval
        
        donnees_visu$statut_doublet <- ifelse(cle_visu %in% cle_apres, "Conservé", "Retiré") # Étiquette chaque évènement échantillonné selon sa présence ou non dans la population post-filtration
      } else { # Si aucune filtration n'a encore été appliquée (ff_apres absent ou vide)
        donnees_visu$statut_doublet <- "Conservé" # Par défaut, tout est considéré comme conservé pour ne pas fausser la lecture avant tri effectif
      } 
      
      donnees_retirees   <- donnees_visu[donnees_visu$statut_doublet == "Retiré", ]   # Isole le sous-ensemble des évènements exclus par le gate de doublets
      donnees_conservees <- donnees_visu[donnees_visu$statut_doublet == "Conservé", ] # Isole le sous-ensemble des évènements retenus (singlets) qui recevront la coloration par densité
      
      couleur_retire <- "darkred" # Teinte neutre et volontairement éloignée de la palette spectrale, pour signaler sans ambiguïté les évènements exclus
      couleur_legende_conserve <- "#2C7FB8" # Nuance bleutée représentative du dégradé de densité, utilisée uniquement comme pastille de légende (n'affecte pas la coloration réelle des points)
      
      # Densité (population conservée) calculée par binning raster (rapide, indépendant du
      # nombre d'événements) plutôt que par ggpointdensity (densité par point, coûteuse).
      df_densite_cons <- calculer_densite_raster(donnees_conservees[[canal_x]], donnees_conservees[[canal_y]], c(0, limite_max), c(0, limite_max))
      
      graphique <- ggplot2::ggplot() # Initialise le biplot morphologique de linéarité en couches séparées afin de distinguer visuellement les évènements retirés des évènements conservés
      if (!is.null(df_densite_cons)) {
        graphique <- graphique +
          ggplot2::geom_raster(data = df_densite_cons, ggplot2::aes(x = X, y = Y, fill = densite), interpolate = TRUE) + # Affiche la densité des évènements conservés sous forme d'image raster (rapide, même avec beaucoup d'événements)
          ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE) # Applique l'échelle pseudo-spectrale standardisée pour cartographier les gradients de population conservée
      }
      if (nrow(donnees_retirees) > 0) { # N'ajoute la couche des évènements retirés que si elle est non vide pour cet échantillon
        graphique <- graphique +
          ggplot2::geom_point(data = donnees_retirees, ggplot2::aes(x = .data[[canal_x]], y = .data[[canal_y]]), color = couleur_retire, size = 0.3, alpha = 0.6) # Superpose en rouge foncé les évènements exclus, par-dessus la densité (raster opaque), pour matérialiser leur rejet
      }
      graphique <- graphique +
        ggplot2::coord_cartesian(xlim = c(0, limite_max), ylim = c(0, limite_max)) + # Force un cadrage isométrique et strict de l'origine jusqu'au maximum pour préserver la diagonale théorique des singlets
        ggplot2::theme_bw() + # Applique un habillage blanc structuré et épuré optimisant l'évaluation des contrastes thermiques de densité
        ggplot2::labs( # Configure l'ensemble des métadonnées, titres et labels scientifiques entourant la figure
          title = paste("Retrait des doublets", type_analyse, ":", nom_echantillon), # Affiche dynamiquement le titre de la manipulation couplé au nom du fichier FCS audité
          subtitle = paste0("Événements conservés : ", format(nb_apres, big.mark=" "), " | ", pourcentage, "%"), # Resitue les métriques d'efficacité du tri (nombre d'événements uniques et rendement)
          x = lbl_x, y = lbl_y # Injecte les libellés biologiques nettoyés aux axes correspondants
        ) + 
        ggplot2::theme(legend.position = "none", plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "darkblue", size = 11)) # Masque la barre d'échelle continue des densités (redondante) et stylise les titres ; la légende conservé/retiré est gérée manuellement ci-dessous
      
      if (infos_gate$type == "poly") { # Couche géométrique adaptative : si le filtre appliqué provient d'un gating polygonal manuel de l'utilisateur
        coordonnees_gate <- as.data.frame(infos_gate$gate@boundaries) # Extrait sous forme de tableau bidimensionnel les sommets géométriques (X, Y) du polygone de tri manuel
        colnames(coordonnees_gate) <- c("x", "y") # Renomme explicitement les colonnes du tableau pour simplifier l'intégration géométrique
        graphique <- graphique + ggplot2::geom_polygon(data = coordonnees_gate, ggplot2::aes(x = x, y = y), fill = NA, color = "darkred", linewidth = 0.6) # Superpose les contours du polygone utilisateur sous forme de ligne rouge continue
      } 
      # NOTE : la barre rouge pointillée (seuil statistique MAD, type == "stat") a été retirée à la demande de l'utilisateur ;
      # la distinction conservé/retiré repose désormais uniquement sur la coloration des points et la légende manuelle ci-dessous
      
      # --- Légende manuelle "Conservé / Retiré", indépendante de toute échelle de couleur pour éviter les conflits avec le dégradé de densité ---
      pas_legende <- limite_max * 0.05 # Définit un espacement vertical proportionnel entre les deux entrées de la légende manuelle
      y_conserve  <- limite_max * 0.97 # Position verticale de la première entrée ("Conservé") en haut du cadre du graphique
      y_retire    <- y_conserve - pas_legende # Position verticale de la seconde entrée ("Retiré"), juste en dessous de la première
      x_pastille  <- limite_max * 0.80 # Position horizontale commune des pastilles de couleur de la légende
      x_texte     <- limite_max * 0.83 # Position horizontale du texte explicatif, légèrement décalée à droite des pastilles
      
      graphique <- graphique +
        ggplot2::annotate("point", x = x_pastille, y = y_conserve, color = couleur_legende_conserve, size = 2.2) + # Pastille représentative de la coloration des évènements conservés (teinte issue du dégradé de densité)
        ggplot2::annotate("text", x = x_texte, y = y_conserve, label = "Conservé", hjust = 0, size = 3.2, color = "black") + # Libellé textuel associé à la pastille des évènements conservés
        ggplot2::annotate("point", x = x_pastille, y = y_retire, color = couleur_retire, size = 2.2) + # Pastille représentative de la couleur grise attribuée aux évènements retirés
        ggplot2::annotate("text", x = x_texte, y = y_retire, label = "Retiré", hjust = 0, size = 3.2, color = "black") # Libellé textuel associé à la pastille des évènements retirés
      
      if (is.null(self$plots_doublets)) self$plots_doublets <- list() # Initialise la structure de liste dédiée au stockage des graphiques de doublets si absente en mémoire de l'objet R6
      self$plots_doublets[[paste0(nom_echantillon, "_", type_analyse)]] <- graphique # Archive l'objet visuel généré au sein du registre d'environnement de la classe R6 en indexant par échantillon et axe
      return(graphique) # Renvoie l'objet graphique ggplot2 complet, prêt pour affichage immédiat ou intégration dans une interface UI Shiny
    },
    
    
    retirer_les_cellules_mortes = function(canal_fsc = "FSC-A", marqueur_viabilite,
                                           points_utilisateur, nom_echantillon = NULL) {
      
      if (is.null(points_utilisateur)) stop("Aucun point fourni pour la viabilité.")
      if (nrow(points_utilisateur) < 3) stop("Un polygone nécessite ≥ 3 points.")
      
      # Source avant viabilité : transformation > doublets > débris > brut
      liste_fcs_source <- if (!is.null(self$post_transformation) && length(self$post_transformation) > 0) {
        self$post_transformation
      } else if (!is.null(self$post_doublets_final) && length(self$post_doublets_final) > 0) {
        self$post_doublets_final
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) {
        self$post_debris
      } else {
        self$get_derniere_source()
      }
      
      if (is.null(self$gate_viabilite)) self$gate_viabilite <- list()
      if (is.null(self$post_viabilite)) self$post_viabilite <- list()
      
      matrice_coordonnees <- as.matrix(points_utilisateur[, 1:2])
      colnames(matrice_coordonnees) <- c(canal_fsc, marqueur_viabilite)
      
      polygone_viabilite <- flowCore::polygonGate(.gate = matrice_coordonnees, filterId = "Gate_Viabilite")
      
      appliquer_le_gate_vivantes = function(nom) {
        ff <- liste_fcs_source[[nom]]
        if (is.null(ff)) return(NULL)
        
        resultat_filtrage <- flowCore::filter(ff, polygone_viabilite)
        
        self$gate_viabilite[[nom]] <- polygone_viabilite
        self$post_viabilite[[nom]] <- ff[resultat_filtrage@subSet, ]
      }
      
      echantillons <- if (is.null(nom_echantillon)) names(liste_fcs_source) else nom_echantillon
      for (nom in echantillons) appliquer_le_gate_vivantes(nom)
      
      self$update_pipeline("viabilite", nom_echantillon)
    },
    
    visualiser_viabilite = function(nom_echantillon) {
      
      if (is.null(self$post_viabilite[[nom_echantillon]])) {
        message("Pas de données Viabilité pour ", nom_echantillon)
        return(NULL)
      }
      
      # Source "avant" gating : si une version transformée existe, on l'utilise
      flowframe_avant <- if (!is.null(self$post_transformation) &&
                             !is.null(self$post_transformation[[nom_echantillon]])) {
        self$post_transformation[[nom_echantillon]]
      } else if (!is.null(self$post_doublets_final[[nom_echantillon]])) {
        self$post_doublets_final[[nom_echantillon]]
      } else if (!is.null(self$post_debris[[nom_echantillon]])) {
        self$post_debris[[nom_echantillon]]
      } else {
        self$get_derniere_source()[[nom_echantillon]]
      }
      
      if (is.null(flowframe_avant)) return(NULL)
      
      donnees_globales <- as.data.frame(flowCore::exprs(flowframe_avant))
      
      total_evenements_avant <- nrow(donnees_globales)
      # Densité calculée par binning raster : on garde systématiquement l'intégralité des
      # événements (pas de sous-échantillonnage), le coût du binning restant indépendant de N.
      donnees_visu <- donnees_globales
      
      gate_polygone <- self$gate_viabilite[[nom_echantillon]]
      
      if (inherits(gate_polygone, "polygonGate")) {
        bound <- gate_polygone@boundaries
        if (is.null(dim(bound)) || ncol(bound) != 2) {
          stop("Le gate polygonal ne contient pas une matrice Nx2 de boundaries.")
        }
        coordonnees_gate <- data.frame(
          x = bound[, 1],
          y = bound[, 2]
        )
        noms_canaux <- colnames(bound)
      } else if (inherits(gate_polygone, "rectangleGate")) {
        mins <- gate_polygone@min
        maxs <- gate_polygone@max
        noms_canaux <- names(mins)
        coordonnees_gate <- data.frame(
          x = c(mins[1], maxs[1], maxs[1], mins[1]),
          y = c(mins[2], mins[2], maxs[2], maxs[2])
        )
      } else if (inherits(gate_polygone, "filterResult")) {
        gate_interne <- gate_polygone@filter
        if (inherits(gate_interne, "polygonGate")) {
          bound <- gate_interne@boundaries
          coordonnees_gate <- data.frame(
            x = bound[, 1],
            y = bound[, 2]
          )
          noms_canaux <- colnames(bound)
        } else if (inherits(gate_interne, "rectangleGate")) {
          mins <- gate_interne@min
          maxs <- gate_interne@max
          noms_canaux <- names(mins)
          coordonnees_gate <- data.frame(
            x = c(mins[1], maxs[1], maxs[1], mins[1]),
            y = c(mins[2], mins[2], maxs[2], maxs[2])
          )
        } else {
          stop("Type de gate non supporté : ", class(gate_interne))
        }
      } else {
        stop("Type de gate non supporté : ", class(gate_polygone))
      }
      
      canal_x <- as.character(noms_canaux[1])
      canal_y <- as.character(noms_canaux[2])
      
      
      if (!(canal_x %in% colnames(donnees_visu))) {
        stop("Canal X introuvable dans les données : ", canal_x,
             "\nColonnes disponibles : ", paste(colnames(donnees_visu), collapse=", "))
      }
      if (!(canal_y %in% colnames(donnees_visu))) {
        stop("Canal Y introuvable dans les données : ", canal_y,
             "\nColonnes disponibles : ", paste(colnames(donnees_visu), collapse=", "))
      }
      
      flowframe_apres <- self$post_viabilite[[nom_echantillon]]
      total_evenements_apres <- nrow(flowCore::exprs(flowframe_apres))
      pourcentage_conservation <- round((total_evenements_apres / total_evenements_avant) * 100, 1)
      
      lbl_x <- if (!is.null(self$get_label)) self$get_label(flowframe_avant, canal_x) else canal_x
      lbl_y <- if (!is.null(self$get_label)) self$get_label(flowframe_avant, canal_y) else canal_y
      
      lim_x_viab <- range(donnees_globales[[canal_x]], na.rm = TRUE)
      lim_y_viab <- range(donnees_globales[[canal_y]], na.rm = TRUE)
      
      # Densité calculée par binning raster (rapide, indépendant du nombre d'événements)
      # plutôt que par ggpointdensity (densité par point, coûteuse avec beaucoup d'événements).
      df_densite_viab <- calculer_densite_raster(donnees_visu[[canal_x]], donnees_visu[[canal_y]], lim_x_viab, lim_y_viab)
      
      graphique_viabilite <- if (is.null(df_densite_viab)) {
        ggplot2::ggplot() + ggplot2::theme_bw() +
          ggplot2::labs(title = paste("Retrait des cellules mortes :", nom_echantillon), subtitle = "Pas assez d'événements pour tracer la densité", x = lbl_x, y = lbl_y)
      } else {
        ggplot2::ggplot(
          df_densite_viab,
          ggplot2::aes(x = X, y = Y, fill = densite)
        ) +
          ggplot2::geom_raster(interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(
            colours = PALETTE_DENSITE
          ) +
          ggplot2::geom_polygon(
            data = coordonnees_gate,
            ggplot2::aes(x = x, y = y),
            fill = NA, color = "grey20", linewidth = 0.6, inherit.aes = FALSE
          ) +
          ggplot2::coord_cartesian(
            xlim = lim_x_viab,
            ylim = lim_y_viab
          ) +
          ggplot2::theme_bw() +
          ggplot2::labs(
            title = paste("Retrait des cellules mortes :", nom_echantillon),
            subtitle = paste0(
              "Cellules vivantes : ", format(total_evenements_apres, big.mark = " "),
              " | ", pourcentage_conservation, "%"
            ),
            x = lbl_x, y = lbl_y
          ) +
          ggplot2::theme(
            legend.position = "none",
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(color = "darkblue", size = 11)
          )
      }
      
      if (is.null(self$plots_viabilite)) self$plots_viabilite <- list()
      self$plots_viabilite[[nom_echantillon]] <- graphique_viabilite
      
      return(graphique_viabilite)
    },
    
    
    transformation_arcsinh = function(canaux = "", echantillon = NULL, cofactor = 400) {
      
      noms_a_traiter <- if (is.null(echantillon)) names(self$echantillons_traites) else echantillon
      
      arc_sinh_transform <- function(x, cf) {
        asinh(x / cf)
      }
      
      self$cofactor_transformation <- cofactor
      
      for (nom in noms_a_traiter) {
        ff <- if (!is.null(self$post_doublets_final[[nom]])) {
          self$post_doublets_final[[nom]]
        } else if (!is.null(self$post_debris[[nom]])) {
          self$post_debris[[nom]]
        } else {
          self$echantillons_traites[[nom]]
        }
        
        if (is.null(ff)) next
        
        matrice_exprs <- flowCore::exprs(ff)
        
        canaux_a_transformer <- if (canaux[1] == "") colnames(matrice_exprs) else canaux
        canaux_existants <- intersect(canaux_a_transformer, colnames(matrice_exprs))
        
        if (length(canaux_existants) == 0) next
        
        matrice_exprs[, canaux_existants] <- arc_sinh_transform(matrice_exprs[, canaux_existants], cofactor)
        
        ff_transforme <- ff
        flowCore::exprs(ff_transforme) <- matrice_exprs
        
        if (is.null(self$post_transformation)) self$post_transformation <- list()
        self$post_transformation[[nom]] <- ff_transforme
      }
      
      invisible(self)
    },
    
    creer_gate = function(nom_gate, type = "polygon", axes = c("FSC-A", "SSC-A"), points = NULL) {
      # Validation de la structure de stockage
      if (is.null(self$gates_personnalisees)) self$gates_personnalisees <- list()
      
      # Création de la porte selon le type
      if (type == "polygon") {
        if (is.null(points) || nrow(points) < 3) stop("Un polygone nécessite au moins 3 points.")
        gate <- flowCore::polygonGate(.gate = as.matrix(points[, 1:2]), filterId = nom_gate)
        
      } else if (type == "rectangle") {
        if (is.null(points) || length(points) != 4) stop("Le rectangle nécessite c(minX, maxX, minY, maxY).")
        gate <- flowCore::rectangleGate(filterId = nom_gate, 
                                        .gate = matrix(c(points[1], points[2], points[3], points[4]), 
                                                       ncol = 2, byrow = TRUE, 
                                                       dimnames = list(NULL, axes)))
      } else {
        stop("Type de gate invalide : utilisez 'polygon' ou 'rectangle'.")
      }
      
      # Sauvegarde
      self$gates_personnalisees[[nom_gate]] <- list(gate = gate, axes = axes, type = type)
      
      message("Succès : Gate '", nom_gate, "' créée sur les axes ", paste(axes, collapse="/"), ".")
      return(invisible(self))
    },
    
    # ===========================================
    # SECTION EXPORT — PRÉTRAITEMENT
    # ===========================================
    
    exporter_fcs_pretraitement = function(noms_echantillons = "all", etapes = c("debris", "doublets", "viabilite"), dossier_export = ".") { # Méthode permettant d'écrire sur le disque les fichiers FCS aux différentes étapes du prétraitement (débris, doublets, cellules mortes), au choix de l'utilisateur
      etapes <- intersect(etapes, c("debris", "doublets", "viabilite")) # Filtre les étapes demandées pour ne garder que les valeurs reconnues
      if (length(etapes) == 0) {
        stop("Aucune étape de prétraitement valide sélectionnée (attendu : 'debris', 'doublets' et/ou 'viabilite').")
      }
      
      liste_sources <- list(debris = self$post_debris, doublets = self$post_doublets_final, viabilite = self$post_viabilite) # Regroupe les trois listes de résultats de prétraitement disponibles pour un accès homogène
      suffixes      <- list(debris = "_post_debris", doublets = "_post_doublets", viabilite = "_post_viabilite") # Associe à chaque étape le suffixe de nommage utilisé pour distinguer les fichiers exportés
      
      noms_disponibles <- unique(unlist(lapply(etapes, function(e) names(liste_sources[[e]])))) # Recense l'ensemble des échantillons disponibles, toutes étapes sélectionnées confondues
      if (length(noms_disponibles) == 0) {
        stop("Aucun résultat de prétraitement disponible : exécutez d'abord le retrait des débris, des doublets et/ou des cellules mortes.")
      }
      
      tubes_a_exporter <- if (length(noms_echantillons) == 1 && noms_echantillons == "all") { # Si l'utilisateur souhaite exporter la totalité des échantillons traités
        noms_disponibles # Sélectionne l'intégralité des échantillons disponibles pour les étapes demandées
      } else { # Sinon, si une liste restreinte de noms a été fournie par l'utilisateur
        intersect(noms_echantillons, noms_disponibles) # Identifie par intersection les échantillons demandés qui existent réellement en mémoire
      }
      if (length(tubes_a_exporter) == 0) {
        stop("Aucun des échantillons spécifiés n'a été trouvé dans les résultats de prétraitement.")
      }
      
      if (!dir.exists(dossier_export)) { # Vérifie si le dossier de destination spécifié n'existe pas encore physiquement sur le disque
        dir.create(dossier_export, recursive = TRUE) # Crée automatiquement l'arborescence des dossiers manquants pour éviter une erreur d'écriture
      }
      
      fichiers_ecrits <- character(0) # Accumule au fil de la boucle les chemins des fichiers FCS effectivement écrits sur le disque
      
      for (nom in tubes_a_exporter) { # Boucle de traitement itérative pour chaque échantillon sélectionné
        for (etape in etapes) { # Boucle interne pour chaque étape de prétraitement demandée (débris, doublets et/ou viabilité)
          fcs_obj <- liste_sources[[etape]][[nom]] # Extrait l'objet flowFrame nettoyé correspondant à cet échantillon et cette étape
          if (is.null(fcs_obj)) next # Sécurité : passe au suivant si cet échantillon n'a pas de résultat pour cette étape précise
          
          nom_fichier_propre <- paste0(gsub("[^a-zA-Z0-9_]", "_", nom), suffixes[[etape]], ".fcs") # Construit un nom de fichier sûr, distinguant l'étape d'origine du nettoyage
          chemin_fcs <- file.path(dossier_export, nom_fichier_propre) # Concatène le chemin du dossier et le nom du fichier pour obtenir l'adresse d'écriture finale
          flowCore::write.FCS(fcs_obj, filename = chemin_fcs) # Enregistre physiquement l'objet flowFrame nettoyé au format binaire FCS standard sur le disque dur
          fichiers_ecrits <- c(fichiers_ecrits, chemin_fcs) # Ajoute le chemin du fichier fraîchement écrit à la liste de suivi
        }
      }
      
      message(paste(length(fichiers_ecrits), "fichier(s) FCS post-prétraitement exporté(s) dans :", dossier_export)) # Affiche un message de confirmation récapitulatif dans la console de commande
      invisible(fichiers_ecrits) # Renvoie de manière invisible la liste des chemins écrits (utile pour zipper ensuite depuis Shiny)
    },
    
    sauvegarder_session_pretraitement_rds = function(nom_fichier = "Pretraitement_Session_Complete.rds") { # Méthode permettant de sauvegarder l'intégralité des paramètres et résultats du prétraitement (bordures, débris, doublets, viabilité) dans un fichier binaire R (.rds)
      sauvegarde <- list( # Initialise une structure de liste imbriquée pour regrouper de manière organisée tous les éléments à archiver
        meta = list( # Sous-liste dédiée aux informations générales et de traçabilité de l'expérience
          date_export  = Sys.time(), # Enregistre l'horodatage exact (date et heure) de la création de la sauvegarde
          echantillons = names(self$echantillons_traites) # Mémorise la liste des échantillons présents dans la cohorte au moment de l'export
        ),
        retrait_bordures = list( # Sous-liste dédiée aux réglages utilisés pour le retrait des événements saturés (Margins)
          canaux_bordures = self$canaux_bordures # Archive les canaux/détecteurs ciblés lors du retrait des bordures
        ),
        debris = list( # Sous-liste dédiée au filtrage des débris
          gates = self$gate_debris # Archive les polygones de sélection des cellules (hors débris) pour chaque échantillon
        ),
        doublets = list( # Sous-liste dédiée au filtrage des doublets
          gates_FSC = self$gate_doublets_FSC, # Archive les seuils/polygones de discrimination des doublets sur l'axe FSC
          gates_SSC = self$gate_doublets_SSC # Archive les seuils/polygones de discrimination des doublets sur l'axe SSC
        ),
        viabilite = list( # Sous-liste dédiée au retrait des cellules mortes
          gates           = self$gate_viabilite, # Archive les polygones de sélection des cellules vivantes
          transformation  = self$cofactor_transformation # Archive le cofacteur Arcsinh utilisé pour la transformation du marqueur de viabilité
        ),
        visualisations = list( # Sous-liste dédiée à l'archivage des rendus graphiques produits pour l'assurance qualité
          plots_debris    = self$plots_debris, # Archive l'historique des figures de densité illustrant le retrait des débris
          plots_doublets  = self$plots_doublets, # Archive l'historique des figures illustrant le retrait des doublets (FSC et SSC)
          plots_viabilite = self$plots_viabilite # Archive l'historique des figures illustrant le retrait des cellules mortes
        )
      )
      
      saveRDS(sauvegarde, file = nom_fichier) # Sérialise et enregistre l'objet liste complet sous forme de fichier binaire compressé (.rds) sur le stockage local
      message(paste("Session de prétraitement et paramètres sauvegardés avec succès dans :", nom_fichier)) # Génère un message de confirmation explicite au sein de la console R de commande
      invisible(nom_fichier) # Renvoie de manière invisible le chemin du fichier généré, prêt à être proposé au téléchargement depuis Shiny
    }
    
    
  ),
  
  private = list(df_control_file = NULL)
)