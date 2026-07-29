library(R6)
library(ggplot2)
library(flowCore)
library(base64enc)

source("R/utils.R")

CARROT <- R6Class(
  classname = "CARROT",
  public = list(
    # Importation des fichiers
    mode = "Conventionnel", # l'utilisateur précise si ses données sont issues d'un cytomètre conventionnel ou spectral
    mapping_canal_fichier = NULL, # correspondance canal -> nom de fichier saisie par l'utilisateur à l'import (non systématiquement utilisée selon le flux d'import choisi)
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
    dossier_monomarques = NULL, # dossier source contenant les fichiers monomarqués/unstained originaux (avant filtrage dans dossier_controles_asp)
    dossier_echantillons = NULL, # dossier source contenant les fichiers échantillons à unmixer
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
    post_flowAI = list(), # Contient les données des fichiers FCS après le contrôle qualité alternatif flowAI (vérification du débit, de la lueur et de la stabilité)
    rapports_flowai = list(), # Contient, par échantillon, le chemin du dossier temporaire regroupant le rapport natif flowAI (HTML, TXT, PNG) généré par flow_auto_qc
    parametres_peacoqc_utilises = NULL, # Mémorise les derniers réglages PeacoQC appliqués à la cohorte (pour affichage dans le résumé PDF)
    parametres_flowai_utilises = NULL, # Mémorise les derniers réglages flowAI appliqués à la cohorte (pour affichage dans le résumé PDF)
    post_retrait_bordures = list(), # Stocke la matrice d'expression des échantillons nettoyée des signaux saturés (valeurs maximales ou minimales des détecteurs)
    gate_debris = list(), # Contient les coordonnées et les structures géométriques des fenêtres (gates) de sélection des cellules (retrait des débris en FSC vs SSC)
    post_debris = list(), # Stocke les données des échantillons filtrées où seuls les événements correspondants aux cellules (hors débris) ont été conservés
    gates_history = list(), # Historique ordonné des gates nommés appliqués : list(nom_gate -> list(nom_echantillon -> list(polygone, cx, cy, post_data, n_avant, n_apres)))
    
    # ── SECTION ANALYSES : gates personnalisés et résultats des analyses avancées ──
    gates_personnalisees = list(), # Gates créés via creer_gate(), un par nom : list(nom_gate -> list(gate, axes, type)). Chaque appel avec un nouveau nom ajoute une entrée sans écraser les précédentes : l'utilisateur peut créer et conserver autant de gates qu'il le souhaite.
    analyses_umap        = list(), # Résultats de projection_UMAP(), indexés par nom de gate : list(nom_gate -> list(embedding, echantillon_origine, canaux, parametres))
    analyses_tsne        = list(), # Résultats de projection_tSNE(), indexés par nom de gate, même structure que analyses_umap
    analyses_pca         = list(), # Résultats de projection_PCA(), indexés par nom de gate : list(nom_gate -> list(embedding, echantillon_origine, canaux, variance_expliquee, rotation, parametres))
    clusters_flowsom     = list(), # Résultats de creer_clusters(), indexés par nom de gate : list(nom_gate -> list(fsom, clusters, metaclusters, echantillon_origine, canaux, parametres))
    groupes_echantillons = list(), # Assignation échantillon -> groupe/condition (definir_groupe_echantillon()), utilisée par comparer_groupes()
    
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
    cofactor_transformation = NULL, # Cofacteur Arcsinh utilisé lors du dernier appel à transformation_arcsinh() (archivé pour le résumé de session RDS)
    canaux_transformes = list(), # Liste des canaux ayant effectivement reçu la transformation, par échantillon
    
    # variable de rappel Shiny (callback) : permet à un module d'être notifié quand le pipeline change sans dépendance circulaire directe
    pipeline_callback = NULL,
    
    # config/dictionnaire marqueurs (pour Shiny)
    config_marqueurs = NULL, # data.frame brut (lignes = tubes/échantillons, colonnes = canaux) des annotations marqueurs saisies par l'utilisateur
    dictionnaire_marqueurs = NULL, # vecteur nommé canal -> marqueur biologique, dérivé de config_marqueurs, utilisé par obtenir_label() pour tous les libellés d'axes
    
    #pour l'import sans fichiers contrôles (échantillons déjà compensés / unmixés)
    sans_controles = FALSE, # TRUE si l'utilisateur n'a fourni aucun tube monomarqué/unstained
    deja_traite = FALSE, # TRUE si l'utilisateur indique que ses échantillons sont déjà compensés (conventionnel) ou unmixés (spectral)
    cofacteur_defaut = 150, # cofacteur Arcsinh utilisé par défaut pour construire trans_list automatiquement en l'absence de contrôles (ajustable ensuite via l'onglet Transformation)
    
    
    #initialisation de la classe
    initialize = function(df_monomarques = NULL, df_echantillons = NULL, chemin_racine = NULL, mode = "Conventionnel", deja_traite = FALSE) { # initialiser la classe R6 
      self$mode <- mode # Enregistre le mode utilisé (Conventionnel ou Spectral)
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
    
    # ============= CHARGER_FCS =====================
    charger_fcs = function() {
      if (!is.null(self$chemins_monomarques) && nrow(self$chemins_monomarques) > 0) { # Si au moins un tube contrôle a été fourni à l'import
        self$tubes_monomarques <- lapply(seq_len(nrow(self$chemins_monomarques)), function(i) { # Lit chaque fichier FCS contrôle un par un, dans l'ordre du tableau d'import
          row <- self$chemins_monomarques[i, ]
          flowCore::read.FCS(row$chemin, transformation = FALSE, truncate_max_range = FALSE) # transformation=FALSE : on garde les valeurs brutes, la transformation (Arcsinh) est appliquée plus tard et séparément ; truncate_max_range=FALSE : ne tronque pas les événements saturés en bord de plage (gérés ensuite par le retrait des bordures)
        })
        noms_tubes <- sapply(seq_len(nrow(self$chemins_monomarques)), function(i) { # Détermine le nom à donner à chaque tube dans self$tubes_monomarques
          row <- self$chemins_monomarques[i, ]
          if (!is.na(row$type) && row$type == "Unstained") return("TUBE_UNSTAINED") # Nom conventionnel fixe pour le tube non marqué, quel que soit son nom de fichier d'origine (repéré ainsi partout ailleurs dans le pipeline)
          return(row$canal) # Pour un tube monomarqué, le nom est directement le canal sur lequel il a été acquis (ex: "PE-A")
        })
        names(self$tubes_monomarques) <- noms_tubes
        self$sans_controles <- FALSE # Au moins un contrôle existe : le pipeline suivra le flux "avec contrôles" (compensation calculée depuis ces tubes)
      } else {
        self$tubes_monomarques <- list()
        self$sans_controles <- TRUE # Aucun contrôle : le pipeline suivra le flux "sans contrôles" (voir plus bas : transformation par défaut, éventuel import de spillover déjà présent dans le FCS)
      }
      
      self$echantillons = lapply(self$chemins_echantillons$chemin, function(f) { # échantillons biologiques
        flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE) # Mêmes réglages que pour les contrôles : lecture brute, sans transformation ni troncature
      }) 
      names(self$echantillons) <- self$chemins_echantillons$tube_name # Les échantillons sont indexés par leur nom de tube partout dans le pipeline (clé de toutes les listes post_*)
      
      if (self$sans_controles && self$mode == "Conventionnel" && self$deja_traite) { # Cas particulier : cytométrie conventionnelle, pas de contrôles fournis, mais données annoncées déjà compensées
        invisible(tryCatch(self$importer_spillover_fcs(), error = function(e) NULL)) # Tente de récupérer la matrice de spillover déjà appliquée, directement depuis les métadonnées du FCS (mot-clé $SPILLOVER ou équivalent) ; échec silencieux si absente (juste informatif, non bloquant)
      }
      
      if (self$sans_controles) { # Sans contrôles (quel que soit le mode) : il faut quand même pouvoir visualiser/gater les données, donc une transformation par défaut est appliquée
        invisible(tryCatch(self$transformer_fcs(cofacteur = self$cofacteur_defaut),
                           error = function(e) NULL)) # Échec silencieux si la transformation par défaut ne peut pas s'appliquer (ex: canaux non standards) ; l'utilisateur pourra la relancer manuellement depuis l'onglet Transformation
      }
      
      if (isTRUE(self$deja_traite)) { # Si l'utilisateur a répondu "Oui" à "Données déjà compensées/démixées ?" , les fichiers importés sont, par définition, déjà dans leur état final
        self$echantillons_traites <- self$echantillons # on copie les échantillons tel quel dans self$echantillons_traites => pour pouvoir passer directement aux étapes de QC et prétraitement
      }
    },
    
    # ============= importer_spillover_fcs =====================
    importer_spillover_fcs = function(nom_echantillon = NULL) { # méthode qui permet de récupérer la matrice de spillover (via métadonnées) d'un échantillon déjà compensé
      if (is.null(self$echantillons) || length(self$echantillons) == 0) {
        stop("Aucun échantillon chargé.")
      }
      nom <- if (!is.null(nom_echantillon)) nom_echantillon else names(self$echantillons)[1] # Par défaut, inspecte le premier échantillon disponible (la matrice de spillover est en principe identique pour toute la cohorte acquise dans la même session)
      fcs <- self$echantillons[[nom]]
      if (is.null(fcs)) stop("Échantillon introuvable : ", nom)
      
      mat <- extraire_spillover_depuis_fcs(fcs) # Fonction utilitaire (utils.R) : cherche et parse tous les formats connus de mot-clé de compensation dans les métadonnées FCS
      
      if (is.null(mat)) { # Aucune matrice exploitable trouvée : construit un message d'erreur informatif plutôt que de simplement échouer
        cles_trouvees <- lister_cles_spillover_fcs(fcs) # Recherche large (toute clé contenant "spill"/"comp") 
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
      
      canaux_valides <- intersect(rownames(mat), flowCore::colnames(fcs)) # # On ne conserve que les canaux réellement présents dans l'échantillon
      if (length(canaux_valides) == 0) {
        return(list(succes = FALSE,
                    message = paste0("Une matrice a été trouvée (canaux : ",
                                     paste(rownames(mat), collapse = ", "),
                                     ") mais aucun de ces canaux ne correspond aux paramètres de l'échantillon (",
                                     paste(flowCore::colnames(fcs), collapse = ", "), ").")))
      }
      
      mat <- mat[canaux_valides, canaux_valides, drop = FALSE] # Restreint la matrice (carrée) aux seuls canaux valides, dans le même ordre pour lignes et colonnes
      self$S_matrix <- mat # Stocke la matrice pour affichage/export (informatif uniquement, cf. en-tête de méthode)
      self$canaux   <- canaux_valides
      return(list(succes = TRUE,
                  message = paste0("Matrice importée pour ", length(canaux_valides), " canaux.")))
    },
    
    # ============= obtenir_label =====================
    obtenir_label = function(fcs, canal) { # Permet d'extraire le nom des marqueurs biologiques à afficher sur les axes des graphiques
      if (is.na(canal) || canal == "") return(canal) # Si aucun canal n'est donné, on retourne tel quel
      
      if (!is.null(self$dictionnaire_marqueurs) && canal %in% names(self$dictionnaire_marqueurs)) {
        marqueur_manuel <- self$dictionnaire_marqueurs[[canal]] # Priorité à l'annotation manuelle de l'utilisateur (onglet "Configuration Marqueurs")
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
    
    # ============= definir_config_marqueurs =====================
    definir_config_marqueurs = function(df_config) { 
      self$config_marqueurs <- df_config  # Enregistre la table de configuration des marqueurs dans l’objet R6
      
      if (is.null(df_config) || nrow(df_config) == 0 || ncol(df_config) == 0) {  # Si la configuration est vide ou inexistante, on réinitialise le dictionnaire
        self$dictionnaire_marqueurs <- NULL
        return(invisible(NULL))  # Sortie silencieuse
      }
      
      dict <- sapply(colnames(df_config), function(cx) { # Pour chaque colonne (chaque marqueur), on extrait les valeurs textuelles non vides
        valeurs <- trimws(as.character(df_config[[cx]]))  # Nettoie les espaces et convertit en texte
        valeurs <- valeurs[!is.na(valeurs) & nchar(valeurs) > 0]  # Filtre les valeurs non vides
        if (length(valeurs) > 0) valeurs[1] else NA_character_  # Retient uniquement la première valeur non vide comme "nom" du marqueur
      })
      
      dict <- dict[!is.na(dict)]  # Supprime les entrées NA
      self$dictionnaire_marqueurs <- as.list(dict)  # Stocke le dictionnaire sous forme de liste nommée
      
      invisible(self$dictionnaire_marqueurs)  # Retour silencieux de la liste
    },
    
    # ============= mettre_a_jour_pipeline =====================
    mettre_a_jour_pipeline = function(etape, nom_echantillon = NULL) { # Méthode permettant d’enregistrer l’avancement du pipeline pour une étape donnée
      
      if (is.null(self$pipeline_status)) self$pipeline_status <- list()  # Initialise la structure de suivi si elle n’existe pas encore
      horodatage <- format(Sys.time(), "%H:%M:%S")   # Capture l’heure courante pour tracer l’exécution
      
      if (is.null(nom_echantillon)) {  # Cas où l’étape concerne l’ensemble des échantillons
        self$pipeline_status[[etape]] <- list(
          statut = "Terminé",
          date = horodatage,
          portee = "Tous les échantillons"
        )
        
      } else { # Cas où l’étape concerne un échantillon spécifique
        if (is.null(self$pipeline_status[[etape]])) self$pipeline_status[[etape]] <- list()
        
        self$pipeline_status[[etape]][[nom_echantillon]] <- list(
          statut = "Terminé",
          date = horodatage
        )
      }
      
      if (is.function(self$pipeline_callback)) { # Si un callback utilisateur est défini, on le déclenche pour notifier l’avancement
        self$pipeline_callback(etape, nom_echantillon)
      }
    },
    
    # ===========================
    # COMPENSATION (Conventionnel)
    # ===========================
    
    # ============= transformer_fcs =====================
    transformer_fcs = function(cofacteur) { # Fonction appliquant une transformation arcsinh aux canaux fluorescents
      
      avec_controles <- !is.null(self$tubes_monomarques) && # Détermine si des tubes monomarqués sont disponibles
        length(self$tubes_monomarques) > 0 # (contrôles nécessaires pour certaines étapes du pipeline)
      
      if (avec_controles) { # Si des contrôles sont présents
        tous_canaux <- flowCore::colnames(self$tubes_monomarques[[1]]) # Récupère les noms de canaux du premier tube monomarqué
      } else if (!is.null(self$echantillons) &&   # Sinon, si des échantillons sont chargés
                 length(self$echantillons) > 0) { # (cas sans contrôles : données déjà compensées)
        tous_canaux <- flowCore::colnames(self$echantillons[[1]])  # Récupère les canaux depuis le premier échantillon
      } else {
        stop("Aucun tube monomarqué ni échantillon chargé.")  # Erreur si aucune source de canaux n’est disponible
      }
      
      canaux_fluo <- tous_canaux[!grepl("fsc|ssc|time", # Filtre les canaux fluorescents (exclut FSC/SSC/Time)
                                        tous_canaux, ignore.case = TRUE)] # via une expression régulière
      
      arsinh_fun <- flowCore::arcsinhTransform(a = 0,   # Crée une transformation arcsinh standard
                                               b = 1/cofacteur,  # Paramètre d’échelle (cofacteur)
                                               c = 0)  # Décalage nul
      
      if (avec_controles) {  # Si des tubes monomarqués sont disponibles
        self$monomarques_trans <- lapply(self$tubes_monomarques, function(fcs) {  # Applique la transformation à chaque tube monomarqué
          
          canaux_presents <- intersect(canaux_fluo, # Identifie les canaux fluorescents présents dans ce FCS
                                       flowCore::colnames(fcs))  # (sécurité si certains canaux manquent)
          
          if (length(canaux_presents) > 0) { # Si au moins un canal fluo est disponible
            funs_locales <- lapply(seq_along(canaux_presents),# Crée une liste de transformations arcsinh
                                   function(x) arsinh_fun)   # (une par canal)
            
            local_trans  <- flowCore::transformList( # Construit un transformList pour ces canaux
              from = canaux_presents,   # Canaux d’entrée
              tfun = funs_locales,  # Transformations arcsinh
              to   = canaux_presents # Canaux de sortie (identiques)
            )
            
            return(flowCore::transform(fcs, local_trans)) # Applique la transformation au FCS
          } else {
            return(fcs) # Si aucun canal fluo : renvoie le FCS inchangé
          }
        })
        
        names(self$monomarques_trans) <- names(self$tubes_monomarques)  # Conserve les noms des tubes monomarqués
      } else {
        self$monomarques_trans <- NULL # Pas de contrôles → pas de transformation monomarquée
      }
      
      funs_globales <- lapply(seq_along(canaux_fluo),  # Crée une liste de transformations arcsinh globales
                              function(x) arsinh_fun) # (pour tous les canaux fluorescents)
      
      self$trans_list <- flowCore::transformList(  # Enregistre la transformation globale
        from = canaux_fluo,  # Canaux fluorescents
        tfun = funs_globales,  # Transformations arcsinh
        to   = canaux_fluo # Canaux de sortie identiques
      )
    }, 
    
    # ============= definir_et_extraire =====================
    definir_et_extraire = function(nom_canal, intervalle_gate_negatif, intervalle_gate_positif, utiliser_unstained = TRUE) {   # Définit les gates négatif/positif et extrait les événements correspondants
      
      if (!nom_canal %in% self$canaux) stop("Le canal spécifié n'existe pas.") # Vérifie que le canal demandé existe dans la configuration
      
      existe_unstained <- "TUBE_UNSTAINED" %in% names(self$monomarques_trans) # Vérifie si un tube unstained transformé est disponible
      nom_tube_neg <- if(utiliser_unstained && existe_unstained) "TUBE_UNSTAINED" else nom_canal # Choisit la source négative : unstained si dispo, sinon tube du canal
      self$source_neg_utilisee[[nom_canal]] <- nom_tube_neg  # Enregistre la source utilisée pour le négatif
      
      source_trans_neg <- self$monomarques_trans[[nom_tube_neg]]   # Récupère le FCS transformé (arcsinh) pour le négatif
      source_brute_neg <- self$tubes_monomarques[[nom_tube_neg]]  # Récupère le FCS brut correspondant
      
      limites_negatif <- setNames(list(intervalle_gate_negatif), nom_canal) # Crée une liste nommée définissant les bornes du gate négatif
      gate_negatif    <- flowCore::rectangleGate(filterId = paste0("Gate_Negatif_", nom_canal), .gate = limites_negatif) # Construit un rectangleGate pour le négatif
      self$bornes_gates_neg[[nom_canal]] <- gate_negatif # Stocke le gate négatif dans l’objet
      
      garde_evts_du_gate_negatif <- flowCore::filter(source_trans_neg, gate_negatif) # Applique le gate négatif sur les données transformées
      self$gates_negatifs[[nom_canal]] <- source_brute_neg[garde_evts_du_gate_negatif@subSet, ]  # Extrait les événements bruts correspondant au gate négatif
      
      if (!is.null(intervalle_gate_positif)) {  # Si un gate positif est fourni
        limites_positif <- setNames(list(intervalle_gate_positif), nom_canal)  # Définition des bornes du gate positif
        gate_positif    <- flowCore::rectangleGate(filterId = paste0("Gate_Positif_", nom_canal), .gate = limites_positif)   # Construction du gate positif
        self$bornes_gates_pos[[nom_canal]] <- gate_positif # Stockage du gate positif
        
        garde_evts_du_gate_positif <- flowCore::filter(self$monomarques_trans[[nom_canal]], gate_positif) # Filtrage du tube monomarqué transformé pour le positif
        self$gates_positifs[[nom_canal]]   <- self$tubes_monomarques[[nom_canal]][garde_evts_du_gate_positif@subSet, ] # Extraction des événements bruts correspondant au gate positif
      }
    },
    
    # ============= graphiques_gates =====================
    graphiques_gates = function(nom_canal = NULL, shiny_neg = NULL, shiny_pos = NULL, afficher_unstained_neg = TRUE) {   # Génère les graphiques de densité pour visualiser les gates négatif/positif
      canaux_a_generer <- if (is.null(nom_canal)) self$canaux else nom_canal # Si aucun canal spécifié → tracer tous les canaux
      
      plots <- lapply(canaux_a_generer, function(canal) {  # Boucle sur chaque canal à tracer
        existe_unstained <- "TUBE_UNSTAINED" %in% names(self$monomarques_trans) # Vérifie si un tube unstained transformé est disponible
        tube_neg_a_tracer <- if(afficher_unstained_neg && existe_unstained) "TUBE_UNSTAINED" else canal  # Choisit la source négative à afficher (unstained ou monomarqué)
        fcs_brut  <- self$tubes_monomarques[[canal]] # Récupère le FCS brut du canal
        nom_bio   <- self$obtenir_label(fcs_brut, canal)  # Récupère le nom biologique du marqueur (ex: CD3, CD19)
        df_trans_pos  <- as.data.frame(flowCore::exprs(self$monomarques_trans[[canal]])) # Données transformées du tube positif
        df_trans_neg  <- as.data.frame(flowCore::exprs(self$monomarques_trans[[tube_neg_a_tracer]]))  # Données transformées du tube négatif
        lim_n <- if (!is.null(shiny_neg)) shiny_neg else c(0, 2)  # Limites du gate négatif (ou valeurs par défaut)
        lim_p <- if (!is.null(shiny_pos)) shiny_pos else c(4, 7) # Limites du gate positif (ou valeurs par défaut)
        
        pct_n <- round(sum(df_trans_neg[[canal]] >= lim_n[1] & df_trans_neg[[canal]] <= lim_n[2], na.rm=TRUE) /  # Pourcentage d’événements dans le gate négatif
                         nrow(df_trans_neg) * 100, 1)
        pct_p <- round(sum(df_trans_pos[[canal]] >= lim_p[1] & df_trans_pos[[canal]] <= lim_p[2], na.rm=TRUE) / # Pourcentage d’événements dans le gate positif
                         nrow(df_trans_pos) * 100, 1)
        
        p <- ggplot() +  # Initialise le graphique
          geom_density(data = df_trans_pos, aes(x = .data[[canal]], y = after_stat(count)),  # Courbe de densité du tube positif
                       fill = "#d90429", alpha = 0.4) +
          geom_vline(xintercept = lim_p, color = "#d90429", linetype = "solid", linewidth = 0.9) # Traits verticaux du gate positif
        + annotate("text", x = mean(lim_p), y = Inf, label = paste0("Pos: ", pct_p, "%"),  # Annotation du pourcentage positif
                   vjust = 1.5, color = "#d90429", fontface = "bold") +
          theme_bw() +
          labs(title = paste("Ajustement des Gates -", canal), x = nom_bio, y = "Nombre d'événements") # Titres et labels
        
        if(tube_neg_a_tracer == "TUBE_UNSTAINED") { # Cas où le négatif affiché est l’unstained
          
          p <- p + geom_density(data = df_trans_neg, aes(x = .data[[canal]], y = after_stat(count)),  # Courbe de densité du tube unstained
                                fill = "#0077b6", alpha = 0.4) +
            geom_vline(xintercept = lim_n, color = "#0077b6", linetype = "dashed", linewidth = 0.9)  # Traits verticaux du gate négatif
          + annotate("text", x = mean(lim_n), y = Inf, label = paste0("Unstained: ", pct_n, "%"),  # Annotation du pourcentage unstained
                     vjust = 3, color = "#0077b6", fontface = "bold") +
            labs(subtitle = "Bleu : Tube Unstained | Rouge : Tube Monomarqué")# Sous-titre explicatif
          
        } else { # Cas où le négatif affiché est le tube monomarqué
          
          p <- p + geom_vline(xintercept = lim_n, color = "#0077b6", linetype = "dashed", linewidth = 0.9) # Traits du gate négatif
          + annotate("text", x = mean(lim_n), y = Inf, label = paste0("Neg: ", pct_n, "%"),   # Annotation du pourcentage négatif interne
                     vjust = 3, color = "#0077b6", fontface = "bold") +
            labs(subtitle = "Distribution du tube monomarqué (Négatif interne)")  # Sous-titre explicatif
        }
        
        return(p)  # Retourne le graphique pour ce canal
      })
      
      names(plots) <- canaux_a_generer  # Nomme les graphiques selon les canaux
      if (!is.null(nom_canal)) return(plots[[nom_canal]])  # Si un canal est demandé → retourne uniquement celui-ci
      return(plots)  # Sinon → retourne la liste complète
    }, 
    
    # ============= calculer_spillover =====================
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
    
    # ============= modifier_spillover =====================
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
    
    # ============= comparer_medianes_spillover =====================
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
    
    # ============= compenser =====================
    compenser = function() { # Applique la compensation spectrale à tous les échantillons
      
      self$echantillons_traites <- lapply(names(self$echantillons), function(nom) { # Boucle sur chaque échantillon à compenser
        fcs <- self$echantillons[[nom]] # Récupère le FCS brut de l’échantillon
        matrice_cible <- if (!is.null(self$S_matrices_par_echantillon[[nom]])) { # Vérifie si une matrice de spillover spécifique existe pour cet échantillon
          self$S_matrices_par_echantillon[[nom]] # Si oui → utiliser la matrice spécifique
        } else {
          self$S_matrix   # Sinon → utiliser la matrice globale calculée précédemment
        }
        
        canaux_communs <- intersect(colnames(matrice_cible), flowCore::colnames(fcs))  # Identifie les canaux présents à la fois dans la matrice et dans le FCS
        sub_matrix <- matrice_cible[canaux_communs, canaux_communs, drop = FALSE]  # Extrait la sous‑matrice correspondant aux canaux disponibles
        return(flowCore::compensate(fcs, sub_matrix))  # Applique la compensation flowCore sur l’échantillon
      })
      
      names(self$echantillons_traites) <- names(self$echantillons)  # Renomme la liste des échantillons traités pour correspondre aux originaux
    },
    
    # ============= visualiser_compensation =====================
    visualiser_compensation = function(nom_echantillon, canal_x, canal_y, affichage = "Both") { # Visualise la compensation entre deux canaux (ou toutes les paires)
      
      # Sécurité : canal_x/canal_y peuvent être NULL ou vides (ex: sélecteurs
      # Shiny pas encore peuplés au tout premier rendu, ou changement
      # d'échantillon en cours). Sans ce garde-fou, la comparaison "canal_x ==
      # 'ALL'" juste en dessous plante avec "argument is of length zero" dès
      # que canal_x est NULL, au lieu de simplement attendre une sélection
      # valide — même pattern défensif que dans controler_monomarques().
      if (is.null(canal_x) || is.null(canal_y) || canal_x == "" || canal_y == "") return(NULL)
      
      # Sécurité : la transformation Arcsinh doit avoir été appliquée au moins
      # une fois (onglet Transformation) avant de pouvoir visualiser quoi que
      # ce soit ici — sans ce contrôle, flowCore::transform() plus bas reçoit
      # un transformList NULL et échoue avec une erreur interne peu explicite.
      if (is.null(self$trans_list)) {
        stop("La transformation (Arcsinh) n'a pas encore été appliquée. Rendez-vous dans l'onglet Transformation avant de visualiser la compensation.")
      }
      
      local_trans <- self$trans_list   # Transformation arcsinh globale à appliquer
      fcs_source  <- self$echantillons[[nom_echantillon]] # FCS brut de l’échantillon
      
      # ───────────────────────────────────────────────
      # Cas "ALL" : toutes les combinaisons fluo
      # ───────────────────────────────────────────────
      if (canal_x == "ALL" || canal_y == "ALL") { # Mode vue d’ensemble : toutes les paires de canaux
        
        if (identical(affichage, "Both")) {  # Interdit d’afficher avant+après en vue d’ensemble
          stop("Le mode 'Vue d'ensemble' ne permet pas l'affichage 'Both' : choisissez 'Before compensation only' ou 'After compensation only'.")
        }
        
        tous_canaux  <- flowCore::colnames(fcs_source)  # Récupère tous les canaux du FCS
        canaux_fluo  <- tous_canaux[!grepl("FSC|SSC|Time", tous_canaux, ignore.case = TRUE)]   # Filtre les canaux fluorescents
        combinaisons <- expand.grid(x = canaux_fluo, y = canaux_fluo, stringsAsFactors = FALSE)  # Toutes les combinaisons possibles
        combinaisons <- combinaisons[combinaisons$x < combinaisons$y, ]  # Garde uniquement x<y pour éviter les doublons
        
        plots_list <- lapply(seq_len(nrow(combinaisons)), function(i) {  # Génère chaque plot via appel récursif
          self$visualiser_compensation(
            nom_echantillon,
            combinaisons$x[i],
            combinaisons$y[i],
            affichage = affichage
          )
        })
        
        plots_list <- Filter(Negate(is.null), plots_list)  # Supprime les NULL
        if (length(plots_list) == 0) return(NULL)  # Aucun plot → NULL
        
        return(structure(plots_list, class = c("carrot_plots_list", "list")))  # Retourne la liste structurée
      }
      
      # ───────────────────────────────────────────────
      # Cas simple : une paire de canaux
      # ───────────────────────────────────────────────
      label_x <- self$obtenir_label(fcs_source, canal_x)  # Nom biologique du canal X
      label_y <- self$obtenir_label(fcs_source, canal_y) # Nom biologique du canal Y
      
      extraire_matrice <- function(fcs, cx, cy) {   # Fonction interne : extrait les deux colonnes transformées
        if (is.null(fcs)) return(NULL)   # Sécurité : FCS absent
        if (!(cx %in% flowCore::colnames(fcs)) || !(cy %in% flowCore::colnames(fcs))) return(NULL) # Sécurité : canaux absents
        fcs_trans <- flowCore::transform(fcs, local_trans)  # Applique la transformation arcsinh
        mat  <- flowCore::exprs(fcs_trans)[, c(cx, cy), drop = FALSE]  # Extrait les deux colonnes
        as.data.frame(mat) # Retourne un data.frame
      }
      
      df_avant <- if (affichage %in% c("Both", "Before compensation only"))  # Données avant compensation
        extraire_matrice(self$echantillons[[nom_echantillon]], canal_x, canal_y) else NULL
      df_apres <- if (affichage %in% c("Both", "After compensation only")) # Données après compensation
        extraire_matrice(self$echantillons_traites[[nom_echantillon]], canal_x, canal_y) else NULL
      
      # Limites communes
      all_data <- rbind(df_avant, df_apres)# Combine pour définir les limites
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
      creer_plot_raster <- function(df, titre) {# Convertit les points en densité rasterisée
        if (is.null(df) || nrow(df) == 0) return(NULL)  
        colnames(df) <- c("X", "Y") 
        resolution <- 400 # Résolution de la grille
        x_breaks <- seq(lim_x[1], lim_x[2], length.out = resolution + 1)                       
        y_breaks <- seq(lim_y[1], lim_y[2], length.out = resolution + 1)                                         
        
        df_binned <- df |>
          dplyr::mutate(
            x_bin = cut(X, breaks = x_breaks, include.lowest = TRUE),                                                
            y_bin = cut(Y, breaks = y_breaks, include.lowest = TRUE)                                                
          ) |>
          dplyr::count(x_bin, y_bin, name = "densite") |>
          tidyr::drop_na()                                                                                         
        
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
            colours = PALETTE_DENSITE,
            values  = PALETTE_DENSITE_STOPS
          ) +
          ggplot2::coord_cartesian(xlim = lim_x, ylim = lim_y) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
          ggplot2::labs(title = titre, x = label_x, y = label_y)
        
      }
      
      plot_avant <- creer_plot_raster(df_avant, paste(nom_echantillon, "- Avant"))  # Plot avant compensation
      plot_apres <- creer_plot_raster(df_apres, paste(nom_echantillon, "- Après"))  # Plot après compensation
      
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
    
    # ============= exporter_fcs_compenses =====================
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
    
    # ============= sauvegarder_session_rds =====================
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
    
    # ============= lancer_asp =====================
    lancer_asp = function(type_cytometre = "aurora") { # Lance la génération du fichier de contrôle AutoSpectral pour le cytomètre choisi
      
      if (self$mode != "Spectral") {# Vérifie que le pipeline est en mode spectral
        stop("Cette méthode nécessite le mode 'Spectral'.") # Stoppe si ce n’est pas le cas
      }
      
      self$asp_config <- AutoSpectral::get.autospectral.param(cytometer = type_cytometre)  # Récupère les paramètres AutoSpectral adaptés au cytomètre (Aurora, etc.)
      
      control_dir <- path.expand(dirname(self$chemins_monomarques$chemin[1]))# Détermine le dossier contenant les tubes monomarqués (répertoire de contrôle)
      fichier_csv <- file.path(path.expand(self$dossier_racine), "fcs_control_file") # Chemin du fichier CSV de contrôle à générer
      
      if (file.exists(paste0(fichier_csv, ".csv"))) { # Si un ancien fichier existe déjà
        file.remove(paste0(fichier_csv, ".csv"))  # On le supprime pour éviter les conflits
      }
      
      AutoSpectral::create.control.file(# Génère le fichier de contrôle AutoSpectral
        control.dir = control_dir,  # Dossier contenant les FCS monomarqués
        asp = self$asp_config, # Paramètres AutoSpectral
        filename = fichier_csv  # Nom du fichier CSV à créer
      )
    },
    
    # ============= verifier_asp =====================
    verifier_asp = function(seuil_warning = 5000, seuil_error = 1000) { # Vérifie la qualité des fichiers de contrôle AutoSpectral (événements min. requis)
      
      if (is.null(self$asp_config)) {  # Vérifie que lancer_asp() a été exécuté
        stop("Erreur : asp_config est NULL. Lancez d'abord lancer_asp().") # Stoppe si la configuration ASP n’est pas disponible
      }
      
      chemin_csv_complet <- file.path(path.expand(self$dossier_racine), "fcs_control_file.csv") # Chemin complet du fichier CSV de contrôle généré par lancer_asp()
      dossier_fcs <- path.expand(dirname(self$chemins_monomarques$chemin[1])) # Déduit le dossier contenant les FCS monomarqués
      
      if (!file.exists(chemin_csv_complet)) {  # Vérifie que le fichier CSV existe
        stop("Fichier de contrôle introuvable à : ", chemin_csv_complet)  # Stoppe si absent
      }
      
      verification <- AutoSpectral::check.control.file( # Appelle la fonction de vérification AutoSpectral
        control.dir       = dossier_fcs,   # Dossier contenant les FCS monomarqués
        control.def.file  = chemin_csv_complet,   # Fichier CSV de définition des contrôles
        asp               = self$asp_config, # Paramètres ASP du cytomètre
        min.event.warning = seuil_warning,  # Seuil d’avertissement (événements insuffisants)
        min.event.error   = seuil_error  # Seuil d’erreur (événements très insuffisants)
      )
      
      if (is.null(verification)) { # Si aucun tableau n’est retourné
      } else {
        print(verification)   # Affiche le tableau de vérification dans la console R
      }
      
      return(verification)  # Retourne le résultat (NULL ou data.frame)
    },
    
    # ============= definir_gates_landmarks (cellules + beads) =====================
    definir_gates_landmarks = function(control_name, n.cells = 2000, percentile = 70, 
                                       grid.n = 100, bandwidth.factor = 1, 
                                       fsc.channel = NULL, ssc.channel = NULL) {
      
      if (is.null(self$asp_config)) {
        stop("La configuration ASP n'est pas initialisée. Appelez d'abord lancer_asp().")
      }   # Vérifie que la configuration AutoSpectral est chargée
      
      old_wd <- getwd()
      setwd(path.expand(self$dossier_racine))
      on.exit(setwd(old_wd))   # Change temporairement le répertoire de travail et garantit le retour à l’ancien
      
      output_dir <- file.path(self$dossier_racine, "figure_gate")
      if (!dir.exists(output_dir)) dir.create(output_dir)   # Crée le dossier de sortie pour les figures si nécessaire
      
      gate_result <- AutoSpectral::define.gate.landmarks(
        control.file = "fcs_control_file.csv",   # Fichier de définition des contrôles ASP
        control.dir = path.expand(dirname(self$chemins_monomarques$chemin[1])),   # Dossier contenant les FCS monomarqués
        asp = self$asp_config,   # Paramètres ASP du cytomètre
        gate.name = control_name,   # Nom du gate à générer
        n.cells = n.cells,   # Nombre d’événements utilisés pour estimer les landmarks
        percentile = percentile,   # Percentile utilisé pour définir les limites du gate
        grid.n = grid.n,   # Taille de la grille pour l’estimation de densité
        bandwidth.factor = bandwidth.factor,   # Facteur d’ajustement du lissage
        fsc.channel = fsc.channel,   # Canal FSC si spécifié
        ssc.channel = ssc.channel,   # Canal SSC si spécifié
        output.dir = output_dir   # Dossier où les figures seront enregistrées
      )   # Appelle AutoSpectral pour calculer les gates landmarks
      
      if (is.null(self$gates)) {
        self$gates <- list()
      }   # Initialise la liste des gates si elle n’existe pas
      
      self$gates[[control_name]] <- gate_result   # Stocke le résultat sous le nom du contrôle
      
      return(gate_result)   # Retourne l’objet gate généré
    },
    
    # ============= definir_gates_density (beads uniquement) =====================
    definir_gates_density = function(control_name, n.cells = 2000, grid.n = 100, 
                                     bandwidth.factor = 1, fsc.channel = NULL, 
                                     ssc.channel = NULL) {
      
      if (is.null(self$asp_config)) {
        stop("La configuration ASP n'est pas initialisée. Appelez d'abord lancer_asp().")
      }   # Vérifie que la configuration AutoSpectral est disponible
      
      old_wd <- getwd()   # Sauvegarde le répertoire de travail actuel
      setwd(path.expand(self$dossier_racine))   # Se déplace dans le dossier racine du projet
      on.exit(setwd(old_wd))   # Garantit le retour au répertoire initial à la fin de la fonction
      
      output_dir <- file.path(self$dossier_racine, "figure_gate")   # Dossier où seront enregistrées les figures ASP
      if (!dir.exists(output_dir)) dir.create(output_dir)   # Crée le dossier s’il n’existe pas
      
      gate_result <- AutoSpectral::define.gate.density(
        control.file = "fcs_control_file.csv",   # Fichier de définition des contrôles ASP
        control.dir = path.expand(dirname(self$chemins_monomarques$chemin[1])),   # Dossier contenant les FCS monomarqués
        asp = self$asp_config,   # Paramètres ASP du cytomètre
        gate.name = control_name,   # Nom du gate à générer
        n.cells = n.cells,   # Nombre d’événements utilisés pour estimer la densité
        grid.n = grid.n,   # Taille de la grille pour la densité 2D
        bandwidth.factor = bandwidth.factor,   # Facteur de lissage pour l’estimation de densité
        fsc.channel = fsc.channel,   # Canal FSC si spécifié
        ssc.channel = ssc.channel,   # Canal SSC si spécifié
        output.dir = output_dir   # Dossier où les figures seront enregistrées
      )   # Appelle AutoSpectral pour calculer le gate basé sur la densité
      
      if (is.null(self$gates)) {
        self$gates <- list()
      }   # Initialise la liste des gates si elle n’existe pas
      
      self$gates[[control_name]] <- gate_result   # Stocke le gate sous le nom du contrôle
      
      return(gate_result)   # Retourne l’objet gate généré
    },
    
    # ============= definir_tune_gates =====================
    definir_tune_gates = function(gate.name, n_cells = 2000, percentile = 70, bandwidth = 1) {
      
      csv_file <- file.path(path.expand(self$dossier_racine), "fcs_control_file.csv") # Construit le chemin complet vers le fichier CSV de contrôle AutoSpectral
      
      if (!file.exists(csv_file)) {
        stop("Fichier CSV introuvable à : ", csv_file, ". Avez-vous bien lancé lancer_asp() ?")
      }  # Vérifie que le fichier de contrôle existe, sinon stoppe l’exécution
      
      dossier_fcs <- path.expand(dirname(self$chemins_monomarques$chemin[1]))  # Déduit le dossier contenant les fichiers FCS monomarqués
      output_dir <- file.path(self$dossier_racine, "figure_gate_tuning") # Dossier où seront enregistrées les figures de tuning
      
      if (!dir.exists(output_dir)) dir.create(output_dir) # Crée le dossier de sortie si nécessaire
      
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
      ) # Appelle AutoSpectral pour ajuster (tuner) le gate existant selon les paramètres fournis
      
      self$gates[[gate.name]] <- gate_tuned # Stocke le gate tuné dans l’objet R6 sous son nom
      
      return(invisible(gate_tuned)) # Retourne le résultat sans l’imprimer dans la console
    },
    
    # ============= charger_et_nettoyer =====================
    charger_et_nettoyer = function() {
      if (is.null(self$gates) || length(self$gates) == 0) {
        message("⚠️ Aucune gate définie. Les fichiers seront traités sans filtrage spatial.")
      } # Avertit l’utilisateur si aucun gate ASP n’a été défini (le nettoyage sera moins précis)
      
      dossier_fcs <- path.expand(dirname(self$chemins_monomarques$chemin[1])) # Déduit le dossier contenant les fichiers FCS monomarqués
      dossier_figures <- path.expand(file.path(self$dossier_racine, "figure_clean_controls")) # Dossier où seront enregistrées les figures de nettoyage
      
      if (!dir.exists(dossier_figures)) {
        dir.create(dossier_figures, recursive = TRUE)
      }  # Crée le dossier de sortie si nécessaire
      
      old_wd <- getwd()
      setwd(dossier_figures)
      on.exit(setwd(old_wd)) # Change temporairement le répertoire de travail pour que AutoSpectral écrive les figures ici
      
      flow_res <- AutoSpectral::define.flow.control(
        control.dir = dossier_fcs, 
        control.def.file = file.path(self$dossier_racine, self$asp_control_file),
        asp = self$asp_config,
        gate.list = if(length(self$gates) > 0) self$gates else NULL
      ) # Définit les contrôles ASP (application des gates, extraction des populations pertinentes)
      
      flow_cleaned <- AutoSpectral::clean.controls(
        flow.control = flow_res,
        asp = self$asp_config,
        main.figures = TRUE 
      ) # Nettoie les contrôles ASP (suppression des événements aberrants, figures de diagnostic)
      
      self$flow.control <- flow_cleaned # Stocke le résultat nettoyé dans l’objet R6
      
      message("✅ Chargement et nettoyage terminés avec succès.") # Message de confirmation
      return(invisible(self$flow.control)) # Retourne l’objet nettoyé sans l’afficher dans la console
    },
    
    # ============= extraire_spectre_fluorophore =====================
    extraire_spectre_fluorophore = function() { # Fonction qui extrait les spectres des fluorophores à partir des contrôles ASP
      if (is.null(self$flow.control)) { # Vérifie que les contrôles ASP ont été chargés et nettoyés
        stop("Erreur : flow.control n'est pas chargé.") # Stoppe si aucune donnée de contrôle n’est disponible
      }
      
      old_wd <- getwd() # Sauvegarde le répertoire de travail actuel
      setwd(path.expand(self$dossier_racine)) # Se place dans le dossier racine pour l’écriture des fichiers ASP
      
      spectra_result <- AutoSpectral::get.fluorophore.spectra( # Appelle AutoSpectral pour extraire les spectres des fluorophores
        flow.control = self$flow.control, # Données de contrôle nettoyées
        asp = self$asp_config # Paramètres ASP du cytomètre
      )
      
      setwd(old_wd) # Restaure le répertoire de travail initial
      self$spectra <- spectra_result # Stocke les spectres dans l’objet R6
      return(invisible(self$spectra)) # Retourne les spectres sans les afficher dans la console
    },
    
    # =========================
    # étapes optionnelles
    # =========================
    
    # ============= extraire_spectre_af =====================
    extraire_spectre_af = function(unstained_fcs_path, tissue_name, refine = TRUE) { # Fonction qui extrait le spectre d’autofluorescence à partir d’un échantillon unstained
      
      if (is.null(self$spectra)) { # Vérifie que les spectres fluorophores ont été extraits auparavant
        stop("Erreur : Les spectres fluorophores n'ont pas été extraits. Lancez extraire_spectre_fluorophore() d'abord.") # Stoppe si les spectres nécessaires ne sont pas disponibles
      }
      
      dossier_figures <- path.expand(file.path(self$dossier_racine, "figure_autofluorescence")) # Dossier où seront enregistrées les figures d’autofluorescence
      dossier_tables <- path.expand(file.path(self$dossier_racine, "table_autofluorescence")) # Dossier où seront enregistrées les tables de résultats
      
      if (!dir.exists(dossier_figures)) dir.create(dossier_figures, recursive = TRUE) # Crée le dossier des figures si nécessaire
      if (!dir.exists(dossier_tables)) dir.create(dossier_tables, recursive = TRUE) # Crée le dossier des tables si nécessaire
      
      af_result <- AutoSpectral::get.af.spectra( # Appelle AutoSpectral pour extraire le spectre d’autofluorescence
        unstained.sample = path.expand(unstained_fcs_path), # Chemin vers le fichier FCS unstained
        asp = self$asp_config, # Paramètres ASP du cytomètre
        spectra = self$spectra, # Spectres fluorophores précédemment extraits
        refine = refine, # Active ou non le raffinage du spectre AF
        figures = TRUE, # Génère les figures automatiquement
        plot.dir = dossier_figures, # Dossier où les figures seront enregistrées
        table.dir = dossier_tables, # Dossier où les tables seront enregistrées
        title = paste("Autofluorescence -", tissue_name) # Titre utilisé pour les sorties graphiques
      )
      
      if (is.null(self$af_spectra)) { # Initialise la liste des spectres AF si elle n’existe pas
        self$af_spectra <- list()
      }
      
      self$af_spectra[[tissue_name]] <- af_result # Stocke le spectre AF sous le nom du tissu
      return(invisible(af_result)) # Retourne le résultat sans l’afficher dans la console
    },
    
    # ============= preparer_variants_spectraux =====================
    preparer_variants_spectraux = function(tissue_af_name = NULL, refine = TRUE) { # Prépare les variantes spectrales en combinant spectres fluorophores et autofluorescence
      
      nom_tissu <- if (!is.null(tissue_af_name)) { # Si un nom de tissu est fourni, on l’utilise
        tissue_af_name
      } else if (length(self$af_spectra) > 0) { # Sinon, si des spectres AF existent, on prend le premier
        names(self$af_spectra)[1]
      } else {
        stop("Erreur : Aucune AF trouvée dans self$af_spectra. Lancez extraire_spectre_af() d'abord.") # Stoppe si aucun spectre AF n’est disponible
      }
      
      if (is.null(self$spectra) || is.null(self$af_spectra[[nom_tissu]])) { # Vérifie que les spectres fluorophores et AF existent pour ce tissu
        stop("Erreur : Spectres ou AF non trouvés pour le tissu : ", nom_tissu)
      }
      
      dossier_variants <- path.expand(file.path(self$dossier_racine, "figure_spectral_variants")) # Dossier où seront enregistrées les figures des variantes spectrales
      if (!dir.exists(dossier_variants)) dir.create(dossier_variants, recursive = TRUE) # Crée le dossier si nécessaire
      
      chemin_dossier_controles <- path.expand(dirname(self$chemins_monomarques$chemin[1])) # Dossier contenant les FCS monomarqués
      chemin_fichier_csv <- path.expand(file.path(self$dossier_racine, self$asp_control_file)) # Fichier CSV de définition des contrôles ASP
      
      variants_result <- AutoSpectral::get.spectral.variants( # Appelle AutoSpectral pour générer les variantes spectrales
        control.dir = chemin_dossier_controles, # Dossier des contrôles
        control.def.file = chemin_fichier_csv, # Fichier de définition ASP
        asp = self$asp_config, # Paramètres ASP du cytomètre
        spectra = self$spectra, # Spectres fluorophores
        af.spectra = self$af_spectra[[nom_tissu]], # Spectre d’autofluorescence du tissu
        refine = refine, # Active ou non le raffinage des variantes
        figures = TRUE, # Génère les figures automatiquement
        output.dir = dossier_variants # Dossier où les figures seront enregistrées
      )
      
      if (is.null(self$variants)) self$variants <- list() # Initialise la liste des variantes si elle n’existe pas
      
      self$variants[[nom_tissu]] <- variants_result # Stocke les variantes spectrales sous le nom du tissu
      return(invisible(variants_result)) # Retourne le résultat sans l’afficher dans la console
    },
    
    # ===========================   
    # unmixing 
    # ===========================
    
    # ============= unmix_fcs =====================
    unmix_fcs = function(fcs_file_path, tissue_name = NULL, method = "AutoSpectral", speed = "slow") { # Fonction qui applique l’unmixing spectral à un fichier FCS donné
      
      dossier_sortie <- file.path(self$dossier_racine, "AutoSpectral_unmixed") # Dossier où seront enregistrés les fichiers FCS unmixed
      if (!dir.exists(dossier_sortie)) dir.create(dossier_sortie) # Crée le dossier de sortie si nécessaire
      
      n_detectors <- ncol(self$spectra[[1]]) # Nombre de détecteurs (colonnes du spectre fluorophore)
      
      af_s <- if (!is.null(tissue_name) && !is.null(self$af_spectra[[tissue_name]])) { # Si un tissu est spécifié et possède un spectre AF, on l’utilise
        self$af_spectra[[tissue_name]]
      } else if (length(self$af_spectra) > 0) { # Sinon, si au moins un spectre AF existe, on prend le premier
        self$af_spectra[[1]]
      } else {
        matrix(nrow = 0, ncol = n_detectors) # Sinon, pas d’AF → matrice vide (AutoSpectral gère ce cas)
      }
      
      var_s <- if (!is.null(tissue_name) && !is.null(self$variants[[tissue_name]])) { # Si des variantes spectrales existent pour ce tissu, on les utilise
        self$variants[[tissue_name]]
      } else {
        NULL # Sinon, pas de variantes spectrales
      }
      
      AutoSpectral::unmix.fcs( # Appelle AutoSpectral pour effectuer l’unmixing spectral
        fcs.file = fcs_file_path, # Chemin du fichier FCS à traiter
        spectra = self$spectra, # Spectres fluorophores
        asp = self$asp_config, # Paramètres ASP du cytomètre
        flow.control = self$flow.control, # Contrôles ASP nettoyés
        method = method, # Méthode d’unmixing (AutoSpectral par défaut)
        af.spectra = af_s, # Spectre d’autofluorescence à utiliser
        spectra.variants = var_s, # Variantes spectrales si disponibles
        speed = speed, # Mode de calcul (slow = plus précis)
        output.dir = dossier_sortie, # Dossier de sortie pour les fichiers unmixed
        parallel = TRUE # Active le traitement parallèle
      )
      
    },
    
    # ============= unmix_folder =====================
    unmix_folder = function(folder_path, tissue_name = NULL, method = "AutoSpectral", speed = "slow") { # Fonction qui applique l’unmixing spectral à tous les fichiers FCS d’un dossier
      
      dossier_sortie <- file.path(self$dossier_racine, "AutoSpectral_unmixed") # Dossier où seront enregistrés les fichiers unmixed
      if (!dir.exists(dossier_sortie)) dir.create(dossier_sortie, recursive = TRUE) # Crée le dossier de sortie si nécessaire
      
      n_detectors <- ncol(self$spectra[[1]]) # Nombre de détecteurs (colonnes du spectre fluorophore)
      
      af_s <- if (!is.null(tissue_name) && !is.null(self$af_spectra[[tissue_name]])) { # Si un tissu est spécifié et possède un spectre AF, on l’utilise
        self$af_spectra[[tissue_name]]
      } else if (length(self$af_spectra) > 0) { # Sinon, si au moins un spectre AF existe, on prend le premier
        self$af_spectra[[1]]
      } else {
        matrix(nrow = 0, ncol = n_detectors) # Sinon, pas d’AF → matrice vide (AutoSpectral gère ce cas)
      }
      
      var_s <- if (!is.null(tissue_name) && !is.null(self$variants[[tissue_name]])) { # Si des variantes spectrales existent pour ce tissu, on les utilise
        self$variants[[tissue_name]]
      } else {
        NULL # Sinon, pas de variantes spectrales
      }
      
      AutoSpectral::unmix.folder( # Appelle AutoSpectral pour effectuer l’unmixing spectral sur tout un dossier
        fcs.dir = folder_path, # Dossier contenant les fichiers FCS à traiter
        spectra = self$spectra, # Spectres fluorophores
        asp = self$asp_config, # Paramètres ASP du cytomètre
        flow.control = self$flow.control, # Contrôles ASP nettoyés
        method = method, # Méthode d’unmixing (AutoSpectral par défaut)
        af.spectra = af_s, # Spectre d’autofluorescence à utiliser
        spectra.variants = var_s, # Variantes spectrales si disponibles
        speed = speed, # Mode de calcul (slow = plus précis)
        output.dir = dossier_sortie, # Dossier de sortie pour les fichiers unmixed
        parallel = TRUE # Active le traitement parallèle
      )
    },
    
    # ============= verifier_qualite_unmix =====================
    verifier_qualite_unmix = function(fluorophore, single_stained_fcs, unstained_fcs, cytometer = "aurora", gate = TRUE) { # Fonction qui vérifie la qualité de l’unmixing pour un fluorophore donné
      
      dossier_gates <- file.path(self$dossier_racine, "figure_gate") # Dossier où sont enregistrées les figures des gates ASP
      if (!dir.exists(dossier_gates)) dir.create(dossier_gates, recursive = TRUE) # Crée le dossier si nécessaire
      
      if (is.null(self$spectra)) { # Vérifie que les spectres fluorophores ont été extraits
        stop("Erreur : Les spectres n'ont pas été extraits. Lancez extraire_spectre_fluorophore() d'abord.") # Stoppe si les spectres sont absents
      }
      
      if (!fluorophore %in% rownames(self$spectra)) { # Vérifie que le fluorophore demandé existe dans les spectres
        stop("Erreur : Fluorophore introuvable dans les spectres.") # Stoppe si le fluorophore n’est pas reconnu
      }
      
      dossier_figures <- path.expand(file.path(self$dossier_racine, "figure_compare_unmix")) # Dossier où seront enregistrées les figures de comparaison unmix
      if (!dir.exists(dossier_figures)) dir.create(dossier_figures, recursive = TRUE) # Crée le dossier si nécessaire
      
      spectre_cible <- self$spectra[fluorophore, ] # Extrait le spectre de référence du fluorophore
      
      qc_result <- AutoSpectral::compare.unmix( # Appelle AutoSpectral pour comparer l’unmixing entre single-stained et unstained
        single.stained.fcs = path.expand(single_stained_fcs), # Fichier FCS du tube monomarqué
        unstained.fcs      = path.expand(unstained_fcs), # Fichier FCS du tube unstained
        fluorophore        = fluorophore, # Nom du fluorophore à tester
        spectra            = self$spectra, # Spectres fluorophores complets
        ref.spectrum       = spectre_cible, # Spectre de référence utilisé pour la comparaison
        test.spectrum      = spectre_cible, # Spectre testé (identique ici, mais structure requise par AutoSpectral)
        cytometer          = cytometer, # Type de cytomètre (Aurora par défaut)
        gate               = gate, # Active ou non le gating automatique pour la comparaison
        plot.dir           = dossier_figures # Dossier où les figures de QC seront enregistrées
      )
      
      return(invisible(qc_result)) # Retourne le résultat sans l’afficher dans la console
    },
    
    # ============= charger_fcs_unmixes =====================
    charger_fcs_unmixes = function(dossier = "AutoSpectral_unmixed") { # Fonction qui charge tous les fichiers FCS déjà unmixed
      
      chemin_complet <- file.path(self$dossier_racine, dossier) # Construit le chemin complet vers le dossier d’unmixing
      if(!dir.exists(chemin_complet)) stop("Dossier introuvable : ", chemin_complet) # Stoppe si le dossier n’existe pas
      
      fichiers <- list.files(chemin_complet, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE) # Liste tous les fichiers FCS du dossier
      
      for (f in fichiers) { # Boucle sur chaque fichier FCS trouvé
        nom_cle <- basename(f) # Utilise le nom du fichier comme clé dans la liste
        self$echantillons_traites[[nom_cle]] <- flowCore::read.FCS(f, truncate_max_range = FALSE) # Charge le FCS unmixed et le stocke
      }
      
      message("Chargement terminé : ", length(fichiers), " échantillons importés depuis ", dossier) # Message récapitulatif
    },
    
    # ============= obtenir_chemins_figures =====================
    obtenir_chemins_figures = function(control_name) { # Fonction qui récupère les chemins des figures ASP pour un contrôle donné
      
      dossier <- file.path(self$dossier_racine, "figure_gate") # Dossier contenant les figures de gates
      if (!dir.exists(dossier)) return(NULL) # Retourne NULL si le dossier n’existe pas
      
      list.files(dossier, pattern=paste0(control_name, ".*\\.png$"), full.names=TRUE) # Retourne tous les PNG correspondant au contrôle
    },
    
    # ============= visualiser_unmixing =====================
    visualiser_unmixing = function(nom_fichier_fcs, canal_x, canal_y, cofacteur = 150) { # Fonction qui visualise le résultat de l’unmixing pour deux canaux donnés
      
      # Sécurité : mêmes précautions que dans visualiser_compensation() —
      # canal_x/canal_y peuvent être NULL ou vides tant que les sélecteurs
      # Shiny ne sont pas encore peuplés.
      if (is.null(canal_x) || is.null(canal_y) || canal_x == "" || canal_y == "") return(NULL)
      
      fcs_unmixed <- self$echantillons_traites[[nom_fichier_fcs]] # Récupère le fichier FCS unmixed en mémoire
      if (is.null(fcs_unmixed)) stop("Fichier introuvable en mémoire.") # Stoppe si le fichier n’a pas été chargé
      
      trans_list <- flowCore::transformList(c(canal_x, canal_y), flowCore::arcsinhTransform(a = 0, b = 1/cofacteur, c = 0)) # Crée une transformation arcsinh pour les deux canaux
      mat <- flowCore::exprs(flowCore::transform(fcs_unmixed, trans_list))[, c(canal_x, canal_y)] # Applique la transformation et extrait les deux colonnes transformées
      df <- as.data.frame(mat) # Convertit en data.frame
      colnames(df) <- c("Axe_X", "Axe_Y") # Renomme les colonnes pour ggplot
      lim_x <- range(df$Axe_X, na.rm = TRUE) # Calcule les limites X pour le plot
      lim_y <- range(df$Axe_Y, na.rm = TRUE) # Calcule les limites Y pour le plot
      df_densite <- calculer_densite_raster(df$Axe_X, df$Axe_Y, lim_x, lim_y) # Calcule la densité rasterisée via ta fonction interne
      
      if (is.null(df_densite)) { # Si la densité n’a pas pu être calculée (trop peu d’événements)
        return(ggplot() + theme_bw() + # Retourne un plot vide avec message
                 labs(title = paste("Résultat après Unmixing :", nom_fichier_fcs),
                      subtitle = "Pas assez d'événements pour tracer la densité",
                      x = self$obtenir_label(fcs_unmixed, canal_x), y = self$obtenir_label(fcs_unmixed, canal_y)))
      }
      
      ggplot(df_densite, aes(x = X, y = Y, fill = densite)) + # Construit le plot rasterisé
        geom_raster(interpolate = TRUE) + 
        scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS) + 
        coord_cartesian(xlim = lim_x, ylim = lim_y) + 
        theme_bw() + 
        theme(legend.position = "none", aspect.ratio = 1) + 
        labs(
          title = paste("Résultat après Unmixing :", nom_fichier_fcs), 
          x = self$obtenir_label(fcs_unmixed, canal_x), 
          y = self$obtenir_label(fcs_unmixed, canal_y)  
        )
    },
    
    # ============= visualiser_figures =====================
    visualiser_figures = function(dossier_nom) { # Fonction qui affiche toutes les images d’un dossier dans une visionneuse HTML
      
      chemin_dossier <- file.path(self$dossier_racine, dossier_nom) # Construit le chemin complet vers le dossier demandé
      fichiers <- list.files(chemin_dossier, pattern="\\.(jpg|jpeg|png)$", full.names=TRUE, ignore.case=TRUE) # Liste toutes les images JPG/PNG du dossier
      
      if (length(fichiers) == 0) return(message("Aucune image.")) # Si aucune image trouvée, message et sortie
      
      html_elements <- sapply(fichiers, function(f) { # Pour chaque image, génère un bloc HTML contenant le nom + l’image encodée en base64
        mime <- ifelse(grepl("\\.(jpg|jpeg)", f, ignore.case=TRUE), "image/jpeg", "image/png") # Détermine le type MIME selon l’extension
        paste0("<div><h3>", basename(f), "</h3><img src='", base64enc::dataURI(file=f, mime=mime), "' style='max-width:100%'></div>") # Génère le bloc HTML
      })
      
      temp_html <- tempfile(fileext=".html") # Crée un fichier HTML temporaire
      writeLines(c("<html><body>", html_elements, "</body></html>"), temp_html) # Écrit le contenu HTML dans le fichier
      
      rstudioapi::viewer(temp_html) # Ouvre la page HTML dans la visionneuse RStudio
    }, 
    
    
    # ============================================================
    #       ️ SECTION PRÉ-TRAITEMENT
    # ============================================================
    
    # ============= obtenir_derniere_source =====================
    obtenir_derniere_source = function() {
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
    
    # ============= appliquer_peacoqc =====================
    appliquer_peacoqc = function(dossier_rapports = NULL, reglages_specifiques = list()) { # Applique PeacoQC sur tous les échantillons unmixed
      
      parametres_par_defaut <- list( # Paramètres QC par défaut basés sur PeacoQC::PeacoQC()
        determine_good_cells = "all",
        min_cells = 150,
        max_bins = 500,
        step = 500,
        MAD = 6,
        IT_limit = 0.6,
        consecutive_bins = 5,
        remove_zeros = FALSE,
        force_IT = 150,
        peak_removal = 1/3,
        min_nr_bins_peakdetection = 10
      )
      
      config_qc <- utils::modifyList(parametres_par_defaut, reglages_specifiques) # Fusionne les réglages spécifiques avec les valeurs par défaut
      
      self$post_PeacoQC <- list() # Initialise la liste des flowFrames filtrés
      self$plots_peacoqc <- list() # Initialise la liste des plots ggplot
      self$plots_peacoqc_natif <- list() # Initialise la liste des PNG natifs
      
      for (nom in names(self$echantillons_traites)) { # Boucle sur chaque échantillon unmixed
        
        message("Exécution PeacoQC sur : ", nom) # Indique l’échantillon en cours
        ff_actuel <- self$echantillons_traites[[nom]] # Récupère le flowFrame
        vrai_canal_temps <- grep("time", colnames(ff_actuel), value = TRUE, ignore.case = TRUE)[1] # Cherche le canal temporel
        if (is.na(vrai_canal_temps)) vrai_canal_temps <- "Time" # Fallback si non trouvé
        tous_canaux <- colnames(ff_actuel) # Liste tous les canaux
        canaux_fluo <- tous_canaux[!grepl("fsc|ssc|time", tous_canaux, ignore.case = TRUE)] # Sélectionne uniquement les canaux fluorescence
        
        dossier_temp_plot <- file.path( # Crée un dossier temporaire unique pour les PNG
          tempdir(),
          paste0("PeacoQC_", gsub("[^a-zA-Z0-9_]", "_", nom), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        )
        dir.create(dossier_temp_plot, recursive = TRUE, showWarnings = FALSE) # Crée le dossier
        
        res <- PeacoQC::PeacoQC( # Exécute PeacoQC sur l’échantillon
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
        
        self$post_PeacoQC[[nom]] <- res$FinalFF # Stocke le flowFrame filtré
        
        self$plots_peacoqc[[nom]] <- res$Plot # Stocke le plot ggplot
        
        pngs <- list.files(dossier_temp_plot, pattern = "\\.png$", full.names = TRUE, recursive = TRUE) # Recherche les PNG générés
        self$plots_peacoqc_natif[[nom]] <- if (length(pngs) > 0) pngs[1] else NULL # Stocke le PNG natif (un seul par échantillon)
      }
      
      self$parametres_peacoqc_utilises <- config_qc # Archive les paramètres utilisés
      
      if (!is.null(self$mettre_a_jour_pipeline)) self$mettre_a_jour_pipeline("PeacoQC") # Met à jour le statut du pipeline
    },
    
    # ============= appliquer_flowai =====================
    appliquer_flowai = function(reglages_specifiques = list()) { # Applique flowAI sur tous les échantillons unmixed
      
      if (length(self$echantillons_traites) == 0) { # Vérifie qu’il existe des échantillons à analyser
        stop("Aucun échantillon traité disponible pour flowAI.") # Stoppe si aucun échantillon n’est chargé
      }
      
      remove_from       <- reglages_specifiques$remove_from       %||% "all" # Paramètre flowAI : quelles étapes de QC appliquer
      timeCh            <- reglages_specifiques$timeCh            %||% NULL # Canal temporel si spécifié
      second_fractionFR <- reglages_specifiques$second_fractionFR %||% 0.1 # Fraction utilisée pour l’algorithme FlowRate
      alphaFR           <- reglages_specifiques$alphaFR           %||% 0.01 # Seuil statistique pour FlowRate
      decompFR          <- if (isTRUE(reglages_specifiques$decompFR)) "cffilter" else FALSE # Active la décomposition FR si demandée
      
      ChExcludeFS <- reglages_specifiques$ChExcludeFS # Canaux exclus pour FlowSignal
      if (is.null(ChExcludeFS) || length(ChExcludeFS) == 0) ChExcludeFS <- c("FSC", "SSC") # Exclut FSC/SSC par défaut
      
      outlier_binsFS <- reglages_specifiques$outlier_binsFS %||% FALSE # Détection des bins aberrants dans FlowSignal
      pen_valueFS    <- reglages_specifiques$pen_valueFS    %||% 500 # Pénalité pour FlowSignal
      max_cptFS      <- reglages_specifiques$max_cptFS      %||% 3 # Nombre max de points de rupture
      
      ChExcludeFM <- reglages_specifiques$ChExcludeFM # Canaux exclus pour FlowMargin
      if (is.null(ChExcludeFM) || length(ChExcludeFM) == 0) ChExcludeFM <- c("FSC", "SSC") # Exclut FSC/SSC par défaut
      
      sideFM       <- reglages_specifiques$sideFM       %||% "both" # Côté du margin à analyser
      neg_valuesFM <- reglages_specifiques$neg_valuesFM %||% 1 # Gestion des valeurs négatives dans FlowMargin
      
      self$post_flowAI <- list() # Initialise la liste des flowFrames QC
      self$rapports_flowai <- list() # Initialise la liste des dossiers de rapports
      
      for (nom in names(self$echantillons_traites)) { # Boucle sur chaque échantillon
        
        fcs_obj <- self$echantillons_traites[[nom]] # Récupère l’échantillon
        dossier_tmp <- file.path(tempdir(), paste0("flowAI_", nom, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))) # Dossier temporaire unique
        dir.create(dossier_tmp, recursive = TRUE, showWarnings = FALSE) # Crée le dossier
        
        res <- tryCatch({ # Exécute flowAI avec gestion d’erreur
          
          flowAI::flow_auto_qc(
            fcsfiles            = fcs_obj, # FlowFrame à analyser
            remove_from         = remove_from, # Étapes QC à appliquer
            output              = 1, # Génère les rapports HTML
            timeCh              = timeCh, # Canal temporel
            second_fractionFR   = second_fractionFR, # Paramètre FlowRate
            alphaFR             = alphaFR, # Seuil FlowRate
            decompFR            = decompFR, # Décomposition FlowRate
            ChExcludeFS         = ChExcludeFS, # Exclusions FlowSignal
            outlier_binsFS      = outlier_binsFS, # Détection outliers FlowSignal
            pen_valueFS         = pen_valueFS, # Pénalité FlowSignal
            max_cptFS           = max_cptFS, # Points de rupture FlowSignal
            ChExcludeFM         = ChExcludeFM, # Exclusions FlowMargin
            sideFM              = sideFM, # Côté FlowMargin
            neg_valuesFM        = neg_valuesFM, # Gestion valeurs négatives FlowMargin
            html_report         = "_QC", # Nom du rapport HTML
            mini_report         = "_QCmini", # Rapport réduit
            fcs_QC              = "_QC", # Fichier FCS QC
            fcs_highQ           = FALSE, # Pas de FCS high quality
            fcs_lowQ            = FALSE, # Pas de FCS low quality
            folder_results      = dossier_tmp # Dossier de sortie
          )
          
        }, error = function(e) { # Capture les erreurs flowAI
          stop(paste("Erreur flowAI pour l'échantillon", nom, ":", conditionMessage(e))) # Stoppe avec message clair
        })
        
        if (inherits(res, "flowSet")) { # flowAI peut retourner un flowSet
          res <- res[[1]] # On extrait le premier élément
        }
        
        self$post_flowAI[[nom]] <- res # Stocke le flowFrame QC
        self$rapports_flowai[[nom]] <- dossier_tmp # Stocke le dossier de rapport
      }
      
      self$mettre_a_jour_pipeline("flowAI") # Met à jour l’état du pipeline
      invisible(TRUE) # Retourne TRUE sans affichage
    },
    
    # ============= exporter_fcs_qc =====================
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
    
    # ============= generer_pdf_resume_qc =====================
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
    
    # ============= retirer_les_bordures =====================
    retirer_les_bordures = function(canal1, canal2, nom_echantillon = NULL) { # Applique RemoveMargins sur un ou plusieurs échantillons
      
      self$canaux_bordures <- c(canal1, canal2) # Enregistre les deux canaux à utiliser pour détecter les marges
      
      liste_source <- self$obtenir_derniere_source() # Récupère la dernière étape QC disponible (PeacoQC, flowAI, ou unmixed)
      
      noms_a_traiter <- if (is.null(nom_echantillon)) names(liste_source) else nom_echantillon # Détermine quels échantillons traiter
      
      if (is.null(noms_a_traiter) || length(noms_a_traiter) == 0) { # Vérifie qu’il y a des échantillons à traiter
        warning("Aucun échantillon trouvé à traiter pour le retrait des bordures.") # Avertit si aucun échantillon
        return(NULL) # Sortie silencieuse
      }
      
      for (nom in noms_a_traiter) { # Boucle sur les échantillons sélectionnés
        if (!is.null(liste_source[[nom]])) { # Vérifie que l’échantillon existe dans la source QC
          message("Retrait des bordures (Margins) sur : ", nom) # Indique l’échantillon en cours
          
          self$post_retrait_bordures[[nom]] <- PeacoQC::RemoveMargins( # Applique l’algorithme RemoveMargins
            ff       = liste_source[[nom]], # FlowFrame à nettoyer
            channels = self$canaux_bordures # Canaux utilisés pour détecter les marges
          )
        }
      }
      
      if (!is.null(self$mettre_a_jour_pipeline)) self$mettre_a_jour_pipeline("bordures", nom_echantillon) # Met à jour l’état du pipeline
    },
    
    # ============= retirer_les_debris =====================
    retirer_les_debris = function(matrice_points, canal_x, canal_y, nom_echantillon = NULL, source_nettoyage = "brutes") { # Applique un polygonGate pour retirer les débris
      
      if (is.null(matrice_points)) stop("Aucune coordonnée de gate fournie.") # Vérifie que la gate est fournie
      if (!is.matrix(matrice_points) && !is.data.frame(matrice_points)) {
        stop("La structure de la gate doit être une matrice ou un data.frame.") # Vérifie le format de la gate
      }
      if (nrow(matrice_points) < 3) stop("Un polygone de filtrage nécessite au moins 3 points.") # Un polygonGate nécessite ≥ 3 points
      
      liste_source <- NULL # Initialise la source de données
      
      if (source_nettoyage == "peacoqc") {
        liste_source <- self$post_PeacoQC # Utilise les données post-PeacoQC
      } else if (source_nettoyage == "flowai") {
        liste_source <- self$post_flowAI # Utilise les données post-flowAI
      }
      
      if (is.null(liste_source) || length(liste_source) == 0) { # Si la source demandée est vide
        if (source_nettoyage %in% c("peacoqc", "flowai")) {
          message("⚠️ [Warning] Source '", source_nettoyage, "' introuvable ou vide. Repli sur les données compensées brutes.") # Avertit du fallback
        }
        liste_source <- self$echantillons_traites # Repli sur les données unmixed brutes
      }
      
      if (is.null(liste_source) || length(liste_source) == 0) {
        stop("Aucune donnée disponible (brute, PeacoQC ou flowAI) pour appliquer le filtre débris.") # Stoppe si aucune donnée exploitable
      }
      
      matrice_points <- as.matrix(matrice_points[, 1:2]) # Convertit la gate en matrice 2 colonnes
      colnames(matrice_points) <- c(canal_x, canal_y) # Associe les noms des canaux
      
      polygone_debris <- flowCore::polygonGate(.gate = matrice_points, filterId = "Gate_Debris") # Crée le polygonGate
      
      if (is.null(self$gate_debris)) self$gate_debris <- list() # Initialise la liste des gates débris
      if (is.null(self$post_debris)) self$post_debris <- list() # Initialise la liste des flowFrames filtrés
      
      appliquer_le_filtrage = function(nom) { # Fonction interne appliquant le filtrage
        flowframe_entree <- liste_source[[nom]] # Récupère l’échantillon
        if (is.null(flowframe_entree)) return(NULL) # Ignore si absent
        
        message("Application du filtre Débris (PolygonGate) via source '", source_nettoyage, "' sur : ", nom) # Indique l’échantillon traité
        
        resultat_filtre <- flowCore::filter(flowframe_entree, polygone_debris) # Applique le polygonGate
        self$gate_debris[[nom]] <- polygone_debris # Stocke la gate utilisée
        self$post_debris[[nom]] <- flowframe_entree[resultat_filtre@subSet, ] # Stocke les événements filtrés
      }
      
      noms_a_traiter <- if (is.null(nom_echantillon)) names(liste_source) else nom_echantillon # Détermine les échantillons à traiter
      
      for (nom in noms_a_traiter) { # Boucle sur les échantillons
        appliquer_le_filtrage(nom) # Applique le filtrage
      }
      
      if (!is.null(self$mettre_a_jour_pipeline)) {
        self$mettre_a_jour_pipeline("debris", nom_echantillon) # Met à jour l’état du pipeline
      }
      
      return(invisible(self)) # Retourne l’objet sans affichage
    },
    
    # ============= appliquer_gate_nomme =====================
    appliquer_gate_nomme = function(nom_gate, matrice_points, canal_x, canal_y,
                                    source_nettoyage = "brutes", nom_echantillon = NULL) { # Applique un gate polygonal nommé et conserve son historique
      
      if (is.null(nom_gate) || nchar(trimws(nom_gate)) == 0) stop("Le gate doit avoir un nom.") # Vérifie que le gate a un nom valide
      if (is.null(matrice_points) || nrow(matrice_points) < 3) stop("Polygone invalide (< 3 points).") # Vérifie que la gate est un polygone valide
      
      if (length(self$gates_history) > 0 && is.null(self$gates_history[[nom_gate]])) { # Cas où un historique existe mais pas pour ce gate
        dernier_gate <- self$gates_history[[length(self$gates_history)]] # Récupère le dernier gate appliqué
        liste_source <- lapply(names(dernier_gate), function(n) dernier_gate[[n]]$post_data) # Utilise les données post-gate précédent
        names(liste_source) <- names(dernier_gate) # Conserve les noms des échantillons
      } else if (!is.null(self$gates_history[[nom_gate]]) && length(self$gates_history) > 1) { # Cas où ce gate existe déjà dans l’historique
        
        idx_gate <- match(nom_gate, names(self$gates_history)) # Trouve la position du gate dans l’historique
        if (idx_gate > 1) {
          gate_precedent <- self$gates_history[[idx_gate - 1]] # Récupère le gate précédent
          liste_source <- lapply(names(gate_precedent), function(n) gate_precedent[[n]]$post_data) # Utilise les données post-gate précédent
          names(liste_source) <- names(gate_precedent)
        } else {
          liste_source <- self$obtenir_derniere_source() # Si premier gate, utilise la dernière source QC disponible
        }
      } else {
        if (source_nettoyage == "peacoqc" && length(self$post_PeacoQC) > 0) {
          liste_source <- self$post_PeacoQC # Utilise les données post-PeacoQC
        } else if (source_nettoyage == "flowai" && length(self$post_flowAI) > 0) {
          liste_source <- self$post_flowAI # Utilise les données post-flowAI
        } else {
          liste_source <- self$obtenir_derniere_source() # Sinon utilise la dernière source QC disponible
        }
      }
      
      if (is.null(liste_source) || length(liste_source) == 0) {
        stop("Aucune donnée disponible pour appliquer le gate.") # Stoppe si aucune donnée exploitable
      }
      
      mat_gate <- as.matrix(matrice_points[, 1:2]) # Convertit la gate en matrice 2 colonnes
      colnames(mat_gate) <- c(canal_x, canal_y) # Associe les noms des canaux
      
      polygone <- flowCore::polygonGate(.gate = mat_gate, filterId = nom_gate) # Crée le polygonGate nommé
      
      deja_initialise  <- !is.null(self$gates_history[[nom_gate]]) && length(self$gates_history[[nom_gate]]) > 0 # Vérifie si ce gate a déjà été appliqué
      noms_a_traiter   <- if (is.null(nom_echantillon) || !deja_initialise) names(liste_source) else nom_echantillon # Détermine les échantillons à traiter
      resultats_gate   <- if (deja_initialise) self$gates_history[[nom_gate]] else list() # Récupère ou initialise l’historique du gate
      
      for (nom in noms_a_traiter) { # Boucle sur les échantillons
        ff_entree <- liste_source[[nom]] # Récupère l’échantillon
        if (is.null(ff_entree)) next # Ignore si absent
        
        n_avant <- nrow(flowCore::exprs(ff_entree)) # Nombre d’événements avant filtrage
        message("Application gate '", nom_gate, "' sur : ", nom) # Indique l’échantillon traité
        
        res_filtre <- flowCore::filter(ff_entree, polygone) # Applique le polygonGate
        ff_apres   <- ff_entree[res_filtre@subSet, ] # Conserve les événements filtrés
        n_apres    <- nrow(flowCore::exprs(ff_apres)) # Nombre d’événements après filtrage
        
        resultats_gate[[nom]] <- list( # Stocke les résultats pour cet échantillon
          polygone  = mat_gate,
          canal_x   = canal_x,
          canal_y   = canal_y,
          post_data = ff_apres,
          n_avant   = n_avant,
          n_apres   = n_apres
        )
      }
      
      self$gates_history[[nom_gate]] <- resultats_gate # Met à jour l’historique du gate
      
      self$gate_debris <- lapply(resultats_gate, function(r) { # Met à jour les gates actifs
        flowCore::polygonGate(.gate = r$polygone, filterId = nom_gate)
      })
      
      self$post_debris <- lapply(resultats_gate, function(r) r$post_data) # Met à jour les données filtrées
      
      return(invisible(self)) # Retourne l’objet sans affichage
    },
    
    # ============= retirer_doublets_FSC =====================
    retirer_doublets_FSC = function(facteur_sensibilite = 4, axe_discrimination = "H_A", nom_echantillon = NULL) { # Détection statistique des doublets via ratios FSC
      
      liste_fcs_source <- if (!is.null(self$post_debris) && length(self$post_debris) > 0) self$post_debris else self$obtenir_derniere_source() # Utilise les données post-débris si disponibles, sinon la dernière source QC
      label_source <- if (!is.null(self$post_debris) && length(self$post_debris) > 0) "post_debris" else "source_parente" # Indique la provenance des données
      
      if (is.null(self$gate_doublets_FSC)) self$gate_doublets_FSC <- list() # Initialise le stockage des seuils FSC
      if (is.null(self$post_doublets_FSC)) self$post_doublets_FSC <- list() # Initialise le stockage des flowFrames filtrés
      if (is.null(self$post_doublets_final)) self$post_doublets_final <- list() # Initialise le stockage final des singlets
      
      calculer_doublets = function(nom, lbl_src) { # Fonction interne appliquant la détection des doublets pour un échantillon
        flowframe_entree <- liste_fcs_source[[nom]] # Récupère l’échantillon
        if (is.null(flowframe_entree)) return(NULL) # Ignore si absent
        
        matrice_exprs <- flowCore::exprs(flowframe_entree) # Extrait la matrice des intensités FSC
        
        if (axe_discrimination == "H_A") { # Choix des axes FSC-H vs FSC-A
          canal_x <- grep("FSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("FSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "W_A") { # Choix des axes FSC-W vs FSC-A
          canal_x <- grep("FSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("FSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "H_W") { # Choix des axes FSC-H vs FSC-W
          canal_x <- grep("FSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("FSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else {
          stop("axe_discrimination invalide. Choisir parmi 'H_A', 'W_A' ou 'H_W'") # Stoppe si axe invalide
        }
        
        if (is.na(canal_x) || is.na(canal_y)) { # Vérifie que les canaux FSC nécessaires existent
          warning("Canaux FSC requis pour l'axe ", axe_discrimination, " introuvables pour ", nom, ". Étape ignorée.") # Avertit si structure FCS incompatible
          self$post_doublets_FSC[[nom]] <- flowframe_entree # Conserve l’échantillon tel quel
          return(NULL) # Ignore cet échantillon
        } 
        
        ratio_Y_X <- matrice_exprs[, canal_y] / (matrice_exprs[, canal_x] + 1e-6) # Ratio FSC pour détecter les doublets (évite division par zéro)
        val_mad <- stats::mad(ratio_Y_X, na.rm = TRUE) # Dispersion robuste des singlets
        if (val_mad == 0) val_mad <- mean(ratio_Y_X, na.rm = TRUE) * 0.05 # Valeur de secours si MAD nul
        
        seuil_statistique <- stats::median(ratio_Y_X, na.rm = TRUE) + (facteur_sensibilite * val_mad) # Seuil d’exclusion des doublets
        
        self$gate_doublets_FSC[[nom]] <- list( # Stocke les paramètres du seuil
          type = "stat", seuil = seuil_statistique, facteur = facteur_sensibilite,
          source = lbl_src, channels = c(canal_x, canal_y)
        ) 
        
        flowframe_filtre <- flowframe_entree[ratio_Y_X < seuil_statistique & is.finite(ratio_Y_X), ] # Retire les doublets et valeurs instables
        self$post_doublets_FSC[[nom]] <- flowframe_filtre # Stocke l’échantillon filtré
        self$post_doublets_final[[nom]] <- flowframe_filtre # Met à jour la version finale des singlets
      } 
      
      noms_a_traiter <- if (is.null(nom_echantillon) || length(self$gate_doublets_FSC) == 0) names(liste_fcs_source) else nom_echantillon # Traite tout au premier passage, puis échantillon ciblé si déjà initialisé
      for (nom in noms_a_traiter) { calculer_doublets(nom, label_source) } # Applique la détection des doublets à chaque échantillon
      
      self$mettre_a_jour_pipeline("doublets_FSC", nom_echantillon) # Met à jour l’état du pipeline
    },
    
    # ============= retirer_doublets_SSC =====================
    retirer_doublets_SSC = function(facteur_sensibilite = 4, axe_discrimination = "H_A", nom_echantillon = NULL) { # Détection statistique des doublets via ratios SSC
      
      liste_fcs_source <- if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) { # Utilise les données post-FSC si disponibles
        self$post_doublets_FSC
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) { # Sinon utilise les données post-débris
        self$post_debris
      } else {
        self$obtenir_derniere_source() # Sinon utilise la dernière source QC disponible
      } 
      
      label_source <- if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) { # Indique la provenance des données
        "post_doublets_FSC"
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) {
        "post_debris"
      } else {
        "source_parente"
      } 
      
      if (is.null(self$gate_doublets_SSC)) self$gate_doublets_SSC <- list() # Initialise le stockage des seuils SSC
      if (is.null(self$post_doublets_SSC)) self$post_doublets_SSC <- list() # Initialise le stockage des flowFrames filtrés
      if (is.null(self$post_doublets_final)) self$post_doublets_final <- list() # Initialise le stockage final des singlets
      
      calculer_doublets = function(nom, lbl_src) { # Fonction interne appliquant la détection des doublets pour un échantillon
        flowframe_entree <- liste_fcs_source[[nom]] # Récupère l’échantillon
        if (is.null(flowframe_entree)) return(NULL) # Ignore si absent
        
        matrice_exprs <- flowCore::exprs(flowframe_entree) # Extrait la matrice des intensités SSC
        
        if (axe_discrimination == "H_A") { # Choix des axes SSC-H vs SSC-A
          canal_x <- grep("SSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("SSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "W_A") { # Choix des axes SSC-W vs SSC-A
          canal_x <- grep("SSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("SSC-A", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else if (axe_discrimination == "H_W") { # Choix des axes SSC-H vs SSC-W
          canal_x <- grep("SSC-H", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
          canal_y <- grep("SSC-W", colnames(matrice_exprs), value = TRUE, ignore.case = TRUE)[1]
        } else {
          stop("axe_discrimination invalide. Choisir parmi 'H_A', 'W_A' ou 'H_W'") # Stoppe si axe invalide
        }
        
        if (is.na(canal_x) || is.na(canal_y)) { # Vérifie que les canaux SSC nécessaires existent
          warning("Canaux SSC requis pour l'axe ", axe_discrimination, " introuvables pour ", nom, ". Étape ignorée.") # Avertit si structure FCS incompatible
          self$post_doublets_SSC[[nom]] <- flowframe_entree # Conserve l’échantillon tel quel
          return(NULL) # Ignore cet échantillon
        } 
        
        ratio_Y_X <- matrice_exprs[, canal_y] / (matrice_exprs[, canal_x] + 1e-6) # Ratio SSC pour détecter les doublets (évite division par zéro)
        val_mad <- stats::mad(ratio_Y_X, na.rm = TRUE) # Dispersion robuste des singlets
        if (val_mad == 0) val_mad <- mean(ratio_Y_X, na.rm = TRUE) * 0.05 # Valeur de secours si MAD nul
        
        seuil_statistique <- stats::median(ratio_Y_X, na.rm = TRUE) + (facteur_sensibilite * val_mad) # Seuil d’exclusion des doublets
        
        self$gate_doublets_SSC[[nom]] <- list( # Stocke les paramètres du seuil
          type = "stat", seuil = seuil_statistique, facteur = facteur_sensibilite,
          source = lbl_src, channels = c(canal_x, canal_y)
        )
        
        flowframe_filtre <- flowframe_entree[ratio_Y_X < seuil_statistique & is.finite(ratio_Y_X), ] # Retire les doublets et valeurs instables
        self$post_doublets_SSC[[nom]] <- flowframe_filtre # Stocke l’échantillon filtré
        self$post_doublets_final[[nom]] <- flowframe_filtre # Met à jour la version finale des singlets
      } 
      
      noms_a_traiter <- if (is.null(nom_echantillon) || length(self$gate_doublets_SSC) == 0) names(liste_fcs_source) else nom_echantillon # Traite tout au premier passage, puis échantillon ciblé si déjà initialisé
      for (nom in noms_a_traiter) { calculer_doublets(nom, label_source) } # Applique la détection des doublets à chaque échantillon
      
      self$mettre_a_jour_pipeline("doublets_SSC", nom_echantillon) # Met à jour l’état du pipeline
    },
    
    # ============= gate_les_doublets_FSC =====================
    gate_les_doublets_FSC = function(points_utilisateur, axe_discrimination = "H_A", nom_echantillon = NULL) { # Gate manuel pour exclure les doublets via FSC
      
      liste_fcs_source <- if (!is.null(self$post_debris) && length(self$post_debris) > 0) self$post_debris else self$obtenir_derniere_source() # Utilise les données post-débris si disponibles, sinon la dernière source QC
      if (length(liste_fcs_source) == 0) stop("Aucune donnée source disponible pour le gating des doublets FSC.") # Stoppe si aucune donnée exploitable
      
      premier_ff  <- liste_fcs_source[[1]] # Récupère un flowFrame pour inspecter les noms de canaux
      tous_canaux <- flowCore::colnames(premier_ff) # Liste des canaux disponibles
      
      if (axe_discrimination == "H_A") { # Choix des axes FSC-H vs FSC-A
        canal_x <- grep("FSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("FSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "W_A") { # Choix des axes FSC-W vs FSC-A
        canal_x <- grep("FSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("FSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "H_W") { # Choix des axes FSC-H vs FSC-W
        canal_x <- grep("FSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("FSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      }
      
      matrice_coords <- as.matrix(points_utilisateur[, 1:2]) # Convertit les points utilisateur en matrice
      colnames(matrice_coords) <- c(canal_x, canal_y) # Associe les noms des canaux aux colonnes du polygone
      
      poly_fsc <- flowCore::polygonGate(.gate = matrice_coords, filterId = "Gate_Doublets_FSC") # Crée le polygonGate
      
      appliquer_fsc = function(nom) { # Fonction interne appliquant le gate à un échantillon
        ff_entree  <- liste_fcs_source[[nom]] # Récupère l’échantillon
        if (is.null(ff_entree)) return(NULL) # Ignore si absent
        
        res_filtre <- flowCore::filter(ff_entree, poly_fsc) # Applique le gate polygonal
        self$gate_doublets_FSC[[nom]] <- list(type = "poly", gate = poly_fsc, channels = c(canal_x, canal_y)) # Stocke le gate et ses métadonnées
        ff_propre <- ff_entree[res_filtre@subSet, ] # Conserve uniquement les événements dans le polygone
        self$post_doublets_FSC[[nom]]    <- ff_propre # Stocke l’échantillon filtré
        self$post_doublets_final[[nom]]  <- ff_propre # Met à jour la version finale des singlets
      } 
      
      noms <- if (is.null(nom_echantillon) || length(self$gate_doublets_FSC) == 0) names(liste_fcs_source) else nom_echantillon # Traite tout au premier passage, puis échantillon ciblé si déjà initialisé
      for (n in noms) { appliquer_fsc(n) } # Applique le gate à chaque échantillon
      
      if (!is.null(self$mettre_a_jour_pipeline)) self$mettre_a_jour_pipeline("doublets_FSC", nom_echantillon) # Met à jour l’état du pipeline
    },
    
    # ============= gate_les_doublets_SSC =====================
    gate_les_doublets_SSC = function(points_utilisateur, axe_discrimination = "H_A", nom_echantillon = NULL) { # Gate manuel pour exclure les doublets via SSC
      
      liste_fcs_source <- if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) { # Utilise les données post-FSC si disponibles
        self$post_doublets_FSC
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) { # Sinon utilise les données post-débris
        self$post_debris
      } else {
        self$obtenir_derniere_source() # Sinon utilise la dernière source QC disponible
      } 
      if (length(liste_fcs_source) == 0) stop("Aucune donnée source disponible pour le gating des doublets SSC.") # Stoppe si aucune donnée exploitable
      
      premier_ff  <- liste_fcs_source[[1]] # Récupère un flowFrame pour inspecter les noms de canaux
      tous_canaux <- flowCore::colnames(premier_ff) # Liste des canaux disponibles
      
      if (axe_discrimination == "H_A") { # Choix des axes SSC-H vs SSC-A
        canal_x <- grep("SSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("SSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "W_A") { # Choix des axes SSC-W vs SSC-A
        canal_x <- grep("SSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("SSC-A", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      } else if (axe_discrimination == "H_W") { # Choix des axes SSC-H vs SSC-W
        canal_x <- grep("SSC-H", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
        canal_y <- grep("SSC-W", tous_canaux, value = TRUE, ignore.case = TRUE)[1]
      }
      
      matrice_coords <- as.matrix(points_utilisateur[, 1:2]) # Convertit les points utilisateur en matrice
      colnames(matrice_coords) <- c(canal_x, canal_y) # Associe les noms des canaux aux colonnes du polygone
      
      poly_ssc <- flowCore::polygonGate(.gate = matrice_coords, filterId = "Gate_Doublets_SSC") # Crée le polygonGate
      
      appliquer_ssc = function(nom) { # Fonction interne appliquant le gate à un échantillon
        ff_entree <- liste_fcs_source[[nom]] # Récupère l’échantillon
        if (is.null(ff_entree)) return(NULL) # Ignore si absent
        
        res_filtre <- flowCore::filter(ff_entree, poly_ssc) # Applique le gate polygonal
        self$gate_doublets_SSC[[nom]] <- list(type = "poly", gate = poly_ssc, channels = c(canal_x, canal_y)) # Stocke le gate et ses métadonnées
        ff_propre <- ff_entree[res_filtre@subSet, ] # Conserve uniquement les événements dans le polygone
        self$post_doublets_SSC[[nom]]   <- ff_propre # Stocke l’échantillon filtré
        self$post_doublets_final[[nom]] <- ff_propre # Met à jour la version finale des singlets
      } 
      
      noms <- if (is.null(nom_echantillon) || length(self$gate_doublets_SSC) == 0) names(liste_fcs_source) else nom_echantillon # Traite tout au premier passage, puis échantillon ciblé si déjà initialisé
      for (n in noms) { appliquer_ssc(n) } # Applique le gate à chaque échantillon
      
      if (!is.null(self$mettre_a_jour_pipeline)) self$mettre_a_jour_pipeline("doublets_SSC", nom_echantillon) # Met à jour l’état du pipeline
    },
    
    # ============= visualiser_peacoqc =====================
    visualiser_peacoqc = function(nom_echantillon) { # Visualisation comparative avant/après PeacoQC
      
      if (is.null(self$post_PeacoQC[[nom_echantillon]])) { # Vérifie que PeacoQC a été exécuté pour cet échantillon
        message("Pas de données PeacoQC pour ", nom_echantillon) # Avertit si aucun résultat QC
        return(NULL) # Sortie propre
      } 
      
      flowframe_initial <- self$echantillons_traites[[nom_echantillon]] # Données brutes avant PeacoQC
      flowframe_nettoye <- self$post_PeacoQC[[nom_echantillon]] # Données filtrées par PeacoQC
      
      donnees_initiales <- as.data.frame(flowCore::exprs(flowframe_initial)) # Convertit les intensités brutes en data.frame
      donnees_nettoyees <- as.data.frame(flowCore::exprs(flowframe_nettoye)) # Convertit les intensités filtrées en data.frame
      
      canal_temps <- grep("time", colnames(donnees_initiales), value = TRUE, ignore.case = TRUE)[1] # Identifie le canal temporel
      canal_fsc   <- grep("FSC",  colnames(donnees_initiales), value = TRUE, ignore.case = TRUE)[1] # Identifie un canal FSC
      
      if (is.na(canal_temps)) canal_temps <- colnames(donnees_initiales)[1] # Fallback si aucun canal time
      if (is.na(canal_fsc))   canal_fsc   <- colnames(donnees_initiales)[2] # Fallback si aucun canal FSC
      
      limites_temps <- range(donnees_initiales[[canal_temps]], na.rm = TRUE) # Limites X
      limites_fsc   <- range(donnees_initiales[[canal_fsc]],   na.rm = TRUE) # Limites Y
      
      total_evenements  <- nrow(donnees_initiales) # Nombre total d’événements bruts
      evenements_gardes <- nrow(donnees_nettoyees) # Nombre d’événements conservés
      pourcentage_conservation <- if(total_evenements > 0) round((evenements_gardes / total_evenements) * 100, 1) else 0 # Rendement QC
      
      colonnes_cle <- intersect(colnames(donnees_initiales), colnames(donnees_nettoyees)) # Colonnes communes pour identifier les exclus
      
      if (evenements_gardes > 0 && length(colonnes_cle) > 0) { # Si comparaison possible
        cle_init <- do.call(paste, c(lapply(colonnes_cle, function(col) donnees_initiales[[col]]), sep = "\r")) # Clé unique brute
        cle_nett <- do.call(paste, c(lapply(colonnes_cle, function(col) donnees_nettoyees[[col]]), sep = "\r")) # Clé unique filtrée
        donnees_exclues <- donnees_initiales[!(cle_init %in% cle_nett), ] # Événements supprimés par PeacoQC
      } else {
        donnees_exclues <- donnees_initiales # Si comparaison impossible, tout est considéré exclu
      }
      
      lbl_x <- if (!is.null(self$obtenir_label)) self$obtenir_label(flowframe_initial, canal_temps) else canal_temps # Label X
      lbl_y <- if (!is.null(self$obtenir_label)) self$obtenir_label(flowframe_initial, canal_fsc)   else canal_fsc   # Label Y
      
      df_densite <- calculer_densite_raster(donnees_nettoyees[[canal_temps]], donnees_nettoyees[[canal_fsc]], limites_temps, limites_fsc) # Densité des événements conservés
      
      graphique_qc <- ggplot2::ggplot() # Initialise le plot
      
      if (!is.null(df_densite)) { # Ajoute la densité si disponible
        graphique_qc <- graphique_qc +
          ggplot2::geom_raster(data = df_densite, ggplot2::aes(x = X, y = Y, fill = densite), interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS)
      }
      
      if (nrow(donnees_exclues) > 0) { # Ajoute les points exclus
        graphique_qc <- graphique_qc +
          ggplot2::geom_point(
            data = donnees_exclues,
            ggplot2::aes(x = .data[[canal_temps]], y = .data[[canal_fsc]]),
            size = 0.3, alpha = 0.6, color = "darkred"
          )
      }
      
      graphique_qc <- graphique_qc +
        ggplot2::coord_cartesian(xlim = limites_temps, ylim = limites_fsc) + # Respecte les limites originales
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "none",
          aspect.ratio = 1,
          plot.title = ggplot2::element_text(face = "bold"),
          plot.subtitle = ggplot2::element_text(color = "darkblue", size = 11)
        ) +
        ggplot2::labs(
          title = paste("Contrôle qualité PeacoQC :", nom_echantillon),
          subtitle = paste0("Événements conservés : ", format(evenements_gardes, big.mark=" "), " | ", pourcentage_conservation, "%"),
          x = lbl_x, y = lbl_y
        )
      
      if (is.null(self$plots_peacoqc)) self$plots_peacoqc <- list() # Initialise le stockage des plots si nécessaire
      self$plots_peacoqc[[nom_echantillon]] <- graphique_qc # Archive le graphique
      
      return(graphique_qc) # Renvoie le plot final
    },
    
    # ============= visualiser_flowai =====================
    visualiser_flowai = function(nom_echantillon) { # Visualisation avant/après flowAI
      
      if (is.null(self$post_flowAI) || is.null(self$post_flowAI[[nom_echantillon]])) { # Vérifie que flowAI a été exécuté
        stop("Aucun résultat flowAI trouvé pour cet échantillon.") # Stoppe si aucun QC flowAI disponible
      } 
      
      ff_nettoye <- self$post_flowAI[[nom_echantillon]] # Données filtrées par flowAI
      fcs_initial <- self$echantillons_traites[[nom_echantillon]] # Données brutes unmixed
      
      exprs_initiales <- as.data.frame(flowCore::exprs(fcs_initial)) # Convertit les intensités brutes en data.frame
      
      canal_temps <- grep("time", colnames(exprs_initiales), ignore.case = TRUE, value = TRUE)[1] # Détecte le canal temporel
      canal_taille <- grep("FSC", colnames(exprs_initiales), ignore.case = TRUE, value = TRUE)[1] # Détecte un canal FSC
      
      if (is.na(canal_temps) || is.na(canal_taille)) { # Fallback si les canaux ne sont pas trouvés
        canal_temps  <- colnames(exprs_initiales)[ncol(exprs_initiales)]
        canal_taille <- colnames(exprs_initiales)[1]
      } 
      
      exprs_nettoyees <- as.data.frame(flowCore::exprs(ff_nettoye)) # Intensités des événements conservés
      
      exprs_initiales$Status <- "Éliminé (flowAI)" # Marque tous les événements comme éliminés par défaut
      
      indices_conserves <- which(exprs_initiales[[canal_temps]] %in% exprs_nettoyees[[canal_temps]]) # Identifie les événements conservés via le canal temporel
      
      exprs_initiales$Status[indices_conserves] <- "Conservé" # Marque les événements conservés
      
      total_pts <- nrow(exprs_initiales) # Nombre total d’événements bruts
      total_conserves <- length(indices_conserves) # Nombre d’événements conservés
      pourcentage_conservation <- round((total_conserves / total_pts) * 100, 1) # Rendement QC
      
      limites_temps  <- range(exprs_initiales[[canal_temps]],  na.rm = TRUE) # Limites X
      limites_taille <- range(exprs_initiales[[canal_taille]], na.rm = TRUE) # Limites Y
      
      donnees_conservees <- exprs_initiales[exprs_initiales$Status == "Conservé", ] # Sous-ensemble conservé
      donnees_eliminees  <- exprs_initiales[exprs_initiales$Status == "Éliminé (flowAI)", ] # Sous-ensemble éliminé
      
      df_densite <- calculer_densite_raster(donnees_conservees[[canal_temps]], donnees_conservees[[canal_taille]], limites_temps, limites_taille) # Densité des événements conservés
      
      graphique_flowai <- ggplot2::ggplot() # Initialise le graphique
      
      if (!is.null(df_densite)) { # Ajoute la densité si disponible
        graphique_flowai <- graphique_flowai +
          ggplot2::geom_raster(data = df_densite, ggplot2::aes(x = X, y = Y, fill = densite), interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS)
      }
      
      if (nrow(donnees_eliminees) > 0) { # Ajoute les points éliminés
        graphique_flowai <- graphique_flowai +
          ggplot2::geom_point(
            data = donnees_eliminees,
            ggplot2::aes(x = .data[[canal_temps]], y = .data[[canal_taille]]),
            size = 0.4, alpha = 0.6, color = "#d62728"
          )
      }
      
      graphique_flowai <- graphique_flowai +
        ggplot2::coord_cartesian(xlim = limites_temps, ylim = limites_taille) + # Respecte les limites originales
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = paste("Contrôle Qualité flowAI :", nom_echantillon),
          subtitle = paste0("Événements conservés : ", format(total_conserves, big.mark=" "), " / ",
                            format(total_pts, big.mark=" "), " (", pourcentage_conservation, "%)"),
          x = paste("Axe du Temps :", canal_temps),
          y = paste("Axe Morphologique :", canal_taille)
        ) +
        ggplot2::theme(
          legend.position = "none",
          aspect.ratio = 1,
          plot.title = ggplot2::element_text(face = "bold", size = 14)
        )
      
      if (is.null(self$plots_flowai)) self$plots_flowai <- list() # Initialise le stockage des graphiques flowAI si nécessaire
      self$plots_flowai[[nom_echantillon]] <- graphique_flowai # Archive le graphique
      
      return(graphique_flowai) # Renvoie le plot final
    },
    
    # ============= visualiser_debris =====================
    visualiser_debris = function(nom_echantillon) { # Visualisation du gate Débris appliqué à un échantillon
      
      if (is.null(self$post_debris[[nom_echantillon]])) { # Vérifie que l’étape Débris a été exécutée
        message("Pas de données Débris pour ", nom_echantillon) # Avertit si aucun résultat
        return(NULL) # Sortie propre
      }
      
      flowframe_avant <- if (!is.null(self$post_retrait_bordures) && length(self$post_retrait_bordures) > 0) { # Priorité : données post-bordures
        self$post_retrait_bordures[[nom_echantillon]]
      } else if (!is.null(self$post_PeacoQC) && length(self$post_PeacoQC) > 0) { # Sinon données post-PeacoQC
        self$post_PeacoQC[[nom_echantillon]]
      } else if (!is.null(self$post_flowAI) && length(self$post_flowAI) > 0) { # Sinon données post-flowAI
        self$post_flowAI[[nom_echantillon]]
      } else {
        self$echantillons_traites[[nom_echantillon]] # Sinon données unmixed brutes
      }
      
      if (is.null(flowframe_avant)) return(NULL) # Stoppe si aucune donnée exploitable
      
      flowframe_apres <- self$post_debris[[nom_echantillon]] # Données filtrées par le gate Débris
      
      donnees_globales <- as.data.frame(flowCore::exprs(flowframe_avant)) # Convertit les intensités avant filtrage
      gate_polygone    <- self$gate_debris[[nom_echantillon]] # Récupère le polygonGate utilisé
      coordonnees_gate <- as.data.frame(gate_polygone@boundaries) # Coordonnées du polygone
      colnames(coordonnees_gate) <- c("x", "y") # Noms des colonnes
      
      params   <- flowCore::parameters(gate_polygone) # Récupère les noms des canaux utilisés
      canal_x  <- params[1] # Canal X du gate
      canal_y  <- params[2] # Canal Y du gate
      
      total_evenements_avant <- nrow(donnees_globales) # Nombre d’événements avant filtrage
      total_evenements_apres <- nrow(flowCore::exprs(flowframe_apres)) # Nombre d’événements après filtrage
      
      pourcentage_conservation <- if (total_evenements_avant > 0)
        round((total_evenements_apres / total_evenements_avant) * 100, 1) else 0 # Rendement du gate
      
      lbl_x <- self$obtenir_label(flowframe_avant, canal_x) # Label biologique X
      lbl_y <- self$obtenir_label(flowframe_avant, canal_y) # Label biologique Y
      
      lim_x_globale <- c(0, max(donnees_globales[[canal_x]], na.rm = TRUE)) # Limites X
      lim_y_globale <- c(0, max(donnees_globales[[canal_y]], na.rm = TRUE)) # Limites Y
      
      df_densite <- calculer_densite_raster(donnees_globales[[canal_x]], donnees_globales[[canal_y]], lim_x_globale, lim_y_globale) # Densité globale
      
      graphique_debris <- if (is.null(df_densite)) { # Cas où la densité ne peut pas être calculée
        ggplot2::ggplot() + ggplot2::theme_bw() +
          ggplot2::labs(
            title = paste("Nettoyage des débris :", nom_echantillon),
            subtitle = "Pas assez d'événements pour tracer la densité",
            x = lbl_x, y = lbl_y
          )
      } else {
        ggplot2::ggplot(df_densite, ggplot2::aes(x = X, y = Y, fill = densite)) + # Plot densité
          ggplot2::geom_raster(interpolate = TRUE) + # Raster haute densité
          ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS) +
          ggplot2::geom_polygon( # Overlay du gate Débris
            data = coordonnees_gate,
            ggplot2::aes(x = x, y = y),
            fill = NA, color = "black", linewidth = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::coord_cartesian(xlim = lim_x_globale, ylim = lim_y_globale) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            legend.position = "none",
            aspect.ratio = 1,
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
      
      if (is.null(self$plots_debris)) self$plots_debris <- list() # Initialise le stockage des plots si nécessaire
      self$plots_debris[[nom_echantillon]] <- graphique_debris # Archive le graphique
      
      return(graphique_debris) # Renvoie le plot final
    },
    
    # ============= visualiser_doublets =====================
    visualiser_doublets = function(nom_echantillon, type_analyse = "FSC") { # Visualisation des doublets FSC ou SSC
      
      infos_gate <- if (type_analyse == "FSC") self$gate_doublets_FSC[[nom_echantillon]] else self$gate_doublets_SSC[[nom_echantillon]] # Récupère les métadonnées du gate utilisé
      if (is.null(infos_gate)) return(NULL) # Stoppe si aucun gate n’a été appliqué
      
      ff_avant <- if (type_analyse == "FSC") { # Source amont pour FSC
        if (!is.null(self$post_debris) && length(self$post_debris) > 0) self$post_debris[[nom_echantillon]] else self$obtenir_derniere_source()[[nom_echantillon]]
      } else { # Source amont pour SSC
        if (!is.null(self$post_doublets_FSC) && length(self$post_doublets_FSC) > 0) self$post_doublets_FSC[[nom_echantillon]]
        else if (!is.null(self$post_debris) && length(self$post_debris) > 0) self$post_debris[[nom_echantillon]]
        else self$obtenir_derniere_source()[[nom_echantillon]]
      }
      if (is.null(ff_avant)) return(NULL) # Stoppe si aucune donnée amont
      
      ff_apres <- if (type_analyse == "FSC") self$post_doublets_FSC[[nom_echantillon]] else self$post_doublets_SSC[[nom_echantillon]] # Données filtrées
      donnees_source <- as.data.frame(flowCore::exprs(ff_avant)) # Données avant filtrage
      nb_avant <- nrow(donnees_source) # Nombre d’événements avant filtrage
      donnees_visu <- donnees_source # Copie pour annotation
      canal_x <- infos_gate$channels[1] # Canal X utilisé pour le gate
      canal_y <- infos_gate$channels[2] # Canal Y utilisé pour le gate
      nb_apres <- if (!is.null(ff_apres)) nrow(flowCore::exprs(ff_apres)) else 0 # Nombre d’événements conservés
      pourcentage <- if (nb_avant > 0) round((nb_apres / nb_avant) * 100, 1) else 0 # Rendement
      limite_max <- max(donnees_source[, c(canal_x, canal_y)], na.rm = TRUE) # Limite max pour cadrage
      lbl_x <- if (!is.null(self$obtenir_label)) self$obtenir_label(ff_avant, canal_x) else canal_x # Label X
      lbl_y <- if (!is.null(self$obtenir_label)) self$obtenir_label(ff_avant, canal_y) else canal_y # Label Y
      
      if (!is.null(ff_apres) && nb_apres > 0) { # Détermination conservé / retiré
        donnees_apres <- as.data.frame(flowCore::exprs(ff_apres)) # Données filtrées
        colonnes_cle <- intersect(colnames(donnees_visu), colnames(donnees_apres)) # Colonnes communes
        cle_visu  <- do.call(paste, c(lapply(colonnes_cle, function(col) donnees_visu[[col]]),  sep = "\r")) # Clé brute
        cle_apres <- do.call(paste, c(lapply(colonnes_cle, function(col) donnees_apres[[col]]), sep = "\r")) # Clé filtrée
        donnees_visu$statut_doublet <- ifelse(cle_visu %in% cle_apres, "Conservé", "Retiré") # Annotation
      } else {
        donnees_visu$statut_doublet <- "Conservé" # Si aucun filtrage
      }
      
      donnees_retirees   <- donnees_visu[donnees_visu$statut_doublet == "Retiré", ] # Points retirés
      donnees_conservees <- donnees_visu[donnees_visu$statut_doublet == "Conservé", ] # Points conservés
      couleur_retire <- "darkred" # Couleur des doublets retirés
      couleur_legende_conserve <- "#2C7FB8" # Couleur de légende pour conservés
      df_densite_cons <- calculer_densite_raster(donnees_conservees[[canal_x]], donnees_conservees[[canal_y]], c(0, limite_max), c(0, limite_max)) # Densité des conservés
      
      graphique <- ggplot2::ggplot() # Plot initial
      
      if (!is.null(df_densite_cons)) {
        graphique <- graphique +
          ggplot2::geom_raster(data = df_densite_cons, ggplot2::aes(x = X, y = Y, fill = densite), interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS)
      }
      
      if (nrow(donnees_retirees) > 0) {
        graphique <- graphique +
          ggplot2::geom_point(data = donnees_retirees, ggplot2::aes(x = .data[[canal_x]], y = .data[[canal_y]]),
                              color = couleur_retire, size = 0.3, alpha = 0.6)
      }
      
      graphique <- graphique +
        ggplot2::coord_cartesian(xlim = c(0, limite_max), ylim = c(0, limite_max)) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = paste("Retrait des doublets", type_analyse, ":", nom_echantillon),
          subtitle = paste0("Événements conservés : ", format(nb_apres, big.mark=" "), " | ", pourcentage, "%"),
          x = lbl_x, y = lbl_y
        ) +
        ggplot2::theme(legend.position = "none", aspect.ratio = 1,
                       plot.title = element_text(face = "bold"),
                       plot.subtitle = element_text(color = "darkblue", size = 11))
      
      if (infos_gate$type == "poly") { # Ajout du polygone si gate manuel
        coordonnees_gate <- as.data.frame(infos_gate$gate@boundaries)
        colnames(coordonnees_gate) <- c("x", "y")
        graphique <- graphique +
          ggplot2::geom_polygon(data = coordonnees_gate, ggplot2::aes(x = x, y = y),
                                fill = NA, color = "darkred", linewidth = 0.6)
      }
      
      pas_legende <- limite_max * 0.05 # Espacement légende
      y_conserve  <- limite_max * 0.97 # Position conservé
      y_retire    <- y_conserve - pas_legende # Position retiré
      x_pastille  <- limite_max * 0.80 # Position pastille
      x_texte     <- limite_max * 0.83 # Position texte
      
      graphique <- graphique +
        ggplot2::annotate("point", x = x_pastille, y = y_conserve, color = couleur_legende_conserve, size = 2.2) +
        ggplot2::annotate("text",  x = x_texte,    y = y_conserve, label = "Conservé", hjust = 0, size = 3.2) +
        ggplot2::annotate("point", x = x_pastille, y = y_retire,   color = couleur_retire, size = 2.2) +
        ggplot2::annotate("text",  x = x_texte,    y = y_retire,   label = "Retiré",   hjust = 0, size = 3.2)
      
      if (is.null(self$plots_doublets)) self$plots_doublets <- list() # Initialise le stockage
      self$plots_doublets[[paste0(nom_echantillon, "_", type_analyse)]] <- graphique # Archive
      
      return(graphique) # Renvoie le plot final
    },
    
    # ============= retirer_les_cellules_mortes =====================
    retirer_les_cellules_mortes = function(canal_fsc = "FSC-A", marqueur_viabilite,
                                           points_utilisateur, nom_echantillon = NULL) { # Gate manuel pour exclure les cellules mortes
      
      if (is.null(points_utilisateur)) stop("Aucun point fourni pour la viabilité.") # Vérifie que des points ont été fournis
      if (nrow(points_utilisateur) < 3) stop("Un polygone nécessite ≥ 3 points.") # Un polygonGate doit avoir au moins 3 sommets
      
      liste_fcs_source <- if (!is.null(self$post_transformation) && length(self$post_transformation) > 0) { # Priorité : données transformées
        self$post_transformation
      } else if (!is.null(self$post_doublets_final) && length(self$post_doublets_final) > 0) { # Sinon données post-doublets
        self$post_doublets_final
      } else if (!is.null(self$post_debris) && length(self$post_debris) > 0) { # Sinon données post-débris
        self$post_debris
      } else {
        self$obtenir_derniere_source() # Sinon dernière source QC disponible
      }
      
      if (is.null(self$gate_viabilite)) self$gate_viabilite <- list() # Initialise le stockage des gates viabilité
      if (is.null(self$post_viabilite)) self$post_viabilite <- list() # Initialise le stockage des flowFrames filtrés
      
      matrice_coordonnees <- as.matrix(points_utilisateur[, 1:2]) # Convertit les points utilisateur en matrice
      colnames(matrice_coordonnees) <- c(canal_fsc, marqueur_viabilite) # Associe les noms des canaux
      polygone_viabilite <- flowCore::polygonGate(.gate = matrice_coordonnees, filterId = "Gate_Viabilite") # Crée le polygonGate
      
      appliquer_le_gate_vivantes = function(nom) { # Applique le gate viabilité à un échantillon
        ff <- liste_fcs_source[[nom]] # Récupère l’échantillon
        if (is.null(ff)) return(NULL) # Ignore si absent
        
        resultat_filtrage <- flowCore::filter(ff, polygone_viabilite) # Applique le gate
        
        self$gate_viabilite[[nom]] <- polygone_viabilite # Stocke le gate utilisé
        self$post_viabilite[[nom]] <- ff[resultat_filtrage@subSet, ] # Conserve uniquement les cellules vivantes
      }
      
      echantillons <- if (is.null(nom_echantillon) || length(self$gate_viabilite) == 0) names(liste_fcs_source) else nom_echantillon # Traite tout au premier passage, puis échantillon ciblé si déjà initialisé
      
      for (nom in echantillons) appliquer_le_gate_vivantes(nom) # Applique le gate à chaque échantillon
      
      self$mettre_a_jour_pipeline("viabilite", nom_echantillon) # Met à jour l’état du pipeline
    },
    
    # ============= visualiser_viabilite =====================
    visualiser_viabilite = function(nom_echantillon) { # Visualisation du gate de viabilité appliqué à un échantillon
      
      if (is.null(self$post_viabilite[[nom_echantillon]])) { # Vérifie que l’étape viabilité a été exécutée
        message("Pas de données Viabilité pour ", nom_echantillon) # Avertit si aucun résultat
        return(NULL) # Sortie propre
      }
      
      flowframe_avant <- if (!is.null(self$post_transformation) &&
                             !is.null(self$post_transformation[[nom_echantillon]])) { # Priorité : données transformées
        self$post_transformation[[nom_echantillon]]
      } else if (!is.null(self$post_doublets_final[[nom_echantillon]])) { # Sinon données post-doublets
        self$post_doublets_final[[nom_echantillon]]
      } else if (!is.null(self$post_debris[[nom_echantillon]])) { # Sinon données post-débris
        self$post_debris[[nom_echantillon]]
      } else {
        self$obtenir_derniere_source()[[nom_echantillon]] # Sinon dernière source QC disponible
      }
      
      if (is.null(flowframe_avant)) return(NULL) # Stoppe si aucune donnée amont
      
      donnees_globales <- as.data.frame(flowCore::exprs(flowframe_avant)) # Convertit les intensités en data.frame
      
      total_evenements_avant <- nrow(donnees_globales) # Nombre d’événements avant filtrage
      donnees_visu <- donnees_globales # Copie pour visualisation
      
      gate_polygone <- self$gate_viabilite[[nom_echantillon]] # Récupère le gate utilisé
      
      if (inherits(gate_polygone, "polygonGate")) { # Cas d’un gate polygonal
        bound <- gate_polygone@boundaries # Coordonnées du polygone
        if (is.null(dim(bound)) || ncol(bound) != 2) stop("Le gate polygonal ne contient pas une matrice Nx2 de boundaries.") # Vérifie la structure
        coordonnees_gate <- data.frame(x = bound[, 1], y = bound[, 2]) # Convertit en data.frame
        noms_canaux <- colnames(bound) # Récupère les noms des canaux
      } else if (inherits(gate_polygone, "rectangleGate")) { # Cas d’un gate rectangulaire
        mins <- gate_polygone@min # Bornes min
        maxs <- gate_polygone@max # Bornes max
        noms_canaux <- names(mins) # Noms des canaux
        coordonnees_gate <- data.frame( # Convertit le rectangle en polygone
          x = c(mins[1], maxs[1], maxs[1], mins[1]),
          y = c(mins[2], mins[2], maxs[2], maxs[2])
        )
      } else if (inherits(gate_polygone, "filterResult")) { # Cas d’un filterResult contenant un gate interne
        gate_interne <- gate_polygone@filter
        if (inherits(gate_interne, "polygonGate")) {
          bound <- gate_interne@boundaries
          coordonnees_gate <- data.frame(x = bound[, 1], y = bound[, 2])
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
          stop("Type de gate non supporté : ", class(gate_interne)) # Stoppe si gate interne non reconnu
        }
      } else {
        stop("Type de gate non supporté : ", class(gate_polygone)) # Stoppe si gate externe non reconnu
      }
      
      canal_x <- as.character(noms_canaux[1]) # Canal X utilisé pour la viabilité
      canal_y <- as.character(noms_canaux[2]) # Canal Y utilisé pour la viabilité
      
      if (!(canal_x %in% colnames(donnees_visu))) stop("Canal X introuvable dans les données : ", canal_x) # Vérifie la présence du canal X
      if (!(canal_y %in% colnames(donnees_visu))) stop("Canal Y introuvable dans les données : ", canal_y) # Vérifie la présence du canal Y
      
      flowframe_apres <- self$post_viabilite[[nom_echantillon]] # Données filtrées
      total_evenements_apres <- nrow(flowCore::exprs(flowframe_apres)) # Nombre d’événements conservés
      pourcentage_conservation <- round((total_evenements_apres / total_evenements_avant) * 100, 1) # Rendement
      
      lbl_x <- if (!is.null(self$obtenir_label)) self$obtenir_label(flowframe_avant, canal_x) else canal_x # Label X
      lbl_y <- if (!is.null(self$obtenir_label)) self$obtenir_label(flowframe_avant, canal_y) else canal_y # Label Y
      
      lim_x_viab <- range(donnees_globales[[canal_x]], na.rm = TRUE) # Limites X
      lim_y_viab <- range(donnees_globales[[canal_y]], na.rm = TRUE) # Limites Y
      
      df_densite_viab <- calculer_densite_raster(donnees_visu[[canal_x]], donnees_visu[[canal_y]], lim_x_viab, lim_y_viab) # Densité
      
      graphique_viabilite <- if (is.null(df_densite_viab)) { # Cas où la densité ne peut pas être calculée
        ggplot2::ggplot() + ggplot2::theme_bw() +
          ggplot2::labs(
            title = paste("Retrait des cellules mortes :", nom_echantillon),
            subtitle = "Pas assez d'événements pour tracer la densité",
            x = lbl_x, y = lbl_y
          )
      } else {
        ggplot2::ggplot(df_densite_viab, ggplot2::aes(x = X, y = Y, fill = densite)) + # Plot densité
          ggplot2::geom_raster(interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS) +
          ggplot2::geom_polygon( # Overlay du gate viabilité
            data = coordonnees_gate,
            ggplot2::aes(x = x, y = y),
            fill = NA, color = "grey20", linewidth = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::coord_cartesian(xlim = lim_x_viab, ylim = lim_y_viab) +
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
            aspect.ratio = 1,
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(color = "darkblue", size = 11)
          )
      }
      
      if (is.null(self$plots_viabilite)) self$plots_viabilite <- list() # Initialise le stockage des plots
      self$plots_viabilite[[nom_echantillon]] <- graphique_viabilite # Archive le graphique
      
      return(graphique_viabilite) # Renvoie le plot final
    },
    
    # ============= transformation_arcsinh =====================
    transformation_arcsinh = function(canaux = "", echantillon = NULL, cofactor = 400) { # Applique une transformation arcsinh sur les canaux choisis
      
      noms_a_traiter <- if (is.null(echantillon)) names(self$echantillons_traites) else echantillon # Détermine quels échantillons transformer
      
      arc_sinh_transform <- function(x, cf) { # Fonction interne appliquant arcsinh(x / cofactor)
        asinh(x / cf)
      }
      
      self$cofactor_transformation <- cofactor # Stocke le cofacteur utilisé pour la transformation
      
      for (nom in noms_a_traiter) { # Boucle sur les échantillons à transformer
        ff <- if (!is.null(self$post_doublets_final[[nom]])) { # Priorité : données post-doublets
          self$post_doublets_final[[nom]]
        } else if (!is.null(self$post_debris[[nom]])) { # Sinon données post-débris
          self$post_debris[[nom]]
        } else {
          self$echantillons_traites[[nom]] # Sinon données unmixed brutes
        }
        
        if (is.null(ff)) next # Ignore si l’échantillon est absent
        
        matrice_exprs <- flowCore::exprs(ff) # Récupère la matrice des intensités
        
        canaux_a_transformer <- if (canaux[1] == "") colnames(matrice_exprs) else canaux # Si aucun canal spécifié → tous les canaux
        canaux_existants <- intersect(canaux_a_transformer, colnames(matrice_exprs)) # Filtre les canaux réellement présents
        
        if (length(canaux_existants) == 0) next # Ignore si aucun canal valide
        
        matrice_exprs[, canaux_existants] <- arc_sinh_transform(matrice_exprs[, canaux_existants], cofactor) # Applique arcsinh sur les canaux sélectionnés
        
        ff_transforme <- ff # Copie du flowFrame
        flowCore::exprs(ff_transforme) <- matrice_exprs # Remplace la matrice d’expression par la version transformée
        
        if (is.null(self$post_transformation)) self$post_transformation <- list() # Initialise la liste si nécessaire
        self$post_transformation[[nom]] <- ff_transforme # Stocke l’échantillon transformé
      }
      
      invisible(self) # Retourne l’objet sans affichage
    },
    
    # ===========================================
    # SECTION ANALYSES — gestion des gates personnalisés (créés dans l'onglet Analyses)
    # ===========================================
    
    # ============= creer_gate =====================
    creer_gate = function(nom_gate, type = "polygon", axes = c("FSC-A", "SSC-A"), points = NULL, gate_parent = NULL, cofacteur = 150, nom_echantillon = NULL) {
      
      if (is.null(self$gates_personnalisees)) self$gates_personnalisees <- list()
      
      deja_initialise <- !is.null(self$gates_personnalisees[[nom_gate]]) && length(self$gates_personnalisees[[nom_gate]]$formes) > 0
      
      # Validation géométrique et mise en forme des points selon le type
      if (type == "polygon") {
        if (is.null(points) || nrow(points) < 3) stop("Un polygone nécessite au moins 3 points.")
        mat_forme <- as.matrix(points[, 1:2])
        colnames(mat_forme) <- NULL # Les noms de colonnes sont réappliqués à la volée dans construire_gate_flowcore()
      } else if (type == "rectangle") {
        if (is.null(points) || length(points) != 4) stop("Le rectangle nécessite c(minX, maxX, minY, maxY).")
        mat_forme <- as.numeric(points)
      } else {
        stop("Type de gate invalide : utilisez 'polygon' ou 'rectangle'.")
      }
      
      if (!deja_initialise) {
        if (!is.null(gate_parent) && nchar(gate_parent) > 0) {
          if (identical(gate_parent, nom_gate)) stop("Un gate ne peut pas être son propre parent.")
          if (is.null(self$gates_personnalisees[[gate_parent]])) stop("Le gate parent '", gate_parent, "' n'existe pas. Créez-le d'abord.")
        } else {
          gate_parent <- NULL
        }
        
        source_par_echantillon <- if (is.null(gate_parent)) self$obtenir_derniere_source() else self$resoudre_population_gate(gate_parent)
        if (is.null(source_par_echantillon) || length(source_par_echantillon) == 0) {
          stop("Aucune donnée disponible pour créer ce gate (vérifiez le prétraitement ou le gate parent choisi).")
        }
        
        formes <- stats::setNames(rep(list(mat_forme), length(source_par_echantillon)), names(source_par_echantillon))
        
        self$gates_personnalisees[[nom_gate]] <- list(axes = axes, type = type, gate_parent = gate_parent, cofacteur = cofacteur, formes = formes)
        
        message("Succès : Gate '", nom_gate, "' créé sur les axes ", paste(axes, collapse = "/"),
                " (cofacteur ", cofacteur, "), appliqué par défaut à ", length(formes), " échantillon(s).",
                if (!is.null(gate_parent)) paste0(" Sous-population de '", gate_parent, "'.") else "")
      } else {
        infos  <- self$gates_personnalisees[[nom_gate]]
        formes <- infos$formes
        cibles <- if (is.null(nom_echantillon)) names(formes) else nom_echantillon
        for (nom in cibles) formes[[nom]] <- mat_forme
        self$gates_personnalisees[[nom_gate]]$formes <- formes
        
        if (is.null(nom_echantillon)) {
          message("Gate '", nom_gate, "' : forme réinitialisée pour tous les échantillons (", length(formes), ").")
        } else {
          message("Gate '", nom_gate, "' : forme mise à jour pour ", nom_echantillon, " uniquement.")
        }
      }
      
      return(invisible(self))
    },
    
    # ============= construire_gate_flowcore =====================
    construire_gate_flowcore = function(type, axes, mat_forme, nom_gate) {
      if (identical(type, "polygon")) {
        mat <- as.matrix(mat_forme)
        colnames(mat) <- axes # OBLIGATOIRE : flowCore::polygonGate() exige des noms de colonnes ("Matrix of gate boundaries must have colnames." sinon)
        flowCore::polygonGate(.gate = mat, filterId = nom_gate)
      } else {
        flowCore::rectangleGate(filterId = nom_gate,
                                .gate = matrix(c(mat_forme[1], mat_forme[2], mat_forme[3], mat_forme[4]),
                                               ncol = 2, byrow = TRUE, dimnames = list(NULL, axes)))
      }
    },
    
    # ============= obtenir_forme_gate =====================
    obtenir_forme_gate = function(nom_gate, nom_echantillon) {
      if (is.null(self$gates_personnalisees) || is.null(self$gates_personnalisees[[nom_gate]])) return(NULL)
      self$gates_personnalisees[[nom_gate]]$formes[[nom_echantillon]]
    },
    
    # ============= resoudre_population_gate =====================
    resoudre_population_gate = function(nom_gate, nom_echantillon = NULL) {
      if (is.null(self$gates_personnalisees) || is.null(self$gates_personnalisees[[nom_gate]])) {
        stop("Le gate '", nom_gate, "' n'existe pas. Créez-le d'abord avec creer_gate().")
      }
      infos_gate <- self$gates_personnalisees[[nom_gate]]
      cofacteur  <- infos_gate$cofacteur %||% 150 # Rétrocompatibilité : les gates créés avant l'ajout de ce réglage utilisent la valeur par défaut
      
      source_par_echantillon <- if (is.null(infos_gate$gate_parent)) {
        self$obtenir_derniere_source() # Pas de parent : part de l'étape la plus avancée disponible du pipeline
      } else {
        self$resoudre_population_gate(infos_gate$gate_parent) # Avec un parent : ne considère que la sous-population déjà sélectionnée par celui-ci (récursif sur toute la chaîne)
      }
      
      if (is.null(source_par_echantillon) || length(source_par_echantillon) == 0) return(list())
      
      trans_axes <- flowCore::transformList(infos_gate$axes, flowCore::arcsinhTransform(a = 0, b = 1 / cofacteur, c = 0)) # Même transformation (et même cofacteur) que celle utilisée pour le tracé interactif du gate
      
      noms_a_traiter <- if (is.null(nom_echantillon)) names(source_par_echantillon) else intersect(nom_echantillon, names(source_par_echantillon))
      
      resultat <- list()
      for (nom in noms_a_traiter) {
        ff <- source_par_echantillon[[nom]]
        if (is.null(ff) || nrow(flowCore::exprs(ff)) == 0) next
        if (!all(infos_gate$axes %in% flowCore::colnames(ff))) next # Sécurité : ignore un échantillon qui ne posséderait pas les 2 canaux du gate
        
        mat_forme <- infos_gate$formes[[nom]]
        if (is.null(mat_forme)) mat_forme <- infos_gate$formes[[1]] # Repli : échantillon ajouté après coup sans forme propre → utilise la première forme disponible comme valeur par défaut
        if (is.null(mat_forme)) next
        
        gate_echantillon <- self$construire_gate_flowcore(infos_gate$type, infos_gate$axes, mat_forme, nom_gate)
        
        ff_transforme   <- flowCore::transform(ff, trans_axes) # Applique la transformation uniquement pour tester l'appartenance au gate
        resultat_filtre <- flowCore::filter(ff_transforme, gate_echantillon)
        ff_gate         <- ff[resultat_filtre@subSet, ] # Sous-échantillonne le flowFrame D'ORIGINE (valeurs non transformées) avec les indices obtenus sur la version transformée
        if (nrow(flowCore::exprs(ff_gate)) > 0) resultat[[nom]] <- ff_gate
      }
      resultat
    },
    
    # ============= resumer_gate =====================
    resumer_gate = function(nom_gate) {
      if (is.null(self$gates_personnalisees) || is.null(self$gates_personnalisees[[nom_gate]])) {
        stop("Le gate '", nom_gate, "' n'existe pas. Créez-le d'abord avec creer_gate().")
      }
      infos_gate <- self$gates_personnalisees[[nom_gate]]
      
      population_gate  <- self$resoudre_population_gate(nom_gate)
      population_amont <- if (is.null(infos_gate$gate_parent)) self$obtenir_derniere_source() else self$resoudre_population_gate(infos_gate$gate_parent)
      
      noms <- union(names(population_gate), names(population_amont))
      if (length(noms) == 0) return(data.frame(echantillon = character(0), n_evenements = integer(0), pct_parent = numeric(0)))
      
      lignes <- lapply(noms, function(nom) {
        n_gate   <- if (!is.null(population_gate[[nom]]))  nrow(flowCore::exprs(population_gate[[nom]]))  else 0L
        n_amont  <- if (!is.null(population_amont[[nom]])) nrow(flowCore::exprs(population_amont[[nom]])) else 0L
        data.frame(
          echantillon  = nom,
          n_evenements = n_gate,
          pct_parent   = if (n_amont > 0) round(100 * n_gate / n_amont, 2) else NA_real_
        )
      })
      do.call(rbind, lignes)
    },
    
    # ============= supprimer_gate =====================
    supprimer_gate = function(nom_gate, cascade = TRUE) {
      if (is.null(self$gates_personnalisees) || is.null(self$gates_personnalisees[[nom_gate]])) {
        message("Le gate '", nom_gate, "' n'existe pas déjà plus.")
        return(invisible(self))
      }
      
      enfants_directs <- names(self$gates_personnalisees)[
        vapply(self$gates_personnalisees, function(g) identical(g$gate_parent, nom_gate), logical(1))
      ]
      
      if (length(enfants_directs) > 0) {
        if (!cascade) stop("Impossible de supprimer '", nom_gate, "' : les gates suivants en dépendent : ", paste(enfants_directs, collapse = ", "), ". Utilisez cascade = TRUE pour les supprimer aussi.")
        for (enfant in enfants_directs) self$supprimer_gate(enfant, cascade = TRUE) # Supprime récursivement toute la descendance avant de supprimer ce gate
      }
      
      self$gates_personnalisees[[nom_gate]] <- NULL
      if (!is.null(self$analyses_umap))     self$analyses_umap[[nom_gate]]     <- NULL
      if (!is.null(self$analyses_tsne))     self$analyses_tsne[[nom_gate]]     <- NULL
      if (!is.null(self$analyses_pca))      self$analyses_pca[[nom_gate]]      <- NULL
      if (!is.null(self$clusters_flowsom))  self$clusters_flowsom[[nom_gate]]  <- NULL
      
      message("Gate '", nom_gate, "' supprimé", if (length(enfants_directs) > 0) paste0(" (ainsi que ", length(enfants_directs), " gate(s) enfant(s))") else "", ".")
      return(invisible(self))
    },
    
    # ===========================================
    # SECTION ANALYSES — projections et clustering sur un gate personnalisé
    # ===========================================
    
    # ============= projection_UMAP =====================
    projection_UMAP = function(nom_gate, canaux = NULL, n_neighbors = 15, min_dist = 0.1,
                               n_components = 2, metric = "euclidean", sous_echantillonnage_max = 50000, cofacteur = 150) { # Réalise une projection UMAP (uwot) de la population sélectionnée par le gate
      if (!requireNamespace("uwot", quietly = TRUE)) stop("Le package 'uwot' est requis pour la projection UMAP (install.packages('uwot')).") # Vérifie la disponibilité du package avant tout calcul
      
      donnees <- private$extraire_donnees_gate(nom_gate, canaux, sous_echantillonnage_max, cofacteur) # Extrait, transforme (Arcsinh) et poole la matrice d'expression des cellules contenues dans le gate, toutes échantillons confondus
      mat <- as.matrix(donnees$expression) # Convertit explicitement en matrice numérique, format attendu par uwot::umap
      
      set.seed(self$seed) # Fixe la graine aléatoire pour la reproductibilité (UMAP est un algorithme stochastique)
      embedding <- uwot::umap( # Calcule la projection non linéaire en dimension réduite
        mat,
        n_neighbors  = n_neighbors, # Nombre de voisins définissant le voisinage local de chaque cellule
        min_dist = min_dist, # Distance minimale autorisée entre points dans l'espace de sortie (compacité des amas)
        n_components = n_components, # Dimensionnalité de sortie souhaitée (2 par défaut, pour visualisation)
        metric  = metric, # Métrique de distance utilisée pour le calcul des plus proches voisins
        n_threads = 1, # Un seul thread pour garantir la reproductibilité malgré la graine fixée
        verbose = TRUE # Affiche la progression du calcul dans la console R
      )
      colnames(embedding) <- paste0("UMAP_", seq_len(ncol(embedding))) # Nomme explicitement les colonnes de sortie (UMAP_1, UMAP_2, ...)
      
      if (is.null(self$analyses_umap)) self$analyses_umap <- list() # Sécurité : initialise la structure de stockage si absente
      self$analyses_umap[[nom_gate]] <- list( # Archive le résultat complet sous le nom du gate, en écrasant un éventuel résultat précédent pour ce même gate
        embedding = embedding, # Coordonnées des cellules dans l'espace UMAP réduit
        echantillon_origine = donnees$echantillon_origine, # Vecteur de traçabilité : échantillon d'origine de chaque ligne de l'embedding
        canaux = donnees$canaux, # Liste des canaux/marqueurs effectivement utilisés pour le calcul
        expression = donnees$expression, # Matrice d'expression brute poolée (mêmes lignes que l'embedding), pour permettre de colorer la projection par marqueur
        parametres = list(n_neighbors = n_neighbors, min_dist = min_dist, n_components = n_components, metric = metric) # Mémorise les réglages utilisés pour traçabilité
      )
      
      message("UMAP calculée pour le gate '", nom_gate, "' sur ", nrow(mat), " cellules (", length(donnees$canaux), " canaux).") # Confirme le succès et résume l'ampleur du calcul effectué
      return(invisible(self$analyses_umap[[nom_gate]])) # Renvoie de manière invisible le résultat complet, prêt à être exploité par l'interface Shiny
    },
    
    # ============= projection_tSNE =====================
    projection_tSNE = function(nom_gate, canaux = NULL, dims = 2, perplexity = 30, theta = 0.5,
                               max_iter = 1000, sous_echantillonnage_max = 20000, cofacteur = 150) { # Réalise une projection t-SNE (Rtsne) de la population sélectionnée par le gate
      if (!requireNamespace("Rtsne", quietly = TRUE)) stop("Le package 'Rtsne' est requis pour la projection t-SNE (install.packages('Rtsne')).") # Vérifie la disponibilité du package avant tout calcul
      
      donnees <- private$extraire_donnees_gate(nom_gate, canaux, sous_echantillonnage_max, cofacteur) # Extrait, transforme (Arcsinh) et poole la matrice d'expression des cellules contenues dans le gate, toutes échantillons confondus
      mat <- as.matrix(donnees$expression) # Convertit explicitement en matrice numérique, format attendu par Rtsne::Rtsne
      
      if (nrow(mat) <= 3 * perplexity) { # Rtsne exige un nombre de cellules suffisant par rapport à la perplexité demandée
        stop("Pas assez de cellules dans le gate (", nrow(mat), ") pour une perplexité de ", perplexity, ". Réduisez la perplexité ou vérifiez le gate.")
      }
      
      set.seed(self$seed) # Fixe la graine aléatoire pour la reproductibilité (t-SNE est un algorithme stochastique)
      resultat_tsne <- Rtsne::Rtsne( # Calcule la projection non linéaire en dimension réduite via l'algorithme de Barnes-Hut
        mat,
        dims= dims, # Dimensionnalité de sortie souhaitée (2 par défaut, pour visualisation)
        perplexity = perplexity, # Contrôle l'équilibre entre structure locale et globale (nombre effectif de voisins considérés)
        theta  = theta, # Paramètre d'approximation de Barnes-Hut (0 = exact mais lent, proche de 1 = rapide mais approximatif)
        max_iter = max_iter, # Nombre d'itérations d'optimisation de la descente de gradient
        check_duplicates = FALSE, # Autorise les lignes dupliquées (fréquentes en cytométrie après arrondis/compensation) sans lever d'erreur
        num_threads = 1, # Un seul thread pour garantir la reproductibilité malgré la graine fixée
        verbose = TRUE # Affiche la progression du calcul dans la console R
      )
      
      embedding <- resultat_tsne$Y # Extrait la matrice de coordonnées finales de l'objet retourné par Rtsne
      colnames(embedding) <- paste0("tSNE_", seq_len(ncol(embedding))) # Nomme explicitement les colonnes de sortie (tSNE_1, tSNE_2, ...)
      
      if (is.null(self$analyses_tsne)) self$analyses_tsne <- list() # Sécurité : initialise la structure de stockage si absente
      self$analyses_tsne[[nom_gate]] <- list( # Archive le résultat complet sous le nom du gate, en écrasant un éventuel résultat précédent pour ce même gate
        embedding = embedding, # Coordonnées des cellules dans l'espace t-SNE réduit
        echantillon_origine = donnees$echantillon_origine, # Vecteur de traçabilité : échantillon d'origine de chaque ligne de l'embedding
        canaux  = donnees$canaux, # Liste des canaux/marqueurs effectivement utilisés pour le calcul
        expression= donnees$expression, # Matrice d'expression brute poolée (mêmes lignes que l'embedding), pour permettre de colorer la projection par marqueur
        parametres = list(dims = dims, perplexity = perplexity, theta = theta, max_iter = max_iter) # Mémorise les réglages utilisés pour traçabilité
      )
      
      message("t-SNE calculée pour le gate '", nom_gate, "' sur ", nrow(mat), " cellules (", length(donnees$canaux), " canaux).") # Confirme le succès et résume l'ampleur du calcul effectué
      return(invisible(self$analyses_tsne[[nom_gate]])) # Renvoie de manière invisible le résultat complet, prêt à être exploité par l'interface Shiny
    },
    
    # ============= projection_PCA =====================
    projection_PCA = function(nom_gate, canaux = NULL, n_components = 2, centrer = TRUE, reduire = TRUE,
                              sous_echantillonnage_max = NULL, cofacteur = 150) { # Réalise une analyse en composantes principales (stats::prcomp) de la population sélectionnée par le gate
      donnees <- private$extraire_donnees_gate(nom_gate, canaux, sous_echantillonnage_max, cofacteur) # Extrait, transforme (Arcsinh) et poole la matrice d'expression des cellules contenues dans le gate, toutes échantillons confondus
      mat <- as.matrix(donnees$expression) # Convertit explicitement en matrice numérique, format attendu par prcomp
      
      resultat_pca <- stats::prcomp(mat, center = centrer, scale. = reduire) # Calcule la décomposition en valeurs propres (SVD) de la matrice d'expression, centrée et/ou réduite selon les réglages
      
      n_disponibles <- ncol(resultat_pca$x) # Nombre total de composantes principales effectivement calculables (borné par le nombre de canaux)
      n_a_garder <- min(n_components, n_disponibles) # Sécurité : ne garde pas plus de composantes que ce qui est mathématiquement disponible
      embedding <- resultat_pca$x[, seq_len(n_a_garder), drop = FALSE] # Extrait les coordonnées des cellules sur les n premières composantes principales
      
      variance_expliquee <- (resultat_pca$sdev^2 / sum(resultat_pca$sdev^2))[seq_len(n_a_garder)] # Calcule la proportion de variance totale expliquée par chacune des composantes conservées
      
      if (is.null(self$analyses_pca)) self$analyses_pca <- list() # Sécurité : initialise la structure de stockage si absente
      self$analyses_pca[[nom_gate]] <- list( # Archive le résultat complet sous le nom du gate, en écrasant un éventuel résultat précédent pour ce même gate
        embedding = embedding, # Coordonnées des cellules sur les composantes principales conservées
        echantillon_origine = donnees$echantillon_origine, # Vecteur de traçabilité : échantillon d'origine de chaque ligne de l'embedding
        canaux = donnees$canaux, # Liste des canaux/marqueurs effectivement utilisés pour le calcul
        expression = donnees$expression, # Matrice d'expression brute poolée (mêmes lignes que l'embedding), pour permettre de colorer la projection par marqueur
        variance_expliquee  = variance_expliquee, # Proportion de variance expliquée par chaque composante conservée
        rotation = resultat_pca$rotation[, seq_len(n_a_garder), drop = FALSE], # Matrice des poids (loadings) de chaque canal sur chaque composante, utile pour interpréter les axes
        parametres = list(n_components = n_a_garder, centrer = centrer, reduire = reduire) # Mémorise les réglages utilisés pour traçabilité
      )
      
      message("PCA calculée pour le gate '", nom_gate, "' sur ", nrow(mat), " cellules (",
              round(sum(variance_expliquee) * 100, 1), "% de variance expliquée sur ", n_a_garder, " composante(s)).") # Confirme le succès et résume la qualité de la réduction de dimension obtenue
      return(invisible(self$analyses_pca[[nom_gate]])) # Renvoie de manière invisible le résultat complet, prêt à être exploité par l'interface Shiny
    },
    
    # ============= creer_clusters =====================
    creer_clusters = function(nom_gate, canaux = NULL, xdim = 10, ydim = 10, n_metaclusters = 10,
                              sous_echantillonnage_max = NULL, reutiliser_donnees_de = NULL, cofacteur = 150) { # Réalise un clustering non supervisé (FlowSOM) de la population sélectionnée par le gate
      if (!requireNamespace("FlowSOM", quietly = TRUE)) stop("Le package 'FlowSOM' est requis pour le clustering (BiocManager::install('FlowSOM')).") # Vérifie la disponibilité du package avant tout calcul
      
      donnees <- if (!is.null(reutiliser_donnees_de)) {
        source_analyse <- if (identical(reutiliser_donnees_de, "umap")) self$analyses_umap[[nom_gate]] else self$analyses_tsne[[nom_gate]]
        if (is.null(source_analyse) || is.null(source_analyse$expression)) {
          stop("Aucune analyse ", reutiliser_donnees_de, " disponible pour le gate '", nom_gate, "' à réutiliser. Calculez-la d'abord, ou laissez reutiliser_donnees_de = NULL pour une extraction indépendante.")
        }
        list(expression = source_analyse$expression, echantillon_origine = source_analyse$echantillon_origine, canaux = source_analyse$canaux)
      } else {
        private$extraire_donnees_gate(nom_gate, canaux, sous_echantillonnage_max, cofacteur) # Extraction, transformation (Arcsinh) et pooling indépendants (comportement historique)
      }
      mat <- as.matrix(donnees$expression) # Convertit explicitement en matrice numérique
      
      ff_pool <- flowCore::flowFrame(mat) # Reconstruit un flowFrame unique à partir de la matrice poolée, format attendu en entrée de FlowSOM::FlowSOM()
      
      set.seed(self$seed) # Fixe la graine aléatoire pour la reproductibilité (l'initialisation de la grille SOM est aléatoire)
      fsom <- FlowSOM::FlowSOM( # Construit la carte auto-organisatrice (SOM) puis réalise la méta-clusterisation en une seule étape
        ff_pool,
        colsToUse = colnames(mat), # Restreint explicitement le clustering aux canaux/marqueurs sélectionnés
        xdim      = xdim, # Largeur de la grille SOM (nombre de nœuds sur l'axe X)
        ydim      = ydim, # Hauteur de la grille SOM (nombre de nœuds sur l'axe Y)
        nClus     = n_metaclusters, # Nombre de métaclusters cibles pour la clusterisation hiérarchique finale des nœuds SOM
        seed      = self$seed # Transmet également la graine en interne à FlowSOM pour une reproductibilité complète de l'algorithme
      )
      
      clusters     <- FlowSOM::GetClusters(fsom) # Récupère l'assignation de chaque cellule à son nœud SOM individuel (cluster fin)
      metaclusters <- FlowSOM::GetMetaclusters(fsom) # Récupère l'assignation de chaque cellule à son métacluster (regroupement de nœuds SOM voisins)
      
      if (is.null(self$clusters_flowsom)) self$clusters_flowsom <- list() # Sécurité : initialise la structure de stockage si absente
      self$clusters_flowsom[[nom_gate]] <- list( # Archive le résultat complet sous le nom du gate, en écrasant un éventuel résultat précédent pour ce même gate
        fsom                = fsom, # Objet FlowSOM complet (utile pour les visualisations natives : arbre SOM, heatmap des marqueurs, etc.)
        clusters            = clusters, # Vecteur d'assignation aux clusters fins (nœuds SOM), un par cellule
        metaclusters        = metaclusters, # Vecteur d'assignation aux métaclusters, un par cellule
        echantillon_origine = donnees$echantillon_origine, # Vecteur de traçabilité : échantillon d'origine de chaque cellule clusterisée
        canaux              = donnees$canaux, # Liste des canaux/marqueurs effectivement utilisés pour le calcul
        aligne_sur          = reutiliser_donnees_de, # Mémorise si ce clustering est aligné sur une projection existante ("umap"/"tsne") ou indépendant (NULL)
        expression          = donnees$expression, # Matrice d'expression brute poolée (mêmes lignes que les clusters), pour permettre les heatmaps de marqueurs par cluster
        parametres          = list(xdim = xdim, ydim = ydim, n_metaclusters = n_metaclusters) # Mémorise les réglages utilisés pour traçabilité
      )
      
      message("Clustering FlowSOM calculé pour le gate '", nom_gate, "' sur ", nrow(mat), " cellules (",
              length(unique(metaclusters)), " métaclusters, grille ", xdim, "x", ydim, ").") # Confirme le succès et résume l'ampleur du clustering obtenu
      return(invisible(self$clusters_flowsom[[nom_gate]])) # Renvoie de manière invisible le résultat complet, prêt à être exploité par l'interface Shiny
    },
    
    # ============= resumer_expression_clusters =====================
    resumer_expression_clusters = function(nom_gate, niveau = "metacluster") { # Résume l’expression médiane par cluster ou metacluster
      
      res <- self$clusters_flowsom[[nom_gate]] # Récupère les résultats FlowSOM associés au gate
      if (is.null(res)) stop("Aucun clustering disponible pour le gate '", nom_gate, "'. Lancez creer_clusters() d'abord.") # Stoppe si aucun clustering n’existe
      
      assignation <- if (identical(niveau, "metacluster")) res$metaclusters else res$clusters # Choisit le niveau d’assignation (clusters ou metaclusters)
      mat <- res$expression # Matrice d’expression normalisée utilisée pour le clustering
      
      ids <- sort(unique(assignation)) # Liste ordonnée des identifiants de clusters/metaclusters
      
      medianes <- do.call(rbind, lapply(ids, function(id) { # Calcule la médiane par canal pour chaque cluster
        sous_mat <- mat[assignation == id, , drop = FALSE] # Sous-matrice des cellules appartenant au cluster
        apply(sous_mat, 2, stats::median, na.rm = TRUE) # Médiane par canal
      }))
      
      effectifs <- vapply(ids, function(id) sum(assignation == id), integer(1)) # Nombre de cellules par cluster
      
      etiquettes <- as.character(ids) # Étiquettes par défaut
      
      if (identical(niveau, "metacluster") && !is.null(res$noms_metaclusters)) { # Ajoute les noms personnalisés si disponibles
        etiquettes <- vapply(as.character(ids), function(id) {
          nom_perso <- res$noms_metaclusters[[id]] # Nom personnalisé du metacluster
          if (!is.null(nom_perso) && nchar(trimws(nom_perso)) > 0) paste0(id, " - ", nom_perso) else id # Combine ID + nom si présent
        }, character(1))
      }
      
      rownames(medianes) <- paste0(etiquettes, " (n=", format(effectifs, big.mark = " "), ")") # Ajoute effectif dans le nom de ligne
      
      medianes # Renvoie la matrice des médianes
    },
    
    # ============= renommer_metacluster =====================
    renommer_metacluster = function(nom_gate, id_metacluster, nouveau_nom) { # Renomme un metacluster FlowSOM pour un gate donné
      
      if (is.null(self$clusters_flowsom) || is.null(self$clusters_flowsom[[nom_gate]])) { # Vérifie que le clustering existe
        stop("Aucun clustering disponible pour le gate '", nom_gate, "'. Lancez creer_clusters() d'abord.") # Stoppe si aucun résultat FlowSOM
      }
      
      if (is.null(self$clusters_flowsom[[nom_gate]]$noms_metaclusters)) self$clusters_flowsom[[nom_gate]]$noms_metaclusters <- list() # Initialise la liste des noms si nécessaire
      
      self$clusters_flowsom[[nom_gate]]$noms_metaclusters[[as.character(id_metacluster)]] <- trimws(nouveau_nom) # Enregistre le nouveau nom (trimé) pour l’ID donné
      
      message("Métacluster ", id_metacluster, " renommé : '", nouveau_nom, "'.") # Message de confirmation
      
      return(invisible(self)) # Retour silencieux de l’objet R6
    },
    
    # ============= definir_groupe_echantillon =====================
    definir_groupe_echantillon = function(nom_echantillon, groupe) { # Assigne un échantillon à un groupe expérimental
      
      if (is.null(self$groupes_echantillons)) self$groupes_echantillons <- list() # Initialise la structure si nécessaire
      groupe <- trimws(groupe) # Nettoie les espaces autour du nom de groupe
      
      if (nchar(groupe) == 0) { # Si le groupe est vide
        self$groupes_echantillons[[nom_echantillon]] <- NULL # Retire l’échantillon de tout groupe
      } else {
        self$groupes_echantillons[[nom_echantillon]] <- groupe # Assigne l’échantillon au groupe spécifié
      }
      
      return(invisible(self)) # Retour silencieux
    },
    
    # ============= comparer_groupes =====================
    comparer_groupes = function(nom_gate, variable = "pourcentage", canal = NULL) { # Compare des groupes expérimentaux sur un gate donné
      
      if (identical(variable, "MFI") && (is.null(canal) || nchar(canal) == 0)) { # Si l’utilisateur demande une MFI, un canal est obligatoire
        stop("Précisez un canal pour comparer la MFI (intensité médiane de fluorescence).")
      }
      
      if (is.null(self$groupes_echantillons) || length(self$groupes_echantillons) == 0) { # Vérifie que des groupes ont été définis
        stop("Aucun groupe défini. Assignez chaque échantillon à un groupe avec definir_groupe_echantillon() avant de comparer.")
      }
      
      resume <- self$resumer_gate(nom_gate) # Récupère le résumé du gate : échantillon / n_evenements / pct_parent
      
      resume$groupe <- vapply(resume$echantillon, function(n) { # Assigne le groupe à chaque échantillon
        g <- self$groupes_echantillons[[n]]
        if (is.null(g)) NA_character_ else g
      }, character(1))
      
      if (identical(variable, "MFI")) { # Cas où l’on compare une intensité médiane
        population_gate <- self$resoudre_population_gate(nom_gate) # Récupère les flowFrames filtrés pour ce gate
        
        resume$valeur <- vapply(resume$echantillon, function(n) { # Calcule la MFI pour chaque échantillon
          ff <- population_gate[[n]]
          if (is.null(ff) || !(canal %in% flowCore::colnames(ff))) return(NA_real_) # Ignore si canal absent
          stats::median(flowCore::exprs(ff)[, canal], na.rm = TRUE) # Médiane du canal
        }, numeric(1))
        
      } else { # Cas où l’on compare le pourcentage parent
        resume$valeur <- resume$pct_parent
      }
      
      resume <- resume[!is.na(resume$groupe) & !is.na(resume$valeur), , drop = FALSE] # Retire les lignes sans groupe ou sans valeur
      
      if (nrow(resume) == 0) { # Vérifie qu’il reste des données exploitables
        stop("Aucun échantillon avec à la fois un groupe assigné et une valeur disponible pour ce gate/canal.")
      }
      
      groupes_presents <- unique(resume$groupe) # Liste des groupes présents
      
      if (length(groupes_presents) < 2) { # Il faut au moins deux groupes pour comparer
        stop("Il faut au moins 2 groupes distincts (avec des échantillons assignés) pour comparer.")
      }
      
      p_value <- NA_real_ # Valeur p par défaut
      methode <- NA_character_ # Nom du test utilisé
      
      if (length(groupes_presents) == 2) { # Comparaison à deux groupes → test de Wilcoxon
        valeurs_a <- resume$valeur[resume$groupe == groupes_presents[1]]
        valeurs_b <- resume$valeur[resume$groupe == groupes_presents[2]]
        
        test <- tryCatch(stats::wilcox.test(valeurs_a, valeurs_b), error = function(e) NULL) # Test non paramétrique
        if (!is.null(test)) { p_value <- test$p.value; methode <- "Wilcoxon-Mann-Whitney" }
        
      } else { # Comparaison multi-groupes → test de Kruskal-Wallis
        test <- tryCatch(stats::kruskal.test(valeur ~ factor(groupe), data = resume), error = function(e) NULL)
        if (!is.null(test)) { p_value <- test$p.value; methode <- "Kruskal-Wallis" }
      }
      
      list( # Retourne un objet structuré contenant les résultats
        donnees  = resume, # Tableau échantillon / n_evenements / pct_parent / groupe / valeur
        variable = variable, # Type de variable comparée (pourcentage ou MFI)
        canal    = canal, # Canal utilisé si MFI
        methode  = methode, # Test statistique utilisé
        p_value  = p_value, # Valeur p
        groupes  = groupes_presents # Groupes inclus dans la comparaison
      )
    },
    
    # ============= exporter_fcs_pretraitement =====================
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
    
    # ============= sauvegarder_session_pretraitement_rds =====================
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
    
    
  ), # fin public list
  
  private = list(
    df_control_file = NULL,
    
    # ============= extraire_donnees_gate =====================
    extraire_donnees_gate = function(nom_gate, canaux = NULL, sous_echantillonnage_max = NULL, cofacteur = 150) {
      if (is.null(self$gates_personnalisees) || is.null(self$gates_personnalisees[[nom_gate]])) { # Vérifie que le gate demandé a bien été créé au préalable
        stop("Le gate '", nom_gate, "' n'existe pas. Créez-le d'abord avec creer_gate().")
      }
      population_gate <- self$resoudre_population_gate(nom_gate)
      if (length(population_gate) == 0) {
        stop("Aucune cellule trouvée dans le gate '", nom_gate, "' pour les échantillons disponibles.")
      }
      
      liste_expr    <- list() # Accumule, pour chaque échantillon, la sous-matrice d'expression des cellules retenues par le gate
      liste_origine <- list() # Accumule, pour chaque échantillon, un vecteur répétant son nom (traçabilité de l'origine de chaque cellule poolée)
      
      for (nom in names(population_gate)) { # Parcourt chaque échantillon où au moins une cellule est retenue par ce gate
        ff_gate <- population_gate[[nom]]
        canaux_disponibles <- flowCore::colnames(ff_gate) # Liste les canaux réellement présents dans ce fichier FCS précis
        
        canaux_a_garder <- if (is.null(canaux)) { # Si l'utilisateur n'a précisé aucun canal explicitement
          canaux_disponibles[!grepl("FSC|SSC|Time", canaux_disponibles, ignore.case = TRUE)] # Par défaut : tous les canaux de fluorescence (exclut les paramètres morphologiques et le temps)
        } else {
          intersect(canaux, canaux_disponibles) # Sinon, ne garde que l'intersection entre les canaux demandés et ceux réellement présents dans ce fichier
        }
        if (length(canaux_a_garder) == 0) next
        trans_analyse <- flowCore::transformList(canaux_a_garder, flowCore::arcsinhTransform(a = 0, b = 1 / cofacteur, c = 0))
        ff_transforme <- flowCore::transform(ff_gate, trans_analyse)
        
        mat <- flowCore::exprs(ff_transforme)[, canaux_a_garder, drop = FALSE] # Extrait la sous-matrice numérique d'expression TRANSFORMÉE pour les canaux retenus
        liste_expr[[nom]]    <- mat # Mémorise cette sous-matrice pour l'échantillon courant
        liste_origine[[nom]] <- rep(nom, nrow(mat)) # Mémorise, pour chaque cellule de cet échantillon, son nom d'origine
      }
      
      if (length(liste_expr) == 0) {
        stop("Aucune cellule trouvée dans le gate '", nom_gate, "' pour les échantillons disponibles.")
      }
      
      canaux_communs <- Reduce(intersect, lapply(liste_expr, colnames)) # Détermine l'ensemble des canaux présents dans TOUS les échantillons retenus, pour un pooling cohérent
      if (length(canaux_communs) == 0) {
        stop("Aucun canal commun entre les échantillons pour ce gate.")
      }
      liste_expr <- lapply(liste_expr, function(m) m[, canaux_communs, drop = FALSE]) # Restreint chaque sous-matrice aux seuls canaux communs, pour garantir des colonnes identiques avant l'empilement
      
      expression_totale <- do.call(rbind, liste_expr) # Empile verticalement les sous-matrices de tous les échantillons en une seule matrice poolée
      origine_totale    <- unlist(liste_origine, use.names = FALSE) # Empile de la même façon les vecteurs de traçabilité, dans le même ordre que les lignes de expression_totale
      
      if (!is.null(sous_echantillonnage_max) && nrow(expression_totale) > sous_echantillonnage_max) { # Si un plafond de sous-échantillonnage est fixé et dépassé (utile pour les algorithmes coûteux comme UMAP/t-SNE)
        set.seed(self$seed) # Fixe la graine pour la reproductibilité du tirage aléatoire
        indices_tires     <- sample(seq_len(nrow(expression_totale)), sous_echantillonnage_max) # Tire aléatoirement le nombre de cellules autorisé, toutes échantillons confondus
        expression_totale <- expression_totale[indices_tires, , drop = FALSE]
        origine_totale    <- origine_totale[indices_tires]
      }
      
      list(expression = expression_totale, echantillon_origine = origine_totale, canaux = canaux_communs) # Renvoie la matrice poolée, sa traçabilité d'origine et la liste des canaux effectivement utilisés
    }
    
  ) # fin private list
) # fin classe R6