# ══════════════════════════════════════════════════════════════════════════════
# utils.R — Fonctions utilitaires partagées entre tous les modules Shiny de CARROT.
# Rien ici ne dépend de l'objet pipeline R6 lui-même de façon impérative (sauf
# get_canaux_filtres, qui lit directement quelques champs) : ce sont des
# fonctions "pures" réutilisables, indépendantes de l'état du pipeline.
# ══════════════════════════════════════════════════════════════════════════════

# Construit, pour un flowFrame donné, la liste des libellés d'axes affichables :
# "canal | marqueur" si un marqueur (description $PnS du FCS) a été renseigné à
# l'acquisition, sinon juste le nom technique du canal. Utilisée notamment pour
# pré-remplir les menus déroulants de canaux dans l'import et la compensation
# (complémentaire de p$obtenir_label(), qui fait la même chose mais pour un seul
# canal à la fois, canal par canal, avec le dictionnaire de marqueurs saisi par
# l'utilisateur en plus des métadonnées du FCS).
get_labels_from_fcs <- function(fcs) {
  if (is.null(fcs)) return(list()) # Sécurité : pas de FCS chargé, rien à décrire
  pd <- tryCatch(flowCore::pData(flowCore::parameters(fcs)), error = function(e) NULL) # Table des paramètres du FCS (nom technique $PnN + description $PnS, une ligne par canal)
  if (is.null(pd)) return(list())
  setNames(lapply(seq_len(nrow(pd)), function(i) {
    d <- pd$desc[i] # Description du canal telle qu'enregistrée dans le FCS (souvent le nom du marqueur, ex: "CD3")
    if (!is.null(d) && !is.na(d) && nchar(trimws(d)) > 0)
      paste0(pd$name[i], " | ", d) # Format standard utilisé partout dans l'app : "canal | marqueur"
    else pd$name[i] # Pas de description disponible : on retombe sur le nom technique seul
  }), pd$name)
}

# Renvoie tous les noms de canaux (colonnes) d'un flowFrame, sans filtrage.
# Contrairement à get_canaux_filtres() ci-dessous, garde FSC/SSC/Time : utilisée
# là où on a besoin de la liste brute complète (ex: sélecteurs d'axes pour le
# gating débris/doublets, qui portent justement sur FSC/SSC).
canaux_fluo_fcs <- function(fcs) {
  colnames(flowCore::exprs(fcs))
}

# Recherche, parmi toutes les métadonnées ($TEXT segment) d'un fichier FCS, les
# clés dont le nom ressemble à une matrice de compensation ("spill" ou "comp",
# insensible à la casse). Ne renvoie PAS la matrice elle-même (voir
# extraire_spillover_depuis_fcs() pour ça) : sert uniquement de repli
# diagnostique quand l'extraction échoue, pour dire à l'utilisateur quelles
# clés existent réellement dans son fichier (utile si le format n'est pas
# reconnu — voir p$importer_spillover_fcs() côté pipeline).
lister_cles_spillover_fcs <- function(fcs) {
  if (is.null(fcs)) return(character(0))
  kw <- tryCatch(flowCore::keyword(fcs), error = function(e) NULL) # Table complète des mots-clés de métadonnées du FCS
  if (is.null(kw) || length(kw) == 0) return(character(0))
  noms <- names(kw)
  noms[grepl("spill|comp", noms, ignore.case = TRUE)]
}

# Extrait et parse la matrice de compensation (spillover) déjà appliquée à un
# fichier FCS, telle qu'enregistrée dans ses métadonnées à l'acquisition ou lors
# d'un export préalable. Essaie plusieurs clés candidates dans l'ordre de
# priorité (les plus standards d'abord), et plusieurs formats possibles pour
# chacune (voir .parser_spillover_valeur() ci-dessous) — les logiciels de
# cytométrie n'utilisent pas tous exactement la même convention de nommage ni
# le même type R renvoyé par flowCore pour cette valeur.
extraire_spillover_depuis_fcs <- function(fcs) {
  if (is.null(fcs)) return(NULL)
  kw <- tryCatch(flowCore::keyword(fcs), error = function(e) NULL)
  if (is.null(kw) || length(kw) == 0) return(NULL)
  
  noms <- names(kw)
  cles_prioritaires <- c("SPILL", "$SPILLOVER", "SPILLOVER", "$SPILL", "$COMP", "COMP") # Clés standards les plus fréquentes, testées en premier
  cles_candidates <- unique(c(
    intersect(cles_prioritaires, noms),
    noms[grepl("spill", noms, ignore.case = TRUE)], # Repli large : toute clé contenant "spill"
    noms[grepl("^\\$?comp$", noms, ignore.case = TRUE)] # Repli large : clé nommée exactement "comp" ou "$comp" (évite de matcher des clés sans rapport contenant juste "comp" ailleurs dans le mot)
  ))
  if (length(cles_candidates) == 0) return(NULL)
  
  for (cle in cles_candidates) { # Essaie chaque clé candidate dans l'ordre jusqu'à en trouver une exploitable
    brut <- kw[[cle]]
    if (is.null(brut)) next
    mat <- .parser_spillover_valeur(brut)
    if (!is.null(mat)) return(mat) # Dès qu'une clé donne une matrice valide, on s'arrête là (pas besoin d'en chercher d'autres)
  }
  NULL
}

# Convertit une valeur brute de mot-clé de compensation (peu importe son type R
# d'origine) en une matrice carrée nommée (lignes = colonnes = canaux). Gère les
# 3 formats rencontrés en pratique selon la version de flowCore et le logiciel
# d'acquisition d'origine : matrice déjà parsée, data.frame, ou chaîne de
# caractères brute au format texte standard FCS ("n,canal1,...,canalN,v11,...").
.parser_spillover_valeur <- function(brut) {
  # Cas 1 : déjà une matrice numérique (flowCore parse parfois directement le mot-clé)
  if (is.matrix(brut)) {
    if (!is.numeric(brut)) storage.mode(brut) <- "numeric" # Force le stockage numérique si flowCore a renvoyé du texte/caractère dans une structure matricielle
    rn <- rownames(brut); cn <- colnames(brut)
    if (!is.null(rn) && !is.null(cn)) return(brut) # Cas idéal : lignes et colonnes déjà nommées par canal
    if (!is.null(cn) && is.null(rn)) { rownames(brut) <- cn; return(brut) } # Seules les colonnes sont nommées : la matrice de spillover étant carrée et symétrique en noms, on réutilise les mêmes noms pour les lignes
    if (!is.null(rn) && is.null(cn)) { colnames(brut) <- rn; return(brut) } # Inverse : seules les lignes sont nommées
    return(NULL) # Ni lignes ni colonnes nommées : impossible de savoir quel canal correspond à quelle ligne/colonne, on abandonne plutôt que de risquer un mauvais mapping
  }
  
  # Cas 2 : data.frame (arrive avec certaines versions de flowCore)
  if (is.data.frame(brut)) {
    mat <- as.matrix(brut)
    return(.parser_spillover_valeur(mat)) # Récursion : une fois converti en matrice, retombe dans le cas 1 ci-dessus
  }
  
  # Cas 3 : chaîne de caractères brute, format FCS standard :
  # "n,canal_1,...,canal_n,v_11,v_12,...,v_nn" (n = nombre de canaux, suivi de
  # leurs noms, puis des n×n valeurs de la matrice à plat, ligne par ligne).
  if (is.character(brut) && length(brut) == 1) {
    champs <- strsplit(trimws(brut), ",")[[1]]
    n <- suppressWarnings(as.integer(champs[1])) # Premier champ = nombre de canaux annoncé
    if (is.na(n) || n <= 0 || length(champs) < (1 + n + n * n)) return(NULL) # Sécurité : format incohérent (n absent/invalide, ou pas assez de champs pour n canaux + n×n valeurs)
    canaux  <- trimws(champs[2:(1 + n)]) # Les n champs suivants sont les noms de canaux, dans l'ordre
    valeurs <- suppressWarnings(as.numeric(champs[(2 + n):(1 + n + n * n)])) # Les n×n champs restants sont les valeurs de la matrice, à plat
    if (any(is.na(valeurs))) return(NULL) # Une valeur non numérique dans le lot rend toute la matrice invalide
    mat <- matrix(valeurs, nrow = n, ncol = n, byrow = TRUE,
                  dimnames = list(canaux, canaux)) # Remet les valeurs à plat sous forme de matrice carrée n×n, remplie ligne par ligne comme l'exige le format FCS
    return(mat)
  }
  NULL # Type de valeur non reconnu (ni matrice, ni data.frame, ni chaîne) : abandon
}

# Renvoie la liste des canaux "utiles" pour la compensation/l'analyse : tous les
# canaux du premier tube monomarqué disponible (ou, à défaut, du premier
# échantillon), en excluant les paramètres morphologiques (FSC/SSC) et le canal
# temporel (Time), qui ne sont jamais des cibles de compensation ou de gating
# par marqueur. Utilisée pour peupler les sélecteurs de canaux dans plusieurs
# modules (compensation, transformation...).
get_canaux_filtres <- function(p) {
  if (!is.null(p$tubes_monomarques) && length(p$tubes_monomarques) > 0) {
    tous_canaux <- flowCore::colnames(p$tubes_monomarques[[1]]) # Priorité aux tubes contrôles : leur panel de canaux est la référence pour la compensation
  } else if (!is.null(p$echantillons) && length(p$echantillons) > 0) {
    tous_canaux <- flowCore::colnames(p$echantillons[[1]]) # Repli : pas de contrôles disponibles (cas "déjà compensé"), on utilise le panel des échantillons eux-mêmes
  } else {
    return(character(0)) # Ni contrôles ni échantillons chargés : rien à proposer
  }
  
  return(tous_canaux[!grepl("fsc|ssc|time", tous_canaux, ignore.case = TRUE)])
}

# ══════════════════════════════════════════════════════════════════════════════
# Densité 2D par binning — cœur visuel partagé par TOUTES les figures de
# l'application (biplots de compensation, gating interactif débris/doublets/
# viabilité, PeacoQC, flowAI, unmixing, projections UMAP/t-SNE...). Toute figure
# de densité de l'app passe par cette même fonction : un changement ici change
# le rendu partout à la fois, donc modifier avec prudence.
# ══════════════════════════════════════════════════════════════════════════════

# Calcule une grille de densité 2D par binning régulier (découpage en cases
# égales via cut() + comptage par case via dplyr), plutôt qu'une estimation de
# densité par noyau (type ggpointdensity) : le coût de calcul ne dépend que de
# la résolution demandée (res × res cases), jamais du nombre d'événements — ce
# qui permet d'afficher des millions de cellules aussi vite que quelques
# milliers, sans sous-échantillonnage.
calculer_densite_raster <- function(x, y, xlim, ylim, res = 400, lissage = FALSE) {
  ok <- is.finite(x) & is.finite(y) # Écarte les valeurs non finies (NA/Inf) qui feraient échouer le découpage en cases
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NULL) # Pas assez de points pour former une densité exploitable
  
  if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5) # Sécurité : évite un intervalle de largeur nulle (tous les points ont la même valeur X) qui ferait planter seq()
  if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
  
  x_breaks <- seq(xlim[1], xlim[2], length.out = res + 1) # res+1 bornes délimitent res cases régulières sur l'axe X
  y_breaks <- seq(ylim[1], ylim[2], length.out = res + 1)
  
  df <- data.frame(X = x, Y = y)
  df_binned <- df |>
    dplyr::mutate(
      x_bin = cut(X, breaks = x_breaks, include.lowest = TRUE), # Assigne chaque point à sa case X (include.lowest : les points exactement sur la borne minimale ne sont pas perdus)
      y_bin = cut(Y, breaks = y_breaks, include.lowest = TRUE)
    ) |>
    dplyr::count(x_bin, y_bin, name = "densite") |> # Compte les événements par case (X,Y) : c'est la densité brute
    tidyr::drop_na() # Retire les combinaisons de cases qui n'existent pas réellement dans les données (évite une grille creuse artificiellement complétée)
  
  if (nrow(df_binned) == 0) return(NULL)
  
  # centres des bins
  x_centers <- (head(x_breaks, -1) + tail(x_breaks, -1)) / 2 # Milieu de chaque case (utilisé comme coordonnée représentative pour l'affichage, plutôt que ses bornes)
  y_centers <- (head(y_breaks, -1) + tail(y_breaks, -1)) / 2
  
  df_binned <- df_binned |>
    dplyr::mutate(
      X = x_centers[as.integer(x_bin)], # as.integer(facteur cut()) donne l'indice de la case → on récupère son centre correspondant
      Y = y_centers[as.integer(y_bin)]
    )
  
  # Compression log1p de la densité avant retour : sans cela, quelques cases
  # extrêmement peuplées (le cœur d'une population) écrasent visuellement toute
  # la variation du reste du nuage vers la même couleur basse de la palette,
  # même avec des millions de cellules au total (voir l'historique de réglage
  # de PALETTE_DENSITE_STOPS ci-dessous pour le calibrage fin des couleurs).
  df_binned$densite <- log1p(df_binned$densite)
  
  df_binned[, c("X", "Y", "densite")]
}

# Palette de densité "pseudo-spectrale" standard de l'application (froid = peu
# dense, chaud = très dense), utilisée par TOUTES les figures de densité,
# statiques (ggplot2::scale_fill_gradientn) comme interactives (colorscale
# plotly ci-dessous).
PALETTE_DENSITE <- c("darkblue", "blue", "cyan", "greenyellow", "yellow", "darkorange", "red")

# Positions (0 à 1) de chaque couleur de PALETTE_DENSITE le long du dégradé.
# Volontairement NON uniformes (contrairement à l'espacement par défaut de 1/6e
# chacune) : les teintes froides sont resserrées sur le début du dégradé et les
# teintes chaudes (orange/rouge) étirées sur une plus grande portion de la fin,
# pour que le rouge apparaisse dès qu'une zone est nettement plus dense que la
# moyenne plutôt que réservé aux tout derniers pixels les plus extrêmes.
PALETTE_DENSITE_STOPS <- c(0, 0.07, 0.16, 0.28, 0.42, 0.55, 1)

# Reconstruit PALETTE_DENSITE/PALETTE_DENSITE_STOPS sous la forme attendue par
# l'argument `colorscale` de plotly (liste de paires [position, couleur]), pour
# que les graphiques interactifs (gating débris/doublets/viabilité, module
# Analyses) utilisent exactement le même dégradé que les figures statiques
# ggplot2 — cohérence visuelle garantie par construction, pas par duplication
# manuelle des couleurs à chaque endroit.
COLORSCALE_DENSITE_PLOTLY <- local({
  lapply(seq_along(PALETTE_DENSITE), function(i) list(PALETTE_DENSITE_STOPS[i], PALETTE_DENSITE[i]))
})

# Génère une image PNG de densité (via calculer_densite_raster() + ggplot2) et
# la renvoie encodée en base64, prête à être intégrée directement dans une
# balise <img> HTML sans passer par un fichier statique servi séparément.
# NOTE : méthode actuellement non appelée par aucun module de l'application
# (utilitaire disponible mais inutilisé pour l'instant) — conservée pour un
# usage futur (ex: export HTML autonome, aperçu embarqué dans un e-mail).
generer_image_densite_base64 <- function(x, y, xlim, ylim, res = 800, largeur_px = 600, hauteur_px = 600) {
  df <- calculer_densite_raster(x, y, xlim, ylim, res = res, lissage = FALSE)
  if (is.null(df)) return(NULL)
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x = X, y = Y, fill = densite)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS) +
    ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0)) + # expand=0 : pas de marge automatique, l'image colle exactement aux limites demandées (important pour un rendu "plein cadre" sans bordure blanche)
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    ggplot2::theme_void() + # Aucun axe, titre ni grille : uniquement l'image de densité, pour une intégration <img> propre
    ggplot2::theme(legend.position = "none",
                   plot.margin = ggplot2::margin(0, 0, 0, 0)) # Supprime toute marge résiduelle autour du graphique
  
  fichier_tmp <- tempfile(fileext = ".png")
  on.exit(unlink(fichier_tmp), add = TRUE) # Nettoie le fichier temporaire dès la fin de la fonction, qu'elle réussisse ou échoue
  ggplot2::ggsave(fichier_tmp, plot = g, width = largeur_px / 96, height = hauteur_px / 96,
                  dpi = 96, bg = "transparent") # Conversion pixels -> pouces à 96 DPI (résolution d'écran standard) ; fond transparent pour une intégration propre sur n'importe quel arrière-plan HTML
  
  base64enc::dataURI(file = fichier_tmp, mime = "image/png")
}