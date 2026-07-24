# R/utils_fcs.R

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

canaux_fluo_fcs <- function(fcs) {
  # On retourne simplement tous les noms de colonnes du fichier
  colnames(flowCore::exprs(fcs))
}

# ─────────────────────────────────────────────────────────────────────────
# Recherche, parmi TOUS les mots-clés des métadonnées FCS, ceux qui ressemblent
# à une matrice de compensation (nom contenant "spill" ou "comp", en ignorant
# la casse et le préfixe "$"). Utile en diagnostic quand l'extraction échoue :
# permet de voir exactement ce que contient le fichier.
# ─────────────────────────────────────────────────────────────────────────
lister_cles_spillover_fcs <- function(fcs) {
  if (is.null(fcs)) return(character(0))
  kw <- tryCatch(flowCore::keyword(fcs), error = function(e) NULL)
  if (is.null(kw) || length(kw) == 0) return(character(0))
  noms <- names(kw)
  noms[grepl("spill|comp", noms, ignore.case = TRUE)]
}

# ─────────────────────────────────────────────────────────────────────────
# Extrait la matrice de compensation (spillover) éventuellement déjà présente
# dans les métadonnées d'un fichier FCS. Recherche insensible à la casse et au
# préfixe "$" (couvre "$SPILLOVER", "SPILL", "spillover", "$SPILL", "$COMP",
# "COMP", "SPILLOVER MATRIX", etc.), et tolère les formats que flowCore peut
# renvoyer : matrice déjà parsée (avec ou sans noms de lignes), data.frame, ou
# chaîne de caractères brute au format FCS standard.
# Retourne NULL si aucune matrice n'est trouvée ou n'est exploitable.
# ─────────────────────────────────────────────────────────────────────────
extraire_spillover_depuis_fcs <- function(fcs) {
  if (is.null(fcs)) return(NULL)
  kw <- tryCatch(flowCore::keyword(fcs), error = function(e) NULL)
  if (is.null(kw) || length(kw) == 0) return(NULL)
  
  noms <- names(kw)
  # Priorité aux clés standards, puis on élargit à toute clé contenant "spill" ou "comp"
  cles_prioritaires <- c("SPILL", "$SPILLOVER", "SPILLOVER", "$SPILL", "$COMP", "COMP")
  cles_candidates <- unique(c(
    intersect(cles_prioritaires, noms),
    noms[grepl("spill", noms, ignore.case = TRUE)],
    noms[grepl("^\\$?comp$", noms, ignore.case = TRUE)]
  ))
  if (length(cles_candidates) == 0) return(NULL)
  
  for (cle in cles_candidates) {
    brut <- kw[[cle]]
    if (is.null(brut)) next
    mat <- .parser_spillover_valeur(brut)
    if (!is.null(mat)) return(mat)
  }
  NULL
}

# Sous-fonction : convertit une valeur brute de mot-clé (matrice, data.frame ou
# chaîne de caractères) en matrice de compensation numérique nommée.
.parser_spillover_valeur <- function(brut) {
  # Cas 1 : déjà une matrice numérique (flowCore parse parfois directement le mot-clé)
  if (is.matrix(brut)) {
    if (!is.numeric(brut)) storage.mode(brut) <- "numeric"
    rn <- rownames(brut); cn <- colnames(brut)
    if (!is.null(rn) && !is.null(cn)) return(brut)
    if (!is.null(cn) && is.null(rn)) { rownames(brut) <- cn; return(brut) }
    if (!is.null(rn) && is.null(cn)) { colnames(brut) <- rn; return(brut) }
    return(NULL) # matrice sans aucun nom de canal exploitable
  }
  
  # Cas 2 : data.frame (arrive avec certaines versions de flowCore)
  if (is.data.frame(brut)) {
    mat <- as.matrix(brut)
    return(.parser_spillover_valeur(mat))
  }
  
  # Cas 3 : chaîne de caractères brute, format FCS standard :
  # "n,chan1,chan2,...,chanN,v11,v12,...,vNN"
  if (is.character(brut) && length(brut) == 1) {
    champs <- strsplit(trimws(brut), ",")[[1]]
    n <- suppressWarnings(as.integer(champs[1]))
    if (is.na(n) || n <= 0 || length(champs) < (1 + n + n * n)) return(NULL)
    canaux  <- trimws(champs[2:(1 + n)])
    valeurs <- suppressWarnings(as.numeric(champs[(2 + n):(1 + n + n * n)]))
    if (any(is.na(valeurs))) return(NULL)
    mat <- matrix(valeurs, nrow = n, ncol = n, byrow = TRUE,
                  dimnames = list(canaux, canaux))
    return(mat)
  }
  NULL
}

get_canaux_filtres <- function(p) {
  # Si des tubes monomarqués sont chargés, on se base dessus (comportement historique).
  # Sinon (import d'échantillons déjà compensés/unmixés, sans contrôles), on se rabat
  # directement sur le premier échantillon biologique disponible.
  if (!is.null(p$tubes_monomarques) && length(p$tubes_monomarques) > 0) {
    tous_canaux <- flowCore::colnames(p$tubes_monomarques[[1]])
  } else if (!is.null(p$echantillons) && length(p$echantillons) > 0) {
    tous_canaux <- flowCore::colnames(p$echantillons[[1]])
  } else {
    return(character(0))
  }
  # Le pattern 'fsc|ssc|time' exclut tout ce qui contient ces mots (insensible à la casse)
  return(tous_canaux[!grepl("fsc|ssc|time", tous_canaux, ignore.case = TRUE)])
}


calculer_densite_raster <- function(x, y, xlim, ylim, res = 400, lissage = FALSE) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NULL)
  
  if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5)
  if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
  
  x_breaks <- seq(xlim[1], xlim[2], length.out = res + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = res + 1)
  
  df <- data.frame(X = x, Y = y)
  df_binned <- df |>
    dplyr::mutate(
      x_bin = cut(X, breaks = x_breaks, include.lowest = TRUE),
      y_bin = cut(Y, breaks = y_breaks, include.lowest = TRUE)
    ) |>
    dplyr::count(x_bin, y_bin, name = "densite") |>
    tidyr::drop_na()
  
  if (nrow(df_binned) == 0) return(NULL)
  
  # centres des bins
  x_centers <- (head(x_breaks, -1) + tail(x_breaks, -1)) / 2
  y_centers <- (head(y_breaks, -1) + tail(y_breaks, -1)) / 2
  
  df_binned <- df_binned |>
    dplyr::mutate(
      X = x_centers[as.integer(x_bin)],
      Y = y_centers[as.integer(y_bin)]
    )
  
  # La colonne "densite" contient à ce stade le nombre BRUT d'événements par
  # bin. Avec de gros échantillons (centaines de milliers à millions
  # d'événements), la population se concentre typiquement dans une poignée de
  # bins qui peuvent contenir des milliers d'événements, alors que le reste de
  # la population (biologiquement tout aussi pertinent) n'en compte que
  # quelques dizaines à centaines. Comme scale_fill_gradientn()/le colorscale
  # plotly font un mappage LINÉAIRE entre le min et le max de "densite", ces
  # quelques bins extrêmes écrasent toute la variation en dessous vers la
  # même couleur basse de la palette (dominante darkblue), même si
  # l'échantillon contient énormément de cellules. On applique donc une
  # compression log1p (standard en cytométrie pour les pseudo-couleurs de
  # densité) : elle conserve l'ordre des bins (les plus denses restent les
  # plus "chauds") mais réduit fortement l'écart entre les valeurs extrêmes et
  # modérées, ce qui répartit les couleurs sur l'ensemble de la population au
  # lieu de les concentrer sur une poignée de pixels. La colonne "densite"
  # n'est utilisée nulle part ailleurs que pour cette coloration (aucun calcul
  # de pourcentage ou de comptage n'en dépend) : ce changement est purement
  # visuel et sans risque, et se propage automatiquement à toutes les figures
  # de l'application (débris, doublets, viabilité, PeacoQC, flowAI, unmixing,
  # compensation) ainsi qu'au gating interactif, puisqu'il est appliqué ici
  # une seule fois dans ce helper partagé.
  df_binned$densite <- log1p(df_binned$densite)
  
  df_binned[, c("X", "Y", "densite")]
}

PALETTE_DENSITE <- c("darkblue", "blue", "cyan", "greenyellow", "yellow", "darkorange", "red")

# Positions (0-1) de chaque couleur de PALETTE_DENSITE le long du dégradé.
# Par défaut (espacement uniforme), chaque couleur occuperait 1/6e de
# l'échelle. Ici, les teintes froides (darkblue/blue/cyan) sont resserrées sur
# le début du dégradé et les teintes chaudes (orange/rouge) sont étirées sur
# une plus grande portion de la fin, afin que le rouge apparaisse dès qu'une
# zone est nettement plus dense que la moyenne, plutôt que réservé aux tout
# derniers pixels les plus extrêmes.
PALETTE_DENSITE_STOPS <- c(0, 0.07, 0.16, 0.28, 0.42, 0.55, 1)

COLORSCALE_DENSITE_PLOTLY <- local({
  lapply(seq_along(PALETTE_DENSITE), function(i) list(PALETTE_DENSITE_STOPS[i], PALETTE_DENSITE[i]))
})


generer_image_densite_base64 <- function(x, y, xlim, ylim, res = 800, largeur_px = 600, hauteur_px = 600) {
  df <- calculer_densite_raster(x, y, xlim, ylim, res = res, lissage = FALSE)
  if (is.null(df)) return(NULL)
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x = X, y = Y, fill = densite)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = PALETTE_DENSITE, values = PALETTE_DENSITE_STOPS) +
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