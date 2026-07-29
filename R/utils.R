# ══════════════════════════════════════════════════════════════════════════════
# utils.R — Fonctions utilitaires partagées entre tous les modules Shiny de CARROT.
# ══════════════════════════════════════════════════════════════════════════════

# ============= get_labels_from_fcs =====================
get_labels_from_fcs <- function(fcs) { # Construit les labels "canal | marqueur" à partir du FCS
  if (is.null(fcs)) return(list()) # Aucun FCS → aucune étiquette
  pd <- tryCatch(flowCore::pData(flowCore::parameters(fcs)), error = function(e) NULL) # Récupère PnN/PnS
  if (is.null(pd)) return(list()) # Si pData absent → rien à décrire
  setNames(lapply(seq_len(nrow(pd)), function(i) { # Pour chaque canal
    d <- pd$desc[i] # Description (marqueur)
    if (!is.null(d) && !is.na(d) && nchar(trimws(d)) > 0)
      paste0(pd$name[i], " | ", d) # Format "canal | marqueur"
    else pd$name[i] # Sinon nom technique seul
  }), pd$name)
}

# ============= canaux_fluo_fcs =====================
canaux_fluo_fcs <- function(fcs) { # Renvoie tous les canaux du flowFrame
  colnames(flowCore::exprs(fcs)) # Pas de filtrage FSC/SSC/Time
}

# ============= lister_cles_spillover_fcs =====================
lister_cles_spillover_fcs <- function(fcs) { # Liste les clés de spillover dans les métadonnées FCS
  if (is.null(fcs)) return(character(0))
  kw <- tryCatch(flowCore::keyword(fcs), error = function(e) NULL) # Récupère toutes les métadonnées
  if (is.null(kw) || length(kw) == 0) return(character(0))
  noms <- names(kw)
  noms[grepl("spill|comp", noms, ignore.case = TRUE)] # Filtre les clés pertinentes
}

# ============= extraire_spillover_depuis_fcs =====================
extraire_spillover_depuis_fcs <- function(fcs) { # Extrait la matrice de compensation depuis les métadonnées
  if (is.null(fcs)) return(NULL)
  kw <- tryCatch(flowCore::keyword(fcs), error = function(e) NULL)
  if (is.null(kw) || length(kw) == 0) return(NULL)
  
  noms <- names(kw)
  cles_prioritaires <- c("SPILL", "$SPILLOVER", "SPILLOVER", "$SPILL", "$COMP", "COMP") # Clés standard
  cles_candidates <- unique(c(
    intersect(cles_prioritaires, noms), # Prioritaires
    noms[grepl("spill", noms, ignore.case = TRUE)], # Repli large
    noms[grepl("^\\$?comp$", noms, ignore.case = TRUE)] # Repli strict sur "comp"
  ))
  if (length(cles_candidates) == 0) return(NULL)
  
  for (cle in cles_candidates) { # Teste chaque clé
    brut <- kw[[cle]]
    if (is.null(brut)) next
    mat <- .parser_spillover_valeur(brut) # Parse la valeur brute
    if (!is.null(mat)) return(mat) # Retourne dès qu’une matrice valide est trouvée
  }
  NULL
}

# ============= .parser_spillover_valeur =====================
.parser_spillover_valeur <- function(brut) { # Convertit une valeur brute en matrice de spillover
  if (is.matrix(brut)) { # Cas matrice déjà parsée
    if (!is.numeric(brut)) storage.mode(brut) <- "numeric" # Force numérique
    rn <- rownames(brut); cn <- colnames(brut)
    if (!is.null(rn) && !is.null(cn)) return(brut) # Lignes/colonnes nommées
    if (!is.null(cn) && is.null(rn)) { rownames(brut) <- cn; return(brut) }
    if (!is.null(rn) && is.null(cn)) { colnames(brut) <- rn; return(brut) }
    return(NULL) # Impossible de mapper les canaux
  }
  
  if (is.data.frame(brut)) { # Cas data.frame
    mat <- as.matrix(brut)
    return(.parser_spillover_valeur(mat)) # Reparsage via matrice
  }
  
  if (is.character(brut) && length(brut) == 1) { # Cas chaîne brute FCS
    champs <- strsplit(trimws(brut), ",")[[1]]
    n <- suppressWarnings(as.integer(champs[1])) # Nombre de canaux
    if (is.na(n) || n <= 0 || length(champs) < (1 + n + n * n)) return(NULL)
    canaux  <- trimws(champs[2:(1 + n)]) # Noms des canaux
    valeurs <- suppressWarnings(as.numeric(champs[(2 + n):(1 + n + n * n)])) # n×n valeurs
    if (any(is.na(valeurs))) return(NULL)
    mat <- matrix(valeurs, nrow = n, ncol = n, byrow = TRUE,
                  dimnames = list(canaux, canaux)) # Matrice carrée nommée
    return(mat)
  }
  NULL # Format non reconnu
}

# ============= get_canaux_filtres =====================
get_canaux_filtres <- function(p) { # Renvoie les canaux utiles (hors FSC/SSC/Time)
  if (!is.null(p$tubes_monomarques) && length(p$tubes_monomarques) > 0) {
    tous_canaux <- flowCore::colnames(p$tubes_monomarques[[1]]) # Priorité aux contrôles
  } else if (!is.null(p$echantillons) && length(p$echantillons) > 0) {
    tous_canaux <- flowCore::colnames(p$echantillons[[1]]) # Repli sur les échantillons
  } else {
    return(character(0))
  }
  tous_canaux[!grepl("fsc|ssc|time", tous_canaux, ignore.case = TRUE)] # Filtre morpho/time
}

# ============= calculer_densite_raster =====================
calculer_densite_raster <- function(x, y, xlim, ylim, res = 400, lissage = FALSE) { # Densité 2D par binning
  ok <- is.finite(x) & is.finite(y) # Retire NA/Inf
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NULL) # Pas assez de points
  
  if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5) # Évite largeur nulle
  if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
  
  x_breaks <- seq(xlim[1], xlim[2], length.out = res + 1) # Bornes X
  y_breaks <- seq(ylim[1], ylim[2], length.out = res + 1) # Bornes Y
  
  df <- data.frame(X = x, Y = y)
  df_binned <- df |>
    dplyr::mutate(
      x_bin = cut(X, breaks = x_breaks, include.lowest = TRUE),
      y_bin = cut(Y, breaks = y_breaks, include.lowest = TRUE)
    ) |>
    dplyr::count(x_bin, y_bin, name = "densite") |>
    tidyr::drop_na()
  
  if (nrow(df_binned) == 0) return(NULL)
  
  x_centers <- (head(x_breaks, -1) + tail(x_breaks, -1)) / 2 # Centres X
  y_centers <- (head(y_breaks, -1) + tail(y_breaks, -1)) / 2 # Centres Y
  
  df_binned <- df_binned |>
    dplyr::mutate(
      X = x_centers[as.integer(x_bin)],
      Y = y_centers[as.integer(y_bin)]
    )
  
  df_binned$densite <- log1p(df_binned$densite) # Compression log1p
  
  df_binned[, c("X", "Y", "densite")]
}

# ============= palettes =====================
PALETTE_DENSITE <- c("darkblue", "blue", "cyan", "greenyellow", "yellow", "darkorange", "red") # Palette standard
PALETTE_DENSITE_STOPS <- c(0, 0.07, 0.16, 0.28, 0.42, 0.55, 1) # Positions non uniformes

# ============= COLORSCALE_DENSITE_PLOTLY =====================
COLORSCALE_DENSITE_PLOTLY <- local({ # Palette plotly cohérente avec ggplot2
  lapply(seq_along(PALETTE_DENSITE), function(i) list(PALETTE_DENSITE_STOPS[i], PALETTE_DENSITE[i]))
})

# ============= generer_image_densite_base64 =====================
generer_image_densite_base64 <- function(x, y, xlim, ylim, res = 800, largeur_px = 600, hauteur_px = 600) { # Génère un PNG base64
  df <- calculer_densite_raster(x, y, xlim, ylim, res = res)
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