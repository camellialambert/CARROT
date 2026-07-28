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
    # ── Première création : réglages globaux fixés, forme appliquée à tous les échantillons disponibles ──
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
    # ── Gate déjà existant : ne touche pas aux réglages globaux, met à jour uniquement la ou les forme(s) ciblée(s) ──
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
}

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
}

# ============= obtenir_forme_gate =====================
obtenir_forme_gate = function(nom_gate, nom_echantillon) {
  if (is.null(self$gates_personnalisees) || is.null(self$gates_personnalisees[[nom_gate]])) return(NULL)
  self$gates_personnalisees[[nom_gate]]$formes[[nom_echantillon]]
}

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
}

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
}

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
}