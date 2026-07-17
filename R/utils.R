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