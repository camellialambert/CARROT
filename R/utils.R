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

get_canaux_filtres <- function(p) {
  tous_canaux <- flowCore::colnames(p$tubes_monomarques[[1]])
  # Le pattern 'fsc|ssc|time' exclut tout ce qui contient ces mots (insensible à la casse)
  return(tous_canaux[!grepl("fsc|ssc|time", tous_canaux, ignore.case = TRUE)])
}

