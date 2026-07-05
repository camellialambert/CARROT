# 1. Chargement du script contenant la classe R6 CARROT
source("~/Desktop/Institut_Cochin/Code/CARROT/pipeline_cytometrie.R")

# 2. Définition des chemins vers vos fichiers FCS déjà compensés
fichier_CARROT <- "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/FCS_Compenses_20260705_102748/ONCOBIO_11_Grandpheno_AllAbs_non_fixe_compense.fcs"
fichier_OMIQ <- "~/Desktop/ONCOBIO_11 Grandpheno_AllAbs_non_fixe.fcs"

mon_pipeline <- CARROT$new()
mon_pipeline$mode <- "Conventionnel" # Ajuster en "Spectral" si nécessaire
nom_tube <- "ONCOBIO_11"

mon_pipeline$echantillons_traites <- list()
mon_pipeline$echantillons_traites[[nom_tube]] <- flowCore::read.FCS(
  fichier_CARROT, 
  transformation = FALSE, 
  truncate_max_range = FALSE
)

mon_pipeline$canaux <- colnames(mon_pipeline$echantillons_traites[["ONCOBIO_11"]])
mon_pipeline$canaux <- setdiff(mon_pipeline$canaux, c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W", "Time", "time"))
mon_pipeline$appliquer_peacoqc()
mon_pipeline$appliquer_flowai(dossier_rapports = NULL, reglages_specifiques = list())

# comparaison 
post_flowAI_CARROT <- mon_pipeline$post_flowAI[["ONCOBIO_11"]]
post_flowAI_OMIQ <- flowCore::read.FCS(fichier_OMIQ,  transformation = FALSE, truncate_max_range = FALSE)


print(paste("Nombre de cellules (OMIQ) :", nrow(post_flowAI_OMIQ)))
print(paste("Nombre de cellules (CARROT) :", nrow(post_flowAI_CARROT)))