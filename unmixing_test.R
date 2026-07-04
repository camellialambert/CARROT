source("~/Desktop/Institut_Cochin/Code/CARROT/unmixing.R")

racine_projet <- "~/Desktop/Institut_Cochin/Jeux de données/Données unmixing"
dossier_controles <- "~/Desktop/Institut_Cochin/Jeux de données/Données unmixing/Reference Group"
dossier_echantillons <- "~/Desktop/Institut_Cochin/Jeux de données/Données unmixing/Echantillons"

fichiers_mono <- c(
  "Unstained (Cells).fcs", "Spark NIR 685 (Cells).fcs",
  "Spark Blue 550 (Cells).fcs", "PerCP-Cy5.5 (Cells).fcs",
  "PE (Cells).fcs", "LIVE DEAD Blue (Cells).fcs",
  "ECD (Cells).fcs", "BV786 (Cells).fcs",
  "BV711 (Cells).fcs", "BV650 (Cells).fcs",
  "BV421 (Cells).fcs", "APC-Vio 770 (Cells).fcs",
  "APC (Cells).fcs", "Alexa Fluor 700 (Cells).fcs",
  "Alexa Fluor 488 (Cells).fcs"
)

df_monomarques <- data.frame(
  chemin = file.path(dossier_controles, fichiers_mono),
  type = c("Unstained", rep("Monomarque", length(fichiers_mono) - 1)),
  canal = c(NA, "Spark NIR 685", "Spark Blue 550", 
            "PerCP-Cy5.5", "PE", "LIVE/DEAD Blue", 
            "PE-Texas Red", "BV786", 
            "BV711", "BV650", "BV421", "APC-Vio770", 
            "APC", "Alexa Fluor 700", "Alexa Fluor 488"),
  stringsAsFactors = FALSE
)

df_echantillons <- data.frame(
  chemin = file.path(dossier_echantillons, "MM.fcs"),
  tube_name = "MM_Sample", 
  stringsAsFactors = FALSE
)


obj <- EchantillonCompense$new(
  df_monomarques = df_monomarques, 
  df_echantillons = df_echantillons, 
  chemin_racine = racine_projet, 
  mode = "Spectral"
)

obj$charger_fcs()
obj$lancer_asp(type_cytometre = "aurora")


# 1. Définir le chemin de ton fichier modifié
chemin_csv <- file.path(obj$dossier_racine, "fcs_control_file.csv")

# 2. Lire le fichier (que ce soit séparé par des virgules ou points-virgules)
df_repare <- read.csv(chemin_csv, sep = ";", stringsAsFactors = FALSE)

# Si jamais il avait quand même été enregistré avec des virgules :
if (ncol(df_repare) == 1) {
  df_repare <- read.csv(chemin_csv, sep = ",", stringsAsFactors = FALSE)
}

# 3. Réenregistrer le fichier avec le séparateur STRICT (la virgule) requis par AutoSpectral
write.csv(df_repare, file = chemin_csv, row.names = FALSE, quote = TRUE)

message("Le fichier fcs_control_file.csv a été nettoyé et réenregistré au format standard.")
obj$verifier_asp(warning = 5000, error = 1000)

cells <- obj$definir_gates_density(control_name = "smallGate_2")
viable_cells <- obj$definir_gates_density(control_name = "viabilityGate_1")
obj$visualiser_figures("figure_gate")
obj$charger_et_nettoyer()
obj$extraire_fluorophore_spectre()
# Extraction de l'AF sur le fichier Unstained pour le tissu "Cells"
chemin_unstained <- df_monomarques$chemin[df_monomarques$type == "Unstained"]

obj$extraire_spectre_af(
  unstained_fcs_path = chemin_unstained, 
  tissue_name = "Cells"
)
obj$preparer_variants_spectraux(tissue_af_name = "Cells")

obj$unmix_folder(
  folder_path = dossier_echantillons, 
  tissue_name = "Cells", 
  method = "WLS"
)

obj$charger_fcs_unmixes(dossier = "AutoSpectral_unmixed")
colnames(obj$echantillons_traites[[1]])
obj$visualiser_unmixing(
  nom_fichier_fcs = "MM AutoSpectral.fcs", 
  canal_x = "BV421-A",  # Mets ici le nom exact trouvé à la première étape
  canal_y = "PE-A"      # Mets ici le nom exact trouvé à la première étape
)

# 1. Vérifier la séparation APC / Alexa Fluor 700 (Deux FAILs proches)
obj$visualiser_unmixing(
  nom_fichier_fcs = "MM AutoSpectral.fcs", 
  canal_x = "APC-A", 
  canal_y = "Alexa Fluor 700-A"
)

# 2. Vérifier le PE (qui était en FAIL) face au BV421 (qui est PASS)
obj$visualiser_unmixing(
  nom_fichier_fcs = "MM AutoSpectral.fcs", 
  canal_x = "BV421-A", 
  canal_y = "PE-A"
)