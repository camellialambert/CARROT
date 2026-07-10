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

### définition des gates

obj$definir_tune_gates(gate.name = "bv421", n_cells = 2000, percentile = 70, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "bv650", n_cells = 2000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "bv711", n_cells = 2000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "A588", n_cells = 2000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "SB550", n_cells = 1000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "PerCP", n_cells = 2000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "PE", n_cells = 1000, percentile = 30, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "SN685", n_cells = 1000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "A700", n_cells = 2000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "ECD", n_cells = 2000, percentile = 50, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "APC770", n_cells = 2000, percentile = 60, bandwidth = 0.5)
obj$definir_tune_gates(gate.name = "viable_cells", n_cells = 1000, percentile = 80, bandwidth = 0.5)


# création des gates réels

obj$definir_gates_landmarks(
  control_name = "bv421",
  n.cells = 2000,
  percentile = 70
)

obj$definir_gates_landmarks(
  control_name = "bv650",
  n.cells = 2000,
  percentile = 50
)

obj$definir_gates_landmarks(
  control_name = "bv711",
  n.cells = 2000,
  percentile = 50
)

obj$definir_gates_landmarks(
  control_name = "A588",
  n.cells = 2000,
  percentile = 50
)

obj$definir_gates_density(
  control_name = "SB550",
  n.cells = 2000,
  grid.n = 100,
  bandwidth.factor = 1
)

obj$definir_gates_landmarks(
  control_name = "PerCP",
  n.cells = 2000,
  percentile = 50
)

obj$definir_gates_landmarks(
  control_name = "PE",
  n.cells = 1000,
  percentile = 30
)

obj$definir_gates_landmarks(
  control_name = "SN685",
  n.cells = 1000,
  percentile = 50
)

obj$definir_gates_landmarks(
  control_name = "A700",
  n.cells = 2000,
  percentile = 50
)

obj$definir_gates_landmarks(
  control_name = "ECD",
  n.cells = 2000,
  percentile = 50
)

obj$definir_gates_landmarks(
  control_name = "APC770",
  n.cells = 2000,
  percentile = 60
)

obj$definir_gates_density(
  control_name = "viable_cells",
  n.cells = 2000,
  grid.n = 100,
  bandwidth.factor = 1
)

# la suite

obj$charger_et_nettoyer()
obj$extraire_fluorophore_spectre()

chemin_unstained <- df_monomarques$chemin[df_monomarques$type == "Unstained"]

obj$extraire_spectre_af(
  unstained_fcs_path = chemin_unstained, 
  tissue_name = "Cells"
)


obj$unmix_folder(
  folder_path = dossier_echantillons, 
  tissue_name = "Cells", 
  method = "WLS"
)

obj$unmix_folder(
  folder_path = dossier_echantillons, 
  tissue_name = "Cells", 
  method = "OLS"
)

obj$unmix_folder(
  folder_path = dossier_echantillons, 
  tissue_name = "Cells", 
  method = "WLS"
)

obj$unmix_folder(
  folder_path = dossier_echantillons, 
  tissue_name = "Cells", 
  method = "AutoSpectral"
)

obj$charger_fcs_unmixes(dossier = "AutoSpectral_unmixed")
colnames(obj$echantillons_traites[[1]])


# 1. Vérifier la séparation APC / Alexa Fluor 700 (Deux FAILs proches)
obj$visualiser_unmixing(
  nom_fichier_fcs = "MM WLS.fcs", 
  canal_x = "BV421-A", 
  canal_y = "BV786-A"
)

obj$visualiser_unmixing(
  nom_fichier_fcs = "MM OLS.fcs", 
  canal_x = "BV421-A", 
  canal_y = "BV786-A"
)

obj$visualiser_unmixing(
  nom_fichier_fcs = "MM AutoSpectral.fcs", 
  canal_x = "SC-A", 
  canal_y = "SSC-A"
)

#######################
# correction

# Charger le fichier en mémoire
fcs_robin <- flowCore::read.FCS(
  "~/Desktop/Institut_Cochin/Jeux de données/Données unmixing/MM_Robin.fcs", 
  truncate_max_range = FALSE
)

# L'ajouter à votre objet
obj$echantillons_traites[["MM_Robin.fcs"]] <- fcs_robin
colnames(fcs_robin)

obj$visualiser_unmixing(
  nom_fichier_fcs = "MM_Robin.fcs", 
  canal_x = "FSC-A",  # Mets ici le nom exact trouvé à la première étape
  canal_y = "SSC-A"      # Mets ici le nom exact trouvé à la première étape
)