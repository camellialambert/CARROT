#test 

source("~/Desktop/Institut_Cochin/Code/CARROT/pipeline_cytometrie.R")

mes_monomarques <- data.frame(
  fichier = c("Compensation_Controls_Alexa.fcs", "Compensation_Controls_FITC.fcs", 
              "Compensation_Controls_PE.fcs", "Compensation_Controls_PerCP.fcs"),
  chemin = c("~/Desktop/Institut_Cochin/Jeux de données/compensation/Compensation_Controls_Alexa.fcs",
             "~/Desktop/Institut_Cochin/Jeux de données/compensation/Compensation_Controls_FITC.fcs",
             "~/Desktop/Institut_Cochin/Jeux de données/compensation/Compensation_Controls_PE.fcs",
             "~/Desktop/Institut_Cochin/Jeux de données/compensation/Compensation_Controls_PerCP.fcs"
             ),
  type = c("Monomarque", "Monomarque", "Monomarque", "Monomarque"),
  canal = c("Alexa Fluor 700-A", "FITC-A", "PE-A", "PerCP-Cy5-5-A"), # NA pour le non-marqué
  stringsAsFactors = FALSE
)

mes_echantillons <- data.frame(
  tube_name = c("MM"),
  chemin = c("~/Desktop/Institut_Cochin/Jeux de données/compensation/Specimen_001_MM.fcs"),
  stringsAsFactors = FALSE
)

# 1. Initialisation
pipeline <- CARROT$new(
  df_monomarques = mes_monomarques, 
  df_echantillons = mes_echantillons, 
  mode = "Conventionnel"
)

# 2. Chargement des fichiers
pipeline$charger_fcs()

# 3. Transformation (indispensable avant gating)
# Choisissez un cofacteur adapté (ex: 150-200 pour du Arcsinh classique)
pipeline$transformer_fcs(cofacteur = 150)

print(pipeline$monomarques_trans[["FITC-A"]])

# 4. Définition des gates (test sur un canal)
# Exemple : canal "FITC-A"
pipeline$definir_et_extraire(
  nom_canal = "FITC-A", 
  intervalle_gate_negatif = c(0, 1.5), 
  intervalle_gate_positif = c(3.5, 6.0)
)

pipeline$definir_et_extraire(
  nom_canal = "PE-A", 
  intervalle_gate_negatif = c(0, 1.5), 
  intervalle_gate_positif = c(3.5, 6.0)
)

pipeline$definir_et_extraire(
  nom_canal = "PerCP-Cy5-5-A", 
  intervalle_gate_negatif = c(0, 1.5), 
  intervalle_gate_positif = c(3.5, 6.0)
)

pipeline$definir_et_extraire(
  nom_canal = "Alexa Fluor 700-A", 
  intervalle_gate_negatif = c(0, 1.5), 
  intervalle_gate_positif = c(3.5, 6.0)
)
# 5. Visualisation des gates (pour valider le positionnement)
# Si vous êtes dans RStudio, cette ligne ouvrira le tracé dans l'onglet Plots
pipeline$graphiques_gates(nom_canal = "FITC-A")

# 6. Calcul de la matrice de spillover
# Cette étape nécessite que vous ayez défini les gates pour TOUS vos canaux
pipeline$calculer_spillover()

# 7. Contrôle visuel de la compensation (biplots)
# Permet de voir si le spillover est bien corrigé
pipeline$controler_monomarques(fichier_monomarque = "Alexa Fluor 700-A", canal_x = "FITC-A", canal_y = "PE-A")

# 8. Application globale
pipeline$compenser()
pipeline$sauvegarder_session_rds(nom_fichier = "Compensation_Session_Complete.rds")
pipeline$visualiser_compensation(nom_echantillon = "MM", canal_x = "FITC-A", canal_y = "PE-A", max_points = 10000, affichage = "Both")



