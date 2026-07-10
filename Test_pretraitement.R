source("~/Desktop/Institut_Cochin/Code/CARROT/pipeline_cytometrie.R")

monomarques <- data.frame(chemin = c("~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/Alexa_Fluor_700_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/APC_H7_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/APC_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/BV_421_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/BV_605_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/BV_650_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/BV_711_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/BV_786_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/DAPI_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/FITC_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/PE_Cy7_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/PE_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/PE_Texas_Red_Stained_Control.fcs",
                                     "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Controles/PerCP_Cy5_5_Stained_Control.fcs"),
                          fichier = c("Alexa_Fluor_700_Stained_Control.fcs",
                                     "APC_H7_Stained_Control.fcs",
                                     "APC_Stained_Control.fcs",
                                     "BV_421_Stained_Control.fcs",
                                     "BV_605_Stained_Control.fcs",
                                     "BV_650_Stained_Control.fcs",
                                     "BV_711_Stained_Control.fcs",
                                     "BV_786_Stained_Control.fcs",
                                     "DAPI_Stained_Control.fcs",
                                     "FITC_Stained_Control.fcs",
                                     "PE_Cy7_Stained_Control.fcs",
                                     "PE_Stained_Control.fcs",
                                     "PE_Texas_Red_Stained_Control.fcs",
                                     "PerCP_Cy5_5_Stained_Control.fcs"), 
                          canal = c("Alexa Fluor 700-A",
                                    "APC-H7-A",
                                    "APC-A",
                                    "BV 421-A",
                                    "BV 605-A",
                                    "BV 650-A",
                                    "BV 711-A",
                                    "BV 786-A",
                                    "DAPI-A",
                                    "FITC-A",
                                    "PE-Cy7-A",
                                    "PE-A",
                                    "PE-Texas Red-A",
                                    "PerCP-Cy5-5-A"
                                    ), 
                          type = c("Monomarque", "Monomarque", "Monomarque", "Monomarque",
                                   "Monomarque", "Monomarque", "Monomarque", "Monomarque",
                                   "Monomarque", "Monomarque", "Monomarque", "Monomarque",
                                   "Monomarque", "Monomarque"),
                          stringsAsFactors = FALSE)

echantillons <- data.frame(chemin = "~/Desktop/Institut_Cochin/Jeux de données/Oncobio/Fichier_compense.fcs", 
                           tube_name = "Echantillon",
                           stringsAsFactors = FALSE)

carrot <- CARROT$new(df_monomarques = monomarques, df_echantillons = echantillons)
carrot$charger_fcs()
carrot$echantillons_traites <- carrot$echantillons
#pretraitement
carrot$appliquer_peacoqc()
#carrot$appliquer_flowai()

## nettoyage des débris à partir des données de peacoqc
pts <- matrix(c(
  66000, 5000,   # Point 1
  160000, 10000,  # Point 2
  210000, 20000,# Point 3
  250000, 30000,# Point 4
  250000, 91000, # Point 5
  240000, 141000, # Point 6
  238000, 200000, # Point 7
  220000, 250000, # Point 8
  200000, 250000, # Point 9
  131000, 246000, # Point 10
  110000, 240000, # Point 11
  58000, 231000,  # Point 12
  8000, 197000,   # Point 13
  5000, 137000,   # Point 14
  5000, 120000,   # Point 15
  5000, 70000,   # Point 16
  12000, 60000,   # Point 17
  60000, 56000   # Point 18
), ncol = 2, byrow = TRUE)


carrot$retirer_les_debris(
  matrice_points = pts,
  canal_x = "FSC-A",
  canal_y = "SSC-A",
  nom_echantillon = "Echantillon",
  source_nettoyage = "peacoqc"
)

#graph <- carrot$visualiser_debris(nom_echantillon = "Echantillon", max_points = 1569700)
#print(graph)

## nettoyage des doublets à partir des données de débris
carrot$retirer_doublets_FSC(
  facteur_sensibilite = 4,   
  axe_discrimination = "W_A", 
  nom_echantillon = "Echantillon"
)

#graph_fsc <- carrot$visualiser_doublets("Echantillon", type_analyse = "FSC", max_points = 554000)
#print(graph_fsc)

carrot$retirer_doublets_SSC(
  facteur_sensibilite = 4.5, 
  axe_discrimination = "W_A",  
  nom_echantillon = "Echantillon"
)

#graph_ssc <- carrot$visualiser_doublets("Echantillon", type_analyse = "SSC", max_points = 554000)
#print(graph_ssc)

## nettoyage des cellules mortes à partir des données de doublets
#d'abord, transformation du canal de viabilité:

carrot$transformation_arcsinh(
  canaux = "AmCyan-A", 
  echantillon = "Echantillon", 
  cofactor = 400
)

pts_viab <- matrix(c(
  -0.9, 0,           # Légèrement en dessous de 0 pour inclure le négatif
  2.6, 0,          # Seuil haut de viabilité sur l'axe X
  2.6, 250000,     # En haut à droite
  -0.9, 250000       # En haut à gauche
), ncol = 2, byrow = TRUE)

carrot$retirer_les_cellules_mortes(
  canal_fsc = "AmCyan-A",
  marqueur_viabilite = "FSC-A", 
  points_utilisateur = pts_viab,
  nom_echantillon = "Echantillon"
)

# 3. Visualisation du résultat
graph_viab <- carrot$visualiser_viabilite("Echantillon", max_points = 540000)
print(graph_viab)


