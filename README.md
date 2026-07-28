# CARROT
**C**ytometry **A**nalysis and **R**eporting with **R**-based **O**pen-source **T**ools

CARROT est un pipeline complet d'analyse de cytométrie en flux. Il permet de traiter des données de cytométrie conventionnelle et spectrale, les étapes pouvant être réalisées sont :  import, compensation, unmixing spectral, contrôle qualité (PeacoQC, flowAI), prétraitement (retrait des débris, doublets et cellules mortes) et analyses avancées (gating, PCA, UMAP, tSNE).
On y retrouve une combinaison de méthodes automatiques et de méthodes manuelles.
L'outil est accompagné d'une interface Shiny pour une utilisation plus simple.

CARROT a été développé à l'**Institut Cochin** (collaboration entre les plateformes BIOINFORMAT'IC et CYBIO), avec le soutien de l'**Association Française de Cytométrie (AFC)**.

---

## Installation

Installation possible uniquement via Docker pour le moment.

* Installez Docker : https://docs.docker.com/get-docker/
* Puis lancez (sous Windows, utilisez l'invite de commandes) :

```bash
docker pull camellialambert/carrot
docker run -dp 80:3838 --rm camellialambert/carrot
```

* Accédez à l'application via votre navigateur habituel, à l'adresse **http://localhost/**
* N'oubliez pas d'arrêter le conteneur une fois terminé (`docker stop <ID_du_conteneur>`, trouvable via `docker ps`)

```bash
docker run -dp 80:3838 --rm -v /chemin/vers/vos/donnees:/data camellialambert/carrot
```

---

## Fonctionnalités

| Étape | Description |
|---|---|
| **Importation** | Chargement de fichiers FCS (cytomètre conventionnel ou spectral), annotation des marqueurs, gestion des groupes/conditions |
| **Compensation** | Gating manuel sur tubes monomarqués, calcul et édition interactive de la matrice de spillover |
| **Unmixing** | Démixage spectral via [AutoSpectral](https://github.com/DrCytometer/AutoSpectral) |
| **Contrôle qualité** | Détection des anomalies d'acquisition avec [PeacoQC](https://github.com/saeyslab/PeacoQC) et [flowAI](https://bioconductor.org/packages/flowAI/) |
| **Prétraitement** | Retrait des débris, doublets et cellules mortes, gating interactif par échantillon |
| **Analyses avancées** | Gating hiérarchique, projections UMAP / t-SNE / PCA, clustering non supervisé ([FlowSOM](https://github.com/saeyslab/FlowSOM)), comparaisons statistiques entre groupes |

Guide d'utilisation de l'interface dans "guide_utilisateur.pdf"

---

## Technologies

R / Shiny · flowCore · AutoSpectral · PeacoQC · flowAI · FlowSOM · uwot · Rtsne · plotly

---

## Auteure

Camellia Lambert — stage M2 Bioinformatique - Ingénierie de Plateformes (Université Paris Cité), Institut Cochin
Encadrée par Benjamin Saintpierre (plateforme BIOINFORMAT'IC) et Muriel Andrieu (plateforme CYBIO)