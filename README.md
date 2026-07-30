# CARROT
**C**ytometry **A**nalysis and **R**eporting with **R**-based **O**pen-source **T**ools

CARROT is a complete flow cytometry analysis pipeline. It supports both conventional and spectral cytometry data, and includes the following steps: import, compensation, spectral unmixing, quality control (PeacoQC, flowAI), preprocessing (removal of debris, doublets, and dead cells), and advanced analyses (gating, PCA, UMAP, t‑SNE).  
It combines automated methods with manual, user‑driven procedures.
The tool comes with a Shiny interface for easier use.

CARROT was developed at Institut Cochin (a collaboration between the BIOINFORMAT’IC and CYBIO platforms), with support from the **Association Française de Cytométrie (AFC)**.

---

## Installation

Installation is currently only available via Docker.

* Install Docker: https://docs.docker.com/get-docker/
* Then, perform on your terminal :

```bash
docker pull --platform linux/amd64 camellialambert/carrot
docker run --platform linux/amd64 -dp 80:3838 --rm camellialambert/carrot
```
**For macOS / Linux (Ubuntu, Debian, etc) Users**
* Download "Lancer_CARROT.sh"
* Go to the path where this file is on your terminal :

```bash
cd /chemin/vers/CARROT
```

* Then run on your terminal : 

```bash
chmod +x Lancer_CARROT.sh
./Lancer_CARROT.sh
```

* On your usual browser access to CARROT at **http://localhost/**
* Once you are finished, close CARROT with :

```bash
docker stop carrot
```

**For Windows Users**

* Download "Lancer_CARROT.bat" and double-click on the file
* This will open a new window and CARROT will automatically open on your usual browser at **http://localhost/**
* Once you are finished, close CARROT by closing the window or with :

```bash
docker stop carrot
```
---

## Features

| Step | Description |
|---|---|
| **Import** | Loading FCS files (conventional or spectral cytometers), marker annotation, management of groups/conditions |
| **Compensation** | Manual gating on single‑stained tubes, computation and interactive editing of the spillover matrix |
| **Unmixing** | Spectral demixing using [AutoSpectral](https://github.com/DrCytometer/AutoSpectral) |
| **Quality control** | Detection of acquisition anomalies with [PeacoQC](https://github.com/saeyslab/PeacoQC) and [flowAI](https://bioconductor.org/packages/flowAI/) |
| **Preprocessing** | Removal of debris, doublets, and dead cells; interactive gating per sample |
| **Analyses avancées** | Hierarchical gating, UMAP / t‑SNE / PCA projections, unsupervised clustering ([FlowSOM](https://github.com/saeyslab/FlowSOM)), statistical comparisons between groups |

A user guide for the interface is available in **guide_utilisateur.pdf**.

---

## Technologies

R / Shiny · flowCore · AutoSpectral · PeacoQC · flowAI · FlowSOM · uwot · Rtsne · plotly

---

## Authors

Camellia Lambert — Master’s internship in Bioinformatics – Platforms Engineering (Université Paris Cité), Institut Cochin

Supervised by Benjamin Saintpierre (BIOINFORMAT’IC platform) and Muriel Andrieu (CYBIO platform)