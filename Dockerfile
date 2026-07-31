# ═══════════════════════════════════════════════════════════════════════════
# Dockerfile pour l'application CARROT (analyse de cytométrie en flux)
# ═══════════════════════════════════════════════════════════════════════════

# Image de base : R + Shiny Server déjà installés et configurés (Ubuntu)
FROM rocker/shiny:4.4.1

# ── 1. Dépendances système (bibliothèques C/C++ dont certains packages R ont besoin) ──
# xml2/curl/ssl      : pour flowCore et le téléchargement de packages
# glpk/gsl           : pour des dépendances de calcul de FlowSOM/flowCore
# png/jpeg/tiff/fontconfig/freetype/harfbuzz/fribidi : pour ggplot2/plotly (rendu graphique)
# cmake/gfortran     : pour compiler certains packages depuis les sources (uwot, Rtsne, AutoSpectralRcpp...)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libglpk-dev \
    libgsl-dev \
    zlib1g-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    gfortran \
    git \
    && rm -rf /var/lib/apt/lists/*

# ── 2. Outils d'installation R (BiocManager pour Bioconductor, remotes pour GitHub) ──
RUN R -e "install.packages(c('BiocManager', 'remotes'), repos = 'https://cloud.r-project.org')"

# ── 3. Packages CRAN ──
# Séparé en plusieurs commandes pour profiter du cache Docker : si l'installation
# échoue ou si on ajoute un package plus tard, on ne refait pas tout depuis zéro.
RUN R -e "install.packages(c( \
    'shiny', 'shinydashboard', 'shinyjs', 'shinyFiles', 'DT', 'plotly', 'ggplot2' \
  ), repos = 'https://cloud.r-project.org')"

RUN R -e "install.packages(c( \
    'R6', 'base64enc', 'dplyr', 'tidyr', 'gridExtra', 'png', 'zip' \
  ), repos = 'https://cloud.r-project.org')"

RUN R -e "install.packages(c( \
    'uwot', 'Rtsne' \
  ), repos = 'https://cloud.r-project.org')"

# ── 4. Packages Bioconductor (bioinformatique / cytométrie) ──
# Installés un par un, chacun vérifié explicitement : si l'un d'eux échoue, le
# build s'arrête net avec un message clair, au lieu de continuer silencieusement
# (BiocManager peut échouer sans faire remonter d'erreur R bloquante sinon).
RUN R -e "BiocManager::install('flowCore', update = FALSE, ask = FALSE); \
          if (!requireNamespace('flowCore', quietly = TRUE)) stop('ÉCHEC : flowCore ne s est pas installé.')"

RUN R -e "BiocManager::install('flowWorkspace', update = FALSE, ask = FALSE); \
          if (!requireNamespace('flowWorkspace', quietly = TRUE)) stop('ÉCHEC : flowWorkspace ne s est pas installé.')"

RUN R -e "BiocManager::install('PeacoQC', update = FALSE, ask = FALSE); \
          if (!requireNamespace('PeacoQC', quietly = TRUE)) stop('ÉCHEC : PeacoQC ne s est pas installé.')"

RUN R -e "BiocManager::install('flowAI', update = FALSE, ask = FALSE); \
          if (!requireNamespace('flowAI', quietly = TRUE)) stop('ÉCHEC : flowAI ne s est pas installé.')"

RUN R -e "BiocManager::install('FlowSOM', update = FALSE, ask = FALSE); \
          if (!requireNamespace('FlowSOM', quietly = TRUE)) stop('ÉCHEC : FlowSOM ne s est pas installé.')"

# ── 5. Packages GitHub (AutoSpectral + son accélérateur C++ AutoSpectralRcpp) ──
RUN R -e "remotes::install_github('DrCytometer/AutoSpectral', upgrade = 'never'); \
          if (!requireNamespace('AutoSpectral', quietly = TRUE)) stop('ÉCHEC : AutoSpectral ne s est pas installé.')"

# AutoSpectralRcpp : package compagnon en C++, qui accélère l'unmixing de 10 à
# 100x selon la documentation officielle (sans lui, AutoSpectral bascule
# silencieusement sur une implémentation "pure R" beaucoup plus lente — c'est
# exactement l'avertissement "Package AutoSpectralRcpp not found" observé,
# cause de la lenteur constatée). Package C++ : nécessite un compilateur
# (déjà présent, voir étape 1) ; peut prendre plusieurs minutes à compiler.
RUN R -e "remotes::install_github('DrCytometer/AutoSpectralRcpp', upgrade = 'never'); \
          if (!requireNamespace('AutoSpectralRcpp', quietly = TRUE)) stop('ÉCHEC : AutoSpectralRcpp ne s est pas installé.')"

# ── 6. Copie de l'application ──
# Retire l'application d'exemple fournie par défaut dans l'image rocker/shiny
RUN rm -rf /srv/shiny-server/*

# Copie directement à la racine de Shiny Server (accessible à http://localhost/)
COPY . /srv/shiny-server/

# Shiny Server doit pouvoir lire les fichiers (utilisateur système "shiny")
RUN chown -R shiny:shiny /srv/shiny-server

# ── 7. Port et lancement ──
EXPOSE 3838

# L'image rocker/shiny lance déjà Shiny Server automatiquement au démarrage du
# conteneur (CMD hérité de l'image de base) : pas besoin de le redéfinir ici.