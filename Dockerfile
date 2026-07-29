FROM rocker/shiny:4.4.1

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

RUN R -e "install.packages(c('BiocManager', 'remotes'), repos = 'https://cloud.r-project.org')"

RUN R -e "install.packages(c( \
    'shiny', 'shinydashboard', 'shinyjs', 'shinyFiles', 'DT', 'plotly', 'ggplot2' \
  ), repos = 'https://cloud.r-project.org')"

RUN R -e "install.packages(c( \
    'R6', 'base64enc', 'dplyr', 'tidyr', 'gridExtra', 'png', 'zip' \
  ), repos = 'https://cloud.r-project.org')"

RUN R -e "install.packages(c( \
    'uwot', 'Rtsne' \
  ), repos = 'https://cloud.r-project.org')"

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

RUN R -e "remotes::install_github('DrCytometer/AutoSpectral', upgrade = 'never'); \
          if (!requireNamespace('AutoSpectral', quietly = TRUE)) stop('ÉCHEC : AutoSpectral ne s est pas installé.')"

RUN rm -rf /srv/shiny-server/*

COPY . /srv/shiny-server/

RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 3838