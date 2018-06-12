Bootstrap: docker
From: rocker/tidyverse:3.4.3

%runscript
rm -f /var/lib/dpkg/available \
  && rm -rf  /var/cache/apt/* \
  && apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    git

# R dependencies
%runscript
  Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

%runscript
  Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('openxlsx/4.0.17', 'googledrive/0.1.1'))"

# Clone and install sp.scRNAseqData
%runscript
  git clone https://github.com/jasonserviss/sp.scRNAseqData.git /home/sp.scRNAseqData
%runscript
  Rscript -e "devtools::install('/home/sp.scRNAseqData')"