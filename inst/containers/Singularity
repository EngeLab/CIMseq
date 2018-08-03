Bootstrap: docker
From: rocker/tidyverse:3.4.3

%post
rm -f /var/lib/dpkg/available \
  && rm -rf  /var/cache/apt/* \
  && apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libudunits2-dev \
    emacs \
    git

# Install CRAN and Bioconductor packages
Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

##sp.scRNAseq imports
Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R'); install_cran(c('mclust/5.4', 'Rtsne/0.13', 'pso/1.0.3', 'matrixStats/0.53.1', 'ggthemes/3.5.0', 'igraph/1.2.1', 'viridis/0.5.1', 'ggraph/1.0.1', 'tidygraph/1.1.0', 'testthat/2.0.0', 'printr/0.1',  'covr/3.1.0', 'Rcpp/0.12.17', 'RcppArmadillo/0.8.500.0', 'RcppEigen/0.3.3.4.0', 'future.apply/0.2.0'))"

Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('S4Vectors', 'BiocStyle'))"

##sp.scRNAseqData imports
Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('openxlsx/4.0.17',  'googledrive/0.1.1'))"

##seqTools imports
Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('e1071/1.6-8', 'doMC/1.3.5'))"

# Clone and install remote R packages
mkdir /home/Github
git clone https://github.com/jasonserviss/sp.scRNAseqData.git /home/Github/sp.scRNAseqData

Rscript -e "devtools::install('/home/Github/sp.scRNAseqData')"
Rscript -e "source('/home/Github/sp.scRNAseqData/inst/rawData/processRaw.R')"
Rscript -e "devtools::install('/home/Github/sp.scRNAseqData')"

git clone https://github.com/jasonserviss/seqTools.git /home/Github/seqTools
Rscript -e "devtools::install('/home/Github/seqTools')"