# =============================================================================
# Dockerfile — reproducible R environment for the HCV PWID model
#
# Optional: CI does not depend on this image. Build and run:
#   docker build -t hepc-model:r .
#   docker run --rm -v "$PWD":/workspace hepc-model:r
#
# Covers the R pipeline (tests, calibration, MCMC, scenarios). Not included:
# cmdstanr/CmdStan (only scripts/model/preprocess.R, which also needs the
# untracked data/sps.csv), ragg (legacy --figures true pies only) and the
# Python figure/manuscript stack (requirements.txt). rocker/r-ver pins CRAN to
# a dated snapshot, so package versions are fixed by the image tag.
# =============================================================================

FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r Rcpp RcppArmadillo jsonlite posterior writexl dplyr lubridate \
    ggplot2 patchwork scales tidyr viridisLite

WORKDIR /workspace
COPY . .

CMD ["Rscript", "scripts/run_tests.R"]
