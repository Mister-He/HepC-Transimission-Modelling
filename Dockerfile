# =============================================================================
# Dockerfile — reproducible R environment for the HCV PWID model
#
# Optional: CI does not depend on this image. Build and run:
#   docker build -t hepc-model:r .
#   docker run --rm -v "$PWD":/workspace hepc-model:r
#
# The Python NPE stack (sbi/torch) is intentionally not included; it is heavy
# and only needed for scripts/run_npe.R. See docker-compose.yml for how to
# run the R pipeline.
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

RUN install2.r Rcpp RcppArmadillo dplyr ggplot2 MASS patchwork

WORKDIR /workspace
COPY . .

CMD ["Rscript", "scripts/run_tests.R"]
