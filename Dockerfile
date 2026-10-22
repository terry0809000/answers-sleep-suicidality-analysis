FROM rocker/r-ver:4.5.1
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
    libglpk-dev libicu-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY renv.lock /app/renv.lock
COPY scripts/setup.R /app/scripts/setup.R
RUN Rscript --vanilla scripts/setup.R
COPY R /app/R
COPY scripts/run_analysis.R /app/scripts/run_analysis.R
COPY tests /app/tests
RUN Rscript --vanilla tests/run_tests.R
ENTRYPOINT ["Rscript", "--vanilla", "/app/scripts/run_analysis.R"]
CMD ["--help"]
