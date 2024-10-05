# Use the official R image as the base image
FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('shiny', 'ggplot2', 'plotly', 'clusterProfiler', 'readxl', 'tidyverse', 'DESeq2', 'biomaRt', 'future', 'tidyr', 'shinyjs', 'org.Hs.eg.db'), repos='http://cran.rstudio.com/')"

# Copy the Shiny app files into the image
COPY . /srv/shiny-server/

# Expose port 3838 for the Shiny app
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server')"]
