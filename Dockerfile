# Start with the official R Shiny image
FROM rocker/shiny:4.0.3

# Install system dependencies (if needed)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install R packages (add more as needed)
RUN R -e "install.packages(c('shiny', 'plumber', 'plotly', 'ggplot2', 'readxl', 'DESeq2', 'org.Hs.eg.db', 'clusterProfiler', 'shinyjs'), repos='https://cran.rstudio.com/')"

# Expose the port for Shiny and Plumber API
EXPOSE 8080 8000

# Set environment variable for Shiny to run on port 8080
ENV PORT=8080

# Copy the Plumber and Shiny files
COPY ./app.R /srv/shiny-server/
COPY ./plumber.R /srv/plumber/

# Run the Shiny app and Plumber API on separate ports
CMD R -e "shiny::runApp('/srv/shiny-server/app.R', port=8080, host='0.0.0.0')" & R -e "plumber::plumb('/srv/plumber/plumber.R')$run(port=8000, host='0.0.0.0')"

