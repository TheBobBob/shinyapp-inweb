library(shiny)
library(ggplot2)
library(plotly)
library(clusterProfiler)
library(readxl)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(future)
library(tidyr)  # Ensure tidyr is loaded for separate_rows
library(shinycssloaders)  # Load the shinycssloaders package
library(shinyjs)

# UI definition
ui <- fluidPage(
  titlePanel("Interactive Volcano Plot with Gene and GO Term Search"),
  useShinyjs(),
  
  # Password protection
  passwordInput("password", "Enter Password:", value = "", placeholder = "Password"),
  actionButton("submit_password", "Submit"),
  
  # Conditional UI for analysis
  uiOutput("main_ui")
)

# Server logic
server <- function(input, output, session) {
  # Password handling
  correct_password <- "my_secret_password"  # Set your password here
  
  observeEvent(input$submit_password, {
    if (input$password == correct_password) {
      showModal(modalDialog(
        title = "Access Granted",
        "Welcome! You can now search for genes and view the volcano plot.",
        easyClose = TRUE,
        footer = NULL
      ))
      
      # Hide the password input and button after validation
      shinyjs::hide("password")
      shinyjs::hide("submit_password")
      
      # Main UI appears after password is correct
      output$main_ui <- renderUI({
        sidebarLayout(
          sidebarPanel(
            fileInput("file", "Choose XLSX File", multiple = FALSE, accept = c(".xlsx", "text/xlsx")),
            actionButton("analyze_button", "Analyze"),
            tags$hr(),
            textInput("gene_search", "Search for a gene or keyword (separate multiple genes with ';'):", ""),
            actionButton("search_gene", "Search Gene"),
            tags$hr(),
            selectInput("GO_search", "Select a GO term:", choices = NULL),
            actionButton("search_GO_term", "Search GO Term"),
            tags$hr(),
            selectInput("description_search", "Search for the name of a pathway:", choices = NULL),
            actionButton("search_description", "Search Description")
          ),
          mainPanel(
            uiOutput("spinner_output")  # Add spinner to the plot output
          )
        )
      })
      
      # Reactive values to store results and search criteria
      searchValues <- reactiveValues(
        gene_search = "",
        GO_search = "All",
        description_search = "",
        df_inverted = NULL  # Store df_inverted here
      )
      
      # Process uploaded file and perform DESeq2 analysis
      observeEvent(input$analyze_button, {
        req(input$file)  # Ensure a file is uploaded
        
        # Start processing and show the spinner
        output$spinner_output <- renderUI({
          withSpinner(plotlyOutput("volcanoPlot"), type = 6, color = "#1E90FF", size = 1)
        })
        
        # Read in gene counts data
        genecounts <- tryCatch({
          read_excel(input$file$datapath, sheet = 1, col_names = TRUE)
        }, error = function(e) {
          showModal(modalDialog(title = "Error", "Could not read the Excel file.", easyClose = TRUE))
          return(NULL)
        })
        
        if (is.null(genecounts)) return(NULL)  # Stop further processing if reading failed
        
        genecounts <- as.data.frame(genecounts)
        rownames(genecounts) <- genecounts[, 1]
        genecounts$Gene_Name <- NULL
        
        num_samples <- ncol(genecounts)
        
        # Check if the number of samples is even
        if (num_samples %% 2 != 0) {
          showModal(modalDialog(
            title = "Error",
            "The number of samples must be even for proper grouping.",
            easyClose = TRUE
          ))
          return(NULL)
        }
        
        # Create the condition data frame
        condition <- data.frame(genotype = rep(c('C', 'R'), each = num_samples / 2), row.names = colnames(genecounts))
        
        # Create DESeq2 dataset
        dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~genotype)
        de <- DESeq(dds)
        res <- results(de)
        
        # Create additional columns for plotting
        res$pvalue_log10 <- -log10(res$pvalue)
        pvalue_threshold <- 0.05
        fold_change_threshold <- 2
        
        res$significance <- ifelse(res$pvalue < pvalue_threshold, "Significant", "Not Significant")
        res$new_column <- rownames(res)
        res$diffexpressed <- ifelse(res$log2FoldChange > 0, "UP", ifelse(res$log2FoldChange < 0, "DOWN", "NO_CHANGE"))
        
        # Generate gene list for GSEA
        organism = "org.Hs.eg.db"
        original_gene_list <- res$log2FoldChange
        names(original_gene_list) <- res$new_column
        gene_list <- na.omit(original_gene_list)
        gene_list = sort(gene_list, decreasing = TRUE)
        
        # Perform GO enrichment analysis
        gse <- gseGO(geneList = gene_list, 
                     ont = "ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = organism, 
                     pAdjustMethod = "none")
        
        # Store inverted results for GO terms in reactive values
        searchValues$df_inverted <- gse@result %>% separate_rows(core_enrichment, sep = "/")
        
        # Update GO term and description choices in UI
        updateSelectInput(session, "GO_search", choices = unique(searchValues$df_inverted$ID))
        updateSelectInput(session, "description_search", choices = unique(searchValues$df_inverted$Description))
        
        # Reactive filtering of results based on user input
        filteredRes <- reactive({
          data <- as.data.frame(res)
          
          # Apply gene search filter
          if (searchValues$gene_search != "") {
            genes <- strsplit(searchValues$gene_search, ";")[[1]]
            genes <- trimws(genes)
            data <- data %>%
              filter(rowSums(sapply(genes, function(gene) grepl(gene, new_column, ignore.case = TRUE))) > 0)
          }
          
          # Apply GO term filter
          if (searchValues$GO_search != "All") {
            selected_genes <- searchValues$df_inverted %>%
              filter(ID == searchValues$GO_search) %>%
              pull(core_enrichment)
            data <- data %>%
              filter(new_column %in% selected_genes)
          }
          
          # Apply description search filter
          if (searchValues$description_search != "") {
            selected_genes <- searchValues$df_inverted %>%
              filter(Description == searchValues$description_search) %>%
              pull(core_enrichment)
            data <- data %>%
              filter(new_column %in% selected_genes)
          }
          
          data
        })
        
        # Render volcano plot based on filtered results
        output$volcanoPlot <- renderPlotly({
          data_res <- filteredRes()
          
          p <- ggplot(data_res, aes(x = log2FoldChange, y = pvalue_log10, 
                                    text = paste("Gene:", new_column, "<br>Log2 Fold Change:", log2FoldChange, 
                                                 "<br>P-value:", pvalue, "<br>Significance:", significance, 
                                                 "<br>Differentially Expressed:", diffexpressed, "<br>-log10 Values:", pvalue_log10))) +
            geom_point(aes(color = log2FoldChange, shape = diffexpressed)) +
            geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dotted", color = "red") +
            geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold), linetype = "dotted", color = "darkblue") +
            xlim(-5, 5) +
            xlab("Log2 Fold Change") + 
            ylab("-log10(P-value)") +
            ggtitle("Volcano Plot")
          
          ggplotly(p, tooltip = "text")
        })
      })
      
      # Update search criteria based on user actions
      observeEvent(input$search_gene, {
        searchValues$gene_search <- input$gene_search
        searchValues$GO_search <- "All"
        searchValues$description_search <- ""
      })
      
      observeEvent(input$search_GO_term, {
        searchValues$GO_search <- input$GO_search
        searchValues$gene_search <- ""
        searchValues$description_search <- ""
      })
      
      observeEvent(input$search_description, {
        searchValues$description_search <- input$description_search
        searchValues$GO_search <- "All"
        searchValues$gene_search <- ""
      })
      
    } else {
      showModal(modalDialog(
        title = "Access Denied",
        "Incorrect password. Please try again.",
        easyClose = TRUE,
        footer = NULL
      ))
    }
  })
}

shinyApp(ui = ui, server = server)
