library(rsconnect)
library(shiny)
library(ggplot2)
library(plotly)
library(clusterProfiler)
library(readxl)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(future)
library(tidyr)
library(shinyjs)
library(org.Hs.eg.db)

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
            actionButton("search_description", "Search Description"),
            textOutput("analysis_status")  # Status text
          ),
          mainPanel(
            plotlyOutput("volcanoPlot")
          )
        )
      })
      
      # Initialize reactive values
      searchValues <- reactiveValues(
        gene_search = "",
        GO_search = "All",
        description_search = "",
        df_inverted = NULL  # Store df_inverted here
      )
      
      # Process uploaded file and perform DESeq2 analysis
      observeEvent(input$analyze_button, {
        req(input$file)  # Ensure a file is uploaded
        
        output$analysis_status <- renderText("Analyzing... Please wait.")  # Inform user
        
        # Read in gene counts data
        genecounts <- tryCatch({
          read_excel(input$file$datapath, sheet = 1, col_names = TRUE)
        }, error = function(e) {
          showModal(modalDialog(title = "Error", "Could not read the Excel file.", easyClose = TRUE))
          output$analysis_status <- renderText("")  # Clear status
          return(NULL)
        })
        
        if (is.null(genecounts)) {
          output$analysis_status <- renderText("")  # Clear status
          return(NULL)
        }  # Stop further processing if reading failed
        
        # Data processing
        genecounts <- as.data.frame(genecounts)
        rownames(genecounts) <- genecounts[, 1]
        genecounts$Gene_Name <- NULL
        
        # Define condition
        condition <- data.frame(genotype = rep(c('C', 'R'), each = 6), row.names = colnames(genecounts))
        
        # Create DESeq2 dataset
        dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~genotype)
        de <- DESeq(dds)
        res <- results(de)
        
        # Generate results and enrichments
        res$pvalue_log10 <- -log10(res$pvalue)
        pvalue_threshold <- 0.05
        fold_change_threshold <- 2
        
        res$significance <- ifelse(res$pvalue < pvalue_threshold, "Significant", "Not Significant")
        res$new_column <- rownames(res)
        res$diffexpressed <- ifelse(res$log2FoldChange > 0, "UP", ifelse(res$log2FoldChange < 0, "DOWN", "NO_CHANGE"))
        
        # Perform GO enrichment analysis
        original_gene_list <- res$log2FoldChange
        names(original_gene_list) <- res$new_column
        gene_list <- na.omit(original_gene_list)
        gene_list = sort(gene_list, decreasing = TRUE)
        
        gse <- gseGO(geneList = gene_list, 
                     ont = "ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = org.Hs.eg.db, 
                     pAdjustMethod = "none")
        
        searchValues$df_inverted <- gse@result %>% separate_rows(core_enrichment, sep = "/")
        
        # Update select input options
        updateSelectInput(session, "GO_search", choices = unique(searchValues$df_inverted$ID))
        updateSelectInput(session, "description_search", choices = unique(searchValues$df_inverted$Description))
        
        # Render volcano plot
        output$volcanoPlot <- renderPlotly({
          data_res <- filteredRes()
          
          # Create volcano plot
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
        
        output$analysis_status <- renderText("Analysis complete!")  # Confirm completion
      })
      
      # Observing search actions
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