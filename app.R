library(shiny)
library(ggplot2)
library(plotly)
library(shinyjs)
library(httr)      # For making API calls
library(jsonlite)  # For processing JSON

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
            textInput("gene_search", "Search for genes (separate multiple genes with ';'):", ""),
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
      
      # Process uploaded file and perform DESeq2 analysis via API
      observeEvent(input$analyze_button, {
        req(input$file)  # Ensure a file is uploaded
        
        output$analysis_status <- renderText("Analyzing... Please wait.")  # Inform user
        
        # Make an API call to the Plumber API
        res <- httr::GET(
          url = "http://localhost:8000/analyze",
          query = list(gene_file = input$file$datapath)
        )
        
        if (res$status_code == 200) {
          # Convert API response to data frame
          data_res <- fromJSON(content(res, "text"))
          
          # Perform GO enrichment analysis using the API results
          original_gene_list <- data_res$log2FoldChange
          names(original_gene_list) <- data_res$new_column
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
          
          # Store inverted results for GO term search
          searchValues$df_inverted <- gse@result %>% separate_rows(core_enrichment, sep = "/")
          
          # Update select input options for GO and pathway searches
          updateSelectInput(session, "GO_search", choices = unique(searchValues$df_inverted$ID))
          updateSelectInput(session, "description_search", choices = unique(searchValues$df_inverted$Description))
          
          # Render volcano plot
          output$volcanoPlot <- renderPlotly({
            p <- ggplot(data_res, aes(x = log2FoldChange, y = pvalue_log10, 
                                      text = paste("Gene:", new_column, "<br>Log2 Fold Change:", log2FoldChange, 
                                                   "<br>P-value:", pvalue, "<br>Significance:", significance, 
                                                   "<br>Differentially Expressed:", diffexpressed, "<br>-log10 Values:", pvalue_log10))) +
              geom_point(aes(color = log2FoldChange, shape = diffexpressed)) +
              geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
              geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "darkblue") +
              xlim(-5, 5) +
              xlab("Log2 Fold Change") + 
              ylab("-log10(P-value)") +
              ggtitle("Volcano Plot")
            
            ggplotly(p, tooltip = "text")
          })
          
          output$analysis_status <- renderText("Analysis complete!")  # Confirm completion
        } else {
          showModal(modalDialog(title = "Error", "API call failed.", easyClose = TRUE))
        }
      })
      
      # Observing gene search for multiple genes
      observeEvent(input$search_gene, {
        gene_search <- strsplit(input$gene_search, ";")[[1]]  # Split gene names by ';'
        searchValues$gene_search <- trimws(gene_search)       # Trim whitespace from each gene
        
        # Highlight selected genes in the volcano plot
        output$volcanoPlot <- renderPlotly({
          data_res <- fromJSON(content(res, "text"))
          
          # Filter the DESeq2 results to highlight searched genes
          filtered_res <- data_res[data_res$new_column %in% searchValues$gene_search, ]
          
          # Create volcano plot with gene highlights
          p <- ggplot(data_res, aes(x = log2FoldChange, y = pvalue_log10, 
                                    text = paste("Gene:", new_column, "<br>Log2 Fold Change:", log2FoldChange, 
                                                 "<br>P-value:", pvalue, "<br>Significance:", significance, 
                                                 "<br>Differentially Expressed:", diffexpressed, "<br>-log10 Values:", pvalue_log10))) +
            geom_point(aes(color = log2FoldChange, shape = diffexpressed)) +
            geom_point(data = filtered_res, color = "red", size = 3) +  # Highlight selected genes in red
            geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
            geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "darkblue") +
            xlim(-5, 5) +
            xlab("Log2 Fold Change") + 
            ylab("-log10(P-value)") +
            ggtitle("Volcano Plot (Gene Search)")
          
          ggplotly(p, tooltip = "text")
        })
      })
      
      # Observing GO term search
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
