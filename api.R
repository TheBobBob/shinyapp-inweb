library(plumber)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(tidyr)

#* @apiTitle Volcano Plot API

#* Perform DESeq2 analysis and return results
#* @param gene_file The gene expression data file
#* @get /analyze
function(gene_file = NULL) {
  if (is.null(gene_file)) {
    return(list(error = "No file provided."))
  }
  
  # Process and analyze gene file
  genecounts <- tryCatch({
    read_excel(gene_file, sheet = 1, col_names = TRUE)
  }, error = function(e) {
    return(list(error = "Could not read the file."))
  })
  
  if (is.null(genecounts)) {
    return(list(error = "File processing failed."))
  }

  genecounts <- as.data.frame(genecounts)
  rownames(genecounts) <- genecounts[, 1]
  genecounts$Gene_Name <- NULL
  
  condition <- data.frame(genotype = rep(c('C', 'R'), each = 6), row.names = colnames(genecounts))
  
  dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~genotype)
  de <- DESeq(dds)
  res <- results(de)
  
  res$pvalue_log10 <- -log10(res$pvalue)
  
  # Prepare response for plotting
  return(as.data.frame(res))
}

# Create plumber router
pr <- plumb("api.R")
pr$run(port = 8000)
