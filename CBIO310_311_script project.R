# Install libraries
install.packages("BiocManager")
install.packages(c("dplyr", "ggplot2", "fgsea", "msigdbr"))
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE", "pathview"))

#Needed Libraries
library(limma)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(pathview)
library(org.Hs.eg.db)  # For human gene annotations
library(DOSE)

----------------------------------------
  ## COVID ##
----------------------------------------
# Data
expr <- read.table("GSE177477_expression.txt",
                   header = TRUE,
                   sep = "\t",
                   row.names = 1)
#Meta Data
meta <- read.table("GSE177477_metadata.txt",
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE)

#Filtration
common_samples <- intersect(colnames(expr), meta$Sample)
expr <- expr[, common_samples]
meta <- meta[meta$Sample %in% common_samples, ]

#Arrangement
meta <- meta[match(colnames(expr), meta$Sample), ]

#Matrix Design
meta$Condition <- factor(meta$Condition,
                         levels = c("Healthy", "Covid"))
design <- model.matrix(~ Condition, data = meta)
design

#Module Building
fit <- lmFit(expr, design)
fit <- eBayes(fit)


#Result
colnames(fit$coefficients)
results <- topTable(fit,
                    coef = "ConditionCovid",
                    number = Inf,
                    adjust.method = "BH")
head(results)

#DEGs filtration
DEGs_covid <- results[
  abs(results$logFC) >= 1 & results$adj.P.Val < 0.05,
]
#COVID signature genes
up_covid   <- DEGs_covid[DEGs_covid$logFC > 0, ]
down_covid <- DEGs_covid[DEGs_covid$logFC < 0, ]

write.csv(DEGs_covid, "COVID_DEGs.csv")

# Up / Down
up_covid   <- DEGs_covid[DEGs_covid$logFC > 0, ]
down_covid <- DEGs_covid[DEGs_covid$logFC < 0, ]

# Save files
write.csv(DEGs_covid, "COVID_DEGs.csv")
write.csv(up_covid,   "COVID_Upregulated.csv")
write.csv(down_covid, "COVID_Downregulated.csv")

----------------------------------------
  ## Psychiatric Disease ##
----------------------------------------
# Load the expression data
expr <- read.table("GSE53987_series_matrix.txt",
                     header = TRUE,
                     sep = "\t",
                     row.names = 1)

# Load the metadata
meta <- read.table("GSE53987_metadata.txt",
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE)

# Ensure matching between the expression data and the metadata
common_samples <- intersect(colnames(expr), meta[,1])
expr <- expr[, common_samples]
meta <- meta[meta[,1] %in% common_samples, ]
meta <- meta[match(colnames(expr), meta[,1]), ]


# Define the condition variable (disease states)
meta$Condition <- factor(meta$Disease.state, 
    levels = c("control", "bipolar disorder", "major depressive disorder", "schizophrenia"))

# Replace spaces with underscores for valid R names
levels(meta$Condition) <- make.names(levels(meta$Condition))

# Create the design matrix for the linear model
design <- model.matrix(~ 0 + Condition, data = meta)
colnames(design) <- levels(meta$Condition)

# Fit the linear model using limma
fit <- lmFit(expr, design)

# Set up the contrast for schizophrenia vs control
contrast.matrix <- makeContrasts(Schizo_vs_Control = schizophrenia - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract the results from the model
results <- topTable(fit2, coef = "Schizo_vs_Control", number = Inf, adjust.method = "BH")

# Filter DEGs based on logFC and adjusted p-value
DEGs_psychiatric <- results[abs(results$logFC) >= 1 & results$adj.P.Val < 0.05, ]

# Separate upregulated and downregulated genes
up_psychiatric <- DEGs_psychiatric[DEGs_psychiatric$logFC > 0, ]
down_psychiatric <- DEGs_psychiatric[DEGs_psychiatric$logFC < 0, ]

# Save files
write.csv(DEGs_psychiatric, "Psychiatric_DEGs.csv")
write.csv(up_psychiatric, "Psychiatric_Upregulated.csv")
write.csv(down_psychiatric, "Psychiatric_Downregulated.csv")

# Handel Missing Data
missing_genes <- gene_df[is.na(gene_df$ENTREZID), ]
write.csv(missing_genes, "missing_genes.csv") 
#Filtration
entrez_genes <- entrez_genes[!is.na(entrez_genes)]
#GSEA Analysis
ego <- enrichGO(
  gene = entrez_genes, 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.1
)
# Checking
if (is.null(ego)) {
  message("No significant pathways found.")
} else {
# Plots
  dotplot(ego, showCategory = 10)
}

----------------------------------------
## Volcano plots ##
----------------------------------------
# Function to create the volcano plot for control vs infected comparison
create_volcano_plot <- function(degs, title) {
    # Create the logP (negative log10 of adjusted p-value)
    degs$logP <- -log10(degs$adj.P.Val)
    
    # Define the categories for points
    degs$Category <- ifelse(degs$logFC > 0 & degs$adj.P.Val < 0.05, "Upregulated", 
                            ifelse(degs$logFC < 0 & degs$adj.P.Val < 0.05, "Downregulated", "Non-significant"))
    
    # Create the volcano plot using ggplot
    ggplot(degs, aes(x = logFC, y = logP)) +
      # Color points based on their significance
      geom_point(aes(color = Category), alpha = 0.6) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray")) +
      labs(title = paste("Volcano Plot - Control vs Infected:", title),
           x = "Log Fold Change (logFC) - Infected vs Control",
           y = "-Log10 Adjusted p-value") +
      theme_minimal() +
      # Add p-value cutoff line
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  
      # Add fold change cutoff lines
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "green") +   
      theme(legend.title = element_blank()) +
      theme(legend.position = "top")
}

# COVID Dataset
create_volcano_plot(Genes, "COVID")

# Psychiatric Disease Dataset
create_volcano_plot(Genes, "Psychiatric Disease")

----------------------------------------
## KEGG PATHWAY ##
----------------------------------------
 
#Read Shared Genes
genes <- read.csv(
    "shared_entrez_numerical_IDs_for_KEGG_GO_Pathway_Enrichment.csv",
    header = TRUE
  )

gene_list <- genes$Symbol
head(gene_list)

#SYMBOL → ENTREZID
gene_df <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_genes <- gene_df$ENTREZID

head(entrez_genes)


#Enrichment
ego <- enrichGO(
  gene          = entrez_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.2
)

#Plots
p1 <- dotplot(ego, showCategory = 10)
print(p1)
