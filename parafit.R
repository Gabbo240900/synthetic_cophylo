# Load required libraries
library(ape)
library(ade4)

# Function to compute ParaFit from file paths
compute_parafit <- function(host_tree_path, parasite_tree_path, association_path, nperm = 500) {
  # Read trees
  host_tree <- read.tree(host_tree_path)
  parasite_tree <- read.tree(parasite_tree_path)
  host_tree$tip.label <- sub(":.*", "", host_tree$tip.label)
  parasite_tree$tip.label <- sub(":.*", "", parasite_tree$tip.label)

  # Read association table: should be two columns (parasite, host)
  assoc <- read.table(association_path, header = FALSE, stringsAsFactors = FALSE)

  # Ensure unique labels
  host_labels <- host_tree$tip.label
  parasite_labels <- parasite_tree$tip.label

  # Filter associations to valid tips only
  assoc <- assoc[assoc[,1] %in% parasite_labels & assoc[,2] %in% host_labels, ]

  # Create distance matrices
  D_host <- cophenetic(host_tree)
  D_parasite <- cophenetic(parasite_tree)

  # Match rows in same order
  D_host <- D_host[sort(rownames(D_host)), sort(colnames(D_host))]
  D_parasite <- D_parasite[sort(rownames(D_parasite)), sort(colnames(D_parasite))]

  # Create binary association matrix
  assoc_matrix <- matrix(0, nrow = length(parasite_labels), ncol = length(host_labels))
  rownames(assoc_matrix) <- parasite_labels
  colnames(assoc_matrix) <- host_labels
  for (i in 1:nrow(assoc)) {
    p <- assoc[i,1]
    h <- assoc[i,2]
    assoc_matrix[p, h] <- 1
  }

    # Remove rows/columns not used in associations
    assoc_matrix <- assoc_matrix[rowSums(assoc_matrix) > 0, colSums(assoc_matrix) > 0]
    D_host <- D_host[colnames(assoc_matrix), colnames(assoc_matrix), drop=FALSE]
    D_parasite <- D_parasite[rownames(assoc_matrix), rownames(assoc_matrix), drop=FALSE]

    # Check if any matrix is empty
    if (nrow(D_host) == 0 || nrow(D_parasite) == 0 || nrow(assoc_matrix) == 0) {
    cat("Skipping: one or more matrices are empty after filtering.\n")
    return(NULL)
    }

    # Run ParaFit with Lingoes correction for negative eigenvalues
    result <- parafit(D_host, D_parasite, assoc_matrix, nperm = nperm, test = TRUE, correction = "lingoes")

    # Print main stats
    cat("ParaFit Global Stat:", result$ParaFitGlobal, "\n")
    cat("ParaFit Global P-value:", result$p.global, "\n")
}


# Automatically run if script is called with arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  host_tree_path <- args[1]
  parasite_tree_path <- args[2]
  association_path <- args[3]
  nperm <- ifelse(length(args) >= 4, as.numeric(args[4]), 500)
  compute_parafit(host_tree_path, parasite_tree_path, association_path, nperm)
}