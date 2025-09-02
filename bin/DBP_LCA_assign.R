#!/usr/bin/env Rscript
# DBP_LCA_assign.R
# Lowest Common Ancestor (LCA) taxonomy assignment for DBP Metabarcoding Pipeline v2

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

# -------- Argument parsing --------
args <- commandArgs(trailingOnly = TRUE)

# Default values
blast_file     <- "result_blastn.txt"       # input from pipeline step 7
tax_file       <- "database/database.txt"   # two-column mapping: accession <TAB> lineage
min_pident     <- 97
max_evalue     <- 1e-20
bitscore_delta <- 5
out_prefix     <- "DBP_results"

# Assign from args if provided
if (length(args) >= 1) min_pident     <- as.numeric(args[1])
if (length(args) >= 2) max_evalue     <- as.numeric(args[2])
if (length(args) >= 3) bitscore_delta <- as.numeric(args[3])
if (length(args) >= 4) out_prefix     <- args[4]

cat("▶ Running LCA with parameters:\n")
cat("   min_pident     =", min_pident, "\n")
cat("   max_evalue     =", max_evalue, "\n")
cat("   bitscore_delta =", bitscore_delta, "\n")
cat("   out_prefix     =", out_prefix, "\n\n")

# -------- Helper functions --------
split_lineage <- function(lin) {
  if (is.na(lin) || lin == "" || toupper(lin) == "NA") return(rep(NA,7))
  parts <- unlist(strsplit(lin, ";"))
  length(parts) <- 7
  return(parts)
}

compute_lca <- function(lineages_list) {
  if (length(lineages_list) == 0) return(list(lineage = NA, rank = NA))
  mat <- do.call(rbind, lineages_list)
  lca <- rep(NA, ncol(mat))
  for (i in seq_len(ncol(mat))) {
    vals <- unique(na.omit(mat[,i]))
    if (length(vals) == 1) lca[i] <- vals else break
  }
  ranks <- c("kingdom","phylum","class","order","family","genus","species")
  deepest_idx <- which(!is.na(lca))
  if (length(deepest_idx) == 0) return(list(lineage = NA, rank = NA))
  deepest <- max(deepest_idx)
  return(list(lineage = paste(lca[1:deepest], collapse=";"), rank = ranks[deepest]))
}

# -------- Load BLAST results --------
cat("Reading BLAST results from:", blast_file, "\n")
blast_raw <- read_tsv(blast_file, col_names = FALSE, comment = "#", progress = FALSE)

if (ncol(blast_raw) >= 12) {
  names(blast_raw)[1:12] <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                              "qstart","qend","sstart","send","evalue","bitscore")
} else {
  names(blast_raw)[1:2] <- c("qseqid","sseqid")
  warning("⚠️ BLAST file has <12 columns, using first two as qseqid/sseqid.")
}

# Ensure required columns exist
for (cname in c("qseqid","sseqid","pident","evalue","bitscore")) {
  if (!cname %in% names(blast_raw)) blast_raw[[cname]] <- NA
}

# Convert numeric safely
blast_raw$pident   <- suppressWarnings(as.numeric(blast_raw$pident))
blast_raw$evalue   <- suppressWarnings(as.numeric(blast_raw$evalue))
blast_raw$bitscore <- suppressWarnings(as.numeric(blast_raw$bitscore))

# -------- Load taxonomy mapping --------
cat("Reading taxonomy mapping from:", tax_file, "\n")
taxonomy <- read_tsv(tax_file, col_names = FALSE, progress = FALSE)
if (ncol(taxonomy) < 2) stop("❌ ERROR: database.txt must be two columns: accession<TAB>lineage")
names(taxonomy)[1:2] <- c("accession","lineage")

# -------- Filter BLAST hits --------
blast <- blast_raw %>%
  filter(!is.na(pident) & pident >= min_pident) %>%
  filter(!is.na(evalue) & evalue <= max_evalue)

if (nrow(blast) == 0) stop("❌ No BLAST hits passed filters. Try adjusting thresholds.")

# -------- LCA assignment --------
queries <- unique(blast$qseqid)
cat("Computing LCA for", length(queries), "queries...\n")

result_list <- vector("list", length(queries))
for (i in seq_along(queries)) {
  q <- queries[i]
  qdf <- blast[blast$qseqid == q, , drop = FALSE]
  result_list[[i]] <- tryCatch({
    top_bits <- max(qdf$bitscore, na.rm = TRUE)
    qdf2 <- qdf %>% filter(bitscore >= (top_bits - bitscore_delta))
    accs <- sapply(str_split(as.character(qdf2$sseqid), "\\s+"), `[`, 1)
    qdf2$accession_clean <- accs
    tax_hits <- left_join(qdf2, taxonomy, by = c("accession_clean" = "accession"))
    lineages <- lapply(tax_hits$lineage, split_lineage)
    lca <- compute_lca(lineages)
    tibble(
      qseqid = q,
      n_hits_considered = nrow(qdf2),
      assigned_lineage  = ifelse(is.na(lca$lineage), "Unassigned", lca$lineage),
      assigned_rank     = ifelse(is.na(lca$rank), "Unassigned", lca$rank),
      supporting_hits   = paste0(qdf2$accession_clean, ":", qdf2$lineage, collapse = ";")
    )
  }, error = function(e) {
    warning("Error with query ", q, ": ", e$message)
    tibble(qseqid = q, n_hits_considered = NA,
           assigned_lineage = "Error", assigned_rank = NA,
           supporting_hits = NA)
  })
}

lca_results <- bind_rows(result_list)

# -------- Save output --------
out_file <- paste0(out_prefix, "_LCA_results.tsv")
write_tsv(lca_results, out_file)
cat("✅ LCA assignments written to:", out_file, "\n")
