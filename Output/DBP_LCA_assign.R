#!/usr/bin/env Rscript
# blast_lca_full_safe.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

## ---------- EDIT these if needed ----------
blast_file <- "result_blastn.txt"
tax_file   <- "database/database.txt"
min_pident <- 97
max_evalue <- 1e-20
bitscore_delta <- 5
top_n <- 0
out_prefix <- "mydataset"
## -----------------------------------------

msg <- function(...) cat(sprintf(...), "\n")

# safe extract size
extract_size <- function(qid) {
  if (is.na(qid) || qid == "") return(1L)
  m <- str_match(qid, "size=(\\d+)")
  if (!is.null(m) && ncol(m) >= 2 && !is.na(m[1,2])) {
    return(as.integer(m[1,2]))
  }
  # try alternative patterns
  m2 <- str_match(qid, "count=(\\d+)")
  if (!is.null(m2) && ncol(m2) >= 2 && !is.na(m2[1,2])) return(as.integer(m2[1,2]))
  return(1L)
}

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

# ---------- Read BLAST robustly ----------
msg("Reading BLAST file: %s", blast_file)
blast_raw <- read_tsv(blast_file, col_names = FALSE, comment = "#", progress = FALSE)

ncol_bl <- ncol(blast_raw)
msg("Detected %d columns in BLAST file (first row shown):", ncol_bl)
print(head(blast_raw, 3))

if (ncol_bl < 2) stop("ERROR: BLAST file appears empty or not tab-delimited.")

# Assign qseqid = first column, sseqid = second column. If at least 12 cols, set names for first 12.
if (ncol_bl >= 12) {
  names(blast_raw)[1:12] <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                              "qstart","qend","sstart","send","evalue","bitscore")
  msg("Assigned standard BLAST column names to first 12 columns.")
} else {
  # Still assign qseqid and sseqid to first two columns
  names(blast_raw)[1:2] <- c("qseqid","sseqid")
  msg("Warning: BLAST has <12 columns (n=%d). I set first two columns to qseqid and sseqid.", ncol_bl)
  msg("If your file uses a different delimiter or has headers, re-check the file format.")
}

# Ensure the columns we need exist (or create as NA if missing)
required_cols <- c("qseqid","sseqid","pident","evalue","bitscore")
for (cname in required_cols) {
  if (!cname %in% names(blast_raw)) {
    blast_raw[[cname]] <- NA
  }
}

# Convert numeric columns safely
blast_raw$pident <- suppressWarnings(as.numeric(blast_raw$pident))
blast_raw$evalue <- suppressWarnings(as.numeric(blast_raw$evalue))
blast_raw$bitscore <- suppressWarnings(as.numeric(blast_raw$bitscore))

msg("Columns after normalization:")
print(colnames(blast_raw)[1:min(12, length(colnames(blast_raw)))])

# Load taxonomy file
msg("Reading taxonomy mapping: %s", tax_file)
taxonomy <- read_tsv(tax_file, col_names = FALSE, progress = FALSE)
if (ncol(taxonomy) < 2) stop("ERROR: taxonomy_file.txt must be two columns: accession[TAB]lineage")
names(taxonomy)[1:2] <- c("accession","lineage")

# ---------- Filter BLAST hits ----------
msg("Filtering hits: pident >= %s, evalue <= %s", min_pident, max_evalue)
blast <- blast_raw %>%
  filter(is.na(pident) == FALSE & pident >= min_pident) %>%
  filter(is.na(evalue) == FALSE & evalue <= max_evalue)

msg("After filtering: %d rows remain", nrow(blast))
if (nrow(blast) == 0) {
  stop("No BLAST hits passed the pident/evalue filters. Try loosening thresholds or inspect your BLAST file.")
}

# ---------- LCA loop (robust) ----------
queries <- unique(blast$qseqid)
msg("Computing LCA for %d queries (this may take a moment)...", length(queries))

result_list <- vector("list", length(queries))
i <- 1L
for (q in queries) {
  qdf <- blast[blast$qseqid == q, , drop = FALSE]
  result_list[[i]] <- tryCatch({
    if (!is.na(top_n) && top_n > 0) {
      qdf2 <- qdf %>% arrange(desc(bitscore)) %>% slice_head(n = top_n)
    } else {
      top_bits <- max(qdf$bitscore, na.rm = TRUE)
      if (is.na(top_bits)) top_bits <- -Inf
      qdf2 <- qdf %>% filter(bitscore >= (top_bits - bitscore_delta))
    }
    accs <- sapply(str_split(as.character(qdf2$sseqid), "\\s+"), `[`, 1)
    qdf2$accession_clean <- accs
    tax_hits <- left_join(qdf2, taxonomy, by = c("accession_clean" = "accession"))
    lineages <- lapply(tax_hits$lineage, split_lineage)
    lca <- compute_lca(lineages)
    tibble::tibble(
      qseqid = q,
      n_hits_considered = nrow(qdf2),
      assigned_lineage = ifelse(is.na(lca$lineage), "Unassigned", lca$lineage),
      assigned_rank = ifelse(is.na(lca$rank), "Unassigned", lca$rank),
      supporting_hits = paste0(
        ifelse(is.null(qdf2$accession_clean), "", qdf2$accession_clean),
        ":",
        ifelse(is.null(qdf2$lineage), "NA", ifelse(is.na(qdf2$lineage), "NA", qdf2$lineage)),
        collapse = ";"
      )
    )
  }, error = function(e) {
    msg("ERROR query %s: %s", q, e$message)
    tibble::tibble(qseqid = q, n_hits_considered = NA_integer_,
                   assigned_lineage = "Error", assigned_rank = NA_character_,
                   supporting_hits = NA_character_)
  })
  i <- i + 1L
}

lca_results <- bind_rows(result_list)
msg("LCA computation done. Rows: %d", nrow(lca_results))
print(head(lca_results, 6))

# ---------- Save outputs ----------
write_tsv(lca_results, paste0(out_prefix, "_lca_assignments.tsv"))
msg("Written: %s", paste0(out_prefix, "_lca_assignments.tsv"))