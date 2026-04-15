# ============================================================================
# Mfuzz c-means clustering on log2 fold-change (2Y / 6M)
# WITH standardise — finds group-differential temporal patterns
# ----------------------------------------------------------------------------
# Produces 3 figures + 3 membership CSVs:
#   FC_GA_std     : GA-stratified  3-dim
#   FC_BPD_std    : BPD-stratified 3-dim
#   FC_Combined_std : 6-dim
# ============================================================================

pkgs <- c("Mfuzz", "Biobase", "openxlsx", "RColorBrewer")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  library(Mfuzz); library(Biobase); library(openxlsx); library(RColorBrewer)
})

set.seed(42)

MIN_N <- 3

# ---- 1. significant metabolite union --------------------------------------
clean <- read.csv("all_tables_clean.csv", check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
sig_mets <- unique(clean$Metabolite[as.numeric(sub("<", "", clean$P)) < 0.05])
cat("Significant metabolites (union):", length(sig_mets), "\n")

# ---- 2. raw Excel ---------------------------------------------------------
# Read with colNames=FALSE to preserve original names (avoid make.names mangling)
raw_all <- read.xlsx("urine_data.xlsx", sheet = 1, startRow = 1, colNames = FALSE)
original_header <- as.character(raw_all[1, ])        # row 1 = metabolite names
colnames(raw_all) <- original_header
raw <- raw_all[-c(1, 2), ]                           # drop header + chem-shift row
raw$Age       <- as.character(raw$Age)
raw$GA3Group  <- as.integer(raw$GA3Group)
raw$BPD3Group <- as.integer(raw$BPD3Group)

met_cols <- colnames(raw)[5:ncol(raw)]
for (c in met_cols) raw[[c]] <- suppressWarnings(as.numeric(raw[[c]]))

mets <- intersect(sig_mets, met_cols)
missing_in_excel <- setdiff(sig_mets, met_cols)
cat("Matched in Excel:", length(mets), "\n")
if (length(missing_in_excel)) cat("  Missing:", paste(missing_in_excel, collapse=", "), "\n")

# ---- 3. helper: group log2 fold-change -------------------------------------
log2fc <- function(met, grpcol, grpval) {
  v6 <- raw[[met]][raw$Age == "6M" & raw[[grpcol]] == grpval]
  v2 <- raw[[met]][raw$Age == "2Y" & raw[[grpcol]] == grpval]
  v6 <- v6[is.finite(v6)]; v2 <- v2[is.finite(v2)]
  if (length(v6) < MIN_N || length(v2) < MIN_N) return(NA_real_)
  m6 <- mean(v6); m2 <- mean(v2)
  if (m6 <= 0 || m2 <= 0) return(NA_real_)
  log2(m2 / m6)
}

# ---- 4. build FC matrices --------------------------------------------------
ga_fc  <- matrix(NA_real_, length(mets), 3,
                 dimnames = list(mets, c("GA_>=37wk","GA_28-32wk","GA_<28wk")))
bpd_fc <- matrix(NA_real_, length(mets), 3,
                 dimnames = list(mets, c("BPD_HC","BPD_No+Mild","BPD_M+S")))
for (m in mets) {
  for (g in 0:2) ga_fc [m, g+1] <- log2fc(m, "GA3Group",  g)
  for (g in 0:2) bpd_fc[m, g+1] <- log2fc(m, "BPD3Group", g)
}
combined_fc <- cbind(ga_fc, bpd_fc)

clean_mat <- function(M) M[complete.cases(M), , drop = FALSE]
mats <- list(
  FC_GA_std       = clean_mat(ga_fc),
  FC_BPD_std      = clean_mat(bpd_fc),
  FC_Combined_std = clean_mat(combined_fc)
)
for (nm in names(mats))
  cat(sprintf("  %-18s: %d metabolites x %d dims\n", nm, nrow(mats[[nm]]), ncol(mats[[nm]])))

# ---- 5. run Mfuzz WITH standardise ----------------------------------------
run_fc_mfuzz_std <- function(M, tag, k) {
  eset   <- ExpressionSet(assayData = M)
  eset_s <- standardise(eset)   # <-- KEY DIFFERENCE: standardise!
  m_val  <- tryCatch(mestimate(eset_s), error = function(e) 1.25)
  if (!is.finite(m_val) || m_val < 1.05) m_val <- 1.25

  k_use <- min(k, nrow(M) - 1)
  cl <- mfuzz(eset_s, c = k_use, m = m_val)

  # --- plot ---
  pdf(paste0("mfuzz_", tag, ".pdf"), width = 11, height = 7)
  mfrow <- c(ceiling(k_use/3), min(k_use, 3))
  mfuzz.plot2(eset_s, cl = cl, mfrow = mfrow,
              time.labels = colnames(M),
              centre = TRUE, x11 = FALSE,
              ylab = "Standardised log2-FC")
  mtext(paste0("Mfuzz log2-FC (standardised) - ", tag,
               "  (k=", k_use, ", m=", round(m_val, 2),
               ", n=", nrow(M), ")"),
        side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
  dev.off()

  # --- CSV ---
  memb <- cl$membership
  colnames(memb) <- paste0("Cluster_", seq_len(ncol(memb)))

  # Also include the raw (unstandardised) FC values for interpretation
  out <- data.frame(Metabolite = rownames(M),
                    HardCluster = cl$cluster,
                    memb, M, check.names = FALSE)
  write.csv(out, paste0("mfuzz_", tag, "_membership.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("  %-18s  k=%d  m=%.2f  n=%d  ->  mfuzz_%s.pdf\n",
              tag, k_use, m_val, nrow(M), tag))

  # --- print cluster members ---
  cat(sprintf("\n  Cluster assignments for %s:\n", tag))
  for (i in 1:k_use) {
    members <- names(cl$cluster[cl$cluster == i])
    cat(sprintf("    Cluster %d (%d metabolites): %s\n",
                i, length(members), paste(members, collapse = ", ")))
  }
  cat("\n")

  invisible(cl)
}

cat("\n=== Running log2-FC Mfuzz (STANDARDISED) ===\n")
run_fc_mfuzz_std(mats$FC_GA_std,       "FC_GA_std",       k = 3)
run_fc_mfuzz_std(mats$FC_BPD_std,      "FC_BPD_std",      k = 3)
run_fc_mfuzz_std(mats$FC_Combined_std, "FC_Combined_std", k = 4)

cat("\nDone. Check mfuzz_FC_*_std.pdf files.\n")
