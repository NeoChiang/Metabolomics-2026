# ============================================================================
# Mfuzz c-means clustering on log2 fold-change (2Y / 6M) trajectories
# ----------------------------------------------------------------------------
# Produces 3 figures + 3 membership CSVs:
#   FC_GA     : GA-stratified  3-dim  [FC_0, FC_1, FC_2]
#   FC_BPD    : BPD-stratified 3-dim  [FC_0, FC_1, FC_2]
#   FC_Combined : 6-dim [GA_FC_0..2, BPD_FC_0..2]
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

MIN_N <- 3   # minimum samples per (age x group) cell

# ---- 1. significant metabolite union --------------------------------------
clean <- read.csv("all_tables_clean.csv", check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
sig_mets <- unique(clean$Metabolite[as.numeric(sub("<", "", clean$P)) < 0.05])
cat("Significant metabolites (union):", length(sig_mets), "\n")

# ---- 2. raw Excel ---------------------------------------------------------
raw <- read.xlsx("urine_data.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
raw <- raw[-1, ]                                      # drop chem-shift row
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

clean_mat    <- function(M) M[complete.cases(M), , drop = FALSE]
dropped_mets <- function(M) rownames(M)[!complete.cases(M)]

raw_mats <- list(FC_GA = ga_fc, FC_BPD = bpd_fc, FC_Combined = combined_fc)
mats     <- lapply(raw_mats, clean_mat)

cat("\n=== Dropped metabolites (any NA in row) per tag ===\n")
for (nm in names(raw_mats)) {
  d <- dropped_mets(raw_mats[[nm]])
  cat(sprintf("  %-12s: kept %d / dropped %d %s\n",
              nm, nrow(mats[[nm]]), length(d),
              if (length(d)) paste0("[", paste(d, collapse=", "), "]") else ""))
}

# ---- 5. run Mfuzz (NO standardisation for FC) -----------------------------
# log2-FC is already on a symmetric, biologically meaningful scale around 0.
# Standardising would wipe out shared-direction signal ("all three groups go up").
run_fc_mfuzz <- function(M, tag, k) {
  eset  <- ExpressionSet(assayData = M)
  m_val <- tryCatch(mestimate(eset), error = function(e) 1.25)
  if (!is.finite(m_val) || m_val < 1.05) m_val <- 1.25

  k_use <- min(k, nrow(M) - 1)
  cl <- mfuzz(eset, c = k_use, m = m_val)

  # --- plot ---
  pdf(paste0("mfuzz_", tag, ".pdf"), width = 11, height = 7)
  mfrow <- c(ceiling(k_use/3), min(k_use, 3))
  mfuzz.plot2(eset, cl = cl, mfrow = mfrow,
              time.labels = colnames(M),
              centre = TRUE, x11 = FALSE,
              ylab = "log2(2Y / 6M)")
  abline(h = 0, lty = 2, col = "grey40")
  mtext(paste0("Mfuzz log2-FC - ", tag,
               "  (k=", k_use, ", m=", round(m_val, 2),
               ", n=", nrow(M), ")"),
        side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
  dev.off()

  # --- CSV ---
  memb <- cl$membership
  colnames(memb) <- paste0("Cluster_", seq_len(ncol(memb)))
  out <- data.frame(Metabolite = rownames(M),
                    HardCluster = cl$cluster,
                    memb, M, check.names = FALSE)
  write.csv(out, paste0("mfuzz_", tag, "_membership.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("  %-12s  k=%d  m=%.2f  n=%d  ->  mfuzz_%s.pdf + _membership.csv\n",
              tag, k_use, m_val, nrow(M), tag))
  invisible(cl)
}

cat("\n=== Running log2-FC Mfuzz ===\n")
run_fc_mfuzz(mats$FC_GA,       "FC_GA",       k = 3)
run_fc_mfuzz(mats$FC_BPD,      "FC_BPD",      k = 3)
run_fc_mfuzz(mats$FC_Combined, "FC_Combined", k = 4)

cat("\nDone.\n")
