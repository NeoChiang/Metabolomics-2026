# ============================================================================
# Mfuzz c-means clustering of 6M-to-2Y urinary metabolite trajectories
# ----------------------------------------------------------------------------
# Input : all_tables_clean.csv       (significant metabolites union)
#         urine_data.xlsx            (raw sample-level intensities)
# Output: 6 PDFs + 6 CSVs of cluster membership
# ============================================================================

# ---- packages --------------------------------------------------------------
pkgs <- c("Mfuzz", "Biobase", "openxlsx", "RColorBrewer")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  library(Mfuzz)
  library(Biobase)
  library(openxlsx)
  library(RColorBrewer)
})

set.seed(42)

# ---- 1. load significant metabolite union ---------------------------------
clean <- read.csv("all_tables_clean.csv", check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
sig_mets <- unique(clean$Metabolite[as.numeric(sub("<", "", clean$P)) < 0.05])
cat("Significant metabolites (union):", length(sig_mets), "\n")

# ---- 2. load raw Excel, compute group-wise means ---------------------------
raw <- read.xlsx("urine_data.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
# Row 1 of the sheet is a chemical-shift annotation row -> drop it
raw <- raw[-1, ]
raw$Age      <- as.character(raw$Age)
raw$GA3Group <- as.integer(raw$GA3Group)
raw$BPD3Group<- as.integer(raw$BPD3Group)

# metabolite columns = 5 : ncol
met_cols <- colnames(raw)[5:ncol(raw)]
for (c in met_cols) raw[[c]] <- suppressWarnings(as.numeric(raw[[c]]))

# keep only the significant metabolites that also exist in the Excel
mets <- intersect(sig_mets, met_cols)
cat("Matched in Excel:", length(mets),
    "  (missing:", setdiff(sig_mets, met_cols), ")\n")

# helper: group mean (NA-safe)
grp_mean <- function(age, grpcol, grpval, met) {
  v <- raw[[met]][raw$Age == age & raw[[grpcol]] == grpval]
  if (length(v) == 0 || all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
}
overall_mean <- function(age, met) {
  v <- raw[[met]][raw$Age == age]
  if (length(v) == 0 || all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
}

# ---- 3. build the 6 input matrices ----------------------------------------
mat_ga6  <- matrix(NA_real_, length(mets), 6,
                   dimnames = list(mets, c("6M_0","6M_1","6M_2","2Y_0","2Y_1","2Y_2")))
mat_bpd6 <- mat_ga6
mat_ga2  <- matrix(NA_real_, length(mets), 2, dimnames = list(mets, c("6M","2Y")))
mat_bpd2 <- mat_ga2
mat_12   <- matrix(NA_real_, length(mets), 12,
                   dimnames = list(mets,
                     c("GA_6M_0","GA_6M_1","GA_6M_2","GA_2Y_0","GA_2Y_1","GA_2Y_2",
                       "BPD_6M_0","BPD_6M_1","BPD_6M_2","BPD_2Y_0","BPD_2Y_1","BPD_2Y_2")))
mat_4    <- matrix(NA_real_, length(mets), 4,
                   dimnames = list(mets, c("GA_6M","GA_2Y","BPD_6M","BPD_2Y")))

for (m in mets) {
  # 6-point GA
  for (g in 0:2) {
    mat_ga6[m, paste0("6M_", g)] <- grp_mean("6M", "GA3Group", g, m)
    mat_ga6[m, paste0("2Y_", g)] <- grp_mean("2Y", "GA3Group", g, m)
  }
  # 6-point BPD
  for (g in 0:2) {
    mat_bpd6[m, paste0("6M_", g)] <- grp_mean("6M", "BPD3Group", g, m)
    mat_bpd6[m, paste0("2Y_", g)] <- grp_mean("2Y", "BPD3Group", g, m)
  }
  # 2-point means
  mat_ga2[m, "6M"] <- overall_mean("6M", m)
  mat_ga2[m, "2Y"] <- overall_mean("2Y", m)
  mat_bpd2[m, ]    <- mat_ga2[m, ]  # identical (same samples)
  # 12-point combined
  mat_12[m, 1:6]   <- mat_ga6[m, ]
  mat_12[m, 7:12]  <- mat_bpd6[m, ]
  # 4-point combined (all equal to overall mean pair, duplicated)
  mat_4[m, ]       <- rep(mat_ga2[m, ], 2)
}

# Drop metabolites with any NA (Mfuzz can't handle)
clean_mat <- function(M) M[complete.cases(M), , drop = FALSE]

mats <- list(
  GA_6pt      = clean_mat(mat_ga6),
  GA_2pt      = clean_mat(mat_ga2),
  BPD_6pt     = clean_mat(mat_bpd6),
  BPD_2pt     = clean_mat(mat_bpd2),
  Combined_12 = clean_mat(mat_12),
  Combined_4  = clean_mat(mat_4)
)
for (nm in names(mats)) cat(sprintf("  %-12s: %d metabolites x %d cols\n",
                                    nm, nrow(mats[[nm]]), ncol(mats[[nm]])))

# ---- 4. Mfuzz clustering ---------------------------------------------------
run_mfuzz <- function(M, tag, k = 4) {
  eset   <- ExpressionSet(assayData = M)
  eset_s <- standardise(eset)
  m_val  <- tryCatch(mestimate(eset_s), error = function(e) 1.25)
  if (!is.finite(m_val) || m_val < 1.05) m_val <- 1.25

  # For 2-point input we can only make k<=2 meaningful; cap k
  k_use <- min(k, nrow(M) - 1, max(2, ncol(M)))

  cl <- mfuzz(eset_s, c = k_use, m = m_val)

  # --- plot ----
  pdf(paste0("mfuzz_", tag, ".pdf"), width = 10, height = 7)
  mfrow <- c(ceiling(k_use/3), min(k_use, 3))
  mfuzz.plot2(eset_s, cl = cl, mfrow = mfrow,
              time.labels = colnames(M),
              centre = TRUE, x11 = FALSE,
              main = paste0("Mfuzz clusters - ", tag,
                            "  (k=", k_use, ", m=", round(m_val, 2), ")"))
  dev.off()

  # --- membership CSV ----
  memb <- cl$membership
  colnames(memb) <- paste0("Cluster_", seq_len(ncol(memb)))
  out <- data.frame(Metabolite = rownames(M),
                    HardCluster = cl$cluster,
                    memb,
                    M, check.names = FALSE)
  write.csv(out, paste0("mfuzz_", tag, "_membership.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("  %-12s done  (k=%d, m=%.2f, n=%d)\n",
              tag, k_use, m_val, nrow(M)))
  invisible(cl)
}

cat("\n=== Running Mfuzz ===\n")
for (nm in names(mats)) {
  # k=4 for 6/12-pt, k=3 for 2/4-pt
  k <- if (ncol(mats[[nm]]) <= 4) 3 else 4
  run_mfuzz(mats[[nm]], nm, k = k)
}

cat("\nAll done. PDFs + CSVs saved in current directory.\n")
