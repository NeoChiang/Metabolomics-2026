# ============================================================================
# Statistical tables: VIP + FC + Wilcoxon for GA / BPD × 6M / 2Y
# ** Creatinine normalisation (no sum normalisation) **
# ============================================================================

pkgs <- c("pls", "openxlsx")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({ library(pls); library(openxlsx) })

set.seed(42)

# ---- 1. load data ----------------------------------------------------------
dat <- read.csv("6M-2Y_Urine_all metabolites_20260330.csv",
                check.names = FALSE, stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
cat("Rows:", nrow(dat), " Cols:", ncol(dat), "\n")

# columns: 1-4 = meta, 5 = Creatinine, 6:(ncol-1) = metabolites, ncol = TSP
crea_col <- "Creatinine"
exclude  <- c(crea_col, "TSP")
met_cols <- setdiff(colnames(dat)[5:ncol(dat)], exclude)
cat("Metabolite columns:", length(met_cols), "\n")

# coerce to numeric
for (mc in c(met_cols, crea_col)) dat[[mc]] <- suppressWarnings(as.numeric(dat[[mc]]))
cat("Creatinine range:", range(dat[[crea_col]], na.rm = TRUE), "\n")

# ---- 2. preprocessing (on ALL samples together) ----------------------------
data_met <- as.matrix(dat[, met_cols])
crea_vec <- dat[[crea_col]]

# Step 1: negative -> 0
data_met[data_met < 0] <- 0

# Step 2: Creatinine normalisation (each sample divided by its own creatinine)
data_norm <- data_met / crea_vec

# Step 3: log2 with half-min imputation for zeros
min_nonzero <- min(data_norm[data_norm > 0], na.rm = TRUE)
data_norm[data_norm == 0] <- min_nonzero / 2
data_log2 <- log2(data_norm)

# Step 4: Pareto scaling
col_means <- colMeans(data_log2)
col_sds   <- apply(data_log2, 2, sd)
col_sds[col_sds == 0] <- 1
data_pareto <- sweep(data_log2,   2, col_means,     "-")
data_pareto <- sweep(data_pareto, 2, sqrt(col_sds), "/")

cat("Preprocessing done (Creatinine normalisation).\n")

# ---- 3. VIP helper (from pls package) -------------------------------------
calc_vip <- function(X, y) {
  if (length(unique(y)) < 2 || nrow(X) < 4) return(rep(NA_real_, ncol(X)))
  y_num <- as.numeric(factor(y)) - 1
  df <- data.frame(y = y_num, X = I(as.matrix(X)))
  fit <- tryCatch(
    plsr(y ~ X, data = df, ncomp = 1, method = "oscorespls", validation = "none"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, ncol(X)))

  W  <- loading.weights(fit)[, 1]
  SS <- drop(fit$Yloadings[1, 1]^2 * sum(scores(fit)[, 1]^2))
  p  <- ncol(X)
  vip <- sqrt(p * (W^2 * SS) / SS)
  vip
}

# ---- 4. format p-value ----------------------------------------------------
fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# ---- 5. build one table ---------------------------------------------------
build_table <- function(age_val, grp_col, comparisons) {
  idx <- dat$Age == age_val
  pareto_sub <- data_pareto[idx, ]
  log2_sub   <- data_log2[idx, ]
  grp        <- dat[[grp_col]][idx]

  results <- data.frame(Metabolite = met_cols, stringsAsFactors = FALSE)

  for (ci in seq_along(comparisons)) {
    comp <- comparisons[[ci]]
    g_num <- comp$num
    g_den <- comp$den
    tag   <- paste0("comp", ci)

    sel <- grp %in% c(g_num, g_den)
    X_par <- pareto_sub[sel, ]
    X_log <- log2_sub[sel, ]
    g_sel <- grp[sel]

    # VIP
    vip_vec <- calc_vip(X_par, g_sel)

    vip_vals <- fc_vals <- p_vals <- numeric(length(met_cols))
    for (j in seq_along(met_cols)) {
      vip_vals[j] <- vip_vec[j]

      # FC: 2^(mean_num - mean_den) using log2 data
      vals_num <- X_log[g_sel == g_num, j]
      vals_den <- X_log[g_sel == g_den, j]
      fc_vals[j] <- 2^(mean(vals_num) - mean(vals_den))

      # Wilcoxon on pareto data
      x_a <- X_par[g_sel == g_num, j]
      x_b <- X_par[g_sel == g_den, j]
      p_vals[j] <- tryCatch(
        wilcox.test(x_a, x_b, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }

    results[[paste0("VIP_", tag)]] <- round(vip_vals, 2)
    results[[paste0("FC_",  tag)]] <- round(fc_vals, 2)
    results[[paste0("P_",   tag)]] <- fmt_p(p_vals)
  }

  # filter: at least one comparison p < 0.05
  p_cols <- grep("^P_comp", names(results), value = TRUE)
  keep <- apply(results[, p_cols, drop = FALSE], 1, function(row) {
    pnums <- suppressWarnings(as.numeric(sub("^<", "", row)))
    any(!is.na(pnums) & pnums < 0.05)
  })
  results <- results[keep, ]

  # sort by comp1 p ascending
  p1_num <- suppressWarnings(as.numeric(sub("^<", "", results$P_comp1)))
  results <- results[order(p1_num), ]
  rownames(results) <- NULL
  results
}

# ---- 6. define comparisons ------------------------------------------------
ga_comps <- list(
  list(num = 2, den = 0),
  list(num = 1, den = 0),
  list(num = 2, den = 1)
)
bpd_comps <- list(
  list(num = 2, den = 0),
  list(num = 1, den = 0),
  list(num = 2, den = 1)
)

cat("\n=== Building tables (Creatinine normalisation) ===\n")
t1 <- build_table("6M", "GA3Group",  ga_comps)
cat("Table1_6M_GA:  ", nrow(t1), "metabolites\n")
t2 <- build_table("2Y", "GA3Group",  ga_comps)
cat("Table2_2Y_GA:  ", nrow(t2), "metabolites\n")
t3 <- build_table("6M", "BPD3Group", bpd_comps)
cat("Table3_6M_BPD: ", nrow(t3), "metabolites\n")
t4 <- build_table("2Y", "BPD3Group", bpd_comps)
cat("Table4_2Y_BPD: ", nrow(t4), "metabolites\n")

# ---- 7. write CSVs --------------------------------------------------------
write.csv(t1, "Table1_6M_GA_Crea.csv",  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(t2, "Table2_2Y_GA_Crea.csv",  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(t3, "Table3_6M_BPD_Crea.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(t4, "Table4_2Y_BPD_Crea.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---- 8. write xlsx --------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "6M_GA");  writeData(wb, "6M_GA",  t1)
addWorksheet(wb, "2Y_GA");  writeData(wb, "2Y_GA",  t2)
addWorksheet(wb, "6M_BPD"); writeData(wb, "6M_BPD", t3)
addWorksheet(wb, "2Y_BPD"); writeData(wb, "2Y_BPD", t4)
saveWorkbook(wb, "Table_summary_Crea.xlsx", overwrite = TRUE)

cat("\nDone. 4 CSVs + Table_summary_Crea.xlsx written.\n")
