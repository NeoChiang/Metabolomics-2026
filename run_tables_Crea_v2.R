# ============================================================================
# Statistical tables with sample-list filtering + P beside each timepoint
# Creatinine normalisation
# ============================================================================

pkgs <- c("pls", "openxlsx")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({ library(pls); library(openxlsx) })

set.seed(42)

# ---- 1. load data + sample lists -------------------------------------------
dat <- read.csv("6M-2Y_Urine_all metabolites_20260330.csv",
                check.names = FALSE, stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
dat$caseno <- sub("_(6M|2Y)$", "", dat[["Sample name"]])

f6 <- read.xlsx("table1/139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx", sheet = 1)
f2 <- read.xlsx("table1/2Y BPD_206 20250721.xlsx", sheet = 1)
ids_6M <- trimws(f6$case.no.)
ids_2Y <- trimws(f2$case.no.)

dat <- dat[(dat$Age == "6M" & dat$caseno %in% ids_6M) |
           (dat$Age == "2Y" & dat$caseno %in% ids_2Y), ]
cat("After sample-list filter: 6M =", sum(dat$Age == "6M"),
    ", 2Y =", sum(dat$Age == "2Y"), "\n")

# ---- 2. preprocessing ------------------------------------------------------
crea_col <- "Creatinine"
exclude  <- c(crea_col, "TSP")
met_cols <- setdiff(colnames(dat)[5:(ncol(dat)-1)], exclude)
# remove caseno helper col from met_cols
met_cols <- setdiff(met_cols, "caseno")

for (mc in c(met_cols, crea_col)) dat[[mc]] <- suppressWarnings(as.numeric(dat[[mc]]))

data_met <- as.matrix(dat[, met_cols])
crea_vec <- dat[[crea_col]]

# Step 1: negative -> 0
data_met[data_met < 0] <- 0
# Step 2: creatinine normalisation
data_norm <- data_met / crea_vec
# Step 3: log2 (half-min for zeros)
min_nz <- min(data_norm[data_norm > 0], na.rm = TRUE)
data_norm[data_norm == 0] <- min_nz / 2
data_log2 <- log2(data_norm)
# Step 4: Pareto scaling
col_means <- colMeans(data_log2)
col_sds   <- apply(data_log2, 2, sd)
col_sds[col_sds == 0] <- 1
data_pareto <- sweep(data_log2,   2, col_means,   "-")
data_pareto <- sweep(data_pareto, 2, sqrt(col_sds), "/")

cat("Preprocessing done.\n")

# ---- 3. helpers -------------------------------------------------------------
calc_vip <- function(X, y) {
  if (length(unique(y)) < 2 || nrow(X) < 4) return(rep(NA_real_, ncol(X)))
  y_num <- as.numeric(factor(y)) - 1
  df <- data.frame(y = y_num, X = I(as.matrix(X)))
  fit <- tryCatch(
    plsr(y ~ X, data = df, ncomp = 1, method = "oscorespls", validation = "none"),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, ncol(X)))
  W  <- loading.weights(fit)[, 1]
  SS <- drop(fit$Yloadings[1, 1]^2 * sum(scores(fit)[, 1]^2))
  sqrt(ncol(X) * (W^2 * SS) / SS)
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# ---- 4. build one table (new layout) ---------------------------------------
# Layout: Metabolite, then for EACH timepoint block:
#   VIP_comp1, FC_comp1, P_comp1, VIP_comp2, FC_comp2, P_comp2, VIP_comp3, FC_comp3, P_comp3
# So columns = Metabolite + 2 timepoints × 3 comparisons × 3 stats = 19 cols

build_table_v2 <- function(grp_col, comparisons, comp_labels) {
  ages <- c("6M", "2Y")
  results <- data.frame(Metabolite = met_cols, stringsAsFactors = FALSE)

  for (a in ages) {
    idx <- dat$Age == a
    pareto_sub <- data_pareto[idx, ]
    log2_sub   <- data_log2[idx, ]
    grp        <- dat[[grp_col]][idx]

    for (ci in seq_along(comparisons)) {
      comp <- comparisons[[ci]]
      g_num <- comp$num; g_den <- comp$den
      prefix <- paste0(a, "_", comp_labels[ci])

      sel <- grp %in% c(g_num, g_den)
      X_par <- pareto_sub[sel, ]
      X_log <- log2_sub[sel, ]
      g_sel <- grp[sel]

      vip_vec <- calc_vip(X_par, g_sel)
      vip_vals <- fc_vals <- p_vals <- numeric(length(met_cols))
      for (j in seq_along(met_cols)) {
        vip_vals[j] <- vip_vec[j]
        vals_num <- X_log[g_sel == g_num, j]
        vals_den <- X_log[g_sel == g_den, j]
        fc_vals[j] <- 2^(mean(vals_num) - mean(vals_den))
        x_a <- X_par[g_sel == g_num, j]
        x_b <- X_par[g_sel == g_den, j]
        p_vals[j] <- tryCatch(
          wilcox.test(x_a, x_b, exact = FALSE)$p.value,
          error = function(e) NA_real_)
      }
      results[[paste0("VIP_", prefix)]] <- round(vip_vals, 2)
      results[[paste0("FC_",  prefix)]] <- round(fc_vals, 2)
      results[[paste0("P_",   prefix)]] <- fmt_p(p_vals)
    }
  }

  # filter: at least one P < 0.05 across all columns
  p_cols <- grep("^P_", names(results), value = TRUE)
  keep <- apply(results[, p_cols, drop = FALSE], 1, function(row) {
    pnums <- suppressWarnings(as.numeric(sub("^<", "", row)))
    any(!is.na(pnums) & pnums < 0.05)
  })
  results <- results[keep, ]

  # sort by first P column (6M comp1) ascending
  first_p <- grep("^P_6M_", names(results), value = TRUE)[1]
  p1_num <- suppressWarnings(as.numeric(sub("^<", "", results[[first_p]])))
  results <- results[order(p1_num), ]
  rownames(results) <- NULL
  results
}

# ---- 5. define comparisons ------------------------------------------------
ga_comps <- list(
  list(num = 2, den = 0),
  list(num = 1, den = 0),
  list(num = 2, den = 1)
)
ga_labels <- c("comp1", "comp2", "comp3")

bpd_comps <- list(
  list(num = 2, den = 0),
  list(num = 1, den = 0),
  list(num = 2, den = 1)
)
bpd_labels <- c("comp1", "comp2", "comp3")

cat("\n=== Building tables ===\n")
t_ga  <- build_table_v2("GA3Group",  ga_comps,  ga_labels)
cat("Table_GA:  ", nrow(t_ga), "metabolites\n")
t_bpd <- build_table_v2("BPD3Group", bpd_comps, bpd_labels)
cat("Table_BPD: ", nrow(t_bpd), "metabolites\n")

# ---- 6. write CSVs --------------------------------------------------------
write.csv(t_ga,  "Table_GA_Crea_v2.csv",  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(t_bpd, "Table_BPD_Crea_v2.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---- 7. write xlsx --------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "GA");  writeData(wb, "GA",  t_ga)
addWorksheet(wb, "BPD"); writeData(wb, "BPD", t_bpd)
saveWorkbook(wb, "Table_summary_Crea_v2.xlsx", overwrite = TRUE)

cat("\nDone. Table_GA_Crea_v2.csv, Table_BPD_Crea_v2.csv, Table_summary_Crea_v2.xlsx written.\n")
