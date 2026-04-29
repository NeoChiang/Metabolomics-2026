# ============================================================================
# Mfuzz c-means clustering on log2 fold-change (2Y / 6M)
# ** Creatinine-normalised version **
# Produces:
#   FC (no standardise):  3 PDFs + 3 CSVs  (FC_GA_Crea, FC_BPD_Crea, FC_Combined_Crea)
#   FC (standardised):    3 PDFs + 3 CSVs  (FC_GA_std_Crea, FC_BPD_std_Crea, FC_Combined_std_Crea)
#   Pretty plots:         6 PDFs           (pretty_FC_*_std_Crea_styleA/B.pdf)
# ============================================================================

pkgs <- c("Mfuzz", "Biobase", "RColorBrewer", "ggplot2", "dplyr", "tidyr", "patchwork", "scales")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(Mfuzz); library(Biobase); library(RColorBrewer)
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)
})

set.seed(42)
MIN_N <- 3

# ---- 1. significant metabolite union ----------------------------------------
clean <- read.csv("all_tables_clean.csv", check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
sig_mets <- unique(clean$Metabolite[as.numeric(sub("<", "", clean$P)) < 0.05])
cat("Significant metabolites (union):", length(sig_mets), "\n")

# ---- 2. load new CSV with Creatinine ----------------------------------------
raw <- read.csv("6M-2Y_Urine_all metabolites_20260330.csv",
                check.names = FALSE, stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
raw$Age       <- as.character(raw$Age)
raw$GA3Group  <- as.integer(raw$GA3Group)
raw$BPD3Group <- as.integer(raw$BPD3Group)

crea_col <- "Creatinine"
exclude  <- c(crea_col, "TSP")
met_cols <- setdiff(colnames(raw)[5:ncol(raw)], exclude)

for (mc in c(met_cols, crea_col)) raw[[mc]] <- suppressWarnings(as.numeric(raw[[mc]]))

# ---- 3. Creatinine normalisation --------------------------------------------
crea_vec <- raw[[crea_col]]
for (mc in met_cols) raw[[mc]] <- raw[[mc]] / crea_vec
cat("Creatinine normalisation done (range:",
    sprintf("%.2e ~ %.2e", min(crea_vec), max(crea_vec)), ")\n")

mets <- intersect(sig_mets, met_cols)
missing_in_csv <- setdiff(sig_mets, met_cols)
cat("Matched in CSV:", length(mets), "\n")
if (length(missing_in_csv))
  cat("  Missing:", paste(missing_in_csv, collapse = ", "), "\n")

# ---- 4. helper: group log2 fold-change --------------------------------------
log2fc <- function(met, grpcol, grpval) {
  v6 <- raw[[met]][raw$Age == "6M" & raw[[grpcol]] == grpval]
  v2 <- raw[[met]][raw$Age == "2Y" & raw[[grpcol]] == grpval]
  v6 <- v6[is.finite(v6)]; v2 <- v2[is.finite(v2)]
  if (length(v6) < MIN_N || length(v2) < MIN_N) return(NA_real_)
  m6 <- mean(v6); m2 <- mean(v2)
  if (m6 <= 0 || m2 <= 0) return(NA_real_)
  log2(m2 / m6)
}

# ---- 5. build FC matrices ---------------------------------------------------
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

raw_mats_fc <- list(FC_GA_Crea = ga_fc, FC_BPD_Crea = bpd_fc, FC_Combined_Crea = combined_fc)
mats_fc     <- lapply(raw_mats_fc, clean_mat)

raw_mats_std <- list(FC_GA_std_Crea = ga_fc, FC_BPD_std_Crea = bpd_fc, FC_Combined_std_Crea = combined_fc)
mats_std     <- lapply(raw_mats_std, clean_mat)

cat("\n=== Dropped metabolites (any NA) per tag ===\n")
for (nm in names(raw_mats_fc)) {
  d <- dropped_mets(raw_mats_fc[[nm]])
  cat(sprintf("  %-22s: kept %d / dropped %d %s\n",
              nm, nrow(mats_fc[[nm]]), length(d),
              if (length(d)) paste0("[", paste(d, collapse=", "), "]") else ""))
}

# ---- 6A. run Mfuzz — NO standardise -----------------------------------------
run_fc_mfuzz <- function(M, tag, k) {
  eset  <- ExpressionSet(assayData = M)
  m_val <- tryCatch(mestimate(eset), error = function(e) 1.25)
  if (!is.finite(m_val) || m_val < 1.05) m_val <- 1.25
  k_use <- min(k, nrow(M) - 1)
  cl <- mfuzz(eset, c = k_use, m = m_val)

  nc <- min(k_use, 3); nr <- ceiling(k_use / nc)
  sq <- 4.5
  pdf(paste0("mfuzz_", tag, ".pdf"), width = nc * sq + 1, height = nr * sq + 0.8)
  par(pty = "s")
  mfuzz.plot2(eset, cl = cl, mfrow = c(nr, nc), time.labels = colnames(M),
              centre = TRUE, x11 = FALSE, ylab = "log2(2Y / 6M)")
  abline(h = 0, lty = 2, col = "grey40")
  mtext(paste0("Mfuzz log2-FC (Crea-norm) - ", tag,
               "  (k=", k_use, ", m=", round(m_val, 2), ", n=", nrow(M), ")"),
        side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
  dev.off()

  memb <- cl$membership; colnames(memb) <- paste0("Cluster_", seq_len(ncol(memb)))
  out <- data.frame(Metabolite = rownames(M), HardCluster = cl$cluster,
                    memb, M, check.names = FALSE)
  write.csv(out, paste0("mfuzz_", tag, "_membership.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("  %-22s  k=%d  m=%.2f  n=%d\n", tag, k_use, m_val, nrow(M)))
  invisible(cl)
}

cat("\n=== Running log2-FC Mfuzz (NO standardise, Crea-norm) ===\n")
run_fc_mfuzz(mats_fc$FC_GA_Crea,       "FC_GA_Crea",       k = 3)
run_fc_mfuzz(mats_fc$FC_BPD_Crea,      "FC_BPD_Crea",      k = 3)
run_fc_mfuzz(mats_fc$FC_Combined_Crea, "FC_Combined_Crea", k = 4)

# ---- 6B. run Mfuzz — WITH standardise ---------------------------------------
run_fc_mfuzz_std <- function(M, tag, k) {
  eset   <- ExpressionSet(assayData = M)
  eset_s <- standardise(eset)
  m_val  <- tryCatch(mestimate(eset_s), error = function(e) 1.25)
  if (!is.finite(m_val) || m_val < 1.05) m_val <- 1.25
  k_use <- min(k, nrow(M) - 1)
  cl <- mfuzz(eset_s, c = k_use, m = m_val)

  nc <- min(k_use, 3); nr <- ceiling(k_use / nc)
  sq <- 4.5
  pdf(paste0("mfuzz_", tag, ".pdf"), width = nc * sq + 1, height = nr * sq + 0.8)
  par(pty = "s")
  mfuzz.plot2(eset_s, cl = cl, mfrow = c(nr, nc), time.labels = colnames(M),
              centre = TRUE, x11 = FALSE, ylab = "Standardised log2-FC")
  mtext(paste0("Mfuzz log2-FC (std, Crea-norm) - ", tag,
               "  (k=", k_use, ", m=", round(m_val, 2), ", n=", nrow(M), ")"),
        side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
  dev.off()

  memb <- cl$membership; colnames(memb) <- paste0("Cluster_", seq_len(ncol(memb)))
  out <- data.frame(Metabolite = rownames(M), HardCluster = cl$cluster,
                    memb, M, check.names = FALSE)
  write.csv(out, paste0("mfuzz_", tag, "_membership.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("  %-22s  k=%d  m=%.2f  n=%d\n", tag, k_use, m_val, nrow(M)))
  cat(sprintf("    Cluster assignments for %s:\n", tag))
  for (i in 1:k_use) {
    members <- names(cl$cluster[cl$cluster == i])
    cat(sprintf("      Cluster %d (%d): %s\n",
                i, length(members), paste(members, collapse = ", ")))
  }
  invisible(cl)
}

cat("\n=== Running log2-FC Mfuzz (STANDARDISED, Crea-norm) ===\n")
run_fc_mfuzz_std(mats_std$FC_GA_std_Crea,       "FC_GA_std_Crea",       k = 3)
run_fc_mfuzz_std(mats_std$FC_BPD_std_Crea,      "FC_BPD_std_Crea",      k = 3)
run_fc_mfuzz_std(mats_std$FC_Combined_std_Crea, "FC_Combined_std_Crea", k = 4)

# ---- 7. Pretty plots (standardised version) ---------------------------------
cat("\n=== Building pretty plots ===\n")

mfuzz_palette <- c("#00CCCC", "#33CC33", "#9966CC", "#CC0066")
cluster_centroid_colour <- "#222222"

load_membership <- function(csv_path, x_levels) {
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8")
  cluster_cols <- grep("^Cluster_", names(df), value = TRUE)
  k <- length(cluster_cols)
  value_cols <- setdiff(names(df), c("Metabolite", "HardCluster", cluster_cols))
  df$MaxMemb <- apply(df[, cluster_cols, drop = FALSE], 1, max)
  val_mat <- as.matrix(df[, value_cols])
  row_means <- rowMeans(val_mat, na.rm = TRUE)
  row_sds   <- apply(val_mat, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0] <- 1
  val_mat_std <- (val_mat - row_means) / row_sds
  df[, value_cols] <- as.data.frame(val_mat_std)
  long <- df %>%
    select(Metabolite, HardCluster, MaxMemb, all_of(value_cols)) %>%
    pivot_longer(cols = all_of(value_cols), names_to = "Group", values_to = "Value") %>%
    mutate(Group = factor(Group, levels = x_levels))
  list(long = long, k = k, value_cols = value_cols)
}

top_per_cluster <- function(csv_path, n = 8) {
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8")
  cluster_cols <- grep("^Cluster_", names(df), value = TRUE)
  out <- list()
  for (i in seq_along(cluster_cols)) {
    sub <- df[df$HardCluster == i, ]
    sub <- sub[order(-sub[[cluster_cols[i]]]), ]
    out[[i]] <- head(sub$Metabolite, n)
  }
  out
}

plot_styleA <- function(long, k, title) {
  centroids <- long %>% group_by(HardCluster, Group) %>%
    summarise(Value = mean(Value), .groups = "drop")
  ggplot(long, aes(x = Group, y = Value, group = Metabolite)) +
    geom_line(aes(colour = MaxMemb), linewidth = 0.7, alpha = 0.85) +
    geom_line(data = centroids, aes(x = Group, y = Value, group = HardCluster),
              colour = "black", linewidth = 1.4, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    facet_wrap(~ HardCluster, ncol = min(k, 3),
               labeller = labeller(HardCluster = function(x) paste("Cluster", x))) +
    scale_colour_gradientn(colours = mfuzz_palette, limits = c(0, 1), name = "Membership") +
    labs(title = title, x = "Group", y = "Standardised log2(2Y / 6M)") +
    theme_minimal(base_size = 11) +
    theme(aspect.ratio = 1,
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right", plot.title = element_text(face = "bold"))
}

plot_styleB <- function(long, k, top_lists, title) {
  centroids <- long %>% group_by(HardCluster, Group) %>%
    summarise(Mean = mean(Value), SD = sd(Value), .groups = "drop")
  p_traj <- ggplot(long, aes(x = Group, y = Value, group = Metabolite)) +
    geom_line(colour = "grey75", linewidth = 0.4, alpha = 0.5) +
    geom_ribbon(data = centroids, aes(x = Group, ymin = Mean - SD, ymax = Mean + SD,
                group = HardCluster), fill = "#E07B6A", alpha = 0.18, inherit.aes = FALSE) +
    geom_line(data = centroids, aes(x = Group, y = Mean, group = HardCluster),
              colour = cluster_centroid_colour, linewidth = 1.4, inherit.aes = FALSE) +
    geom_point(data = centroids, aes(x = Group, y = Mean),
               colour = cluster_centroid_colour, size = 2.4, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    facet_wrap(~ HardCluster, ncol = min(k, 3),
               labeller = labeller(HardCluster = function(x) paste("Cluster", x))) +
    labs(title = title, x = "Group", y = "Standardised log2(2Y / 6M)") +
    theme_minimal(base_size = 11) +
    theme(aspect.ratio = 1,
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"))
  list_df <- do.call(rbind, lapply(seq_along(top_lists), function(i) {
    data.frame(Cluster = paste("Cluster", i), Rank = seq_along(top_lists[[i]]),
               Met = top_lists[[i]], stringsAsFactors = FALSE)
  }))
  list_df$Cluster <- factor(list_df$Cluster, levels = unique(list_df$Cluster))
  p_list <- ggplot(list_df, aes(x = 1, y = -Rank, label = Met)) +
    geom_text(size = 3.2, hjust = 0) +
    facet_wrap(~ Cluster, ncol = min(k, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.95, 5)) +
    labs(x = NULL, y = NULL, caption = "Top metabolites by cluster membership") +
    theme_void(base_size = 10) +
    theme(strip.text = element_blank(),
          plot.caption = element_text(hjust = 0, face = "italic", size = 9))
  p_traj / p_list + plot_layout(heights = c(3.2, 1.6))
}

pretty_configs <- list(
  list(tag = "FC_GA_std_Crea",
       file = "mfuzz_FC_GA_std_Crea_membership.csv",
       x = c("GA_>=37wk", "GA_28-32wk", "GA_<28wk"),
       title = "log2-FC clusters (Crea-norm) by GA"),
  list(tag = "FC_BPD_std_Crea",
       file = "mfuzz_FC_BPD_std_Crea_membership.csv",
       x = c("BPD_HC", "BPD_No+Mild", "BPD_M+S"),
       title = "log2-FC clusters (Crea-norm) by BPD"),
  list(tag = "FC_Combined_std_Crea",
       file = "mfuzz_FC_Combined_std_Crea_membership.csv",
       x = c("GA_>=37wk", "GA_28-32wk", "GA_<28wk", "BPD_HC", "BPD_No+Mild", "BPD_M+S"),
       title = "log2-FC clusters (Crea-norm) combined GA + BPD")
)

for (cfg in pretty_configs) {
  loaded <- load_membership(cfg$file, cfg$x)
  top_lists <- top_per_cluster(cfg$file, n = 8)
  pA <- plot_styleA(loaded$long, loaded$k, cfg$title)
  pB <- plot_styleB(loaded$long, loaded$k, top_lists, cfg$title)
  nc <- min(loaded$k, 3); nr <- ceiling(loaded$k / nc)
  w_A <- nc * 4.5 + 2; h_A <- nr * 4.5 + 1.2
  w_B <- nc * 4.5 + 1; h_B <- nr * 4.5 + nr * 2.5 + 1
  ggsave(paste0("pretty_", cfg$tag, "_styleA.pdf"), plot = pA, width = w_A, height = h_A)
  ggsave(paste0("pretty_", cfg$tag, "_styleB.pdf"), plot = pB, width = w_B, height = h_B)
  cat("  -> pretty_", cfg$tag, "_styleA.pdf + _styleB.pdf\n", sep = "")
}

cat("\nAll done.\n")
