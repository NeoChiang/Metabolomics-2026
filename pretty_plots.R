# ============================================================================
# Pretty plotting of Mfuzz FC-standardised clustering results
# ----------------------------------------------------------------------------
# Reads:  mfuzz_FC_GA_std_membership.csv
#         mfuzz_FC_BPD_std_membership.csv
#         mfuzz_FC_Combined_std_membership.csv
# Writes: pretty_<tag>_styleA.pdf  -- Mfuzz-style with metabolite labels
#         pretty_<tag>_styleB.pdf  -- centroid + metabolite list on side
# ============================================================================

pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)
})

# ---- config ---------------------------------------------------------------
# Mfuzz-style colour ramp: low membership = cyan-ish, high = magenta/red
mfuzz_palette <- c("#00CCCC", "#33CC33", "#9966CC", "#CC0066")

# Cluster colours for centroid (Style B): use heatmap-friendly tones
cluster_centroid_colour <- "#222222"

# ---- helper ---------------------------------------------------------------
# Convert wide membership/value CSV to long format
load_membership <- function(csv_path, x_levels) {
  df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8")
  # value cols = everything after Cluster_* cols
  cluster_cols <- grep("^Cluster_", names(df), value = TRUE)
  k <- length(cluster_cols)
  value_cols <- setdiff(names(df),
                        c("Metabolite", "HardCluster", cluster_cols))

  # max membership per metabolite (used as line alpha)
  df$MaxMemb <- apply(df[, cluster_cols, drop = FALSE], 1, max)

  # ── Row-wise standardisation (same as Mfuzz standardise()) ──
  # Each metabolite: (value - row_mean) / row_sd
  val_mat <- as.matrix(df[, value_cols])
  row_means <- rowMeans(val_mat, na.rm = TRUE)
  row_sds   <- apply(val_mat, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0] <- 1  # prevent division by zero
  val_mat_std <- (val_mat - row_means) / row_sds
  df[, value_cols] <- as.data.frame(val_mat_std)

  long <- df %>%
    select(Metabolite, HardCluster, MaxMemb, all_of(value_cols)) %>%
    pivot_longer(cols = all_of(value_cols),
                 names_to = "Group", values_to = "Value") %>%
    mutate(Group = factor(Group, levels = x_levels))

  list(long = long, k = k, value_cols = value_cols)
}

# Top-N metabolites per cluster (by membership in their own cluster)
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

# ---- STYLE A: Mfuzz-style with metabolite labels --------------------------
plot_styleA <- function(long, k, title, x_label = "Group",
                        y_label = "Standardised log2(2Y / 6M)") {
  centroids <- long %>%
    group_by(HardCluster, Group) %>%
    summarise(Value = mean(Value), .groups = "drop")

  # Pick top-3 metabolites per cluster to label at the rightmost x point
  labels_df <- long %>%
    filter(Group == levels(Group)[length(levels(Group))]) %>%
    group_by(HardCluster) %>%
    slice_max(MaxMemb, n = 3) %>%
    ungroup()

  ggplot(long, aes(x = Group, y = Value, group = Metabolite)) +
    geom_line(aes(colour = MaxMemb), linewidth = 0.7, alpha = 0.85) +
    geom_line(data = centroids,
              aes(x = Group, y = Value, group = HardCluster),
              colour = "black", linewidth = 1.4, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    facet_wrap(~ HardCluster, ncol = min(k, 3),
               labeller = labeller(HardCluster = function(x) paste("Cluster", x))) +
    scale_colour_gradientn(colours = mfuzz_palette,
                           limits = c(0, 1), name = "Membership") +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right",
          plot.title = element_text(face = "bold"))
}

# ---- STYLE B: Centroid + metabolite list on side --------------------------
plot_styleB <- function(long, k, top_lists, title, x_label = "Group",
                        y_label = "Standardised log2(2Y / 6M)") {
  centroids <- long %>%
    group_by(HardCluster, Group) %>%
    summarise(Mean = mean(Value),
              SD   = sd(Value),
              .groups = "drop")

  # Build the trajectory panel
  p_traj <- ggplot(long, aes(x = Group, y = Value, group = Metabolite)) +
    geom_line(colour = "grey75", linewidth = 0.4, alpha = 0.5) +
    geom_ribbon(data = centroids,
                aes(x = Group, ymin = Mean - SD, ymax = Mean + SD,
                    group = HardCluster),
                fill = "#E07B6A", alpha = 0.18, inherit.aes = FALSE) +
    geom_line(data = centroids,
              aes(x = Group, y = Mean, group = HardCluster),
              colour = cluster_centroid_colour, linewidth = 1.4,
              inherit.aes = FALSE) +
    geom_point(data = centroids,
               aes(x = Group, y = Mean),
               colour = cluster_centroid_colour, size = 2.4,
               inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    facet_wrap(~ HardCluster, ncol = min(k, 3),
               labeller = labeller(HardCluster = function(x) paste("Cluster", x))) +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"))

  # Build the metabolite-list table panel
  list_df <- do.call(rbind, lapply(seq_along(top_lists), function(i) {
    data.frame(Cluster = paste("Cluster", i),
               Rank    = seq_along(top_lists[[i]]),
               Met     = top_lists[[i]],
               stringsAsFactors = FALSE)
  }))
  list_df$Cluster <- factor(list_df$Cluster, levels = unique(list_df$Cluster))

  p_list <- ggplot(list_df,
                   aes(x = 1, y = -Rank, label = Met)) +
    geom_text(size = 3.2, hjust = 0) +
    facet_wrap(~ Cluster, ncol = min(k, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.95, 5)) +
    labs(x = NULL, y = NULL,
         caption = "Top metabolites by cluster membership") +
    theme_void(base_size = 10) +
    theme(strip.text = element_blank(),
          plot.caption = element_text(hjust = 0, face = "italic", size = 9))

  p_traj / p_list + plot_layout(heights = c(3.2, 1.6))
}

# ---- run for each tag -----------------------------------------------------
configs <- list(
  list(tag = "FC_GA_std",
       file = "mfuzz_FC_GA_std_membership.csv",
       x = c("GA_>=37wk", "GA_28-32wk", "GA_<28wk"),
       title = "log2-FC clusters stratified by gestational age"),
  list(tag = "FC_BPD_std",
       file = "mfuzz_FC_BPD_std_membership.csv",
       x = c("BPD_HC", "BPD_No+Mild", "BPD_M+S"),
       title = "log2-FC clusters stratified by BPD severity"),
  list(tag = "FC_Combined_std",
       file = "mfuzz_FC_Combined_std_membership.csv",
       x = c("GA_>=37wk", "GA_28-32wk", "GA_<28wk",
             "BPD_HC", "BPD_No+Mild", "BPD_M+S"),
       title = "log2-FC clusters: combined GA + BPD stratification")
)

for (cfg in configs) {
  cat("\n=== Building plots for", cfg$tag, "===\n")
  loaded <- load_membership(cfg$file, cfg$x)
  top_lists <- top_per_cluster(cfg$file, n = 8)

  pA <- plot_styleA(loaded$long, loaded$k, cfg$title)
  pB <- plot_styleB(loaded$long, loaded$k, top_lists, cfg$title)

  ggsave(paste0("pretty_", cfg$tag, "_styleA.pdf"),
         plot = pA, width = 11, height = 5.5)
  ggsave(paste0("pretty_", cfg$tag, "_styleB.pdf"),
         plot = pB, width = 11, height = 7.5)
  cat("  -> pretty_", cfg$tag, "_styleA.pdf  (Mfuzz-style with labels)\n", sep = "")
  cat("  -> pretty_", cfg$tag, "_styleB.pdf  (centroid + metabolite list)\n", sep = "")
}

cat("\nAll done. 6 PDFs created.\n")
