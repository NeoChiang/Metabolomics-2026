# =============================================================================
#  generate_table1_6M_2Y.R
#
#  Produces two demographic tables comparing BPD cohort characteristics at the
#  6-month and 2-year corrected-age follow-ups:
#
#      Table_1_GA_6M_2Y.docx   — grouped by GA3Group
#      Table_1_BPD_6M_2Y.docx  — grouped by BPD3Group
#
#  Cohort
#  ------
#      Each timepoint uses its OWN sample list independently:
#        6M  -> all cases in Blood_Urine sheet  (n ~ 139)
#        2Y  -> all cases in 収案名單 sheet     (n ~ 206)
#      Group sizes therefore differ between 6M and 2Y columns.
#
#  Data sources
#  ------------
#      6M  : 139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx
#            sheet "Blood_Urine".  GA3Group / BPD3Group from this sheet.
#      2Y  : 2Y BPD_206 20250721.xlsx , sheet "收案名單".
#            GA3Group / BPD3Group from this same sheet (verified 100 %
#            agreement with the 6M file on the 139 overlapping cases).
#
#  Time-invariant fields (sex, GA, birth body weight) in the 6M data are
#  backfilled from 収案名單 where the 6M Blood_Urine row is empty.
#
#  Statistical methods
#  -------------------
#      Continuous  : mean +/- SD , p-value = one-way ANOVA (three groups)
#      Categorical : n (%)       , p-value = chi-square, Fisher's exact when
#                                  any expected count < 5
#      P < 0.05 is rendered in bold in the output docx.
#
#  Dependencies (Windows-ready):
#      readxl , dplyr , officer , flextable
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(officer)
  library(flextable)
})

try(Sys.setlocale("LC_ALL", "C.UTF-8"), silent = TRUE)
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------ file paths
f6m <- "139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx"
f2y <- "2Y BPD_206 20250721.xlsx"
stopifnot(file.exists(f6m), file.exists(f2y))

# ------------------------------------------------------------------ load raw
raw6 <- read_excel(f6m, sheet = "Blood_Urine", .name_repair = "minimal")
raw6 <- raw6[!is.na(raw6$`case no.`), ]

rawS <- read_excel(f2y, sheet = "收案名單",   .name_repair = "minimal")
rawS <- rawS[!is.na(rawS$`case no.`), ]

# ------------------------------------------------------------------ 6M tidy
d6 <- data.frame(
  case_no    = as.character(raw6$`case no.`),
  GA3Group   = as.numeric(raw6$GA3Group),
  BPD3Group  = as.numeric(raw6$BPD3Group),
  sex_male   = as.numeric(raw6$`男1/女0`),
  GA_wk      = suppressWarnings(as.numeric(raw6$GA)),
  BWt_g      = suppressWarnings(as.numeric(raw6$BWt)),
  age_mo     = suppressWarnings(as.numeric(raw6$`右邊BW/BH/BMI 測量CA(月份)`)),
  Wt_kg      = suppressWarnings(as.numeric(raw6$`Wt (kg)`)),
  Ht_cm      = suppressWarnings(as.numeric(raw6$`Ht (cm)`)),
  BMI        = suppressWarnings(as.numeric(raw6[[32]])),
  BF_raw     = suppressWarnings(as.numeric(raw6[[54]])),
  Sepsis_raw = suppressWarnings(as.numeric(raw6[[41]])),
  check.names = FALSE
)
d6$BF_ge6 <- ifelse(is.na(d6$BF_raw), NA_integer_,
                    ifelse(d6$BF_raw == 2, 1L, 0L))
d6$Sepsis_ever <- ifelse(is.na(d6$Sepsis_raw), NA_integer_,
                         ifelse(d6$Sepsis_raw >= 1, 1L, 0L))
d6$BF_raw <- NULL; d6$Sepsis_raw <- NULL

# ------------------------------------------------------------------ 2Y tidy
sepsis_raw <- as.character(rawS$`Sepsis, ever`)
sepsis_num <- suppressWarnings(as.numeric(sepsis_raw))
sepsis_num[!is.na(sepsis_raw) & !sepsis_raw %in% c("0","1")] <- NA

ca_years <- suppressWarnings(as.numeric(rawS[[21]]))
ca_years[ca_years < 0 | ca_years > 5] <- NA

parse_ga_str <- function(s) {
  s <- as.character(s)
  m <- regmatches(s, regexec("^\\s*([0-9]+)\\s*\\+\\s*([0-9]*)\\s*$", s))
  sapply(seq_along(s), function(i) {
    mi <- m[[i]]
    if (length(mi) < 3 || mi[1] == "") {
      as.numeric(s[i])
    } else {
      wk <- as.numeric(mi[2])
      d  <- suppressWarnings(as.numeric(mi[3]))
      if (is.na(d)) d <- 0
      wk + d / 7
    }
  })
}

d2 <- data.frame(
  case_no    = as.character(rawS$`case no.`),
  GA3Group   = as.numeric(rawS$GA3Group),
  BPD3Group  = as.numeric(rawS$BPD3Group),
  sex_male   = as.numeric(rawS$sex),
  GA_wk      = parse_ga_str(rawS$GA),
  BWt_g      = suppressWarnings(as.numeric(rawS$BWt)),
  age_mo     = ca_years * 12,
  Wt_kg      = suppressWarnings(as.numeric(rawS$`2Y weight`)),
  Ht_cm      = suppressWarnings(as.numeric(rawS$`2Y height`)),
  BMI        = suppressWarnings(as.numeric(rawS[[32]])),
  BF_raw     = suppressWarnings(as.numeric(rawS[[grep("Breastfeeding", names(rawS))[1]]])),
  Sepsis_ever= sepsis_num,
  check.names = FALSE
)
d2$BF_ge6 <- ifelse(is.na(d2$BF_raw), NA_integer_,
                    ifelse(d2$BF_raw == 2, 1L, 0L))
d2$BF_raw <- NULL

# ------------------------------------------------------------------ cohort
# Each timepoint uses its own sample list — NO intersection filtering.
# Backfill 6M time-invariant fields (sex, GA, BWt) from 収案名單 for cases
# whose 6M Blood_Urine row had empty demographics.
m_6toS <- match(d6$case_no, d2$case_no)
coalesce2 <- function(a, b) ifelse(!is.na(a), a, b)
for (v in c("sex_male", "GA_wk", "BWt_g")) {
  d6[[v]] <- coalesce2(d6[[v]], d2[[v]][m_6toS])
}

cat("\n=============================================\n")
cat(" Sample sizes (each timepoint independent)\n")
cat("=============================================\n")
cat("  6M (Blood_Urine)   n =", nrow(d6), "\n")
cat("  2Y (收案名單)      n =", nrow(d2), "\n")
cat("  overlap (same case) =", length(intersect(d6$case_no, d2$case_no)), "\n")

cat("\n  6M GA3Group  0/1/2:", sum(d6$GA3Group==0,na.rm=TRUE), "/",
    sum(d6$GA3Group==1,na.rm=TRUE), "/", sum(d6$GA3Group==2,na.rm=TRUE), "\n")
cat("  6M BPD3Group 0/1/2:", sum(d6$BPD3Group==0,na.rm=TRUE), "/",
    sum(d6$BPD3Group==1,na.rm=TRUE), "/", sum(d6$BPD3Group==2,na.rm=TRUE), "\n")
cat("  2Y GA3Group  0/1/2:", sum(d2$GA3Group==0,na.rm=TRUE), "/",
    sum(d2$GA3Group==1,na.rm=TRUE), "/", sum(d2$GA3Group==2,na.rm=TRUE), "\n")
cat("  2Y BPD3Group 0/1/2:", sum(d2$BPD3Group==0,na.rm=TRUE), "/",
    sum(d2$BPD3Group==1,na.rm=TRUE), "/", sum(d2$BPD3Group==2,na.rm=TRUE), "\n")

cat("\n  Missing counts per variable:\n")
varlist <- c("sex_male","GA_wk","BWt_g","age_mo","Wt_kg","Ht_cm","BMI",
             "BF_ge6","Sepsis_ever")
for (v in varlist) {
  cat(sprintf("    %-14s  6M missing = %3d / %3d   2Y missing = %3d / %3d\n",
              v, sum(is.na(d6[[v]])), nrow(d6),
              sum(is.na(d2[[v]])), nrow(d2)))
}

# =============================================================================
#  Helpers
# =============================================================================

fmt_p <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}

fmt_meanSD <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("-")
  sprintf("%s ± %s",
          formatC(mean(x), format="f", digits=2),
          formatC(sd(x),   format="f", digits=2))
}

fmt_npct <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  n <- sum(x == 1)
  sprintf("%d (%s%%)", n, formatC(100 * n / length(x), format="f", digits=1))
}

p_continuous <- function(value, group) {
  ok <- is.finite(value) & !is.na(group)
  value <- value[ok]; group <- factor(group[ok])
  if (length(unique(group)) < 2 || length(value) < 3) return(NA_real_)
  tryCatch(summary(aov(value ~ group))[[1]][["Pr(>F)"]][1],
           error = function(e) NA_real_)
}

p_categorical <- function(value, group) {
  ok <- !is.na(value) & !is.na(group)
  tab <- table(group[ok], value[ok])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  expc <- tryCatch(suppressWarnings(chisq.test(tab)$expected),
                   error = function(e) NULL)
  if (!is.null(expc) && any(expc < 5)) {
    tryCatch(fisher.test(tab, workspace = 2e7)$p.value,
             error = function(e) NA_real_)
  } else {
    tryCatch(suppressWarnings(chisq.test(tab)$p.value),
             error = function(e) NA_real_)
  }
}

row_spec <- list(
  list(label = "Sex (male, n(%))",                    var = "sex_male",    type = "cat"),
  list(label = "Gestational age (wk)",                var = "GA_wk",       type = "cont"),
  list(label = "Birth body weight (g)",               var = "BWt_g",       type = "cont"),
  list(label = "Age corrected (month)",               var = "age_mo",      type = "cont"),
  list(label = "Body weight (kg)",                    var = "Wt_kg",       type = "cont"),
  list(label = "Body height (cm)",                    var = "Ht_cm",       type = "cont"),
  list(label = "BMI (kg/m²)",                    var = "BMI",         type = "cont"),
  list(label = "Breastfeeding ≥6 months, n(%)",  var = "BF_ge6",      type = "cat"),
  list(label = "Sepsis ever, n(%)",                   var = "Sepsis_ever", type = "cat")
)

# =============================================================================
#  Build one docx given a grouping variable + labels
# =============================================================================

build_table_docx <- function(group_var, group_levels, group_labels,
                             caption, out_file) {

  data_by <- function(df, lvl) df[df[[group_var]] == lvl, , drop = FALSE]
  K <- length(group_levels)

  # Column layout: Char | 6M_g1..6M_gK | P(6M) | 2Y_g1..2Y_gK | P(2Y)
  ncols <- 1 + K + 1 + K + 1   # = 2K + 3
  body      <- matrix("",    nrow = length(row_spec), ncol = ncols)
  bold_flag <- matrix(FALSE, nrow = length(row_spec), ncol = ncols)

  col_6M_start <- 2
  col_P6       <- 2 + K          # P(6M) right after 6M groups
  col_2Y_start <- 2 + K + 1
  col_P2       <- ncols           # P(2Y) right after 2Y groups

  for (i in seq_along(row_spec)) {
    spec <- row_spec[[i]]
    body[i, 1] <- spec$label
    for (k in seq_along(group_levels)) {
      lvl <- group_levels[k]
      s6 <- data_by(d6, lvl)[[spec$var]]
      s2 <- data_by(d2, lvl)[[spec$var]]
      c6 <- col_6M_start + (k - 1)
      c2 <- col_2Y_start + (k - 1)
      if (spec$type == "cont") {
        body[i, c6] <- fmt_meanSD(s6)
        body[i, c2] <- fmt_meanSD(s2)
      } else {
        body[i, c6] <- fmt_npct(s6)
        body[i, c2] <- fmt_npct(s2)
      }
    }
    g6 <- d6[[group_var]]; g2 <- d2[[group_var]]
    v6 <- d6[[spec$var]];  v2 <- d2[[spec$var]]
    if (spec$type == "cont") {
      p6 <- p_continuous(v6, g6); p2 <- p_continuous(v2, g2)
    } else {
      p6 <- p_categorical(v6, g6); p2 <- p_categorical(v2, g2)
    }
    body[i, col_P6] <- fmt_p(p6)
    body[i, col_P2] <- fmt_p(p2)
    if (!is.na(p6) && p6 < 0.05) bold_flag[i, col_P6] <- TRUE
    if (!is.na(p2) && p2 < 0.05) bold_flag[i, col_P2] <- TRUE
  }

  n6 <- sapply(group_levels, function(lvl) sum(d6[[group_var]] == lvl, na.rm=TRUE))
  n2 <- sapply(group_levels, function(lvl) sum(d2[[group_var]] == lvl, na.rm=TRUE))

  # Internal col_keys (must be unique)
  col_keys <- c("char",
                sprintf("c6M_%d", seq_len(K)), "pval6",
                sprintf("c2Y_%d", seq_len(K)), "pval2")

  # Display labels for level-2 (sub) header
  sub_hdr <- c("Characteristics",
               sprintf("%s (n=%d)", group_labels, n6), "P",
               sprintf("%s (n=%d)", group_labels, n2), "P")

  # Display labels for level-1 (top) header — 6M spans K+1 (groups + P), same for 2Y
  top_hdr <- c("Characteristics",
               rep("6M", K + 1),
               rep("2Y", K + 1))

  df <- as.data.frame(body, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- col_keys

  # ---------- flextable ----------
  ft <- flextable(df, col_keys = col_keys)

  hdr_df <- data.frame(
    col_keys = col_keys,
    top      = top_hdr,
    sub      = sub_hdr,
    stringsAsFactors = FALSE
  )
  ft <- set_header_df(ft, mapping = hdr_df, key = "col_keys")
  ft <- merge_h(ft, part = "header")
  ft <- merge_v(ft, j = 1, part = "header")

  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "body")
  ft <- bold(ft, part = "header")
  ft <- fontsize(ft, size = 9,  part = "body")
  ft <- fontsize(ft, size = 10, part = "header")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2, part = "all")
  ft <- border_outer(ft, border = fp_border(color = "black", width = 1))
  ft <- border_inner_h(ft, border = fp_border(color = "grey60", width = 0.5))
  ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1), part = "header")

  sep <- fp_border(color = "black", width = 0.8)
  ft <- vline(ft, j = 1,          border = sep, part = "all")   # after Characteristics
  ft <- vline(ft, j = col_P6,     border = sep, part = "all")   # after P(6M), before 2Y
  ft <- vline(ft, j = col_P2 - 1, border = sep, part = "body")  # before P(2Y)
  ft <- vline(ft, j = col_P6 - 1, border = sep, part = "body")  # before P(6M)

  for (i in seq_len(nrow(bold_flag))) {
    for (j in seq_len(ncol(bold_flag))) {
      if (bold_flag[i, j]) ft <- bold(ft, i = i, j = j, part = "body")
    }
  }

  ft <- set_table_properties(ft, layout = "autofit")

  # ---------- docx ----------
  footnote <- paste(
    "Data shown are mean ± SD or number (%) of patients as appropriate.",
    "GA, gestational age; BPD, bronchopulmonary dysplasia; HC, healthy",
    "controls; No+Mild BPD, no and mild BPD; M+S BPD, moderate and severe BPD;",
    "wk, week; g, gram; cm, centimeter; BMI, body mass index.",
    "All P-values < 0.05, which is in bold, are significant."
  )

  doc <- read_docx()
  doc <- body_add_par(doc, caption, style = "Normal")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")
  doc <- body_add_par(doc, footnote, style = "Normal")
  print(doc, target = out_file)

  list(n6 = n6, n2 = n2)
}

# =============================================================================
#  Generate
# =============================================================================

cat("\n=============================================\n")
cat(" Generating docx files\n")
cat("=============================================\n")

out1 <- "Table_1_GA_6M_2Y.docx"
info1 <- build_table_docx(
  group_var    = "GA3Group",
  group_levels = c(0, 1, 2),
  group_labels = c("≥37 weeks", "28-32 weeks", "<28 weeks"),
  caption      = paste("Table 1. Comparisons of the demographic",
                       "characteristics among full-term and preterm",
                       "infants less than 32 weeks of gestational age",
                       "categorized by different GA at the corrected",
                       "ages of 6 months and 2 years."),
  out_file     = out1)
cat(" [GA ] 6M group n:", paste(info1$n6, collapse="/"),
    "  2Y group n:", paste(info1$n2, collapse="/"), "\n")
cat(" Wrote:", normalizePath(out1, mustWork = TRUE), "\n")

out2 <- "Table_1_BPD_6M_2Y.docx"
info2 <- build_table_docx(
  group_var    = "BPD3Group",
  group_levels = c(0, 1, 2),
  group_labels = c("HC", "No+Mild BPD", "M+S BPD"),
  caption      = paste("Table 1. Comparisons of the demographic",
                       "characteristics among full-term and preterm",
                       "infants less than 32 weeks of gestational age",
                       "categorized by BPD severity at the corrected",
                       "ages of 6 months and 2 years."),
  out_file     = out2)
cat(" [BPD] 6M group n:", paste(info2$n6, collapse="/"),
    "  2Y group n:", paste(info2$n2, collapse="/"), "\n")
cat(" Wrote:", normalizePath(out2, mustWork = TRUE), "\n")

cat("\nDone.\n")
