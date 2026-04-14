# =============================================================================
#  generate_table1_6M_2Y.R
#
#  Produces two demographic tables comparing BPD cohort characteristics at the
#  6-month and 2-year corrected-age follow-ups:
#
#      Table_1_GA_6M_2Y.docx   — grouped by GA3Group
#                                (0 = >=37 weeks, 1 = 28-32 weeks, 2 = <28 weeks)
#      Table_1_BPD_6M_2Y.docx  — grouped by BPD3Group
#                                (0 = HC, 1 = No+Mild BPD, 2 = Moderate+Severe)
#
#  Cohort
#  ------
#      Only subjects present in BOTH the 6M (Blood_Urine sheet) and the 2Y
#      (collected in 2Y file) datasets are retained.  Grouping (GA3Group /
#      BPD3Group) is taken from the 6M file.  Expected n = 140.
#
#  Data sources
#  ------------
#      6M  : 139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx , sheet
#            "Blood_Urine".
#      2Y  : 2Y BPD_206 20250721.xlsx , sheet "收案名單".
#            The sheet "收案名單" is used instead of "2Y問卷" because
#            "2Y問卷" is sparsely populated (~36/140 cases) and does not
#            contain a sex column, whereas "收案名單" contains cleanly
#            curated sex / 2Y height / 2Y weight / BMI / Breastfeeding>=6mo /
#            Sepsis,ever and its GA3Group / BPD3Group agree 100% with the
#            6M file on the 139 overlapping cases.  "收案名單" effectively
#            IS the cleaned 2Y demographic table.
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
#
#  Run from the folder containing this script and the two xlsx files.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(officer)
  library(flextable)
})

# UTF-8 locale (harmless on Windows where UTF-8 is default R 4.2+)
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
  sex_male   = as.numeric(raw6$`男1/女0`),            # 1 = male
  GA_wk      = suppressWarnings(as.numeric(raw6$GA)),
  BWt_g      = suppressWarnings(as.numeric(raw6$BWt)),
  age_mo     = suppressWarnings(as.numeric(raw6$`右邊BW/BH/BMI 測量CA(月份)`)),
  Wt_kg      = suppressWarnings(as.numeric(raw6$`Wt (kg)`)),
  Ht_cm      = suppressWarnings(as.numeric(raw6$`Ht (cm)`)),
  BMI        = suppressWarnings(as.numeric(raw6[[32]])),   # second "BMI" col
  BF_raw     = suppressWarnings(as.numeric(raw6[[54]])),   # 母奶(...>6mo=2)
  Sepsis_raw = suppressWarnings(as.numeric(raw6[[41]])),   # Ever Sepsis/ABC/Culture
  check.names = FALSE
)
# Breastfeeding >=6 months: coded 0 or 2 -> recode to 0/1
d6$BF_ge6 <- ifelse(is.na(d6$BF_raw), NA_integer_,
                    ifelse(d6$BF_raw == 2, 1L, 0L))
d6$Sepsis_ever <- ifelse(is.na(d6$Sepsis_raw), NA_integer_,
                         ifelse(d6$Sepsis_raw >= 1, 1L, 0L))
d6$BF_raw <- NULL; d6$Sepsis_raw <- NULL

# ------------------------------------------------------------------ 2Y tidy
#  "收案名單" column positions of interest:
#     1  case no.       5  sex     9  GA (text "wk+day")   10 BWt
#    20  blood date_2Y  21  corrected age at blood_2Y  (years)
#    28  GA3Group       29  BPD3Group
#    30  2Y height      31  2Y weight       32  BMI
#    38  Breastfeeding>=6 months             40  Sepsis, ever
sepsis_raw <- as.character(rawS$`Sepsis, ever`)
sepsis_num <- suppressWarnings(as.numeric(sepsis_raw))
# Non-clean entries ("*", "1-0") -> missing
sepsis_num[!is.na(sepsis_raw) & !sepsis_raw %in% c("0","1")] <- NA

ca_years <- suppressWarnings(as.numeric(rawS[[21]]))
# Sentinel negatives (~-122) mean "not collected" -> drop
ca_years[ca_years < 0 | ca_years > 5] <- NA

# 収案名單's GA is text like "39+3" (weeks+days).  Parse -> decimal weeks.
parse_ga_str <- function(s) {
  s <- as.character(s)
  m <- regmatches(s, regexec("^\\s*([0-9]+)\\s*\\+\\s*([0-9]*)\\s*$", s))
  sapply(seq_along(s), function(i) {
    mi <- m[[i]]
    if (length(mi) < 3 || mi[1] == "") {
      as.numeric(s[i])    # already plain numeric? else NA
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
  GA3Group_S = as.numeric(rawS$GA3Group),
  BPD3Group_S= as.numeric(rawS$BPD3Group),
  sex_male   = as.numeric(rawS$sex),                # 1 = male (verified vs 6M)
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
overlap_cases <- intersect(d6$case_no, d2$case_no)
cat("\n=============================================\n")
cat(" Cohort: subjects in BOTH 6M and 2Y\n")
cat("=============================================\n")
cat("  n(6M file)           =", length(unique(d6$case_no)), "\n")
cat("  n(2Y file 收案名單)  =", length(unique(d2$case_no)), "\n")
cat("  n(intersection)      =", length(overlap_cases), "\n")

d6 <- d6[d6$case_no %in% overlap_cases, ]
d2 <- d2[d2$case_no %in% overlap_cases, ]

# apply 6M grouping universally (already true for overlap cases, but be explicit)
grp6 <- d6[, c("case_no","GA3Group","BPD3Group")]
d2 <- merge(d2[, setdiff(names(d2), c("GA3Group_S","BPD3Group_S"))],
            grp6, by = "case_no", all.x = TRUE)

# -------------------------------------------------------------------
# Time-invariant demographic attributes (sex, GA at birth, birth body weight)
# do not change between the 6M and 2Y visits.  A number of rows in the 6M
# Blood_Urine sheet lack these values (empty demographic block for cases
# whose 6M blood/urine was not collected), but they are fully populated in
# the 2Y 収案名單 sheet.  We coalesce the two sources so each subject has a
# single authoritative value; that value is then shown in BOTH the 6M and
# the 2Y columns (which are therefore identical for these three rows —
# as they must be, since sex / GA / birth weight cannot change with age).
# -------------------------------------------------------------------
m_6to2 <- match(d6$case_no, d2$case_no)
m_2to6 <- match(d2$case_no, d6$case_no)
coalesce2 <- function(a, b) ifelse(!is.na(a), a, b)
for (v in c("sex_male","GA_wk","BWt_g")) {
  d6_v <- coalesce2(d6[[v]], d2[[v]][m_6to2])
  d2_v <- coalesce2(d2[[v]], d6[[v]][m_2to6])
  d6[[v]] <- d6_v
  d2[[v]] <- d2_v
}
cat("\n  Coalesced time-invariant fields (sex, GA, BWt): 6M uses 6M's value,\n",
    "  backfilled from 2Y where 6M was missing (and vice versa).\n", sep = "")

# Confirm the two files agree (sanity)
chk <- merge(
  d6[, c("case_no","GA3Group","BPD3Group")],
  merge(d2[, "case_no", drop=FALSE],
        rawS[, c("case no.","GA3Group","BPD3Group")] %>%
          rename(case_no = `case no.`, GA3_S = GA3Group, BPD3_S = BPD3Group),
        by = "case_no"),
  by = "case_no")
cat("\n  Group agreement check (6M vs 收案名單 within cohort):\n")
cat("    GA3Group mismatches :", sum(chk$GA3Group  != chk$GA3_S,  na.rm=TRUE), "\n")
cat("    BPD3Group mismatches:", sum(chk$BPD3Group != chk$BPD3_S, na.rm=TRUE), "\n")

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
  sprintf("%s \u00b1 %s",
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
  out <- tryCatch(summary(aov(value ~ group))[[1]][["Pr(>F)"]][1],
                  error = function(e) NA_real_)
  out
}

p_categorical <- function(value, group) {
  ok <- !is.na(value) & !is.na(group)
  tab <- table(group[ok], value[ok])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  expc <- tryCatch(suppressWarnings(chisq.test(tab)$expected),
                   error = function(e) NULL)
  if (!is.null(expc) && any(expc < 5)) {
    out <- tryCatch(fisher.test(tab, workspace = 2e7)$p.value,
                    error = function(e) NA_real_)
  } else {
    out <- tryCatch(suppressWarnings(chisq.test(tab)$p.value),
                    error = function(e) NA_real_)
  }
  out
}

# Row specification: (display label, variable name, type)
#   type: "cont" or "cat"
row_spec <- list(
  list(label = "Sex (male, n(%))",          var = "sex_male",    type = "cat"),
  list(label = "Gestational age (wk)",      var = "GA_wk",       type = "cont"),
  list(label = "Birth body weight (g)",     var = "BWt_g",       type = "cont"),
  list(label = "Age corrected (month)",     var = "age_mo",      type = "cont"),
  list(label = "Body weight (kg)",          var = "Wt_kg",       type = "cont"),
  list(label = "Body height (cm)",          var = "Ht_cm",       type = "cont"),
  list(label = "BMI (kg/m\u00b2)",          var = "BMI",         type = "cont"),
  list(label = "Breastfeeding \u22656 months, n(%)",
                                            var = "BF_ge6",      type = "cat"),
  list(label = "Sepsis ever, n(%)",         var = "Sepsis_ever", type = "cat")
)

# =============================================================================
#  Build one docx given a grouping variable + labels
# =============================================================================

build_table_docx <- function(group_var, group_levels, group_labels,
                             caption, out_file) {

  # split data by group and timepoint
  data_by <- function(df, lvl) df[df[[group_var]] == lvl, , drop = FALSE]

  # Column layout:
  #   1                                           Characteristics
  #   2..(1+K)                                    6M : one sub-col per group
  #   (2+K)..(1+2K)                               2Y : one sub-col per group
  #   (2+2K)                                      P (6M)
  #   (3+2K)                                      P (2Y)
  K     <- length(group_levels)
  ncols <- 1 + 2 * K + 2
  body      <- matrix("",    nrow = length(row_spec), ncol = ncols)
  bold_flag <- matrix(FALSE, nrow = length(row_spec), ncol = ncols)

  col_6M_start <- 2
  col_2Y_start <- 2 + K
  col_P6       <- ncols - 1
  col_P2       <- ncols

  for (i in seq_along(row_spec)) {
    spec <- row_spec[[i]]
    body[i, 1] <- spec$label
    # per group summaries — 6M columns first, then 2Y columns
    for (k in seq_along(group_levels)) {
      lvl <- group_levels[k]
      s6 <- data_by(d6, lvl)[[spec$var]]
      s2 <- data_by(d2, lvl)[[spec$var]]
      col6 <- col_6M_start + (k - 1)
      col2 <- col_2Y_start + (k - 1)
      if (spec$type == "cont") {
        body[i, col6] <- fmt_meanSD(s6)
        body[i, col2] <- fmt_meanSD(s2)
      } else {
        body[i, col6] <- fmt_npct(s6)
        body[i, col2] <- fmt_npct(s2)
      }
    }
    # p-values across the three groups
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

  # group n at each timepoint
  n6 <- sapply(group_levels, function(lvl) sum(d6[[group_var]] == lvl, na.rm = TRUE))
  n2 <- sapply(group_levels, function(lvl) sum(d2[[group_var]] == lvl, na.rm = TRUE))

  # Level-2 header (subgroup labels, with n)
  sub_hdr <- c("Characteristics",
               sprintf("%s (n=%d)", group_labels, n6),
               sprintf("%s (n=%d)", group_labels, n2),
               "P (6M)", "P (2Y)")

  # Unique internal col_keys for flextable (display labels set via set_header_df)
  col_keys <- c("char",
                sprintf("c6M_%d", seq_len(K)),
                sprintf("c2Y_%d", seq_len(K)),
                "pval6", "pval2")

  df <- as.data.frame(body, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- col_keys

  # ---------- flextable ----------
  ft <- flextable(df, col_keys = col_keys)

  # Two-row header:
  #   row 1 (top): Characteristics | 6M (span K) | 2Y (span K) | P-value | P-value
  #   row 2      :                 | subgroup labels x K       | subgroup x K | P (6M) | P (2Y)
  hdr_df <- data.frame(
    col_keys = col_keys,
    top = c("Characteristics",
            rep("6M", K), rep("2Y", K),
            "P-value", "P-value"),
    sub = sub_hdr,
    stringsAsFactors = FALSE
  )
  ft <- set_header_df(ft, mapping = hdr_df, key = "col_keys")
  ft <- merge_h(ft, part = "header")
  ft <- merge_v(ft, j = 1, part = "header")   # first column label spans both rows

  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "body")
  ft <- bold(ft, part = "header")
  ft <- fontsize(ft, size = 9,  part = "body")
  ft <- fontsize(ft, size = 10, part = "header")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2, part = "all")
  ft <- border_outer(ft, border = fp_border(color = "black", width = 1))
  ft <- border_inner_h(ft, border = fp_border(color = "grey60", width = 0.5))
  ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1), part = "header")
  # Vertical separator between 6M block, 2Y block, and P-value block
  sep <- fp_border(color = "black", width = 0.8)
  ft <- vline(ft, j = col_6M_start - 1, border = sep, part = "all")   # after Characteristics
  ft <- vline(ft, j = col_2Y_start - 1, border = sep, part = "all")   # between 6M and 2Y
  ft <- vline(ft, j = col_2Y_start + K - 1, border = sep, part = "all")  # before P-value

  # Bold cells where p < 0.05
  for (i in seq_len(nrow(bold_flag))) {
    for (j in seq_len(ncol(bold_flag))) {
      if (bold_flag[i, j]) {
        ft <- bold(ft, i = i, j = j, part = "body")
      }
    }
  }

  ft <- set_table_properties(ft, layout = "autofit")

  # ---------- docx ----------
  footnote <- paste(
    "Data shown are mean \u00b1 SD or number (%) of patients as appropriate.",
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
#  Pre-flight diagnostics
# =============================================================================

cat("\n=============================================\n")
cat(" Pre-flight diagnostics (cohort n = ", nrow(d6), ")\n", sep = "")
cat("=============================================\n")

cat("\n  Group sizes (6M file):\n")
cat("    GA3Group  0 (>=37wk):", sum(d6$GA3Group == 0, na.rm=TRUE),
    "  1 (28-32wk):", sum(d6$GA3Group == 1, na.rm=TRUE),
    "  2 (<28wk):",  sum(d6$GA3Group == 2, na.rm=TRUE), "\n")
cat("    BPD3Group 0 (HC):",    sum(d6$BPD3Group == 0, na.rm=TRUE),
    "  1 (No+Mild):",           sum(d6$BPD3Group == 1, na.rm=TRUE),
    "  2 (M+S):",                sum(d6$BPD3Group == 2, na.rm=TRUE), "\n")

cat("\n  Missing counts per variable within the n=", nrow(d6),
    " cohort:\n", sep = "")
varlist <- c("sex_male","GA_wk","BWt_g","age_mo","Wt_kg","Ht_cm","BMI",
             "BF_ge6","Sepsis_ever")
for (v in varlist) {
  cat(sprintf("    %-14s  6M missing = %3d   2Y missing = %3d\n", v,
              sum(is.na(d6[[v]])), sum(is.na(d2[[v]]))))
}

# =============================================================================
#  Generate docx files
# =============================================================================

cat("\n=============================================\n")
cat(" Generating docx files\n")
cat("=============================================\n")

out1 <- "Table_1_GA_6M_2Y.docx"
info1 <- build_table_docx(
  group_var    = "GA3Group",
  group_levels = c(0, 1, 2),
  group_labels = c("\u226537 weeks", "28-32 weeks", "<28 weeks"),
  caption      = paste("Table 1. Comparisons of the demographic",
                       "characteristics among full-term and preterm",
                       "infants less than 32 weeks of gestational age",
                       "categorized by different GA at the corrected",
                       "ages of 6 months and 2 years."),
  out_file     = out1)
cat(" [GA ] group n at 6M:", paste(info1$n6, collapse="/"),
    "  at 2Y:", paste(info1$n2, collapse="/"), "\n")
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
cat(" [BPD] group n at 6M:", paste(info2$n6, collapse="/"),
    "  at 2Y:", paste(info2$n2, collapse="/"), "\n")
cat(" Wrote:", normalizePath(out2, mustWork = TRUE), "\n")

cat("\nDone.\n")
