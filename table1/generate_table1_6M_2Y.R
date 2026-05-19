# =============================================================================
#  generate_table1_6M_2Y.R
#
#  Produces two demographic tables (landscape .docx):
#      Table_1_GA_6M_2Y.docx   — grouped by GA3Group
#      Table_1_BPD_6M_2Y.docx  — grouped by BPD3Group
#
#  Cohort: filtered by "sample list.xlsx" (6M = 139, 2Y = 115).
#
#  Data sources:
#      6M : Blood_Urine sheet.  2Y : 收案名單 sheet.
#      Body weight / height / BMI use PERCENTILE columns.
#      6M BMI percentile computed from WHO LMS parameters (0-24 months).
#
#  Header layout:
#      Char | --- 6M subgroups --- | P | --- 2Y subgroups --- | P
#
#  Dependencies: readxl, dplyr, officer, flextable
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
f6m <- "table1/139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx"
f2y <- "table1/2Y BPD_206 20250721.xlsx"
fsl <- "table1/sample list.xlsx"
if (!file.exists(f6m)) { f6m <- sub("table1/","",f6m); f2y <- sub("table1/","",f2y); fsl <- sub("table1/","",fsl) }
stopifnot(file.exists(f6m), file.exists(f2y), file.exists(fsl))

# ------------------------------------------------------------------ load raw
raw6 <- read_excel(f6m, sheet = "Blood_Urine", .name_repair = "minimal")
raw6 <- raw6[!is.na(raw6$`case no.`), ]

rawS <- read_excel(f2y, sheet = "收案名單", .name_repair = "minimal")
rawS <- rawS[!is.na(rawS$`case no.`), ]

# ------------------------------------------------------------------ sample list
sl <- read_excel(fsl, sheet = 1, .name_repair = "minimal")
sl_6m <- na.omit(as.character(sl$`6m`))
sl_2y <- na.omit(as.character(sl$`2y`))

# ------------------------------------------------------------------ WHO LMS
# BMI-for-age LMS parameters (WHO Child Growth Standards, 0-24 months).
# Used to compute 6M BMI percentile which is absent from the source file.
who_bmi_lms <- data.frame(
  month = rep(0:24, 2),
  sex   = rep(c(1L, 0L), each = 25),   # 1 = male, 0 = female
  L = c(# boys
        -0.3053, 0.2560, 0.2313,-0.0590,-0.4019,-0.6996,-0.9200,-1.0698,
        -1.1661,-1.2236,-1.2541,-1.2667,-1.2667,-1.2581,-1.2430,-1.2231,
        -1.1998,-1.1738,-1.1458,-1.1161,-1.0850,-1.0529,-1.0199,-0.9863,-0.9524,
        # girls
        -0.0631, 0.3711, 0.2024,-0.1530,-0.4840,-0.7341,-0.9009,-1.0046,
        -1.0639,-1.0935,-1.1017,-1.0947,-1.0773,-1.0531,-1.0247,-0.9937,
        -0.9612,-0.9280,-0.8942,-0.8604,-0.8266,-0.7930,-0.7599,-0.7274,-0.6955),
  M = c(# boys
        13.4069,14.8590,16.4459,17.2485,17.5625,17.5949,17.4375,17.1863,
        16.9131,16.6461,16.4000,16.1778,15.9756,15.7896,15.6185,15.4621,
        15.3199,15.1903,15.0713,14.9641,14.8684,14.7835,14.7073,14.6370,14.5780,
        # girls
        13.3363,14.5679,15.9859,16.6853,16.9529,16.9491,16.7982,16.5860,
        16.3574,16.1335,15.9229,15.7275,15.5471,15.3800,15.2253,15.0828,
        14.9519,14.8322,14.7234,14.6241,14.5358,14.4558,14.3820,14.3139,14.2504),
  S = c(# boys
        0.09530,0.08573,0.08251,0.08051,0.07902,0.07813,0.07768,0.07756,
        0.07768,0.07797,0.07838,0.07890,0.07949,0.08016,0.08090,0.08170,
        0.08257,0.08349,0.08447,0.08548,0.08654,0.08762,0.08873,0.08987,0.09103,
        # girls
        0.09300,0.08597,0.08372,0.08218,0.08106,0.08042,0.08016,0.08019,
        0.08042,0.08082,0.08133,0.08195,0.08268,0.08350,0.08442,0.08542,
        0.08649,0.08764,0.08884,0.09010,0.09140,0.09275,0.09413,0.09555,0.09699)
)

bmi_to_percentile <- function(bmi, age_months, sex_male) {
  n <- length(bmi)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    b <- bmi[i]; a <- age_months[i]; s <- sex_male[i]
    if (any(is.na(c(b, a, s))) || a < 0 || a > 24) next
    sub <- who_bmi_lms[who_bmi_lms$sex == s, ]
    # linear interpolation of L, M, S
    L <- approx(sub$month, sub$L, xout = a, rule = 2)$y
    M <- approx(sub$month, sub$M, xout = a, rule = 2)$y
    S <- approx(sub$month, sub$S, xout = a, rule = 2)$y
    z <- if (abs(L) > 0.001) ((b / M)^L - 1) / (L * S) else log(b / M) / S
    out[i] <- pnorm(z) * 100
  }
  out
}

# ------------------------------------------------------------------ 6M tidy
d6 <- data.frame(
  case_no      = as.character(raw6$`case no.`),
  GA3Group     = as.numeric(raw6$GA3Group),
  BPD3Group    = as.numeric(raw6$BPD3Group),
  sex_male     = as.numeric(raw6$`男1/女0`),
  GA_wk        = suppressWarnings(as.numeric(raw6$GA)),
  BWt_g        = suppressWarnings(as.numeric(raw6$BWt)),
  age_mo       = suppressWarnings(as.numeric(raw6$`右邊BW/BH/BMI 測量CA(月份)`)),
  Wt_pct       = suppressWarnings(as.numeric(raw6[[29]])),   # Wt (kg) percentile
  Ht_pct       = suppressWarnings(as.numeric(raw6[[31]])),   # Ht (cm) percentile
  BMI_raw      = suppressWarnings(as.numeric(raw6[[32]])),   # raw BMI for percentile calc
  BF_raw       = suppressWarnings(as.numeric(raw6[[54]])),
  Sepsis_raw   = suppressWarnings(as.numeric(raw6[[41]])),
  check.names  = FALSE
)
d6$BF_ge6      <- ifelse(is.na(d6$BF_raw), NA_integer_, ifelse(d6$BF_raw == 2, 1L, 0L))
d6$Sepsis_ever <- ifelse(is.na(d6$Sepsis_raw), NA_integer_, ifelse(d6$Sepsis_raw >= 1, 1L, 0L))
d6$BMI_pct     <- round(bmi_to_percentile(d6$BMI_raw, d6$age_mo, d6$sex_male), 1)
d6$BF_raw <- NULL; d6$Sepsis_raw <- NULL; d6$BMI_raw <- NULL

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
    if (length(mi) < 3 || mi[1] == "") as.numeric(s[i])
    else {
      wk <- as.numeric(mi[2])
      d  <- suppressWarnings(as.numeric(mi[3]))
      if (is.na(d)) d <- 0
      wk + d / 7
    }
  })
}

d2 <- data.frame(
  case_no      = as.character(rawS$`case no.`),
  GA3Group     = as.numeric(rawS$GA3Group),
  BPD3Group    = as.numeric(rawS$BPD3Group),
  sex_male     = as.numeric(rawS$sex),
  GA_wk        = parse_ga_str(rawS$GA),
  BWt_g        = suppressWarnings(as.numeric(rawS$BWt)),
  age_mo       = ca_years * 12,
  Wt_pct       = suppressWarnings(as.numeric(rawS[[33]])),   # 2Y weight, percentile
  Ht_pct       = suppressWarnings(as.numeric(rawS[[34]])),   # 2Y height, percentile
  BMI_pct      = suppressWarnings(as.numeric(rawS[[35]])),   # 2Y BMI, percentile
  BF_raw       = suppressWarnings(as.numeric(rawS[[grep("Breastfeeding", names(rawS))[1]]])),
  Sepsis_ever  = sepsis_num,
  check.names  = FALSE
)
d2$BF_ge6 <- ifelse(is.na(d2$BF_raw), NA_integer_, ifelse(d2$BF_raw == 2, 1L, 0L))
d2$BF_raw <- NULL

# ------------------------------------------------------------------ filter
d6 <- d6[d6$case_no %in% sl_6m, ]
d2 <- d2[d2$case_no %in% sl_2y, ]

# backfill 6M time-invariant (sex, GA, BWt) from 収案名單
m_6toS <- match(d6$case_no, d2$case_no)
coal <- function(a, b) ifelse(!is.na(a), a, b)
for (v in c("sex_male", "GA_wk", "BWt_g")) d6[[v]] <- coal(d6[[v]], d2[[v]][m_6toS])

cat("\n=============================================\n")
cat(" Sample sizes (filtered by sample list)\n")
cat("=============================================\n")
cat("  6M n =", nrow(d6), "  2Y n =", nrow(d2), "\n")
cat("  overlap =", length(intersect(d6$case_no, d2$case_no)), "\n")
cat("  6M GA3Group  0/1/2:", sum(d6$GA3Group==0,na.rm=T),"/",sum(d6$GA3Group==1,na.rm=T),"/",sum(d6$GA3Group==2,na.rm=T),"\n")
cat("  6M BPD3Group 0/1/2:", sum(d6$BPD3Group==0,na.rm=T),"/",sum(d6$BPD3Group==1,na.rm=T),"/",sum(d6$BPD3Group==2,na.rm=T),"\n")
cat("  2Y GA3Group  0/1/2:", sum(d2$GA3Group==0,na.rm=T),"/",sum(d2$GA3Group==1,na.rm=T),"/",sum(d2$GA3Group==2,na.rm=T),"\n")
cat("  2Y BPD3Group 0/1/2:", sum(d2$BPD3Group==0,na.rm=T),"/",sum(d2$BPD3Group==1,na.rm=T),"/",sum(d2$BPD3Group==2,na.rm=T),"\n")
varlist <- c("sex_male","GA_wk","BWt_g","age_mo","Wt_pct","Ht_pct","BMI_pct","BF_ge6","Sepsis_ever")
cat("\n  Missing counts:\n")
for (v in varlist) cat(sprintf("    %-12s  6M %3d/%3d   2Y %3d/%3d\n", v,
                               sum(is.na(d6[[v]])), nrow(d6), sum(is.na(d2[[v]])), nrow(d2)))

# =============================================================================
#  Helpers
# =============================================================================

fmt_p <- function(p) { if (is.na(p)) "-" else if (p < 0.001) "<0.001" else formatC(p, format="f", digits=3) }

fmt_meanSD <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("-")
  sprintf("%s ± %s", formatC(mean(x), format="f", digits=2), formatC(sd(x), format="f", digits=2))
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
  tryCatch(summary(aov(value ~ group))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
}

p_categorical <- function(value, group) {
  ok <- !is.na(value) & !is.na(group)
  tab <- table(group[ok], value[ok])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  expc <- tryCatch(suppressWarnings(chisq.test(tab)$expected), error = function(e) NULL)
  if (!is.null(expc) && any(expc < 5))
    tryCatch(fisher.test(tab, workspace=2e7)$p.value, error = function(e) NA_real_)
  else
    tryCatch(suppressWarnings(chisq.test(tab)$p.value), error = function(e) NA_real_)
}

row_spec <- list(
  list(label = "Sex (male, n(%))",                      var = "sex_male",    type = "cat"),
  list(label = "Gestational age (wk)",                  var = "GA_wk",       type = "cont"),
  list(label = "Birth body weight (g)",                 var = "BWt_g",       type = "cont"),
  list(label = "Age corrected (month)",                 var = "age_mo",      type = "cont"),
  list(label = "Body weight (percentile)",              var = "Wt_pct",      type = "cont"),
  list(label = "Body height (percentile)",              var = "Ht_pct",      type = "cont"),
  list(label = "BMI (percentile)",                      var = "BMI_pct",     type = "cont"),
  list(label = "Breastfeeding ≥6 months, n(%)",    var = "BF_ge6",      type = "cat"),
  list(label = "Sepsis ever, n(%)",                     var = "Sepsis_ever", type = "cat")
)

# =============================================================================
#  Build one docx (landscape)
# =============================================================================

build_table_docx <- function(group_var, group_levels, group_labels,
                             caption, out_file) {

  data_by <- function(df, lvl) df[df[[group_var]] == lvl, , drop = FALSE]
  K <- length(group_levels)

  # Char | 6M_g1..K | P(6M) | 2Y_g1..K | P(2Y)
  ncols <- 1 + K + 1 + K + 1
  body      <- matrix("",    nrow = length(row_spec), ncol = ncols)
  bold_flag <- matrix(FALSE, nrow = length(row_spec), ncol = ncols)

  col_P6 <- 2 + K
  col_2Y <- 2 + K + 1
  col_P2 <- ncols

  for (i in seq_along(row_spec)) {
    spec <- row_spec[[i]]
    body[i, 1] <- spec$label
    for (k in seq_along(group_levels)) {
      lvl <- group_levels[k]
      s6 <- data_by(d6, lvl)[[spec$var]]
      s2 <- data_by(d2, lvl)[[spec$var]]
      c6 <- 1 + k
      c2 <- col_2Y + (k - 1)
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
    if (spec$type == "cont") { p6 <- p_continuous(v6,g6); p2 <- p_continuous(v2,g2) }
    else                     { p6 <- p_categorical(v6,g6); p2 <- p_categorical(v2,g2) }
    body[i, col_P6] <- fmt_p(p6)
    body[i, col_P2] <- fmt_p(p2)
    if (!is.na(p6) && p6 < 0.05) bold_flag[i, col_P6] <- TRUE
    if (!is.na(p2) && p2 < 0.05) bold_flag[i, col_P2] <- TRUE
  }

  n6 <- sapply(group_levels, function(l) sum(d6[[group_var]]==l, na.rm=TRUE))
  n2 <- sapply(group_levels, function(l) sum(d2[[group_var]]==l, na.rm=TRUE))

  col_keys <- c("char", sprintf("c6_%d",1:K), "pval6", sprintf("c2_%d",1:K), "pval2")
  sub_hdr  <- c("Characteristics",
                sprintf("%s (n=%d)", group_labels, n6), "P",
                sprintf("%s (n=%d)", group_labels, n2), "P")
  top_hdr  <- c("Characteristics", rep("6M", K+1), rep("2Y", K+1))

  df <- as.data.frame(body, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- col_keys

  ft <- flextable(df, col_keys = col_keys)
  hdr <- data.frame(col_keys = col_keys, top = top_hdr, sub = sub_hdr, stringsAsFactors = FALSE)
  ft <- set_header_df(ft, mapping = hdr, key = "col_keys")
  ft <- merge_h(ft, part = "header")
  ft <- merge_v(ft, j = 1, part = "header")

  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "body")
  ft <- bold(ft, part = "header")
  ft <- fontsize(ft, size = 9,  part = "body")
  ft <- fontsize(ft, size = 10, part = "header")
  ft <- padding(ft, padding.top = 2, padding.bottom = 2, part = "all")
  ft <- border_outer(ft, border = fp_border(color="black", width=1))
  ft <- border_inner_h(ft, border = fp_border(color="grey60", width=0.5))
  ft <- hline_bottom(ft, border = fp_border(color="black", width=1), part = "header")
  sep <- fp_border(color="black", width=0.8)
  ft <- vline(ft, j = 1,        border = sep, part = "all")
  ft <- vline(ft, j = col_P6,   border = sep, part = "all")
  ft <- vline(ft, j = col_P6-1, border = sep, part = "body")
  ft <- vline(ft, j = col_P2-1, border = sep, part = "body")

  for (i in seq_len(nrow(bold_flag)))
    for (j in seq_len(ncol(bold_flag)))
      if (bold_flag[i,j]) ft <- bold(ft, i=i, j=j, part="body")

  ft <- set_table_properties(ft, layout = "autofit")

  footnote <- paste(
    "Data shown are mean ± SD or number (%) of patients as appropriate.",
    "Body weight, body height and BMI are presented as percentiles.",
    "GA, gestational age; BPD, bronchopulmonary dysplasia; HC, healthy controls;",
    "No+Mild BPD, no and mild BPD; M+S BPD, moderate and severe BPD;",
    "wk, week; g, gram; BMI, body mass index.",
    "6M BMI percentile was computed from WHO Child Growth Standards LMS parameters.",
    "All P-values < 0.05, which is in bold, are significant.")

  # landscape page
  doc <- read_docx()
  doc <- body_set_default_section(doc,
           prop_section(page_size = page_size(width = 11, height = 8.5, orient = "landscape")))
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
cat(" Generating docx files (landscape)\n")
cat("=============================================\n")

out1 <- "table1/Table_1_GA_6M_2Y.docx"
info1 <- build_table_docx(
  group_var = "GA3Group", group_levels = c(0,1,2),
  group_labels = c("≥37 weeks", "28-32 weeks", "<28 weeks"),
  caption = paste("Table 1. Comparisons of the demographic characteristics",
                  "among full-term and preterm infants categorized by",
                  "different GA at the corrected ages of 6 months and 2 years."),
  out_file = out1)
cat(" [GA ] 6M:", paste(info1$n6, collapse="/"), " 2Y:", paste(info1$n2, collapse="/"), "\n")
cat(" Wrote:", normalizePath(out1, mustWork=TRUE), "\n")

out2 <- "table1/Table_1_BPD_6M_2Y.docx"
info2 <- build_table_docx(
  group_var = "BPD3Group", group_levels = c(0,1,2),
  group_labels = c("HC", "No+Mild BPD", "M+S BPD"),
  caption = paste("Table 1. Comparisons of the demographic characteristics",
                  "among full-term and preterm infants categorized by",
                  "BPD severity at the corrected ages of 6 months and 2 years."),
  out_file = out2)
cat(" [BPD] 6M:", paste(info2$n6, collapse="/"), " 2Y:", paste(info2$n2, collapse="/"), "\n")
cat(" Wrote:", normalizePath(out2, mustWork=TRUE), "\n")

cat("\nDone.\n")
