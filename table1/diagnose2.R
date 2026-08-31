suppressPackageStartupMessages({ library(readxl); library(dplyr) })
Sys.setlocale("LC_ALL","C.UTF-8")

f6m <- "139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx"
f2y <- "2Y BPD_206 20250721.xlsx"

# ---- 6M ----
d6 <- read_excel(f6m, sheet="Blood_Urine", .name_repair="minimal")
cat("6M nrow=", nrow(d6), " unique case no.=", length(unique(d6$`case no.`)), "\n")
cat("6M missing case no.:\n")
print(d6[is.na(d6$`case no.`) | d6$`case no.` == "", c("case no.","name","GA3Group","BPD3Group")])

# ---- 2Y questionnaire ----
d2 <- read_excel(f2y, sheet="2Y問卷", .name_repair="minimal")
# strip empty rows
d2 <- d2[!is.na(d2$`case no.`) & nzchar(trimws(as.character(d2$`case no.`))), ]
cat("\n2Y問卷 non-empty rows=", nrow(d2), "\n")
cat("2Y問卷 columns of interest:\n")
for (col in c("case no.","name","Ht (cm)","Wt (kg)","BMI","GA week","GA day",
              "Sepsis","once Sepsis","Breast feeding","Breast feeding duration",
              "B Wt (g)","Term/Preterm")) {
  if (col %in% names(d2)) {
    cat(sprintf("  %-30s | example: %s\n", col,
                paste(utils::head(unique(na.omit(d2[[col]])), 5), collapse=" / ")))
  } else {
    cat(sprintf("  %-30s | MISSING\n", col))
  }
}

# search for sex column
sex_cols <- grep("sex|男|gender|female|male|性別", names(d2), ignore.case=TRUE, value=TRUE)
cat("\nPossible sex columns in 2Y問卷:\n"); print(sex_cols)

# Other sheets - check 收案名單 for sex / GA info
cat("\n=== 收案名單 ===\n")
ds <- read_excel(f2y, sheet="收案名單", .name_repair="minimal")
ds <- ds[!is.na(ds[[1]]) & nzchar(trimws(as.character(ds[[1]]))), ]
cat("rows=", nrow(ds), " cols=", ncol(ds), "\n")
cat("columns:\n"); print(names(ds))

cat("\n=== 6M問卷 (might give sex for 2Y-only cases) ===\n")
d6q <- read_excel(f2y, sheet="6M問卷", .name_repair="minimal")
d6q <- d6q[!is.na(d6q[[1]]) & nzchar(trimws(as.character(d6q[[1]]))), ]
cat("rows=", nrow(d6q), " cols=", ncol(d6q), "\n")
sex_cols2 <- grep("sex|男|gender|female|male|性別", names(d6q), ignore.case=TRUE, value=TRUE)
cat("sex-like columns in 6M問卷:\n"); print(sex_cols2)
ga_cols <- grep("GA|wk|week|gestat", names(d6q), ignore.case=TRUE, value=TRUE)
cat("GA-like columns in 6M問卷:\n"); print(ga_cols)

# show ranges
cat("\n2Y Ht/Wt/BMI summary:\n")
for (col in c("Ht (cm)","Wt (kg)","BMI")) {
  v <- suppressWarnings(as.numeric(d2[[col]]))
  cat(sprintf("  %-10s n=%d mean=%.2f sd=%.2f range=[%.2f, %.2f]\n",
              col, sum(!is.na(v)), mean(v,na.rm=TRUE), sd(v,na.rm=TRUE),
              min(v,na.rm=TRUE), max(v,na.rm=TRUE)))
}

# 2Y has no explicit corrected-age column; ~24 mo by definition.
# Show GA week
gw <- suppressWarnings(as.numeric(d2$`GA week`))
gd <- suppressWarnings(as.numeric(d2$`GA day`))
cat("\n2Y GA week summary: n=", sum(!is.na(gw)), " range=[", min(gw,na.rm=TRUE), ",", max(gw,na.rm=TRUE), "]\n")

cat("\n2Y Sepsis / once Sepsis values:\n")
print(table(d2$Sepsis, useNA="ifany"))
print(table(d2$`once Sepsis`, useNA="ifany"))

cat("\n2Y Breast feeding duration head:\n")
print(head(d2$`Breast feeding duration`, 30))
