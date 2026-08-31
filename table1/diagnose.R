# Diagnostic script: inspect 6M / 2Y files, print column mapping,
# group counts, and case no. overlap between the two files.
# Run from the table1/ folder.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

f6m <- "139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx"
f2y <- "2Y BPD_206 20250721.xlsx"

stopifnot(file.exists(f6m), file.exists(f2y))

cat("=========================================\n")
cat(" 6M file:", f6m, "\n")
cat("=========================================\n")
d6 <- read_excel(f6m, sheet = "Blood_Urine", .name_repair = "minimal")
cat("rows x cols:", nrow(d6), "x", ncol(d6), "\n")
cat("Columns:\n")
for (i in seq_along(names(d6))) cat(sprintf("  [%2d] %s\n", i, names(d6)[i]))

cat("\n2Y file sheets:\n")
print(excel_sheets(f2y))

cat("\n=========================================\n")
cat(" 2Y file (sheet '2Y問卷'):\n")
cat("=========================================\n")
d2 <- read_excel(f2y, sheet = "2Y問卷", .name_repair = "minimal")
cat("rows x cols:", nrow(d2), "x", ncol(d2), "\n")
cat("Columns:\n")
for (i in seq_along(names(d2))) cat(sprintf("  [%2d] %s\n", i, names(d2)[i]))

cat("\n=========================================\n")
cat(" 6M GA3Group / BPD3Group counts\n")
cat("=========================================\n")
cat("GA3Group  (0=>=37wk, 1=28-32wk, 2=<28wk):\n"); print(table(d6$GA3Group, useNA="ifany"))
cat("BPD3Group (0=HC, 1=No+Mild, 2=M+S):\n");      print(table(d6$BPD3Group, useNA="ifany"))

cat("\n=========================================\n")
cat(" case no. overlap (6M vs 2Y)\n")
cat("=========================================\n")
case6 <- as.character(d6[["case no."]])

# Detect case no. column in 2Y
cand <- intersect(names(d2), c("case no.","case no","CaseNo","caseno","case_no","Case no."))
if (length(cand) == 0) {
  # fallback: first column that looks like an id
  cand <- names(d2)[grepl("case", names(d2), ignore.case = TRUE)]
}
cat("2Y case-id candidate columns:", paste(cand, collapse=" | "), "\n")
case2 <- as.character(d2[[cand[1]]])

cat("n unique 6M case no.:", length(unique(case6)), "\n")
cat("n unique 2Y case no.:", length(unique(case2)), "\n")
both <- intersect(case6, case2)
only6 <- setdiff(case6, case2)
only2 <- setdiff(case2, case6)
cat("intersect (in both):", length(both), "\n")
cat("only in 6M         :", length(only6), "\n")
cat("only in 2Y         :", length(only2), "\n")

cat("\nFirst 20 only-in-2Y (would need new GA/BPD assignment):\n")
print(head(only2, 20))
cat("\nFirst 20 only-in-6M (no 2Y follow-up):\n")
print(head(only6, 20))
