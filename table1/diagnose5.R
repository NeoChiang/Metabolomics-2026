suppressPackageStartupMessages({ library(readxl); library(dplyr) })
Sys.setlocale("LC_ALL","C.UTF-8")

ds <- read_excel("2Y BPD_206 20250721.xlsx", sheet="收案名單", .name_repair="minimal")
ds <- ds[!is.na(ds$`case no.`), ]

# corrected age at 2Y blood draw is column 21 (right after "blood date_2Y" at col 20)
ca2y <- suppressWarnings(as.numeric(ds[[21]]))
cat("== 収案名單 col 21 (corrected age at blood_2Y) ==\n")
cat("all values sorted (NA removed):\n")
print(sort(ca2y[!is.na(ca2y)]))
cat("\nvalues in (0, 4]:\n")
print(sort(ca2y[!is.na(ca2y) & ca2y > 0 & ca2y <= 4]))
cat("\nvalues > 10 (suspect months not years):\n")
print(sort(ca2y[!is.na(ca2y) & ca2y > 10]))
cat("\nHistogram-ish by bucket:\n")
br <- c(-Inf, -10, 0, 1.5, 2.5, 3.5, 5, 10, 30, Inf)
print(table(cut(ca2y, breaks=br), useNA="ifany"))

# check blood date_2Y (col 20) for a few rows — is it a date?
cat("\n== blood date_2Y (col 20) sample ==\n")
print(head(ds[[20]], 20))
cat("class:", class(ds[[20]]), "\n")

# Sepsis ever clean
cat("\n== Sepsis, ever raw values ==\n")
print(table(ds$`Sepsis, ever`, useNA="ifany"))

# Breastfeeding≧6 months
cat("\n== Breastfeeding≧6 months raw values ==\n")
bcol <- names(ds)[grepl("Breastfeeding", names(ds))][1]
print(table(ds[[bcol]], useNA="ifany"))
