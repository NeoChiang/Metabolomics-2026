suppressPackageStartupMessages({ library(readxl); library(dplyr) })
Sys.setlocale("LC_ALL","C.UTF-8")

f6m <- "139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx"
f2y <- "2Y BPD_206 20250721.xlsx"

d6 <- read_excel(f6m, sheet="Blood_Urine", .name_repair="minimal")
d6 <- d6[!is.na(d6$`case no.`), ]
ds <- read_excel(f2y, sheet="收案名單", .name_repair="minimal")
ds <- ds[!is.na(ds$`case no.`), ]
d2 <- read_excel(f2y, sheet="2Y問卷", .name_repair="minimal")
d2 <- d2[!is.na(d2$`case no.`), ]

cat("== 收案名單 GA3Group / BPD3Group ==\n")
print(table(ds$GA3Group, useNA="ifany"))
print(table(ds$BPD3Group, useNA="ifany"))

# verify overlap matches 6M values
cat("\n== Compare 6M GA3Group vs 收案名單 GA3Group on overlap ==\n")
m <- inner_join(
  d6 %>% select(`case no.`, GA3Group_6M=GA3Group, BPD3Group_6M=BPD3Group),
  ds %>% select(`case no.`, GA3Group_S=GA3Group, BPD3Group_S=BPD3Group),
  by="case no."
)
cat("nrow overlap=", nrow(m), "\n")
cat("GA3Group disagreement:\n")
print(table(m$GA3Group_6M, m$GA3Group_S, useNA="ifany"))
cat("BPD3Group disagreement:\n")
print(table(m$BPD3Group_6M, m$BPD3Group_S, useNA="ifany"))

cat("\n== 收案名單 sex / corrected age / 2Y measurements ==\n")
cat("sex values (M/F count):\n"); print(table(ds$sex, useNA="ifany"))
ca_cols <- which(names(ds)=="corrected age")
cat("corrected age columns at positions:", ca_cols, "\n")
# the first 'corrected age' (col 8) corresponds to original baseline measurement?
# col 21 corresponds to blood date_2Y
for (k in ca_cols) {
  v <- suppressWarnings(as.numeric(ds[[k]]))
  cat(sprintf("  pos %2d  preceded-by '%s'  n=%d range=[%.2f, %.2f]\n",
              k, names(ds)[k-1], sum(!is.na(v)), min(v,na.rm=TRUE), max(v,na.rm=TRUE)))
}

cat("\n== 收案名單 2Y height / 2Y weight / BMI completeness ==\n")
for (col in c("2Y height","2Y weight","BMI","Breastfeeding\u2267 months","Sepsis, ever")) {
  if (col %in% names(ds)) {
    v <- ds[[col]]
    cat(sprintf("  %-30s | n_nonNA=%d / %d\n", col, sum(!is.na(v)), nrow(ds)))
  }
}
print(names(ds)[grepl("Breastfeeding|Sepsis", names(ds))])
bcol <- names(ds)[grepl("Breastfeeding", names(ds))][1]
scol <- names(ds)[grepl("Sepsis", names(ds))][1]
cat(sprintf("\nBreastfeeding col: %s -- values:\n", bcol)); print(table(ds[[bcol]], useNA="ifany"))
cat(sprintf("\nSepsis col: %s -- values:\n", scol));     print(table(ds[[scol]], useNA="ifany"))

# Also see how many of 收案名單 cases actually have 2Y measurements
v2h <- suppressWarnings(as.numeric(ds$`2Y height`))
v2w <- suppressWarnings(as.numeric(ds$`2Y weight`))
cat("\n2Y height n_nonNA=", sum(!is.na(v2h)), " 2Y weight n_nonNA=", sum(!is.na(v2w)), "\n")

# corrected age at 2Y blood date (column index 21)
ca_2y <- suppressWarnings(as.numeric(ds[[21]]))
cat("\ncorrected age at blood date_2Y: n=", sum(!is.na(ca_2y)),
    " mean=", round(mean(ca_2y, na.rm=TRUE),2), " sd=", round(sd(ca_2y, na.rm=TRUE),2),
    " range=[", round(min(ca_2y,na.rm=TRUE),2), ",", round(max(ca_2y,na.rm=TRUE),2), "]\n")

# Sex for the 67 only-in-2Y cases - source from 收案名單
only2 <- setdiff(ds$`case no.`, d6$`case no.`)
cat("\n67 only-in-2Y: 收案名單 has sex for", sum(!is.na(ds$sex[ds$`case no.` %in% only2])), "of", length(only2), "\n")
cat("\n67 only-in-2Y: 收案名單 has GA3Group for", sum(!is.na(ds$GA3Group[ds$`case no.` %in% only2])), "of", length(only2), "\n")
cat("67 only-in-2Y: 收案名單 has BPD3Group for", sum(!is.na(ds$BPD3Group[ds$`case no.` %in% only2])), "of", length(only2), "\n")
