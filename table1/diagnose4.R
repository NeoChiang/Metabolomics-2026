suppressPackageStartupMessages({ library(readxl); library(dplyr) })
Sys.setlocale("LC_ALL","C.UTF-8")

f2y <- "2Y BPD_206 20250721.xlsx"
d2 <- read_excel(f2y, sheet="2Y問卷", .name_repair="minimal", col_types="text")
d2 <- d2[!is.na(d2$`case no.`), ]
cat("2Y問卷 rows=", nrow(d2), "\n")

# raw vs numeric counts
for (col in c("Ht (cm)","Wt (kg)","BMI","B Wt (g)","GA week","Sepsis","once Sepsis","Breast feeding","Breast feeding duration")) {
  raw_n <- sum(!is.na(d2[[col]]) & nzchar(trimws(d2[[col]])))
  num_n <- sum(!is.na(suppressWarnings(as.numeric(d2[[col]]))))
  cat(sprintf("  %-30s raw_nonempty=%3d  numeric=%3d\n", col, raw_n, num_n))
}

cat("\nSample of Ht (cm) and Wt (kg) raw values (first 25 non-empty):\n")
nz <- d2[!is.na(d2$`Ht (cm)`) & nzchar(trimws(d2$`Ht (cm)`)), c("case no.","Ht (cm)","Wt (kg)","BMI","BMI")[1:4]]
print(head(nz, 25))

# Breast feeding duration parsing examples
cat("\nUnique-ish Breast feeding duration patterns:\n")
print(sort(table(d2$`Breast feeding duration`, useNA="ifany"), decreasing=TRUE)[1:25])

# Confirm Sex coding using 6M file: 男1/女0
d6 <- read_excel("139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx",
                 sheet="Blood_Urine", .name_repair="minimal")
ds <- read_excel(f2y, sheet="收案名單", .name_repair="minimal")
mtab <- inner_join(
  d6 %>% select(`case no.`, sex_6M=`男1/女0`),
  ds %>% select(`case no.`, sex_S=sex),
  by="case no."
)
cat("\n6M sex (1=male, 0=female) vs 收案名單 sex on overlap:\n")
print(table(mtab$sex_6M, mtab$sex_S, useNA="ifany"))
