suppressPackageStartupMessages({ library(readxl); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")

raw6 <- read_excel("139_BPD_NMR_Urine_6M_list_孟翰_Table 1_20230303.xlsx",
                   sheet="Blood_Urine", .name_repair="minimal")
cat("raw6 class of column 5 (男1/女0):", class(raw6[[5]]), "\n")
cat("raw6 nrow:", nrow(raw6), "\n")
cat("# distinct column names:", length(unique(names(raw6))), "\n")
cat("# total columns:", ncol(raw6), "\n")

# how does readxl expose the duplicate-named columns?
dups <- names(raw6)[duplicated(names(raw6))]
cat("duplicate names (repeated):", paste(unique(dups), collapse=", "), "\n")

# extract sex by position vs by name
cat("\nBy position [[5]] first 10:\n"); print(head(raw6[[5]], 10))
cat("\nBy name $`男1/女0` first 10:\n"); print(head(raw6$`男1/女0`, 10))

# The missing 17 rows: which case_nos?
v5 <- raw6[[5]]
idx <- which(is.na(v5))
cat("\nRows with NA at col 5 (index):", idx, "\n")
cat("Their case no. and name:\n")
print(raw6[idx, c(1,2,3)])

# Full table of case no. vs col 5 availability
cat("\nNumber of rows with non-NA col 5 and non-NA col 17 (BWt):\n")
cat("  col5 non-NA:", sum(!is.na(raw6[[5]])), "\n")
cat("  col17 non-NA:", sum(!is.na(raw6[[17]])), "\n")
cat("  col32 (BMI) non-NA:", sum(!is.na(raw6[[32]])), "\n")
cat("  col28 (Wt)  non-NA:", sum(!is.na(raw6[[28]])), "\n")

# Check which cases in the 139-overlap have NA sex
rawS <- read_excel("2Y BPD_206 20250721.xlsx",
                   sheet="收案名單", .name_repair="minimal")
rawS <- rawS[!is.na(rawS$`case no.`), ]
raw6f <- raw6[!is.na(raw6$`case no.`), ]
ov <- intersect(raw6f$`case no.`, rawS$`case no.`)
cat("\nCases in 6M but not in 収案名單:\n")
print(setdiff(raw6f$`case no.`, ov))

sub6 <- raw6f[raw6f$`case no.` %in% ov, ]
cat("\nIn the 139 overlap, col5 NAs:",  sum(is.na(sub6[[5]])),  "\n")
cat("In the 139 overlap, col17 NAs:", sum(is.na(sub6[[17]])), "\n")
cat("In the 139 overlap, col28 NAs:", sum(is.na(sub6[[28]])), "\n")
cat("In the 139 overlap, col32 NAs:", sum(is.na(sub6[[32]])), "\n")
cat("\nCase numbers in overlap with col28 (Wt) NA:\n")
print(sub6$`case no.`[is.na(sub6[[28]])])
