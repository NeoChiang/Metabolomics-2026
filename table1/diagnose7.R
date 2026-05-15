suppressPackageStartupMessages({ library(readxl); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")

rawS <- read_excel("2Y BPD_206 20250721.xlsx",
                   sheet="收案名單", .name_repair="minimal")
rawS <- rawS[!is.na(rawS$`case no.`), ]

cat("収案名單 GA column — class:", class(rawS$GA), "\n")
cat("non-NA count:", sum(!is.na(rawS$GA)), " / ", nrow(rawS), "\n")
cat("first 30 values:\n"); print(head(rawS$GA, 30))
cat("non-NA first 10:\n"); print(head(rawS$GA[!is.na(rawS$GA)], 10))

cat("\n収案名單 BWt — non-NA:", sum(!is.na(rawS$BWt)), " / ", nrow(rawS), "\n")
cat("first 20:\n"); print(head(rawS$BWt, 20))
