library(tableone)
library(flextable)
library(officer)

data <- read.csv("/home/user/Metabolomics-2026/table1.csv",
                 stringsAsFactors = FALSE, strip.white = TRUE,
                 check.names = FALSE)

colnames(data) <- gsub("[\\r\\n\\s]+", "", colnames(data), perl = TRUE)
colnames(data)[grepl("cell", colnames(data), ignore.case = TRUE)] <- "Analysis_cell_type"

data$group <- factor(data$group, levels = c(12, 3),
                     labels = c("Group 1+2", "Group 3"))

data$Recur <- factor(data$Recur, levels = c(0, 1), labels = c("No", "Yes"))
data$Death <- factor(data$Death, levels = c(0, 1), labels = c("No", "Yes"))
data$Analysis_cell_type <- factor(data$Analysis_cell_type,
                                  levels = c(1, 2),
                                  labels = c("Type 1", "Type 2"))
data$cx <- factor(data$cx, levels = c(0, 1), labels = c("No", "Yes"))
data$adnexa <- factor(data$adnexa, levels = c(0, 1), labels = c("No", "Yes"))
data$PLN <- factor(data$PLN, levels = c(0, 1), labels = c("No", "Yes"))
data$PALN <- factor(data$PALN, levels = c(0, 1), labels = c("No", "Yes"))
data$omentum <- factor(data$omentum, levels = c(0, 1), labels = c("No", "Yes"))

data$MM <- as.numeric(data$MM)

vars <- c("Analysis_cell_type", "age", "Recur", "Death",
          "BH", "BW", "BMI", "CA-125", "WBC", "MM",
          "cx", "adnexa", "PLN", "PALN", "omentum")

catVars <- c("Analysis_cell_type", "Recur", "Death",
             "cx", "adnexa", "PLN", "PALN", "omentum")

nonnormal_vars <- c("CA-125", "WBC", "MM")

tab <- CreateTableOne(vars = vars,
                      strata = "group",
                      data = data,
                      factorVars = catVars,
                      addOverall = TRUE)

tab_mat <- print(tab,
                 nonnormal = nonnormal_vars,
                 exact = catVars,
                 showAllLevels = TRUE,
                 printToggle = FALSE,
                 noSpaces = TRUE)

cat("\n=== Table 1 Preview ===\n")
print(tab_mat)

tab_df <- as.data.frame(tab_mat, stringsAsFactors = FALSE)
tab_df <- cbind(Variable = rownames(tab_df), tab_df)
rownames(tab_df) <- NULL

tab_df$test <- NULL

var_labels <- c(
  "n" = "n",
  "Analysis_cell_type (%)" = "Cell type, n (%)",
  "age (mean (SD))" = "Age, years",
  "Recur (%)" = "Recurrence, n (%)",
  "Death (%)" = "Death, n (%)",
  "BH (mean (SD))" = "Body height, cm",
  "BW (mean (SD))" = "Body weight, kg",
  "BMI (mean (SD))" = "BMI, kg/m²",
  "CA-125 (median [IQR])" = "CA-125, U/mL",
  "WBC (median [IQR])" = "WBC, 10³/µL",
  "MM (median [IQR])" = "MI, %",
  "cx (%)" = "Cervical invasion, n (%)",
  "adnexa (%)" = "Adnexal involvement, n (%)",
  "PLN (%)" = "Pelvic LN metastasis, n (%)",
  "PALN (%)" = "Para-aortic LN metastasis, n (%)",
  "omentum (%)" = "Omentum metastasis, n (%)"
)

for (i in seq_len(nrow(tab_df))) {
  v <- trimws(tab_df$Variable[i])
  if (v %in% names(var_labels)) {
    tab_df$Variable[i] <- var_labels[v]
  }
}

ft <- flextable(tab_df)
ft <- set_header_labels(ft, Variable = "", level = "")
ft <- theme_booktabs(ft)

ft <- fontsize(ft, size = 9, part = "all")
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- padding(ft, padding.top = 1, padding.bottom = 1,
              padding.left = 3, padding.right = 3, part = "all")

ft <- width(ft, j = 1, width = 1.7)
ft <- width(ft, j = 2, width = 0.5)
ft <- width(ft, j = 3, width = 1.25)
ft <- width(ft, j = 4, width = 1.25)
ft <- width(ft, j = 5, width = 1.25)
ft <- width(ft, j = 6, width = 0.55)

ft <- set_table_properties(ft, layout = "fixed", width = 1)

ft <- bold(ft, part = "header")
ft <- align(ft, align = "center", part = "header")
ft <- align(ft, j = 1, align = "left", part = "body")
ft <- align(ft, j = 2:ncol(tab_df), align = "center", part = "body")

p_col <- which(colnames(tab_df) == "p")
if (length(p_col) > 0) {
  for (i in seq_len(nrow(tab_df))) {
    pval <- trimws(tab_df[i, p_col])
    if (!is.na(pval) && pval != "" && pval != "p") {
      pnum <- suppressWarnings(as.numeric(gsub("<", "", pval)))
      if (!is.na(pnum) && pnum < 0.05) {
        ft <- bold(ft, i = i, j = p_col, part = "body")
      }
      if (grepl("^<", pval)) {
        ft <- bold(ft, i = i, j = p_col, part = "body")
      }
    }
  }
}

ft <- border_inner_h(ft, border = fp_border(color = "gray70", width = 0.5))
ft <- hline_top(ft, border = fp_border(color = "black", width = 1.5), part = "header")
ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1.5), part = "header")
ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1.5), part = "body")

output_path <- "/home/user/Metabolomics-2026/Table_1_output.docx"

doc <- read_docx()

doc <- body_add_par(doc,
  "Table 1. Comparisons of the demographic and clinical characteristics between groups.",
  style = "Normal")

doc <- body_add_flextable(doc, ft)

doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc,
  paste0(
    "Data shown are mean ± SD, median [IQR], or number (%) of patients as appropriate. ",
    "BMI, body mass index; CA-125, cancer antigen 125; WBC, white blood cell count; ",
    "MI, myometrial invasion; LN, lymph node. ",
    "P-values < 0.05 (shown in bold) are considered statistically significant. ",
    "Fisher’s exact test was used for categorical variables."
  ),
  style = "Normal")

print(doc, target = output_path)
cat("\nTable 1 saved to:", output_path, "\n")
