suppressPackageStartupMessages({ library(Mfuzz); library(Biobase) })
set.seed(42); MIN_N <- 3

clean_t <- read.csv("all_tables_clean.csv", check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
sig_mets <- unique(clean_t$Metabolite[as.numeric(sub("<","",clean_t$P)) < 0.05])

dat <- read.csv("6M-2Y_Urine_all metabolites_20260330.csv", check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
dat$Age <- as.character(dat$Age); dat$GA3Group <- as.integer(dat$GA3Group); dat$BPD3Group <- as.integer(dat$BPD3Group)
exclude <- c("Creatinine","TSP")
met_cols <- setdiff(colnames(dat)[5:ncol(dat)], exclude)
for (mc in c(met_cols,"Creatinine")) dat[[mc]] <- suppressWarnings(as.numeric(dat[[mc]]))
crea <- dat$Creatinine
for (mc in met_cols) dat[[mc]] <- dat[[mc]] / crea
mets <- intersect(sig_mets, met_cols)

log2fc <- function(met, grpcol, grpval) {
  v6 <- dat[[met]][dat$Age=="6M" & dat[[grpcol]]==grpval]
  v2 <- dat[[met]][dat$Age=="2Y" & dat[[grpcol]]==grpval]
  v6 <- v6[is.finite(v6)]; v2 <- v2[is.finite(v2)]
  if (length(v6)<MIN_N || length(v2)<MIN_N) return(NA_real_)
  m6 <- mean(v6); m2 <- mean(v2)
  if (m6<=0 || m2<=0) return(NA_real_)
  log2(m2/m6)
}

bpd_fc <- matrix(NA_real_, length(mets), 3, dimnames=list(mets, c("BPD_HC","BPD_No+Mild","BPD_M+S")))
for (m in mets) for (g in 0:2) bpd_fc[m,g+1] <- log2fc(m,"BPD3Group",g)
M <- bpd_fc[complete.cases(bpd_fc),,drop=FALSE]

eset   <- ExpressionSet(assayData = M)
eset_s <- standardise(eset)
m_val  <- tryCatch(mestimate(eset_s), error = function(e) 1.25)
if (!is.finite(m_val) || m_val < 1.05) m_val <- 1.25
k_use <- 3
cl <- mfuzz(eset_s, c = k_use, m = m_val)

cat("Cluster assignments (original, no remap):\n")
for (i in 1:3) {
  members <- names(cl$cluster[cl$cluster == i])
  cat(sprintf("  Cluster %d (%d): %s\n", i, length(members), paste(members, collapse=", ")))
}

nc <- 3; nr <- 1; sq <- 4.5
pdf("mfuzz_FC_BPD_std_Crea_red.pdf", width = nc*sq+1, height = nr*sq+0.8)
par(pty = "s")
mfuzz.plot2(eset_s, cl = cl, mfrow = c(nr, nc), time.labels = colnames(M),
            centre = TRUE, centre.col = "darkred", centre.lwd = 2.5,
            x11 = FALSE, ylab = "Standardised log2-FC")
mtext(sprintf("Mfuzz log2-FC (std, Crea-norm) - FC_BPD_std_Crea  (k=3, m=%.2f, n=%d)", m_val, nrow(M)),
      side = 3, line = -1.5, outer = TRUE, cex = 1.1, font = 2)
dev.off()

memb <- cl$membership
colnames(memb) <- paste0("Cluster_", seq_len(ncol(memb)))
out <- data.frame(Metabolite = rownames(M), HardCluster = cl$cluster, memb, M, check.names = FALSE)
write.csv(out, "mfuzz_FC_BPD_std_Crea_membership.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("\nDone -> mfuzz_FC_BPD_std_Crea_red.pdf + mfuzz_FC_BPD_std_Crea_membership.csv\n")
