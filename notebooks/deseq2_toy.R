# deseq2_toy.R
# Minimal, self-contained DESeq2 demo (2 control, 2 treated)

# 1) Load DESeq2 (install if missing)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)

# 2) Where to put outputs (this script lives in notebooks/, so results/ is ..\results)
results_dir <- file.path("..", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# 3) Tiny toy count matrix (rows = genes, columns = samples)
#    Columns: ctrl1, ctrl2, trt1, trt2  (2 per group so DESeq2 can estimate dispersion)
genes <- paste0("gene", 1:6)
cts <- matrix(c(
  # ctrl1
  100,  50, 300,  80, 500, 20,
  # ctrl2
  120,  60, 310,  90, 480, 25,
  # trt1
  400,  55, 200,  85, 700, 30,
  # trt2
  420,  65, 190,  95, 680, 35
), nrow = 6, ncol = 4, byrow = FALSE,
dimnames = list(genes, c("ctrl1","ctrl2","trt1","trt2")))

# 4) Sample info (design: treated vs control)
colData <- data.frame(
  condition = factor(c("control","control","treated","treated"),
                     levels = c("control","treated")),
  row.names = colnames(cts)
)

# 5) Build DESeq2 object & run
dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ condition)
dds <- DESeq(dds)

# 6) Results: treated vs control
res <- results(dds, contrast = c("condition","treated","control"))
res <- res[order(res$padj), ]

# 7) Also write normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# 8) Save outputs (CSV + MA plot PDF)
res_csv  <- file.path(results_dir, "deseq2_toy_results.csv")
nrm_csv  <- file.path(results_dir, "deseq2_toy_normalized_counts.csv")
ma_pdf   <- file.path(results_dir, "deseq2_toy_MAplot.pdf")

write.csv(as.data.frame(res),        res_csv, row.names = TRUE)
write.csv(as.data.frame(norm_counts), nrm_csv, row.names = TRUE)

pdf(ma_pdf, width = 6, height = 5)
plotMA(res, main = "DESeq2 MA plot (toy)")
abline(h = 0, col = "gray60", lty = 2)
dev.off()

message("Done! Files written:")
message("  - ", res_csv)
message("  - ", nrm_csv)
message("  - ", ma_pdf)