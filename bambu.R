# ============================================================
# Bambu + DESeq2 DRS数据分析完整流程
# 功能：为每个样本单独生成转录本/基因counts，然后进行差异分析
# ============================================================

# 加载必要的包
library(bambu)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# ============================================================
# 第一部分：设置文件路径和参数
# ============================================================
# 工作目录
work_dir <- "/data2/users/zhouyiming/DT4_ROOT/DEG"
setwd(work_dir)

# BAM文件路径
bam_files <- c(
  "/data2/users/zhouyiming/DT4_ROOT/analysis/aligned/ctrlR_1/ctrlR_1_genome_sort.bam",
  "/data2/users/zhouyiming/DT4_ROOT/analysis/aligned/ctrlR_2/ctrlR_2_genome_sort.bam",
  "/data2/users/zhouyiming/DT4_ROOT/analysis/aligned/dt4R_1/dt4R_1_genome_sort.bam",
  "/data2/users/zhouyiming/DT4_ROOT/analysis/aligned/dt4R_2/dt4R_2_genome_sort.bam"
)

# 样本名称
sample_names <- c("ctrlL_1", "ctrlL_2", "dt4L_1", "dt4L_2")
names(bam_files) <- sample_names

# 参考文件路径
ref_gtf <- "/data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/annotation/osa_IRGSP_1.annotation.gtf"
ref_genome <- "/data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/genome/osa_IRGSP_1.fa"


output_dir <- file.path(work_dir, "bambu_deseq2_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for(sample in sample_names) {
  dir.create(file.path(output_dir, sample), showWarnings = FALSE, recursive = TRUE)
}

# ============================================================
# 第二部分：准备参考注释
# ============================================================

message("\n=== 读取参考注释文件 ===")
annotations <- prepareAnnotations(ref_gtf)
message("注释文件读取完成!")

# ============================================================
# 第三部分：定义安全的Bambu运行函数（处理warnings bug）
# ============================================================

# 方法1：使用 tryCatch 包装 bambu 函数
run_bambu_safe <- function(reads, annotations, genome, ncore = 4, discovery = TRUE, quant = TRUE) {
  
  result <- tryCatch({
    # 尝试正常运行
    se <- bambu(
      reads = reads,
      annotations = annotations,
      genome = genome,
      ncore = ncore,
      discovery = discovery,
      quant = quant,
      verbose = TRUE
    )
    return(se)
  }, error = function(e) {
    # 如果是 warnings 相关的错误，尝试用不同方式运行
    if(grepl("names|warnings", e$message)) {
      message("检测到warnings处理错误，尝试替代方法...")
      
      # 方法2：分步运行
      tryCatch({
        # 先生成read class
        rcFile <- bambu.processReads(
          reads = reads,
          annotations = annotations,
          genome = genome,
          ncore = ncore
        )
        
        # 再进行定量
        se <- bambu(
          reads = rcFile,
          annotations = annotations,
          genome = genome,
          ncore = ncore,
          discovery = discovery,
          quant = quant,
          rcOutDir = NULL,
          verbose = TRUE
        )
        return(se)
      }, error = function(e2) {
        message("替代方法也失败，尝试最小化运行...")
        
        # 方法3：关闭discovery，只做定量
        se <- bambu(
          reads = reads,
          annotations = annotations,
          genome = genome,
          ncore = ncore,
          discovery = FALSE,  # 关闭新转录本发现
          quant = TRUE,
          verbose = TRUE
        )
        return(se)
      })
    } else {
      stop(e)
    }
  })
  
  return(result)
}

# ============================================================
# 第四部分：对每个样本单独运行Bambu
# ============================================================

message("\n=== 对每个样本单独运行Bambu分析 ===")

for(i in seq_along(sample_names)) {
  
  sample <- sample_names[i]
  bam <- bam_files[i]
  sample_dir <- file.path(output_dir, sample)
  
  message(paste0("\n--- 处理样本: ", sample, " (", i, "/", length(sample_names), ") ---"))
  message(paste0("时间: ", Sys.time()))
  
  tryCatch({
    # 使用安全的bambu运行函数
    se <- run_bambu_safe(
      reads = bam,
      annotations = annotations,
      genome = ref_genome,
      ncore = 4,  # 减少核心数可能有助于稳定性
      discovery = TRUE,
      quant = TRUE
    )
    
    # 保存Bambu对象
    saveRDS(se, file = file.path(sample_dir, paste0(sample, "_bambu_se.rds")))
    
    # 提取转录本counts
    tx_counts <- assays(se)$counts
    colnames(tx_counts) <- sample
    tx_counts_df <- data.frame(
      transcript_id = rownames(tx_counts),
      counts = as.vector(tx_counts)
    )
    colnames(tx_counts_df)[2] <- sample
    
    # 提取基因counts
    gene_counts <- rowsum(as.matrix(tx_counts), rowData(se)$GENEID)
    gene_counts_df <- data.frame(
      gene_id = rownames(gene_counts),
      counts = as.vector(gene_counts)
    )
    colnames(gene_counts_df)[2] <- sample
    
    # 提取CPM
    if("CPM" %in% names(assays(se))) {
      tx_cpm <- assays(se)$CPM
      tx_cpm_df <- data.frame(
        transcript_id = rownames(tx_cpm),
        CPM = as.vector(tx_cpm)
      )
      colnames(tx_cpm_df)[2] <- paste0(sample, "_CPM")
      write.csv(tx_cpm_df, file.path(sample_dir, paste0(sample, "_transcript_CPM.csv")), row.names = FALSE)
    }
    
    # 保存counts文件
    write.csv(tx_counts_df, file.path(sample_dir, paste0(sample, "_transcript_counts.csv")), row.names = FALSE)
    write.csv(gene_counts_df, file.path(sample_dir, paste0(sample, "_gene_counts.csv")), row.names = FALSE)
    
    message(paste0("样本 ", sample, " 完成!"))
    message(paste0("  - 转录本数量: ", nrow(tx_counts)))
    message(paste0("  - 基因数量: ", nrow(gene_counts)))
    
  }, error = function(e) {
    message(paste0("样本 ", sample, " 处理出错: ", e$message))
    message("继续处理下一个样本...")
  })
}

# ============================================================
# 第五部分：合并所有样本运行Bambu
# ============================================================

message("\n=== 合并所有样本进行Bambu分析 ===")
message(paste0("时间: ", Sys.time()))

# 使用安全的bambu运行函数
se_all <- tryCatch({
  run_bambu_safe(
    reads = bam_files,
    annotations = annotations,
    genome = ref_genome,
    ncore = 4,
    discovery = TRUE,
    quant = TRUE
  )
}, error = function(e) {
  message("合并分析出错，尝试关闭discovery模式...")
  bambu(
    reads = bam_files,
    annotations = annotations,
    genome = ref_genome,
    ncore = 4,
    discovery = FALSE,
    quant = TRUE,
    verbose = TRUE
  )
})

# 保存合并的Bambu对象
saveRDS(se_all, file = file.path(output_dir, "bambu_all_samples_se.rds"))

# 提取counts矩阵
transcript_counts_all <- assays(se_all)$counts
colnames(transcript_counts_all) <- sample_names
gene_counts_all <- rowsum(transcript_counts_all, rowData(se_all)$GENEID)

# 保存合并的counts矩阵
write.csv(transcript_counts_all, file.path(output_dir, "all_samples_transcript_counts.csv"), row.names = TRUE)
write.csv(gene_counts_all, file.path(output_dir, "all_samples_gene_counts.csv"), row.names = TRUE)

# 保存转录本注释信息
tx_info <- as.data.frame(rowData(se_all))
write.csv(tx_info, file.path(output_dir, "transcript_annotation_info.csv"), row.names = TRUE)

message("合并完成!")
message(paste0("  - 总转录本数量: ", nrow(transcript_counts_all)))
message(paste0("  - 总基因数量: ", nrow(gene_counts_all)))

novel_tx <- sum(grepl("^tx\\.", rownames(transcript_counts_all)))
message(paste0("  - 新发现转录本数量: ", novel_tx))

# ============================================================
# 第六部分：DESeq2 差异表达分析
# ============================================================

message("\n=== 开始DESeq2差异表达分析 ===")

sample_info <- data.frame(
  sample = sample_names,
  condition = factor(c("ctrl", "ctrl", "dt4", "dt4"), levels = c("ctrl", "dt4")),
  row.names = sample_names
)

message("样本信息:")
print(sample_info)

# ----- 基因水平分析 -----
message("\n--- 基因水平差异分析 ---")

gene_counts_int <- round(gene_counts_all)
mode(gene_counts_int) <- "integer"

keep_gene <- rowSums(gene_counts_int >= 10) >= 2
gene_counts_filtered <- gene_counts_int[keep_gene, ]
message(paste0("过滤后基因数量: ", nrow(gene_counts_filtered)))

dds_gene <- DESeqDataSetFromMatrix(
  countData = gene_counts_filtered,
  colData = sample_info,
  design = ~ condition
)
dds_gene <- DESeq(dds_gene)
saveRDS(dds_gene, file = file.path(output_dir, "dds_gene.rds"))

res_gene <- results(dds_gene, contrast = c("condition", "dt4", "ctrl"))
res_gene_df <- as.data.frame(res_gene[order(res_gene$padj), ])
res_gene_df$gene_id <- rownames(res_gene_df)
res_gene_df$significant <- ifelse(
  !is.na(res_gene_df$padj) & res_gene_df$padj < 0.05 & abs(res_gene_df$log2FoldChange) > 1,
  ifelse(res_gene_df$log2FoldChange > 1, "Up", "Down"),
  "NS"
)

write.csv(res_gene_df, file.path(output_dir, "DESeq2_gene_level_results.csv"), row.names = FALSE)

sig_gene_up <- sum(res_gene_df$significant == "Up", na.rm = TRUE)
sig_gene_down <- sum(res_gene_df$significant == "Down", na.rm = TRUE)
message(paste0("显著上调基因: ", sig_gene_up))
message(paste0("显著下调基因: ", sig_gene_down))

# ----- 转录本水平分析 -----
message("\n--- 转录本水平差异分析 ---")

transcript_counts_int <- round(transcript_counts_all)
mode(transcript_counts_int) <- "integer"

keep_tx <- rowSums(transcript_counts_int >= 10) >= 2
transcript_counts_filtered <- transcript_counts_int[keep_tx, ]
message(paste0("过滤后转录本数量: ", nrow(transcript_counts_filtered)))

dds_tx <- DESeqDataSetFromMatrix(
  countData = transcript_counts_filtered,
  colData = sample_info,
  design = ~ condition
)
dds_tx <- DESeq(dds_tx)
saveRDS(dds_tx, file = file.path(output_dir, "dds_transcript.rds"))

res_tx <- results(dds_tx, contrast = c("condition", "dt4", "ctrl"))
res_tx_df <- as.data.frame(res_tx[order(res_tx$padj), ])
res_tx_df$transcript_id <- rownames(res_tx_df)
res_tx_df$significant <- ifelse(
  !is.na(res_tx_df$padj) & res_tx_df$padj < 0.05 & abs(res_tx_df$log2FoldChange) > 1,
  ifelse(res_tx_df$log2FoldChange > 1, "Up", "Down"),
  "NS"
)

write.csv(res_tx_df, file.path(output_dir, "DESeq2_transcript_level_results.csv"), row.names = FALSE)

sig_tx_up <- sum(res_tx_df$significant == "Up", na.rm = TRUE)
sig_tx_down <- sum(res_tx_df$significant == "Down", na.rm = TRUE)
message(paste0("显著上调转录本: ", sig_tx_up))
message(paste0("显著下调转录本: ", sig_tx_down))

# ============================================================
# 第七部分：可视化
# ============================================================

message("\n=== 生成可视化图表 ===")

# PCA图
vsd <- vst(dds_gene, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 5) +
  geom_text(vjust = -1.5, size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  ggtitle("PCA Plot - Gene Level") +
  scale_color_manual(values = c("ctrl" = "#2166AC", "dt4" = "#B2182B"))

ggsave(file.path(output_dir, "PCA_plot.pdf"), pca_plot, width = 8, height = 6)
ggsave(file.path(output_dir, "PCA_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

# 火山图 - 基因
volcano_gene <- ggplot(res_gene_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC", "NS" = "gray70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot - Gene Level (dt4 vs ctrl)",
    subtitle = paste0("Up: ", sig_gene_up, " | Down: ", sig_gene_down),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Volcano_plot_gene.pdf"), volcano_gene, width = 8, height = 6)
ggsave(file.path(output_dir, "Volcano_plot_gene.png"), volcano_gene, width = 8, height = 6, dpi = 300)

# 火山图 - 转录本
volcano_tx <- ggplot(res_tx_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC", "NS" = "gray70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot - Transcript Level (dt4 vs ctrl)",
    subtitle = paste0("Up: ", sig_tx_up, " | Down: ", sig_tx_down),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Volcano_plot_transcript.pdf"), volcano_tx, width = 8, height = 6)
ggsave(file.path(output_dir, "Volcano_plot_transcript.png"), volcano_tx, width = 8, height = 6, dpi = 300)

# 热图
sig_genes <- res_gene_df[res_gene_df$significant != "NS" & !is.na(res_gene_df$padj), ]
if(nrow(sig_genes) >= 2) {
  top_n <- min(50, nrow(sig_genes))
  top_genes <- head(sig_genes[order(sig_genes$padj), "gene_id"], top_n)
  mat <- assay(vsd)[top_genes[top_genes %in% rownames(assay(vsd))], ]
  
  if(nrow(mat) >= 2) {
    mat_scaled <- t(scale(t(mat)))
    annotation_col <- data.frame(Condition = sample_info$condition, row.names = sample_names)
    
    pdf(file.path(output_dir, "Heatmap_top_genes.pdf"), width = 8, height = 12)
    pheatmap(mat_scaled, annotation_col = annotation_col, 
             show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE,
             main = paste0("Top ", nrow(mat), " DEGs"))
    dev.off()
  }
}

# ============================================================
# 第八部分：生成报告
# ============================================================

message("\n=== 生成分析报告 ===")

summary_report <- paste0(
  "========================================\n",
  "DRS Bambu + DESeq2 分析报告\n",
  "========================================\n",
  "完成时间: ", Sys.time(), "\n\n",
  "【结果汇总】\n",
  "总转录本: ", nrow(transcript_counts_all), "\n",
  "总基因: ", nrow(gene_counts_all), "\n",
  "新转录本: ", novel_tx, "\n\n",
  "【差异表达 - 基因水平】\n",
  "上调: ", sig_gene_up, "\n",
  "下调: ", sig_gene_down, "\n\n",
  "【差异表达 - 转录本水平】\n",
  "上调: ", sig_tx_up, "\n",
  "下调: ", sig_tx_down, "\n",
  "========================================\n"
)

cat(summary_report)
writeLines(summary_report, file.path(output_dir, "analysis_summary.txt"))

message("\n=== 分析完成! ===")
message(paste0("结果保存在: ", output_dir))
