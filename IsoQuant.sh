#!/bin/bash
# ============================================================
# IsoQuant 转录本定量分析流程
# ============================================================

# ==================== 配置参数 ====================
GENOME="/data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/genome/osa_IRGSP_1.fa"
GTF="/data2/users/zhouyiming/DT4_ROOT/analysis/isoquant_results_1/all_samples/osa_IRGSP_1.annotation.corrected.gtf"
WORKDIR="/data2/users/zhouyiming/DT4_ROOT/analysis"
THREADS=16

# 输出目录
OUTDIR="${WORKDIR}/isoquant_results"
mkdir -p ${OUTDIR}

cd ${WORKDIR}

# ==================== BAM文件和样本名称 ====================
declare -A BAM_FILES
BAM_FILES=(
    ["ctrlR_1"]="${WORKDIR}/aligned/ctrlR_1/ctrlR_1_genome_sort.bam"
    ["ctrlR_2"]="${WORKDIR}/aligned/ctrlR_2/ctrlR_2_genome_sort.bam"
    ["dt4R_1"]="${WORKDIR}/aligned/dt4R_1/dt4R_1_genome_sort.bam"
    ["dt4R_2"]="${WORKDIR}/aligned/dt4R_2/dt4R_2_genome_sort.bam"
)

# ==================== Step 1: 对每个样本单独运行IsoQuant ====================
echo "============================================================"
echo "Step 1: 对每个样本运行IsoQuant"
echo "时间: $(date)"
echo "============================================================"

for sample in ctrlR_1 ctrlR_2 dt4R_1 dt4R_2; do
    echo ""
    echo "--- 处理样本: $sample ---"
    echo "时间: $(date)"
    
    BAM="${BAM_FILES[$sample]}"
    SAMPLE_OUT="${OUTDIR}/${sample}"
    
    # 检查BAM文件是否存在
    if [ ! -f "$BAM" ]; then
        echo "错误: BAM文件不存在 - $BAM"
        continue
    fi
    
    # 检查是否已经运行过
    if [ -f "${SAMPLE_OUT}/${sample}/${sample}.transcript_counts.tsv" ]; then
        echo "样本 $sample 已处理，跳过..."
        continue
    fi
    
    echo "输入BAM: $BAM"
    echo "输出目录: $SAMPLE_OUT"
    
    # 运行IsoQuant
    isoquant.py \
        --reference ${GENOME} \
        --genedb ${GTF} \
        --bam ${BAM} \
        --data_type nanopore \
        --stranded forward \
        --threads ${THREADS} \
        --labels ${sample} \
        --prefix ${sample} \
        -o ${SAMPLE_OUT} \
        2>&1 | tee ${SAMPLE_OUT}_isoquant.log
    
    # 检查是否成功
    if [ $? -eq 0 ]; then
        echo "✓ 样本 $sample 完成!"
    else
        echo "✗ 样本 $sample 失败!"
    fi
done

# ==================== Step 2: 合并所有样本运行IsoQuant（可选）====================
echo ""
echo "============================================================"
echo "Step 2: 合并所有样本运行IsoQuant"
echo "时间: $(date)"
echo "============================================================"

# 准备BAM文件列表
BAM_LIST=""
LABEL_LIST=""
for sample in ctrlR_1 ctrlR_2 dt4R_1 dt4R_2; do
    BAM="${BAM_FILES[$sample]}"
    if [ -f "$BAM" ]; then
        BAM_LIST="${BAM_LIST} ${BAM}"
        LABEL_LIST="${LABEL_LIST} ${sample}"
    fi
done

# 去除开头的空格
BAM_LIST=$(echo $BAM_LIST | xargs)
LABEL_LIST=$(echo $LABEL_LIST | xargs)

echo "BAM文件: $BAM_LIST"
echo "样本标签: $LABEL_LIST"

# 合并运行
MERGED_OUT="${OUTDIR}/all_samples"
mkdir -p ${MERGED_OUT}

isoquant.py \
    --reference ${GENOME} \
    --genedb ${GTF} \
    --bam ${BAM_LIST} \
    --data_type nanopore \
    --stranded forward \
    --threads ${THREADS} \
    --labels ${LABEL_LIST} \
    --prefix all_samples \
    -o ${MERGED_OUT} \
    2>&1 | tee ${MERGED_OUT}_isoquant.log

# ==================== Step 3: 整理输出结果 ====================
echo ""
echo "============================================================"
echo "Step 3: 整理输出结果"
echo "============================================================"

# 创建汇总目录
SUMMARY_DIR="${OUTDIR}/summary"
mkdir -p ${SUMMARY_DIR}

# 复制各样本的counts文件到汇总目录
for sample in ctrlR_1 ctrlR_2 dt4R_1 dt4R_2; do
    SAMPLE_OUT="${OUTDIR}/${sample}/${sample}"
    
    if [ -d "$SAMPLE_OUT" ]; then
        # 复制转录本counts
        if [ -f "${SAMPLE_OUT}/${sample}.transcript_counts.tsv" ]; then
            cp ${SAMPLE_OUT}/${sample}.transcript_counts.tsv ${SUMMARY_DIR}/${sample}_transcript_counts.tsv
        fi
        
        # 复制基因counts
        if [ -f "${SAMPLE_OUT}/${sample}.gene_counts.tsv" ]; then
            cp ${SAMPLE_OUT}/${sample}.gene_counts.tsv ${SUMMARY_DIR}/${sample}_gene_counts.tsv
        fi
        
        # 复制TPM文件
        if [ -f "${SAMPLE_OUT}/${sample}.transcript_tpm.tsv" ]; then
            cp ${SAMPLE_OUT}/${sample}.transcript_tpm.tsv ${SUMMARY_DIR}/${sample}_transcript_tpm.tsv
        fi
        
        echo "已复制 $sample 的结果文件"
    fi
done

# 复制合并运行的结果
if [ -d "${MERGED_OUT}/all_samples" ]; then
    cp ${MERGED_OUT}/all_samples/*.tsv ${SUMMARY_DIR}/ 2>/dev/null
    echo "已复制合并分析的结果文件"
fi

# ==================== Step 4: 生成文件列表 ====================
echo ""
echo "============================================================"
echo "输出文件列表"
echo "============================================================"

echo ""
echo "各样本单独结果:"
ls -la ${OUTDIR}/*/

echo ""
echo "汇总目录:"
ls -la ${SUMMARY_DIR}/

echo ""
echo "============================================================"
echo "IsoQuant分析完成!"
echo "时间: $(date)"
echo "============================================================"
echo ""
echo "下一步: 运行DESeq2差异分析"
echo "Rscript ${WORKDIR}/isoquant_deseq2_analysis.R"
