#!/bin/bash

# ============================================
# RNA 修饰位点和 PolyA 批量分析流程
# ============================================

# ---------- 参数设置 ----------
THREADS=8
REF="/data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/mRNA/osa_IRGSP_1.txdb.fa"
INPUT_DIR="/data1/liuqi/software/dorado/test/leaf"  # 包含BAM文件的目录
MAIN_OUTDIR="/data2/users/zhouyiming/DT4_min10/leaf/test"

# ---------- 创建主输出目录 ----------
mkdir -p ${MAIN_OUTDIR}
cd ${MAIN_OUTDIR}

# ============================================
# 批量处理所有BAM文件
# ============================================
for BAM_FILE in ${INPUT_DIR}/*.bam; do
    # 提取不带路径和扩展名的纯文件名
    SAMPLE=$(basename "${BAM_FILE}" .bam)
    
    # 创建样本专属输出目录
    SAMPLE_OUTDIR="${MAIN_OUTDIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUTDIR}"
    
    echo "========================================"
    echo "开始处理样本: ${SAMPLE}"
    echo "输出目录: ${SAMPLE_OUTDIR}"
    echo "========================================"
    
    cd "${SAMPLE_OUTDIR}"
    
    # ============================================
    # Step 1: BAM文件预处理 (排序和索引)
    # ============================================
    echo "[${SAMPLE}] 步骤1: 排序和索引BAM文件..."
    
    # 排序BAM文件
    samtools sort -@ ${THREADS} -o "${SAMPLE}_sort.bam" "${BAM_FILE}"
    
    # 创建索引
    samtools index "${SAMPLE}_sort.bam"
    
    # 验证修饰标签是否保留
    echo "[${SAMPLE}] 检查修饰标签..."
    samtools view "${SAMPLE}_sort.bam" | head -1 | grep -oE "(MM|ML|pt):[^[:space:]]+" | head -3
    
    # ============================================
    # Step 2: 提取 PolyA 长度
    # ============================================
    echo "[${SAMPLE}] 步骤2: 提取polyA长度..."
    
    samtools view "${SAMPLE}_sort.bam" | awk '{
        for(i=12; i<=NF; i++) {
            if($i ~ /^pt:i:/) {
                split($i, a, ":");
                print $1"\t"a[3]
            }
        }
    }' > "${SAMPLE}_polya_lengths.txt"
    
    # PolyA 统计
    echo "[${SAMPLE}] PolyA 长度统计："
    awk '{sum+=$2; count++} END {
        print "  Read数量: " count;
        print "  平均长度: " (count ? sum/count : 0);
    }' "${SAMPLE}_polya_lengths.txt"
    
    # ============================================
    # Step 3: 提取修饰位点
    # ============================================
    echo "[${SAMPLE}] 步骤3: 提取修饰位点..."
    
    modkit pileup \
        --threads ${THREADS} \
        --ref ${REF} \
        --filter-threshold C:0.99 \
        --filter-threshold A:0.99 \
        --filter-threshold T:0.99 \
        --filter-threshold G:0.99 \
        --log-filepath "${SAMPLE}_modkit_pileup.log" \
        --with-header \
        "${SAMPLE}_sort.bam" \
        "${SAMPLE}_modifications.bed"
    
    # ============================================
    # Step 4: 结果统计
    # ============================================
    echo "[${SAMPLE}] 步骤4: 生成统计摘要..."
    
    echo "========== ${SAMPLE} 分析结果统计 ==========" | tee "${SAMPLE}_summary.txt"
    echo "总 reads 数: $(samtools view -c "${SAMPLE}_sort.bam")" | tee -a "${SAMPLE}_summary.txt"
    echo "比对 reads 数: $(samtools view -c -F 4 "${SAMPLE}_sort.bam")" | tee -a "${SAMPLE}_summary.txt"
    echo "PolyA reads 数: $(wc -l < "${SAMPLE}_polya_lengths.txt")" | tee -a "${SAMPLE}_summary.txt"
    echo "总修饰位点数: $(grep -v "^#" "${SAMPLE}_modifications.bed" | wc -l)" | tee -a "${SAMPLE}_summary.txt"
    
    echo "[${SAMPLE}] 处理完成！输出文件在: ${SAMPLE_OUTDIR}"
    cd "${MAIN_OUTDIR}"
done

echo "所有样本处理完成！"
