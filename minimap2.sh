#!/bin/bash
# ============================================================
# DRS数据分析流程 - Part 1: minimap2比对到基因组
# ============================================================

# ==================== 配置参数 ====================
GENOME="/data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/genome/osa_IRGSP_1.fa"
GTF="/data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/annotation/osa_IRGSP_1.annotation.gtf"
WORKDIR="/data2/users/zhouyiming/DT4_ROOT/analysis"
THREADS=16

# 创建工作目录
mkdir -p ${WORKDIR}
cd ${WORKDIR}

# ==================== FASTQ文件路径 ====================
declare -A FASTQ_FILES
FASTQ_FILES=(
    ["ctrlR_1"]="/data2/users/zhouyiming/DT4_ROOT/root_ct1/drs-ct-R-1.pass.fq"
    ["ctrlR_2"]="/data2/users/zhouyiming/DT4_ROOT/root_ct2/drs-ct-R-2.pass.fq"
    ["dt4R_1"]="/data2/users/zhouyiming/DT4_ROOT/root_dt4_1/drs-dt4-R-1.pass.fq"
    ["dt4R_2"]="/data2/users/zhouyiming/DT4_ROOT/root_dt4-2/drs-dt4-R-2.pass.fq"
)

# ==================== Step 1: 创建基因组索引 ====================
echo "============================================================"
echo "Step 1: 创建minimap2基因组索引"
echo "时间: $(date)"
echo "============================================================"

if [ ! -f "${WORKDIR}/osa_genome.mmi" ]; then
    echo "创建索引中..."
    minimap2 -d ${WORKDIR}/osa_genome.mmi ${GENOME}
    echo "索引创建完成!"
else
    echo "索引已存在，跳过创建"
fi

# ==================== Step 2: 比对每个样本 ====================
echo ""
echo "============================================================"
echo "Step 2: 比对reads到基因组"
echo "============================================================"

for sample in ctrlR_1 ctrlR_2 dt4R_1 dt4R_2; do
    echo ""
    echo "--- 处理样本: $sample ---"
    echo "时间: $(date)"
    
    FASTQ="${FASTQ_FILES[$sample]}"
    OUTDIR="${WORKDIR}/aligned/${sample}"
    mkdir -p ${OUTDIR}
    
    # 检查FASTQ文件是否存在
    if [ ! -f "$FASTQ" ]; then
        echo "错误: 文件不存在 - $FASTQ"
        continue
    fi
    
    echo "输入文件: $FASTQ"
    echo "输出目录: $OUTDIR"
    
    # 统计输入reads数量
    echo "统计输入reads数..."
    INPUT_READS=$(wc -l < "$FASTQ")
    INPUT_READS=$((INPUT_READS / 4))
    echo "输入reads数: ${INPUT_READS}"
    
    # minimap2比对
    # -ax splice: RNA-seq模式，考虑剪接
    # -uf: 针对Nanopore direct RNA (正链)
    # -k14: k-mer大小
    # --secondary=no: 不输出次优比对
    echo "开始比对..."
    minimap2 -ax splice \
        -uf \
        -k14 \
        --secondary=no \
        -t ${THREADS} \
        ${WORKDIR}/osa_genome.mmi \
        ${FASTQ} \
        2> ${OUTDIR}/${sample}_minimap2.log | \
    samtools view -bS -@ ${THREADS} -F 2308 - | \
    samtools sort -@ ${THREADS} -o ${OUTDIR}/${sample}_genome_sort.bam
    
    # 创建索引
    echo "创建BAM索引..."
    samtools index ${OUTDIR}/${sample}_genome_sort.bam
    
    # 比对统计
    echo "比对统计:"
    samtools flagstat ${OUTDIR}/${sample}_genome_sort.bam > ${OUTDIR}/${sample}_flagstat.txt
    cat ${OUTDIR}/${sample}_flagstat.txt
    
    # 检查染色体命名
    echo "染色体命名检查:"
    samtools view -H ${OUTDIR}/${sample}_genome_sort.bam | grep "^@SQ" | head -5
    
    echo "样本 $sample 比对完成!"
done

# ==================== Step 3: 生成比对汇总报告 ====================
echo ""
echo "============================================================"
echo "Step 3: 生成比对汇总报告"
echo "============================================================"

REPORT="${WORKDIR}/alignment_summary.txt"
echo "========================================" > ${REPORT}
echo "Minimap2 比对汇总报告" >> ${REPORT}
echo "生成时间: $(date)" >> ${REPORT}
echo "========================================" >> ${REPORT}
echo "" >> ${REPORT}

for sample in ctrlR_1 ctrlR_2 dt4R_1 dt4R_2; do
    OUTDIR="${WORKDIR}/aligned/${sample}"
    if [ -f "${OUTDIR}/${sample}_flagstat.txt" ]; then
        echo "--- ${sample} ---" >> ${REPORT}
        TOTAL=$(grep "in total" ${OUTDIR}/${sample}_flagstat.txt | cut -d' ' -f1)
        MAPPED=$(grep "mapped (" ${OUTDIR}/${sample}_flagstat.txt | head -1 | cut -d' ' -f1)
        RATE=$(grep "mapped (" ${OUTDIR}/${sample}_flagstat.txt | head -1 | grep -oP '\(\K[^%]+')
        echo "总reads: ${TOTAL}" >> ${REPORT}
        echo "比对reads: ${MAPPED}" >> ${REPORT}
        echo "比对率: ${RATE}%" >> ${REPORT}
        echo "" >> ${REPORT}
    fi
done

cat ${REPORT}

echo ""
echo "============================================================"
echo "比对完成!"
echo "BAM文件位置: ${WORKDIR}/aligned/"
echo "============================================================"
echo ""
echo "下一步: 运行R脚本进行Bambu分析"
echo "Rscript ${WORKDIR}/bambu_deseq2_analysis.R"
