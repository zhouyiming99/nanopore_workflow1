#!/bin/bash
# ============================================
# 批量转录本级别 PolyA 长度统计脚本
# ============================================

# ---------- 参数设置 ----------
MAIN_OUTDIR="/data2/users/zhouyiming/DT4_min10/leaf/test"  # 主输出目录

# ============================================
# 函数：处理单个样本
# ============================================
process_sample() {
    SAMPLE_DIR="$1"
    SAMPLE=$(basename "${SAMPLE_DIR}")
    BAM="${SAMPLE_DIR}/${SAMPLE}_sort.bam"
    OUTPUT_PREFIX="${SAMPLE_DIR}/${SAMPLE}"
    
    echo "========================================"
    echo "正在处理样本: ${SAMPLE}"
    echo "BAM 文件: ${BAM}"
    echo "输出前缀: ${OUTPUT_PREFIX}"
    echo "========================================"
    
    # 检查输入文件
    if [ ! -f "${BAM}" ]; then
        echo "错误: BAM文件不存在: ${BAM}"
        return 1
    fi

    # ============================================
    # Step 1: 从 BAM 提取 read 级别 polyA 信息
    # ============================================
    echo "Step 1: 提取 PolyA 信息..."
    
    samtools view ${BAM} | awk '{
        transcript_id = $3
        for(i=12; i<=NF; i++) {
            if($i ~ /^pt:i:/) {
                split($i, a, ":")
                if(transcript_id != "*") {
                    print transcript_id"\t"a[3]
                }
            }
        }
    }' > ${OUTPUT_PREFIX}_polya_tmp.txt

    echo "  提取完成，共 $(wc -l < ${OUTPUT_PREFIX}_polya_tmp.txt) 条记录"

    # ============================================
    # Step 2: 按转录本汇总统计
    # ============================================
    echo "Step 2: 按转录本汇总统计..."
    
    echo -e "transcript_id\tread_count\tmean_polya\tmedian_polya\tmin_polya\tmax_polya\tstd_polya" > ${OUTPUT_PREFIX}_polya_per_transcript.txt

    # 使用 awk 计算统计量
    sort -k1,1 ${OUTPUT_PREFIX}_polya_tmp.txt | awk '
    {
        if($1 != prev && NR > 1) {
            # 输出上一个转录本的统计
            n = count
            mean = sum / n
            
            # 计算标准差
            variance = (sumsq / n) - (mean * mean)
            std = sqrt(variance > 0 ? variance : 0)
            
            # 排序计算中位数
            asort(values)
            if(n % 2 == 1) {
                median = values[int(n/2) + 1]
            } else {
                median = (values[n/2] + values[n/2 + 1]) / 2
            }
            
            printf "%s\t%d\t%.2f\t%.2f\t%d\t%d\t%.2f\n", prev, n, mean, median, min_val, max_val, std
            
            # 重置
            delete values
            count = 0
            sum = 0
            sumsq = 0
            min_val = ""
            max_val = ""
        }
        
        prev = $1
        count++
        sum += $2
        sumsq += $2 * $2
        values[count] = $2
        
        if(min_val == "" || $2 < min_val) min_val = $2
        if(max_val == "" || $2 > max_val) max_val = $2
    }
    END {
        # 输出最后一个转录本
        n = count
        mean = sum / n
        variance = (sumsq / n) - (mean * mean)
        std = sqrt(variance > 0 ? variance : 0)
        
        asort(values)
        if(n % 2 == 1) {
            median = values[int(n/2) + 1]
        } else {
            median = (values[n/2] + values[n/2 + 1]) / 2
        }
        
        printf "%s\t%d\t%.2f\t%.2f\t%d\t%d\t%.2f\n", prev, n, mean, median, min_val, max_val, std
    }' >> ${OUTPUT_PREFIX}_polya_per_transcript.txt

    # 按 read_count 降序排序（保留表头）
    head -1 ${OUTPUT_PREFIX}_polya_per_transcript.txt > ${OUTPUT_PREFIX}_polya_per_transcript_sorted.txt
    tail -n +2 ${OUTPUT_PREFIX}_polya_per_transcript.txt | sort -t$'\t' -k2 -rn >> ${OUTPUT_PREFIX}_polya_per_transcript_sorted.txt
    mv ${OUTPUT_PREFIX}_polya_per_transcript_sorted.txt ${OUTPUT_PREFIX}_polya_per_transcript.txt

    # 清理临时文件
    rm -f ${OUTPUT_PREFIX}_polya_tmp.txt

    # ============================================
    # Step 3: 统计汇总
    # ============================================
    echo ""
    echo "========== ${SAMPLE} 统计汇总 =========="

    TOTAL_TRANSCRIPTS=$(tail -n +2 ${OUTPUT_PREFIX}_polya_per_transcript.txt | wc -l)
    TOTAL_READS=$(awk 'NR>1 {sum+=$2} END {print sum}' ${OUTPUT_PREFIX}_polya_per_transcript.txt)
    AVG_POLYA=$(awk 'NR>1 {sum+=$3*$2; total+=$2} END {printf "%.2f", sum/total}' ${OUTPUT_PREFIX}_polya_per_transcript.txt)

    echo "  转录本数量: ${TOTAL_TRANSCRIPTS}"
    echo "  总 reads 数: ${TOTAL_READS}"
    echo "  加权平均 PolyA 长度: ${AVG_POLYA}"
    echo ""

    # ============================================
    # Step 4: 输出预览
    # ============================================
    echo "========== ${SAMPLE} 输出文件预览 =========="
    echo "文件: ${OUTPUT_PREFIX}_polya_per_transcript.txt"
    echo ""
    echo "Top 15 转录本（按 read 数排序）:"
    head -16 ${OUTPUT_PREFIX}_polya_per_transcript.txt | column -t -s $'\t'

    echo ""
    echo "完成样本 ${SAMPLE} 的处理！"
}

# ============================================
# 主循环：遍历所有样本目录
# ============================================
cd ${MAIN_OUTDIR}

# 查找所有样本目录（名称如 ctrlL_1, ctrlL_2 等）
for SAMPLE_DIR in ${MAIN_OUTDIR}/*; do
    if [ -d "${SAMPLE_DIR}" ]; then
        process_sample "${SAMPLE_DIR}"
    fi
done

echo "所有样本处理完成！"
