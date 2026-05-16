#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Y轴范围设置示例
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

# ============================================
# 1. 数据准备（同前）
# ============================================
sns.set_style("whitegrid")
plt.rcParams.update({
    'font.family': 'Arial',
    'axes.linewidth': 1.5,
    'axes.edgecolor': 'black',
})

# 设置文件路径
files = {
    'C1': "/data2/users/zhouyiming/DT4_min10/leaf/test/ctrlL_1/c1_m6A_strict.bed",
    'C2': "/data2/users/zhouyiming/DT4_min10/leaf/test/ctrlL_2/c2_m6A_strict.bed",
    'T1': "/data2/users/zhouyiming/DT4_min10/leaf/test/dt4L_1/t1_m6A_strict.bed",
    'T2': "/data2/users/zhouyiming/DT4_min10/leaf/test/dt4L_2/t2_m6A_strict.bed"
}

column_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", 
                "thickStart", "thickEnd", "color", "valid_coverage", 
                "percent_modified", "count_modified", "count_canonical", 
                "count_other_mode", "count_delete", "count_fail", 
                "count_diff", "count_nocall"]

def read_data(file_path, sample):
    try:
        data = pd.read_csv(file_path, sep='\t', header=None, names=column_names,
                          on_bad_lines='skip')
    except TypeError:
        data = pd.read_csv(file_path, sep='\t', header=None, names=column_names,
                          error_bad_lines=False, warn_bad_lines=True)
    
    data['percent_modified'] = pd.to_numeric(data['percent_modified'], errors='coerce')
    data['valid_coverage'] = pd.to_numeric(data['valid_coverage'], errors='coerce')
    data = data[data['valid_coverage'] >= 10].copy()
    data['Sample'] = sample
    return data

# 读取数据
c1_df = read_data(files['C1'], 'C1')
c2_df = read_data(files['C2'], 'C2')
t1_df = read_data(files['T1'], 'T1')
t2_df = read_data(files['T2'], 'T2')

combined_data = pd.concat([c1_df, c2_df, t1_df, t2_df], ignore_index=True)

c1_data = c1_df['percent_modified'].dropna()
c2_data = c2_df['percent_modified'].dropna()
t1_data = t1_df['percent_modified'].dropna()
t2_data = t2_df['percent_modified'].dropna()

# ============================================
# 2. Y轴范围设置参数 ⭐⭐⭐
# ============================================
# 使用数据的实际范围，加上一定的边距
Y_MIN = 0          # Y轴最小值
Y_MAX = 150        # Y轴最大值（百分比数据通常0-100）
# ========== 方法1: 手动固定范围 ==========
# Y_MIN = 0          # Y轴最小值
# Y_MAX = 100        # Y轴最大值（百分比数据通常0-100）

# ========== 方法2: 基于数据自动计算 ==========
# 使用数据的实际范围，加上一定的边距
# y_min_auto = combined_data['percent_modified'].min()
# y_max_auto = combined_data['percent_modified'].max()
# Y_MIN = max(0, y_min_auto - (y_max_auto - y_min_auto) * 0.05)
# Y_MAX = y_max_auto * 1.15

# ========== 方法3: 基于分位数（去除极端值） ==========
# Y_MIN = combined_data['percent_modified'].quantile(0.01)  # 1%分位数
# Y_MAX = combined_data['percent_modified'].quantile(0.99)  # 99%分位数

# ========== 方法4: 基于IQR（四分位距） ==========
# Q1 = combined_data['percent_modified'].quantile(0.25)
# Q3 = combined_data['percent_modified'].quantile(0.75)
# IQR = Q3 - Q1
# Y_MIN = max(0, Q1 - 1.5 * IQR)
# Y_MAX = Q3 + 1.5 * IQR

# 选择使用哪种方法（这里使用方法1：固定范围）
print(f"Y轴范围设置: {Y_MIN} - {Y_MAX}")

# ============================================
# 3. 统计检验
# ============================================

def format_pval(p):
    if p < 0.0001:
        return "p < 0.0001"
    elif p < 0.001:
        return f"p = {p:.2e}"
    else:
        return f"p = {p:.4f}"

def get_sig_symbol(p):
    if p < 0.0001: return "****"
    elif p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return "ns"

pvals = {
    'C1 vs T1': stats.mannwhitneyu(c1_data, t1_data, alternative='two-sided')[1],
    'C1 vs T2': stats.mannwhitneyu(c1_data, t2_data, alternative='two-sided')[1],
    'C2 vs T1': stats.mannwhitneyu(c2_data, t1_data, alternative='two-sided')[1],
    'C2 vs T2': stats.mannwhitneyu(c2_data, t2_data, alternative='two-sided')[1],
}

# ============================================
# 4. 绑制箱线图
# ============================================

plot_data = combined_data[['percent_modified', 'Sample']].dropna()
plot_data.columns = ['Methylation', 'Sample']

fig, ax = plt.subplots(figsize=(10, 8))

palette = {
    'C1': '#2c7bb6', 'C2': '#74a9cf',
    'T1': '#d7191c', 'T2': '#f4a582'
}

sns.boxplot(
    x='Sample', y='Methylation', data=plot_data,
    order=['C1', 'C2', 'T1', 'T2'],
    palette=palette, width=0.65, showfliers=False,
    boxprops={'linewidth': 2, 'edgecolor': 'black'},
    whiskerprops={'linewidth': 1.5, 'linestyle': '--'},
    medianprops={'linewidth': 2.5, 'color': '#ffd700'},
    ax=ax
)

ax.axvline(1.5, color='gray', linestyle='--', linewidth=2, alpha=0.7)

# ============================================
# 5. 基于固定Y轴范围计算p值标注位置 ⭐⭐⭐
# ============================================

y_range = Y_MAX - Y_MIN

# p值标注的基础高度（从Y_MAX往上留出空间）
# 这里使用Y_MAX的85%作为第一条线的位置
base_y = Y_MIN + y_range * 0.72   # 第一条线位置
step_y = y_range * 0.07           # 每条线间距

# 4组比较参数
comparisons = [
    (0, 2, 'C1-T1', 'C1 vs T1'),
    (1, 3, 'C2-T2', 'C2 vs T2'),
    (0, 3, 'C1-T2', 'C1 vs T2'),
    (1, 2, 'C2-T1', 'C2 vs T1'),
]

for i, (x1, x2, label, key) in enumerate(comparisons):
    y_line = base_y + i * step_y
    p = pvals[key]
    sig = get_sig_symbol(p)
    
    text_color = '#c00000' if p < 0.05 else '#666666'
    lw = 1.5 if p < 0.05 else 1.0
    
    bracket_drop = y_range * 0.012
    ax.plot([x1, x1, x2, x2], 
            [y_line - bracket_drop, y_line, y_line, y_line - bracket_drop], 
            color='black', lw=lw, solid_capstyle='round')
    
    mid_x = (x1 + x2) / 2
    pval_text = f"{label}: {format_pval(p)} {sig}"
    ax.text(mid_x, y_line + y_range * 0.012, pval_text,
            ha='center', va='bottom', fontsize=9, fontweight='bold',
            color=text_color)

# ============================================
# 6. 设置Y轴范围 ⭐⭐⭐
# ============================================

# 设置固定的Y轴范围（留出顶部空间给p值标注）
ax.set_ylim(Y_MIN, Y_MAX)

# 设置Y轴刻度（可选：自定义刻度间隔）
# ax.set_yticks(np.arange(Y_MIN, Y_MAX + 1, 10))  # 每10个单位一个刻度

# ============================================
# 7. 美化图表
# ============================================

ax.text(0.5, Y_MIN - y_range * 0.06, 'Control', ha='center', fontsize=13, 
        fontweight='bold', color='#2c7bb6')
ax.text(2.5, Y_MIN - y_range * 0.06, 'Treatment', ha='center', fontsize=13, 
        fontweight='bold', color='#d7191c')

ax.set_xlabel('')
ax.set_ylabel('m6A Methylation Level (%)', fontsize=13, fontweight='bold')
ax.set_title('m6A Methylation: 4 Pairwise Comparisons (Mann-Whitney U Test)', 
             fontsize=14, fontweight='bold', pad=10)
ax.tick_params(axis='both', labelsize=11)
ax.set_xticklabels(['Control 1', 'Control 2', 'Treatment 1', 'Treatment 2'])
ax.grid(True, axis='y', linestyle='--', alpha=0.3)
ax.set_facecolor('#f8f8f8')

control_patch = mpatches.Patch(color='#2c7bb6', label='Control')
treatment_patch = mpatches.Patch(color='#d7191c', label='Treatment')
ax.legend(handles=[control_patch, treatment_patch], loc='lower right', fontsize=10)

plt.tight_layout()
plt.savefig('m6A_methylation_4comparisons.png', dpi=600, bbox_inches='tight')
plt.savefig('m6A_methylation_4comparisons.pdf', bbox_inches='tight')
plt.close()

print("图片已保存!")
