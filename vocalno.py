#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import functools
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 1. 参数设置 (宽松模式)
# ============================================================

FILES = {
    'c1': '/data2/users/zhouyiming/DT4_min10/leaf/test/ctrlL_1/c1_m6A_strict.bed',  # 请修改路径
    'c2': '/data2/users/zhouyiming/DT4_min10/leaf/test/ctrlL_2/c2_m6A_strict.bed',
    't1': '/data2/users/zhouyiming/DT4_min10/leaf/test/dt4L_1/t1_m6A_strict.bed',
    't2': '/data2/users/zhouyiming/DT4_min10/leaf/test/dt4L_2/t2_m6A_strict.bed'
}

# --- 降低阈值以获得更多位点 ---
DIFF_CUTOFF = 5        # 差异阈值降为 5% (原为10)
PVALUE_CUTOFF = 0.05   # p < 0.05
MIN_GROUP_COV = 10     # 组内总覆盖度要求 (例如 c1+c2 >= 10 即可，不再要求单样本)
Y_AXIS_MAX = 50        # 火山图Y轴最大值

# ============================================================
# 2. 数据读取
# ============================================================

def load_m5c_data(name, filepath):
    print(f"正在读取 {name}: {filepath}")
    try:
        df = pd.read_csv(filepath, sep='\t')
        
        # 构建ID
        df['site_id'] = df['chrom'] + ":" + df['chromStart'].astype(str)
        
        # 只提取 Count 数据 (百分比我们后面自己重算，这样更准确)
        out_cols = {
            'site_id': 'site_id',
            'count_modified': f'{name}_m',
            'count_canonical': f'{name}_u'
        }
        
        # 这里不做 MIN_COVERAGE 过滤，推迟到合并后按组过滤
        df = df[list(out_cols.keys())].rename(columns=out_cols)
        return df
        
    except Exception as e:
        print(f"  读取失败: {e}")
        return None

print("=" * 60)
print("1. 读取并合并数据 (Outer Merge)")
print("=" * 60)

dfs = []
for name in ['c1', 'c2', 't1', 't2']:
    if name in FILES:
        res = load_m5c_data(name, FILES[name])
        if res is not None:
            dfs.append(res)

if len(dfs) == 0:
    print("错误: 未读取到数据")
    exit()

# 使用 outer merge (并集)，保留只要在任意样本中出现的位点
merged = functools.reduce(lambda left, right: pd.merge(left, right, on='site_id', how='outer'), dfs)

# 将缺失值填充为 0 (表示该样本在该位点没有read覆盖)
merged = merged.fillna(0)

print(f"原始合并位点总数 (含低覆盖度): {len(merged)}")

# ============================================================
# 3. 数据聚合与过滤
# ============================================================

print("\n" + "=" * 60)
print("2. 聚合 Counts 并过滤")
print("=" * 60)

# 聚合 Control 组 (c1 + c2)
merged['ctrl_m_total'] = merged['c1_m'] + merged['c2_m']
merged['ctrl_u_total'] = merged['c1_u'] + merged['c2_u']
merged['ctrl_cov_total'] = merged['ctrl_m_total'] + merged['ctrl_u_total']

# 聚合 Treatment 组 (t1 + t2)
merged['treat_m_total'] = merged['t1_m'] + merged['t2_m']
merged['treat_u_total'] = merged['t1_u'] + merged['t2_u']
merged['treat_cov_total'] = merged['treat_m_total'] + merged['treat_u_total']

# --- 关键过滤步骤 ---
# 只要 "对照组总和" 和 "处理组总和" 够多就算有效，不再要求每个重复都存在
valid_mask = (merged['ctrl_cov_total'] >= MIN_GROUP_COV) & \
             (merged['treat_cov_total'] >= MIN_GROUP_COV)

filtered = merged[valid_mask].copy()

print(f"过滤后用于统计的位点数: {len(filtered)} (Group Coverage >= {MIN_GROUP_COV})")

if len(filtered) == 0:
    print("错误: 过滤后没有剩余位点，请检查数据或降低覆盖度阈值。")
    exit()

# ============================================================
# 4. 计算差异与统计检验
# ============================================================

# 计算加权平均修饰率 (Weighted Percentage)
# 这种算法对覆盖度不均更鲁棒
filtered['ctrl_pct'] = (filtered['ctrl_m_total'] / filtered['ctrl_cov_total']) * 100
filtered['treat_pct'] = (filtered['treat_m_total'] / filtered['treat_cov_total']) * 100

filtered['diff_pct'] = filtered['treat_pct'] - filtered['ctrl_pct']

# Fisher's Exact Test
print("正在进行统计检验...")
pvalues = []
for idx, row in filtered.iterrows():
    # 2x2 列联表
    table = [
        [int(row['ctrl_m_total']), int(row['ctrl_u_total'])],
        [int(row['treat_m_total']), int(row['treat_u_total'])]
    ]
    try:
        _, p = stats.fisher_exact(table, alternative='two-sided')
        pvalues.append(p)
    except:
        pvalues.append(1.0)

filtered['pvalue'] = pvalues
filtered['neg_log10_pvalue'] = -np.log10(filtered['pvalue'])

# 处理 p=0 的无穷大问题
max_val = filtered.loc[filtered['neg_log10_pvalue'] != np.inf, 'neg_log10_pvalue'].max()
if pd.isna(max_val): max_val = 50
filtered['neg_log10_pvalue'] = filtered['neg_log10_pvalue'].replace(np.inf, max_val * 1.1)

# ============================================================
# 5. 结果分类与保存
# ============================================================

def classify(row):
    if pd.isna(row['pvalue']): return 'NS'
    if row['pvalue'] < PVALUE_CUTOFF:
        if row['diff_pct'] > DIFF_CUTOFF: return 'Up'
        if row['diff_pct'] < -DIFF_CUTOFF: return 'Down'
    return 'NS'

filtered['change'] = filtered.apply(classify, axis=1)

# 保存
out_cols = ['site_id', 'ctrl_pct', 'treat_pct', 'diff_pct', 'pvalue', 'change', 
            'ctrl_cov_total', 'treat_cov_total']
filtered.sort_values('pvalue').to_csv('m5C_diff_results.txt', sep='\t', index=False, float_format='%.4f')

sig_df = filtered[filtered['change'] != 'NS']
sig_df[out_cols].to_csv('m5C_diff_significant.txt', sep='\t', index=False, float_format='%.4f')

print("\n结果统计:")
print(filtered['change'].value_counts())
print(f"显著差异位点已保存至: m5C_diff_significant.txt")

# ============================================================
# 6. 绘制火山图
# ============================================================

plot_df = filtered.dropna(subset=['diff_pct', 'neg_log10_pvalue']).copy()

if len(plot_df) > 0:
    plt.figure(figsize=(8, 7))
    
    # Y轴截断
    plot_df['plot_y'] = plot_df['neg_log10_pvalue'].clip(upper=Y_AXIS_MAX)
    
    colors = {'Down': '#377eb8', 'NS': '#bdbdbd', 'Up': '#e41a1c'}
    alphas = {'Down': 0.6, 'NS': 0.15, 'Up': 0.6}
    
    for cat in ['NS', 'Down', 'Up']:
        subset = plot_df[plot_df['change'] == cat]
        if len(subset) > 0:
            plt.scatter(
                subset['diff_pct'],
                subset['plot_y'],
                c=colors[cat],
                label=f"{cat} ({len(subset)})",
                alpha=alphas[cat],
                s=10 if cat == 'NS' else 25, # 稍微加大显著点的尺寸
                edgecolors='none'
            )
            
    # 阈值线
    plt.axhline(-np.log10(PVALUE_CUTOFF), color='gray', linestyle='--', linewidth=1)
    plt.axvline(DIFF_CUTOFF, color='gray', linestyle='--', linewidth=1)
    plt.axvline(-DIFF_CUTOFF, color='gray', linestyle='--', linewidth=1)
    
    # 标注阈值
    plt.text(x=DIFF_CUTOFF, y=Y_AXIS_MAX*0.98, s=f' Diff={DIFF_CUTOFF}%', 
             color='gray', ha='left', va='top', rotation=90, fontsize=8)
    
    # 坐标轴范围
    plt.ylim(0, Y_AXIS_MAX)
    # X轴范围自适应，但至少显示+/-20
    x_limit = max(20, abs(plot_df['diff_pct'].max()), abs(plot_df['diff_pct'].min())) * 1.1
    plt.xlim(-x_limit, x_limit)
    
    plt.xlabel('Difference in Methylation Level (%)', fontsize=12)
    plt.ylabel(r'$-Log_{10}(P$-value$)$', fontsize=12)
    plt.title('m5C Methylation Volcano Plot (Relaxed)', fontsize=14)
    plt.legend(loc='upper right', frameon=False)
    
    plt.tight_layout()
    plt.savefig('m5C_volcano.pdf')
    plt.savefig('m5C_volcano.png', dpi=300)
    print("火山图已保存: m5C_volcano.png")

else:
    print("无数据绘图")
