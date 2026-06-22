import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt

# =============================================================================
# 任务目标：对比 CellChat 和 LIANA（你的方法）在 LR 对、source-target 分布、
#           以及 prob 与 spatial_communication_score 的相关性，验证方法稳健性
# =============================================================================

print("=== 开始加载数据 ===")

# 步骤1：读取三个 CSV 文件
cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
liana_df    = pd.read_csv("liana_consensus_lr_pairs.csv")
adaptive_df = pd.read_csv("adaptive_distance_spatial_scores.csv")

print(f"CellChat 行数: {len(cellchat_df)}")
print(f"LIANA 行数: {len(liana_df)}")
print(f"Adaptive 行数: {len(adaptive_df)}")

# 步骤2：统一 LR 对标识（ligand_receptor），方便后续比较
cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']
liana_df['lr_pair']    = liana_df['ligand'] + '_' + liana_df['receptor']
adaptive_df['lr_pair'] = adaptive_df['ligand'] + '_' + adaptive_df['receptor']

print("\n=== 步骤3：计算 LR 对集合的重叠情况 ===")

cellchat_lrs = set(cellchat_df['lr_pair'].unique())
liana_lrs    = set(liana_df['lr_pair'].unique())
adaptive_lrs = set(adaptive_df['lr_pair'].unique())

intersection = cellchat_lrs.intersection(liana_lrs)
union = cellchat_lrs.union(liana_lrs)
jaccard = len(intersection) / len(union) if len(union) > 0 else 0
print(f"CellChat 与 LIANA 的 LR 对 Jaccard 相似度: {jaccard:.4f}")

print(f"CellChat 唯一 LR 对数: {len(cellchat_lrs)}")
print(f"LIANA 唯一 LR 对数: {len(liana_lrs)}")
print(f"Adaptive 唯一 LR 对数: {len(adaptive_lrs)}")
print(f"CellChat 与 LIANA 共同 LR 对数: {len(intersection)}")

cellchat_unique = cellchat_lrs - liana_lrs
liana_unique    = liana_lrs - cellchat_lrs
print("\nCellChat 独有 LR 对（前10）：", list(cellchat_unique)[:10])
print("LIANA 独有 LR 对（前10）：", list(liana_unique)[:10])

print("\n生物学解读示例（可手动补充）：")
print("- CellChat 独有示例：Glu-SLC17A6_GLS_Grm5 → 谷氨酸信号，兴奋性神经传递")
 
# 步骤4：对比 source → target 分布（修复：先转字符串）
print("\n=== 步骤4：通信方向分布对比 ===")

# 修复关键：把 source 和 target 转为字符串后再拼接
cellchat_df['st'] = cellchat_df['source'].astype(str) + ' → ' + cellchat_df['target'].astype(str)
liana_df['st']    = liana_df['source'].astype(str) + ' → ' + liana_df['target'].astype(str)

cellchat_st_counts = Counter(cellchat_df['st'])
liana_st_counts    = Counter(liana_df['st'])

print("CellChat 前5活跃通信方向：")
for st, cnt in cellchat_st_counts.most_common(5):
    print(f"{st}: {cnt} 次")

print("\nLIANA 前5活跃通信方向：")
for st, cnt in liana_st_counts.most_common(5):
    print(f"{st}: {cnt} 次")

# 步骤5：相关性分析（prob vs spatial_score）
print("\n=== 步骤5：相关性分析 ===")

cellchat_agg = cellchat_df.groupby('lr_pair')['prob'].mean().reset_index()
adaptive_agg = adaptive_df.groupby('lr_pair')['spatial_communication_score'].mean().reset_index()

common_df = pd.merge(cellchat_agg, adaptive_agg, on='lr_pair', how='inner')

if len(common_df) > 0:
    corr, pval = spearmanr(common_df['prob'], common_df['spatial_communication_score'])
    print(f"CellChat prob 与 Adaptive spatial_score 的 Spearman 相关系数: {corr:.4f}")
    print(f"（p-value: {pval:.4e}）")
else:
    print("没有共同的 LR 对，无法计算相关性")

# LIANA rank vs spatial_score
liana_agg = liana_df.groupby('lr_pair')[['magnitude_rank', 'specificity_rank']].mean().reset_index()
liana_agg['liana_score'] = (liana_agg['magnitude_rank'] + liana_agg['specificity_rank']) / 2

common_liana_adapt = pd.merge(liana_agg, adaptive_agg, on='lr_pair', how='inner')
if len(common_liana_adapt) > 0:
    corr_liana, pval_liana = spearmanr(common_liana_adapt['liana_score'], common_liana_adapt['spatial_communication_score'])
    print(f"LIANA 平均 rank 与 Adaptive spatial_score 的 Spearman 相关系数: {corr_liana:.4f}")
    print(f"（p-value: {pval_liana:.4e}）")

# 步骤6：CellChat 高置信对保留率
print("\n=== 步骤6：高置信对保留率 ===")

median_prob = cellchat_agg['prob'].median()
top_cellchat = set(cellchat_agg[cellchat_agg['prob'] > median_prob]['lr_pair'])

retained_liana    = top_cellchat.intersection(liana_lrs)
retained_adaptive = top_cellchat.intersection(adaptive_lrs)

print(f"CellChat 高 prob 对（>中位数）总数: {len(top_cellchat)}")
print(f"其中被 LIANA 保留的比例: {len(retained_liana)/len(top_cellchat):.2%}")
print(f"其中被 Adaptive 模块保留的比例: {len(retained_adaptive)/len(top_cellchat):.2%}")

# 步骤7：脑相关通路偏好对比
print("\n=== 步骤7：脑相关通路偏好 ===")

brain_keywords = ['Glu', 'GABR', 'Wnt', 'Fgf', 'Nrg', 'Cadm', 'Apoe', 'App', 'Ptn', 'Sst', 'Nlgn', 'Nrxn']
cellchat_brain = sum(any(kw in lr for kw in brain_keywords) for lr in cellchat_lrs)
liana_brain    = sum(any(kw in lr for kw in brain_keywords) for lr in liana_lrs)

print(f"脑相关信号通路占比：")
print(f"  CellChat: {cellchat_brain / len(cellchat_lrs):.2%} ({cellchat_brain}/{len(cellchat_lrs)})")
print(f"  LIANA:    {liana_brain / len(liana_lrs):.2%} ({liana_brain}/{len(liana_lrs)})")

print("\n=== 分析完成 ===")
print("建议：把 Jaccard、低重叠率、保留率、相关系数、脑通路占比这些数字直接填进论文表格或 PPT。")
print("下一步可做：画相关性散点图、热图，或继续任务二（空间敏感性验证）。")
