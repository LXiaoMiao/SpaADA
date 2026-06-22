# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# from scipy.stats import spearmanr
# import matplotlib.ticker as ticker
#
# # 支持中文和负号（虽然现在全英文，但保留以防万一）
# # plt.rcParams['font.sans-serif'] = ['SimHei']
# # plt.rcParams['axes.unicode_minus'] = False
#
# print("Loading data...")
# cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
# adaptive_df = pd.read_csv("adaptive_distance_spatial_scores.csv")
#
# # Create lr_pair column
# cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']
# adaptive_df['lr_pair'] = adaptive_df['ligand'] + '_' + adaptive_df['receptor']
#
# # Aggregate mean values per LR pair
# cellchat_agg = cellchat_df.groupby('lr_pair')['prob'].mean().reset_index()
# adaptive_agg = adaptive_df.groupby('lr_pair')['spatial_communication_score'].mean().reset_index()
#
# # Merge common LR pairs
# common_df = pd.merge(cellchat_agg, adaptive_agg, on='lr_pair', how='inner')
# print(f"Number of common LR pairs: {len(common_df)}")
#
# corr, pval = spearmanr(common_df['prob'], common_df['spatial_communication_score'])
# print(f"Spearman ρ = {corr:.4f} (p = {pval:.2e})")
#
# fig, ax = plt.subplots(figsize=(10, 8))
#
# sns.scatterplot(
#     ax=ax,
#     data=common_df,
#     x='prob',
#     y='spatial_communication_score',
#     alpha=0.7,
#     s=80,
#     color='royalblue'
# )
#
# ax.set_xscale('log')
# ax.set_yscale('log')
#
# # Force compact scientific notation (1e-7 style)
# def sci_notation(x, pos):
#     if x == 0:
#         return '0'
#     return f'{x:.0e}'
#
# ax.xaxis.set_major_formatter(ticker.FuncFormatter(sci_notation))
# ax.yaxis.set_major_formatter(ticker.FuncFormatter(sci_notation))
#
# # Annotation with correlation
# ax.text(
#     0.05, 0.95,
#     f'Spearman ρ = {corr:.4f}\np = {pval:.2e}',
#     transform=ax.transAxes,
#     fontsize=14,
#     verticalalignment='top',
#     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
# )
#
# # English title and labels (integrated SpaADA method name)
# ax.set_title('Correlation between CellChat Probability and\nSpaADA Spatial Communication Score', fontsize=16, pad=15)
# ax.set_xlabel('CellChat Communication Probability', fontsize=14)
# ax.set_ylabel('SpaADA Spatial Communication Score', fontsize=14)
#
# ax.grid(True, linestyle='--', alpha=0.5)
#
# plt.tight_layout()
# plt.savefig('correlation_prob_vs_score_SpaADA.png', dpi=300, bbox_inches='tight')
# print("Figure X saved as correlation_prob_vs_score_SpaADA.png")
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import matplotlib.ticker as ticker

# ===============================
# ⭐ 设置中文字体（宋体）
# ===============================
plt.rcParams['font.sans-serif'] = ['SimSun']
plt.rcParams['axes.unicode_minus'] = False

print("正在加载数据...")
cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
adaptive_df = pd.read_csv("adaptive_distance_spatial_scores.csv")

# ===============================
# 构建 lr_pair
# ===============================
cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']
adaptive_df['lr_pair'] = adaptive_df['ligand'] + '_' + adaptive_df['receptor']

# ===============================
# 按 LR 对取均值
# ===============================
cellchat_agg = cellchat_df.groupby('lr_pair')['prob'].mean().reset_index()
adaptive_agg = adaptive_df.groupby('lr_pair')['spatial_communication_score'].mean().reset_index()

# 合并共同 LR 对
common_df = pd.merge(cellchat_agg, adaptive_agg, on='lr_pair', how='inner')
print(f"共同 LR 对数量: {len(common_df)}")

# ===============================
# 计算 Spearman 相关性
# ===============================
corr, pval = spearmanr(common_df['prob'], common_df['spatial_communication_score'])
print(f"Spearman ρ = {corr:.4f} (p = {pval:.2e})")

# ===============================
# 绘图
# ===============================
fig, ax = plt.subplots(figsize=(10, 8))

sns.scatterplot(
    ax=ax,
    data=common_df,
    x='prob',
    y='spatial_communication_score',
    alpha=0.7,
    s=80,
    color='royalblue'
)

# 对数坐标
ax.set_xscale('log')
ax.set_yscale('log')

# 科学计数法格式
def sci_notation(x, pos):
    if x == 0:
        return '0'
    return f'{x:.0e}'

ax.xaxis.set_major_formatter(ticker.FuncFormatter(sci_notation))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(sci_notation))

# ===============================
# ⭐ 中文注释（相关性）
# ===============================
ax.text(
    0.05, 0.95,
    f'Spearman ρ = {corr:.4f}\np = {pval:.2e}',
    transform=ax.transAxes,
    fontsize=14,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)

# ===============================
# ⭐ 中文标题与坐标轴
# ===============================
ax.set_title(
    'CellChat 通信概率与\nSpaADA 空间通信评分的相关性',
    fontsize=16,
    pad=15
)

ax.set_xlabel('CellChat 通信概率', fontsize=14)
ax.set_ylabel('SpaADA 空间通信评分', fontsize=14)

# 网格
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()

# ===============================
# 保存
# ===============================
plt.savefig(
    '相关性_通信概率_vs_空间评分_SpaADA.png',
    dpi=300,
    bbox_inches='tight'
)

plt.close()

print("✅ 图已生成：相关性_通信概率_vs_空间评分_SpaADA.png")
