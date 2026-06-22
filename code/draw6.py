# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# print("正在读取文件...")
#
# # 读取 CellChat 和 LIANA 文件
# cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
# liana_df = pd.read_csv("liana_consensus_lr_pairs.csv")
#
# # 创建 lr_pair
# cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']
# liana_df['lr_pair'] = liana_df['ligand'] + '_' + liana_df['receptor']
#
# # 步骤1: 计算 CellChat prob 中位数
# cellchat_agg = cellchat_df.groupby('lr_pair')['prob'].mean().reset_index()
# median_prob = cellchat_agg['prob'].median()
# print(f"CellChat prob 中位数: {median_prob}")
#
# # 步骤2: CellChat 高 prob 对（> 中位数）
# high_prob_lrs = set(cellchat_agg[cellchat_agg['prob'] > median_prob]['lr_pair'])
# print(f"CellChat 高 prob 对总数: {len(high_prob_lrs)}")
#
# # 步骤3: SpaADA (LIANA) 保留的数量和比例
# retained_lrs = high_prob_lrs & set(liana_df['lr_pair'].unique())
# retention_rate = len(retained_lrs) / len(high_prob_lrs) * 100 if len(high_prob_lrs) > 0 else 0
# print(f"SpaADA 保留数量: {len(retained_lrs)}")
# print(f"保留比例: {retention_rate:.2f}%")
#
# # 步骤4: 准备柱状图数据
# data = pd.DataFrame({
#     'Category': ['CellChat High-Prob Pairs', 'Retained by SpaADA'],
#     'Count': [len(high_prob_lrs), len(retained_lrs)]
# })
#
# # 绘制柱状图
# plt.figure(figsize=(6, 4))
# sns.barplot(x='Category', y='Count', data=data, palette=['#1f77b4', '#ff7f0e'])
# plt.title('Retention of CellChat High-Confidence LR Pairs by SpaADA', fontsize=10)
# plt.ylabel('Number of LR Pairs')
# plt.xlabel('')
# plt.ylim(0, len(high_prob_lrs) * 1.2)  # 留出空间显示百分比
#
# # 在柱子上标注数量和百分比
# for i, v in enumerate(data['Count']):
#     plt.text(i, v + 5, f'{v}\n({v/len(high_prob_lrs)*100:.2f}%)' if i == 1 else f'{v}',
#              ha='center', fontsize=8, fontweight='bold')
#
# plt.tight_layout()
# plt.savefig('Fig6_retention.png', dpi=300, bbox_inches='tight')
# plt.show()
#
# print("柱状图已保存为 Fig6_retention.png")
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ===============================
# ⭐ 设置中文字体（宋体）
# ===============================
plt.rcParams['font.sans-serif'] = ['SimSun']
plt.rcParams['axes.unicode_minus'] = False

print("正在读取文件...")

# 读取数据
cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
liana_df = pd.read_csv("liana_consensus_lr_pairs.csv")

# 构建 lr_pair
cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']
liana_df['lr_pair'] = liana_df['ligand'] + '_' + liana_df['receptor']

# ===============================
# 步骤1：计算 CellChat 中位数
# ===============================
cellchat_agg = cellchat_df.groupby('lr_pair')['prob'].mean().reset_index()
median_prob = cellchat_agg['prob'].median()
print(f"CellChat prob 中位数: {median_prob}")

# ===============================
# 步骤2：筛选高概率 LR 对
# ===============================
high_prob_lrs = set(cellchat_agg[cellchat_agg['prob'] > median_prob]['lr_pair'])
print(f"CellChat 高概率 LR 对数量: {len(high_prob_lrs)}")

# ===============================
# 步骤3：SpaADA 保留情况
# ===============================
retained_lrs = high_prob_lrs & set(liana_df['lr_pair'].unique())
retention_rate = len(retained_lrs) / len(high_prob_lrs) * 100 if len(high_prob_lrs) > 0 else 0

print(f"SpaADA 保留数量: {len(retained_lrs)}")
print(f"保留比例: {retention_rate:.2f}%")

# ===============================
# 步骤4：柱状图数据（中文）
# ===============================
data = pd.DataFrame({
    '类别': ['CellChat高概率LR对', 'SpaADA保留的LR对'],
    '数量': [len(high_prob_lrs), len(retained_lrs)]
})

# ===============================
# 绘图
# ===============================
plt.figure(figsize=(6, 4))

sns.barplot(
    x='类别',
    y='数量',
    data=data,
    palette=['#1f77b4', '#ff7f0e']
)

plt.title('SpaADA对CellChat高置信度LR对的保留情况', fontsize=10)
plt.ylabel('LR对数量')
plt.xlabel('')

plt.ylim(0, len(high_prob_lrs) * 1.2)

# ===============================
# ⭐ 标注数量 & 百分比
# ===============================
for i, v in enumerate(data['数量']):
    if i == 1:
        plt.text(
            i, v + len(high_prob_lrs) * 0.02,
            f'{v}\n({retention_rate:.2f}%)',
            ha='center',
            fontsize=8,
            fontweight='bold'
        )
    else:
        plt.text(
            i, v + len(high_prob_lrs) * 0.02,
            f'{v}',
            ha='center',
            fontsize=8,
            fontweight='bold'
        )

plt.tight_layout()

# ===============================
# 保存
# ===============================
plt.savefig(
    'Fig6_保留分析_中文.png',
    dpi=300,
    bbox_inches='tight'
)

plt.close()

print("✅ 柱状图已保存为 Fig6_保留分析_中文.png")
