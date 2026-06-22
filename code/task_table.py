import pandas as pd

print("正在读取文件...")

# 读取 CellChat 文件
cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']

# 读取 LIANA / SpaADA 文件
liana_df = pd.read_csv("liana_consensus_lr_pairs.csv")
liana_df['lr_pair'] = liana_df['ligand'] + '_' + liana_df['receptor']

# 读取 Adaptive 文件
adaptive_df = pd.read_csv("adaptive_distance_spatial_scores.csv")
adaptive_df['lr_pair'] = adaptive_df['ligand'] + '_' + adaptive_df['receptor']

# 步骤1: 提取唯一 LR 对集合
cellchat_lrs = set(cellchat_df['lr_pair'].unique())
liana_lrs = set(liana_df['lr_pair'].unique())

# 步骤2: 计算三种 LR 对（无重复）
common_lrs = cellchat_lrs & liana_lrs
cellchat_only = cellchat_lrs - liana_lrs
spaeda_only = liana_lrs - cellchat_lrs

print(f"CellChat 唯一 LR 对数: {len(cellchat_lrs)}")
print(f"SpaADA (LIANA) 唯一 LR 对数: {len(liana_lrs)}")
print(f"共同 LR 对数: {len(common_lrs)}")
print(f"CellChat 独有 LR 对数: {len(cellchat_only)}")
print(f"SpaADA 独有 LR 对数: {len(spaeda_only)}")

# 步骤3: 共同 LR 对（全部 24 个，按 CellChat prob 降序）
common_with_prob = cellchat_df[cellchat_df['lr_pair'].isin(common_lrs)].groupby('lr_pair')['prob'].mean().reset_index()
common_with_score = adaptive_df[adaptive_df['lr_pair'].isin(common_lrs)].groupby('lr_pair')['spatial_communication_score'].mean().reset_index()
common_df = pd.merge(common_with_prob, common_with_score, on='lr_pair')
common_df = common_df.sort_values(by='prob', ascending=False)

print("\n=== 共同 LR 对 (总共 {} 个，按 CellChat prob 降序) ===".format(len(common_df)))
print(common_df.to_string(index=False))

# 步骤4: CellChat 独有 LR 对中的高 prob 前 30
cellchat_only_df = cellchat_df[cellchat_df['lr_pair'].isin(cellchat_only)]
cellchat_only_top = cellchat_only_df.groupby('lr_pair')['prob'].mean().reset_index()
cellchat_only_top = cellchat_only_top.sort_values(by='prob', ascending=False).head(30)

print("\n=== CellChat 独有 LR 对 高 prob 前 30（按 prob 降序）===")
print(cellchat_only_top.to_string(index=False))

# 步骤5: SpaADA 独有 LR 对中的高 score 前 30
spaeda_only_df = adaptive_df[adaptive_df['lr_pair'].isin(spaeda_only)]
spaeda_only_top = spaeda_only_df.groupby('lr_pair')['spatial_communication_score'].mean().reset_index()
spaeda_only_top = spaeda_only_top.sort_values(by='spatial_communication_score', ascending=False).head(30)

print("\n=== SpaADA 独有 LR 对 高 score 前 30（按 score 降序）===")
print(spaeda_only_top.to_string(index=False))
