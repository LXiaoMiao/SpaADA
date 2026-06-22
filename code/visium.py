# ====================== SpaADA + COMMOT 正式实验版 ======================
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import commot as ct
import os

# ==================== 1. 设置路径 ====================
data_path = r"F:\第二个点\code\data\visium_mouse_brain"
FOUND_ERROR = FileNotFoundError(f"Visium 数据路径不存在：{data_path}")
lr_file = r"F:\第二个点\code\liana_consensus_lr_pairs.csv"   # 请确认这个路径正确

if not os.path.exists(data_path):
    raise FOUND_ERROR
if not os.path.exists(lr_file):
    raise FileNotFoundError(f"SpaADA LR 文件不存在：{lr_file}")

print("数据路径和 LR 文件确认完成")

# ==================== 2. 加载 Visium 数据 ====================
print("正在加载 Visium 数据...")
adata = sc.read_visium(
    path=data_path,
    count_file='filtered_feature_bc_matrix.h5',
    load_images=True
)

adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['counts'] = adata.X.copy()

print(f"Visium 数据加载完成！形状: {adata.shape}")

# ==================== 3. 读取你的 SpaADA LR 表格 ====================
print("正在读取 SpaADA LR 表格...")
lr_df = pd.read_csv(lr_file)

print(f"SpaADA LR 表格读取完成，共 {len(lr_df)} 行")
print("列名：", lr_df.columns.tolist())
print(lr_df.head())

# ==================== 4. 构建 COMMOT 需要的 df_ligrec ====================
# 使用 ligand 和 receptor 单基因列
df_ligrec = lr_df[['ligand', 'receptor']].copy()
df_ligrec = df_ligrec.rename(columns={'ligand': 'ligand', 'receptor': 'receptor'})

# 为重要通路添加 pathway 标签（你可以根据需要调整）
def assign_pathway(row):
    ligand = str(row['ligand']).lower()
    receptor = str(row['receptor']).lower()
    if 'psap' in ligand or 'gpr37' in receptor:
        return 'PSAP_signaling'
    elif 'app' in ligand or 'hspa8' in ligand or 'aplp' in receptor:
        return 'APP_signaling'
    elif 'apoe' in ligand and 'lrp' in receptor:
        return 'APOE_LRP_signaling'
    else:
        return 'Other_signaling'

df_ligrec['pathway'] = df_ligrec.apply(assign_pathway, axis=1)

print(f"自定义 df_ligrec 构建完成，共 {len(df_ligrec)} 个 LR 对")
print("Pathway 分布：\n", df_ligrec['pathway'].value_counts())

# ==================== 5. 计算空间通讯（使用你的 SpaADA LR） ====================
print("开始使用 SpaADA LR 计算空间通讯...")

ct.tl.spatial_communication(
    adata,
    database_name='user_database',
    df_ligrec=df_ligrec,
    dis_thr=500,
    heteromeric=True,
    pathway_sum=True
)

print("基于 SpaADA LR 的通讯计算完成！")

# ==================== 6. 把结果转到 .obs 并可视化 ====================
print("正在生成接收信号强度热图...")

# 使用整体 sum-receiver（因为单个 pathway 键可能不稳定）
receiver_key = 'commot-user_database-sum-receiver'

if receiver_key in adata.obsm:
    adata.obs['spaada_receiver_strength'] = adata.obsm[receiver_key].sum(axis=1)
    color_key = 'spaada_receiver_strength'
else:
    # fallback
    color_key = 'commot-cellchat-sum-receiver' if 'commot-cellchat-sum-receiver' in adata.obsm else None

fig, ax = plt.subplots(figsize=(10, 8))

sc.pl.spatial(
    adata,
    color=color_key,
    cmap='viridis',
    alpha=0.85,
    size=1.8,
    ax=ax,
    title='SpaADA-based Overall Received Signal Strength',
    show=False
)

plt.tight_layout()

# 保存图片
output_dir = os.path.join(data_path, "commot_output")
os.makedirs(output_dir, exist_ok=True)
save_path = os.path.join(output_dir, "SpaADA_received_signal.png")
fig.savefig(save_path, dpi=300, bbox_inches='tight')

print(f"热图已成功保存到：{save_path}")
print("你可以直接打开这个 png 文件查看基于 SpaADA LR 的结果。")

# 保存处理后的 adata
adata.write_h5ad(os.path.join(output_dir, "visium_spaada_commot_processed.h5ad"))
print("处理后的 adata 已保存！")
