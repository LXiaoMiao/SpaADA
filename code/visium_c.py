# ====================== COMMOT Visium Mouse Brain - PyCharm 兼容最终版 ======================
import scanpy as sc
import matplotlib.pyplot as plt
import commot as ct
import os

# ==================== 1. 设置数据路径 ====================
data_path = r"F:\第二个点\code\data\visium_mouse_brain"

if not os.path.exists(data_path):
    raise FileNotFoundError(f"路径不存在，请检查：{data_path}")

print("数据路径：", data_path)

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

print(f"数据加载完成！形状: {adata.shape}")

# ==================== 3. 加载 CellChat 数据库 ====================
print("正在加载 CellChat 配体-受体数据库...")
df_ligrec = ct.pp.ligand_receptor_database(database='CellChat', species='mouse')
print(f"数据库加载完成，共有 {len(df_ligrec)} 个 LR 对")

# ==================== 4. 计算空间通讯 ====================
print("开始计算空间通讯（这可能需要 2-5 分钟，请耐心等待）...")

ct.tl.spatial_communication(
    adata,
    database_name='cellchat',
    df_ligrec=df_ligrec,
    dis_thr=500,
    heteromeric=True,
    pathway_sum=True
)

print("通讯计算完成！")

# ==================== 5. 把结果转到 .obs ====================
receiver_key = 'commot-cellchat-sum-receiver'
if receiver_key in adata.obsm:
    adata.obs['receiver_strength'] = adata.obsm[receiver_key].sum(axis=1)
    print(f"已创建 obs['receiver_strength']，用于绘图")
else:
    print(f"未找到 {receiver_key}")

# ==================== 6. 可视化并保存图片（避开 PyCharm bug） ====================
print("正在生成并保存接收信号强度热图...")

fig, ax = plt.subplots(figsize=(10, 8))

sc.pl.spatial(
    adata,
    color='receiver_strength',
    cmap='viridis',
    alpha=0.85,
    size=1.8,
    ax=ax,
    title='Overall Received Signal Strength (CellChat sum)',
    show=False
)

plt.tight_layout()

# 保存图片（不使用 plt.show()，避免 PyCharm 后端错误）
output_dir = os.path.join(data_path, "commot_output")
os.makedirs(output_dir, exist_ok=True)
save_path = os.path.join(output_dir, "overall_received_signal.png")
fig.savefig(save_path, dpi=300, bbox_inches='tight')

print(f"热图已成功保存到：{save_path}")
print("你可以直接打开这个 png 文件查看结果。")

# 可选：也保存处理后的 adata
adata.write_h5ad(os.path.join(output_dir, "visium_mouse_brain_commot_processed.h5ad"))
print("处理后的 adata 已保存！")
