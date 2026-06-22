import scanpy as sc
import liana as li
from liana.method import rank_aggregate
import pandas as pd
import numpy as np

# =====================================================
# 步骤①：加载 Visium Adult Mouse Brain 空间转录组数据
# =====================================================

print("步骤①：加载 Visium Mouse Brain 数据 ...")

# ⚠ 注意：这里填的是“目录”，不是 h5 文件
visium_path = "data/visium_mouse_brain"

# 读取 Visium 数据
adata = sc.read_visium(
    path=visium_path,
    load_images=True
)

print(f"Spot 数量: {adata.n_obs}")
print(f"基因数量: {adata.n_vars}")
print("空间坐标是否存在：", "spatial" in adata.obsm)

# =====================================================
# 步骤②：基础质量控制（不破坏空间结构）
# =====================================================

print("\n步骤②：基础 QC + 归一化 ...")

# 去除极低表达基因
sc.pp.filter_genes(adata, min_cells=3)

# 仅做总量归一化（⚠ 不做 log1p）
sc.pp.normalize_total(adata, target_sum=1e4)

# =====================================================
# 步骤③：对 spot 进行聚类（用于 sender / receiver 定义）
# =====================================================

print("\n步骤③：spot 聚类 ...")

# PCA
sc.pp.pca(adata, n_comps=50)

# 基于表达的邻接图（此处不引入空间）
sc.pp.neighbors(adata)

# Leiden 聚类
sc.tl.leiden(
    adata,
    resolution=0.5,
    flavor="igraph",
    n_iterations=2
)

# 将聚类结果作为“通信群体”
adata.obs["cell_type"] = adata.obs["leiden"].astype(str)

print("spot 群体：", adata.obs["cell_type"].unique().tolist())

# =====================================================
# 步骤④：修复 NaN / Inf（数值安全）
# =====================================================

print("\n步骤④：数值检查 ...")

if hasattr(adata.X, "data"):  # 稀疏矩阵
    mask = np.isnan(adata.X.data) | np.isinf(adata.X.data)
    adata.X.data[mask] = 0
else:
    adata.X[np.isnan(adata.X) | np.isinf(adata.X)] = 0

print("数值检查完成。")

# =====================================================
# 步骤⑤：LIANA 共识通信推断
# =====================================================

print("\n步骤⑤：运行 LIANA 共识推断 ...")

rank_aggregate(
    adata=adata,
    groupby="cell_type",
    resource_name="mouseconsensus",  # ⭐ 关键：小鼠 LR 资源
    use_raw=False,
    expr_prop=0.1,
    n_perms=100,
    verbose=True,
    return_all_lrs=True
)

# =====================================================
# 步骤⑥：筛选高置信 LR 对
# =====================================================

print("\n步骤⑥：筛选高置信 LR 对 ...")

res = adata.uns["liana_res"]

high_conf = res[
    (res["magnitude_rank"] <= 0.05) &
    (res["specificity_rank"] <= 0.05)
].copy().reset_index(drop=True)

print(f"高置信 LR 对数量: {len(high_conf)}")

# =====================================================
# 步骤⑦：拆分复合物为单基因
# =====================================================

def explode_complex(df, col):
    return df[col].str.split("_").explode().reset_index(drop=True)

high_conf["ligand"] = explode_complex(high_conf, "ligand_complex")
high_conf["receptor"] = explode_complex(high_conf, "receptor_complex")

high_conf_single = high_conf[
    [
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "ligand",
        "receptor",
        "magnitude_rank",
        "specificity_rank"
    ]
]

print(high_conf_single.head())

# =====================================================
# 步骤⑧：导出 CSV（供 Adaptive Distance-aware 模块使用）
# =====================================================

output_file = "liana_consensus_lr_pairs.csv"
high_conf_single.to_csv(output_file, index=False)

print(f"\n✅ 已导出：{output_file}")
print("下一步可直接结合 adata.obsm['spatial'] 进行距离建模。")
