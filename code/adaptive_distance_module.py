import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
from scipy.spatial.distance import cdist
from scipy.sparse import issparse

# ===============================
# Step 1. 加载 Visium 数据
# ===============================
print("Step 1: 加载 Visium 数据")

adata = sc.read_visium("data/visium_mouse_brain")
adata.var_names_make_unique()

X = adata.X.toarray() if issparse(adata.X) else adata.X
coords = adata.obsm["spatial"]
N, G = X.shape

print(f"Spot 数: {N}")
print(f"基因数: {G}")

# ===============================
# Step 2. 读取 LR 对
# ===============================
lr_df = pd.read_csv("liana_consensus_lr_pairs.csv")
gene_names = adata.var_names.tolist()
gene_to_idx = {g: i for i, g in enumerate(gene_names)}

valid_lrs = [
    (r["ligand_complex"], r["receptor"])
    for _, r in lr_df.iterrows()
    if r["ligand_complex"] in gene_to_idx and r["receptor"] in gene_to_idx
]

print(f"有效 LR 对数量: {len(valid_lrs)}")

# ===============================
# Step 3. 构建 KNN 邻域
# ===============================
print("Step 3: 构建空间 KNN 邻域")

dist = cdist(coords, coords)
K = 20
knn_idx = np.argsort(dist, axis=1)[:, 1:K+1]
knn_dist = np.take_along_axis(dist, knn_idx, axis=1)
knn_dist = knn_dist / knn_dist.max()  # 归一化到 0–1

# ===============================
# Step 4. 构建 LR 特征
# ===============================
print("Step 4: 计算 μ_l, μ_r, d_lr")

features = []
C_cache = []

for l, r in valid_lrs:
    li, ri = gene_to_idx[l], gene_to_idx[r]

    mu_l = np.log(X[:, li].mean() + 1e-6)
    mu_r = np.log(X[:, ri].mean() + 1e-6)
    d_lr = knn_dist.mean()

    features.append([mu_l, mu_r, d_lr])

    # 表达驱动 C_ij（只取 KNN）
    C = X[:, li][:, None] * X[knn_idx[:, :], ri]
    C_cache.append(torch.tensor(C, dtype=torch.float32))

features = torch.tensor(features, dtype=torch.float32)

# ===============================
# Step 5. σ 模型（0–1）
# ===============================
class SigmaModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(3, 16),
            nn.ReLU(),
            nn.Linear(16, 1)
        )

    def forward(self, x):
        return torch.sigmoid(self.net(x)).squeeze()

model = SigmaModel()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)

knn_dist_t = torch.tensor(knn_dist, dtype=torch.float32)

# ===============================
# Step 6. 训练
# ===============================
print("Step 6: 训练 Adaptive Distance 模型")

for epoch in range(100):
    optimizer.zero_grad()

    sigma = model(features) + 1e-6  # (L,)
    loss = 0.0

    for i, C in enumerate(C_cache):
        s = sigma[i]
        Kij = torch.exp(- (knn_dist_t ** 2) / (2 * s ** 2))
        alpha = Kij / (Kij.sum(dim=1, keepdim=True) + 1e-8)
        loss -= (C * alpha).mean()

    loss = loss / len(C_cache)
    loss.backward()
    optimizer.step()

    if epoch % 10 == 0:
        print(f"Epoch {epoch:03d} | loss = {loss.item():.4f}")

# ===============================
# Step 7. 推断 & 归一化 score
# ===============================
print("Step 7: 推断 spatial_communication_score")

with torch.no_grad():
    sigma = model(features).cpu().numpy()

scores = []

for (l, r), s, C in zip(valid_lrs, sigma, C_cache):
    Kij = np.exp(- (knn_dist ** 2) / (2 * s ** 2 + 1e-8))
    alpha = Kij / (Kij.sum(axis=1, keepdims=True) + 1e-8)
    scores.append((C.numpy() * alpha).mean())

scores = np.array(scores)
scores = (scores - scores.min()) / (scores.max() - scores.min() + 1e-8)

# ===============================
# 保存
# ===============================
pd.DataFrame({
    "ligand": [x[0] for x in valid_lrs],
    "receptor": [x[1] for x in valid_lrs],
    "sigma_lr": sigma,
    "spatial_communication_score": scores
}).to_csv("adaptive_distance_spatial_scores.csv", index=False)

print("✅ 完成：adaptive_distance_spatial_scores.csv 已保存")
