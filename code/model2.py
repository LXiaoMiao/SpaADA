import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.neighbors import NearestNeighbors

# =========================
# Step 1. 读取数据
# =========================

print("Step 1: 加载 Visium 数据")

adata = sc.read_visium("data/visium_mouse_brain")
adata.var_names_make_unique()

coords = adata.obsm["spatial"].astype(np.float32)
X = adata.X.toarray().astype(np.float32)
gene_names = adata.var_names.tolist()

N, G = X.shape
print("Spot 数:", N)
print("基因数:", G)

# =========================
# Step 2. 读取 LIANA LR 对
# =========================

lr_df = pd.read_csv("liana_consensus_lr_pairs.csv")

gene_to_idx = {g: i for i, g in enumerate(gene_names)}

valid_rows = []
for _, r in lr_df.iterrows():
    if r["ligand_complex"] in gene_to_idx and r["receptor"] in gene_to_idx:
        valid_rows.append(r)

lr_df = pd.DataFrame(valid_rows).reset_index(drop=True)
print("有效 LR 对数量:", len(lr_df))

# =========================
# Step 3. 构建空间 KNN 邻域
# =========================

print("Step 3: 构建空间 KNN 邻域")

K = 30  # 局部邻域，避免 N^2 爆炸
nbrs = NearestNeighbors(n_neighbors=K, metric="euclidean").fit(coords)
distances, neighbors = nbrs.kneighbors(coords)

distances = torch.tensor(distances, dtype=torch.float32)
neighbors = torch.tensor(neighbors, dtype=torch.long)

# =========================
# Step 4. 预计算 μ_l, μ_r, d_lr
# =========================

print("Step 4: 计算 μ_l, μ_r, d_lr")

mu_l_list = []
mu_r_list = []
d_lr_list = []

for _, row in lr_df.iterrows():
    l = gene_to_idx[row["ligand_complex"]]
    r = gene_to_idx[row["receptor"]]

    xl = X[:, l]
    xr = X[:, r]

    mu_l = xl.mean()
    mu_r = xr.mean()

    mask_i = xl > 0
    mask_j = xr > 0

    if mask_i.sum() == 0 or mask_j.sum() == 0:
        d_lr = distances.mean().item()
    else:
        ii = np.where(mask_i)[0]
        jj = np.where(mask_j)[0]
        d_lr = np.mean(
            np.linalg.norm(
                coords[ii][:, None, :] - coords[jj][None, :, :],
                axis=-1
            )
        )

    mu_l_list.append(mu_l)
    mu_r_list.append(mu_r)
    d_lr_list.append(d_lr)

mu_l = torch.tensor(mu_l_list, dtype=torch.float32)
mu_r = torch.tensor(mu_r_list, dtype=torch.float32)
d_lr = torch.tensor(d_lr_list, dtype=torch.float32)

# =========================
# Step 5. 可训练模型定义
# =========================

class AdaptiveDistanceModel(nn.Module):
    def __init__(self, n_lr, sigma_min=20.0, sigma_max=500.0):
        super().__init__()

        self.sigma_min = sigma_min
        self.sigma_max = sigma_max

        # 可学习参数
        self.w1 = nn.Parameter(torch.tensor(0.1))
        self.w2 = nn.Parameter(torch.tensor(0.1))
        self.w3 = nn.Parameter(torch.tensor(0.1))
        self.b = nn.Parameter(torch.tensor(0.0))

    def forward(self, mu_l, mu_r, d_lr, dist_knn):
        """
        dist_knn: (N, K)
        """
        eps = 1e-6

        z = (
            self.w1 * torch.log(mu_l + eps)
            + self.w2 * torch.log(mu_r + eps)
            + self.w3 * d_lr
            + self.b
        )

        # 有界 σ（关键）
        sigma = self.sigma_min + (self.sigma_max - self.sigma_min) * torch.sigmoid(z)

        # 计算空间注意力（log-softmax 稳定）
        # dist_knn: [N, K]
        dist2 = dist_knn ** 2

        sigma_expand = sigma.view(-1, 1, 1)  # [LR, 1, 1]
        dist2_expand = dist2.unsqueeze(0)    # [1, N, K]

        logits = -dist2_expand / (2 * sigma_expand ** 2)

        logits = logits - logits.max(dim=-1, keepdim=True).values
        alpha = torch.softmax(logits, dim=-1)

        return sigma, alpha


# =========================
# Step 6. 构建训练目标
# =========================

model = AdaptiveDistanceModel(n_lr=len(lr_df))
optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)

X_tensor = torch.tensor(X, dtype=torch.float32)

lig_idx = torch.tensor([gene_to_idx[g] for g in lr_df["ligand_complex"]])
rec_idx = torch.tensor([gene_to_idx[g] for g in lr_df["receptor"]])

# =========================
# Step 7. 训练
# =========================

print("Step 6: 训练 Adaptive Distance 模型")

EPOCHS = 50
lambda_reg = 1e-4

for epoch in range(EPOCHS):
    optimizer.zero_grad()

    sigma, alpha = model(mu_l, mu_r, d_lr, distances)

    # 表达强度 C_ij
    X_l = X_tensor[:, lig_idx]       # [N, LR]
    X_r = X_tensor[:, rec_idx]       # [N, LR]

    # 取邻域
    X_r_knn = X_r[neighbors]         # [N, K, LR]
    X_l_exp = X_l.unsqueeze(1)       # [N, 1, LR]

    C = X_l_exp * X_r_knn            # [N, K, LR]

    # 空间通信强度
    S = C * alpha.permute(1, 2, 0)

    score = S.mean()

    # σ 正则（防止坍缩）
    reg = ((sigma - d_lr.mean()) ** 2).mean()

    loss = -score + lambda_reg * reg

    loss.backward()
    optimizer.step()

    if epoch % 10 == 0:
        print(f"Epoch {epoch:03d} | loss = {loss.item():.4f}")

# =========================
# Step 8. 保存结果
# =========================

with torch.no_grad():
    sigma_final, _ = model(mu_l, mu_r, d_lr, distances)

result_df = pd.DataFrame({
    "ligand": lr_df["ligand_complex"],
    "receptor": lr_df["receptor"],
    "sigma_lr": sigma_final.cpu().numpy()
})

result_df.to_csv("adaptive_distance_spatial_scores.csv", index=False)

print("✅ 训练完成，已保存 adaptive_distance_spatial_scores.csv")
