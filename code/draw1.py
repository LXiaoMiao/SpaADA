# import matplotlib
# matplotlib.use("Agg")
#
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
#
#
# top = 50
# # ===============================
# # 步骤①：读取 LIANA 输出
# # ===============================
# df = pd.read_csv("liana_consensus_lr_pairs.csv")
#
# # ===============================
# # 步骤②：构建绘图变量
# # ===============================
# # ligand → receptor
# df["lr_pair"] = df["ligand"] + " → " + df["receptor"]
#
# # source → target
# df["source_target"] = df["source"].astype(str) + " → " + df["target"].astype(str)
#
# # 将 rank 转为可视化得分
# df["magnitude_plot"] = -np.log10(df["magnitude_rank"] + 1e-20)
# df["specificity_plot"] = -np.log10(df["specificity_rank"] + 1e-20)
#
# # ===============================
# # 步骤③：筛选 Top 50 LR 对（按最大 magnitude）
# # ===============================
# top_lr = (
#     df.groupby("lr_pair")["magnitude_plot"]
#     .max()
#     .sort_values(ascending=False)
#     .head(top)
#     .index
# )
#
# df_plot = df[df["lr_pair"].isin(top_lr)].copy()
#
# # Y 轴顺序（强度高的在上）
# lr_order = (
#     df_plot.groupby("lr_pair")["magnitude_plot"]
#     .max()
#     .sort_values(ascending=False)
#     .index
# )
#
# # ===============================
# # 步骤④：绘制 Dotplot
# # ===============================
# plt.figure(figsize=(28, 16))  # 稍微加大画布，避免 50 行拥挤
#
# sns.scatterplot(
#     data=df_plot,
#     x="source_target",
#     y="lr_pair",
#     size="specificity_plot",     # 点面积 ← specificity
#     hue="magnitude_plot",        # 颜色 ← magnitude
#     sizes=(30, 300),
#     palette="Reds",
#     alpha=0.85,
#     edgecolor=None
# )
#
# plt.yticks(range(len(lr_order)), lr_order)
# plt.xlabel("Source → Target", fontsize=16)
# plt.ylabel("Ligand → Receptor", fontsize=16)
# plt.title(
#     f"Top {top} High-Confidence Ligand–Receptor Communications",
#     fontsize=18,
#     pad=10
# )
#
# plt.xticks(rotation=45, ha="right")
#
# plt.legend(
#     title="Magnitude / Specificity",
#     bbox_to_anchor=(1, 1),
#     loc="upper left"
# )
#
# plt.tight_layout()
# plt.savefig(
#     f"lr_communication_dotplot_top{top}.png",
#     dpi=300,
#     bbox_inches="tight"
# )
# plt.close()

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ===============================
# ⭐ 设置中文字体（宋体）
# ===============================
plt.rcParams['font.sans-serif'] = ['SimSun']
plt.rcParams['axes.unicode_minus'] = False

top = 50

# ===============================
# 步骤①：读取 LIANA 输出
# ===============================
df = pd.read_csv("liana_consensus_lr_pairs.csv")

# ===============================
# 步骤②：构建绘图变量
# ===============================
# 配体 → 受体
df["lr_pair"] = df["ligand"] + " → " + df["receptor"]

# 发送细胞 → 接收细胞
df["source_target"] = df["source"].astype(str) + " → " + df["target"].astype(str)

# 将 rank 转为可视化得分
df["magnitude_plot"] = -np.log10(df["magnitude_rank"] + 1e-20)
df["specificity_plot"] = -np.log10(df["specificity_rank"] + 1e-20)

# ===============================
# 步骤③：筛选 Top 50 LR 对
# ===============================
top_lr = (
    df.groupby("lr_pair")["magnitude_plot"]
    .max()
    .sort_values(ascending=False)
    .head(top)
    .index
)

df_plot = df[df["lr_pair"].isin(top_lr)].copy()

# Y 轴顺序（强度高的在上）
lr_order = (
    df_plot.groupby("lr_pair")["magnitude_plot"]
    .max()
    .sort_values(ascending=False)
    .index
)

# ===============================
# ⭐ 步骤④：重命名列（用于图例中文）
# ===============================
df_plot = df_plot.rename(columns={
    "magnitude_plot": "强度",
    "specificity_plot": "特异性"
})

# ===============================
# 步骤⑤：绘制 Dotplot
# ===============================
plt.figure(figsize=(28, 16))

sns.scatterplot(
    data=df_plot,
    x="source_target",
    y="lr_pair",
    size="特异性",     # 点大小 ← 特异性
    hue="强度",        # 颜色 ← 强度
    sizes=(30, 300),
    palette="Reds",
    alpha=0.85,
    edgecolor=None
)

plt.yticks(range(len(lr_order)), lr_order)

# 中文坐标轴
plt.xlabel("发送细胞 → 接收细胞", fontsize=16)
plt.ylabel("配体 → 受体", fontsize=16)

# 中文标题
plt.title(
    f"前{top}个高置信度LR通讯点图",
    fontsize=18,
    pad=10
)

plt.xticks(rotation=45, ha="right")

# ⭐ 图例（自动使用中文变量名）
plt.legend(
    bbox_to_anchor=(1, 1),
    loc="upper left"
)

plt.tight_layout()

plt.savefig(
    f"配体受体通讯_点图_Top{top}.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()
