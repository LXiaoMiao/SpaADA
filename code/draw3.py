# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# # ===============================
# # Step 1: 读取结果数据
# # ===============================
# df = pd.read_csv("adaptive_distance_spatial_scores.csv")
#
# # 如果你的 sigma 或 score 已经是 0–1，这一步可省略
# # df["sigma_lr"] = df["sigma_lr"].clip(0, 1)
# # df["spatial_communication_score"] = df["spatial_communication_score"].clip(0, 1)
#
# # ===============================
# # Step 2: 画图
# # ===============================
# plt.figure(figsize=(6, 5))
#
# sns.scatterplot(
#     data=df,
#     x="sigma_lr",
#     y="spatial_communication_score",
#     marker="o",          # 圆形点
#     s=8,                # 点大小（推荐 8–15）
#     alpha=0.6,           # 半透明，防止重叠
#     edgecolor=None
# )
#
# # ===============================
# # Step 3: 美化（论文级）
# # ===============================
# plt.xlabel(r"$\sigma_{l,r}$", fontsize=9)
# plt.ylabel("Spatial communication score", fontsize=9)
#
# plt.xticks(fontsize=7)
# plt.yticks(fontsize=7)
#
# plt.title(
#     r"Relationship between $\sigma_{l,r}$ and spatial communication score",
#     fontsize=10
# )
#
# plt.tight_layout()
#
# # ===============================
# # Step 4: 保存
# # ===============================
# plt.savefig(
#     "Fig3_sigma_vs_spatial_score.png",
#     dpi=300,
#     bbox_inches="tight"
# )
# plt.close()
#
# print("✅ 图 2 已生成：Fig2_sigma_vs_spatial_score.png")
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

# ===============================
# 步骤①：读取结果数据
# ===============================
df = pd.read_csv("adaptive_distance_spatial_scores.csv")

# ===============================
# 步骤②：绘制散点图
# ===============================
plt.figure(figsize=(6, 5))

sns.scatterplot(
    data=df,
    x="sigma_lr",
    y="spatial_communication_score",
    marker="o",
    s=8,
    alpha=0.6,
    edgecolor=None
)

# ===============================
# 步骤③：中文标注
# ===============================
plt.xlabel(r"$\sigma_{l,r}$", fontsize=9)
plt.ylabel("空间通信评分", fontsize=9)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.title(
    r"$\sigma_{l,r}$与空间通信评分的关系",
    fontsize=10
)

plt.tight_layout()

# ===============================
# 步骤④：保存图片
# ===============================
plt.savefig(
    "Fig3_sigma_vs_spatial_score_中文.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("✅ 图 3 已生成：Fig3_sigma_vs_spatial_score_中文.png")
