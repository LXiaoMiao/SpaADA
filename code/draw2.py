# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# # 读取结果
# df = pd.read_csv("adaptive_distance_spatial_scores.csv")
#
# plt.figure(figsize=(6, 4))
#
# sns.histplot(
#     df["sigma_lr"],
#     bins=40,
#     kde=True,
#     color="#4C72B0",
#     edgecolor="black",
# )
#
# plt.xticks(fontsize=7)
# plt.yticks(fontsize=7)
#
# plt.xlabel(r"$\sigma_{l,r}$", fontsize=9)
# plt.ylabel("Number of LR pairs", fontsize=9)
# plt.title(r"Distribution of $\sigma_{l,r}$", fontsize=10)
#
# plt.tight_layout()
# plt.savefig("Fig2_sigma_distribution.png", dpi=300)
# plt.show()
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ===============================
# ⭐ 设置中文字体（宋体）
# ===============================
plt.rcParams['font.sans-serif'] = ['SimSun']
plt.rcParams['axes.unicode_minus'] = False

# 读取结果
df = pd.read_csv("adaptive_distance_spatial_scores.csv")

plt.figure(figsize=(6, 4))

sns.histplot(
    df["sigma_lr"],
    bins=40,
    kde=True,
    color="#4C72B0",
    edgecolor="black",
)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

# ===============================
# ⭐ 中文标签
# ===============================
plt.xlabel(r"$\sigma_{l,r}$", fontsize=9)
plt.ylabel("LR对数量", fontsize=9)
plt.title(r"$\sigma_{l,r}$的分布", fontsize=10)

plt.tight_layout()
plt.savefig("Fig2_sigma_distribution_中文.png", dpi=300)
