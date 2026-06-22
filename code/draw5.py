# import pandas as pd
# import matplotlib.pyplot as plt
# from matplotlib_venn import venn2
#
# # 读取 CellChat 文件
# cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
# # 创建 LR 对唯一标识
# cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']
# # 提取唯一 LR 对集合
# cellchat_lrs = set(cellchat_df['lr_pair'].unique())
#
# # 读取 LIANA 文件（SpaADA 基于 LIANA）
# liana_df = pd.read_csv("liana_consensus_lr_pairs.csv")
# # 创建 LR 对唯一标识
# liana_df['lr_pair'] = liana_df['ligand'] + '_' + liana_df['receptor']
# # 提取唯一 LR 对集合
# liana_lrs = set(liana_df['lr_pair'].unique())
#
# print(f"CellChat 唯一 LR 对数: {len(cellchat_lrs)}")
# print(f"SpaADA (LIANA) 唯一 LR 对数: {len(liana_lrs)}")
# print(f"共同 LR 对数: {len(cellchat_lrs & liana_lrs)}")
#
# # 绘制 Venn 图
# plt.figure(figsize=(5, 4))  # 调整画布大小，更紧凑
# venn = venn2(subsets=[cellchat_lrs, liana_lrs], set_labels=('CellChat', 'SpaADA'))
#
# # 设置区域标签（显示数量），字体缩小
# if venn.get_label_by_id('10'):
#     venn.get_label_by_id('10').set_text(f'CellChat only\n{len(cellchat_lrs - liana_lrs)}')
#     venn.get_label_by_id('10').set_fontsize(8)
# if venn.get_label_by_id('01'):
#     venn.get_label_by_id('01').set_text(f'SpaADA only\n{len(liana_lrs - cellchat_lrs)}')
#     venn.get_label_by_id('01').set_fontsize(9)
# if venn.get_label_by_id('11'):
#     venn.get_label_by_id('11').set_text(f'Common\n{len(cellchat_lrs & liana_lrs)}')
#     venn.get_label_by_id('11').set_fontsize(7)
#
# # 设置 set_labels 字体大小（CellChat 和 SpaADA 标签）
# for label in venn.set_labels:
#     label.set_fontsize(9)
#
# # 标题字体稍大，但比默认小
# plt.title('Overlap of Unique LR Pairs between CellChat and SpaADA', fontsize=10)
#
# # 保存高清图
# plt.tight_layout()
# plt.savefig('Fig5_venn_lr_overlap.png', dpi=300, bbox_inches='tight')
# plt.show()
#
# print("Venn 图已保存为 Fig5_venn_lr_overlap.png")
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# ===============================
# ⭐ 设置中文字体（宋体）
# ===============================
plt.rcParams['font.sans-serif'] = ['SimSun']
plt.rcParams['axes.unicode_minus'] = False

# ===============================
# 读取 CellChat 数据
# ===============================
cellchat_df = pd.read_csv("cellchat_lr_pairs.csv")
cellchat_df['lr_pair'] = cellchat_df['ligand'] + '_' + cellchat_df['receptor']
cellchat_lrs = set(cellchat_df['lr_pair'].unique())

# ===============================
# 读取 LIANA（SpaADA）
# ===============================
liana_df = pd.read_csv("liana_consensus_lr_pairs.csv")
liana_df['lr_pair'] = liana_df['ligand'] + '_' + liana_df['receptor']
liana_lrs = set(liana_df['lr_pair'].unique())

print(f"CellChat 唯一 LR 对数: {len(cellchat_lrs)}")
print(f"SpaADA（LIANA）唯一 LR 对数: {len(liana_lrs)}")
print(f"共同 LR 对数: {len(cellchat_lrs & liana_lrs)}")

# ===============================
# 绘制 Venn 图
# ===============================
plt.figure(figsize=(5, 4))

venn = venn2(
    subsets=[cellchat_lrs, liana_lrs],
    set_labels=('CellChat', 'SpaADA')  # 方法名保留英文（论文规范）
)

# ===============================
# ⭐ 区域文字改中文
# ===============================
if venn.get_label_by_id('10'):
    venn.get_label_by_id('10').set_text(
        f'CellChat独有\n{len(cellchat_lrs - liana_lrs)}'
    )
    venn.get_label_by_id('10').set_fontsize(8)

if venn.get_label_by_id('01'):
    venn.get_label_by_id('01').set_text(
        f'SpaADA独有\n{len(liana_lrs - cellchat_lrs)}'
    )
    venn.get_label_by_id('01').set_fontsize(8)

if venn.get_label_by_id('11'):
    venn.get_label_by_id('11').set_text(
        f'共同\n{len(cellchat_lrs & liana_lrs)}'
    )
    venn.get_label_by_id('11').set_fontsize(8)

# ===============================
# 方法标签字体
# ===============================
for label in venn.set_labels:
    label.set_fontsize(9)

# ===============================
# ⭐ 中文标题
# ===============================
plt.title(
    'CellChat与SpaADA的LR对重叠情况',
    fontsize=10
)

plt.tight_layout()

# ===============================
# 保存
# ===============================
plt.savefig(
    'Fig5_venn_LR重叠_中文.png',
    dpi=300,
    bbox_inches='tight'
)

plt.close()

print("✅ Venn 图已保存为 Fig5_venn_LR重叠_中文.png")
