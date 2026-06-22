
print(sort(cluster_sizes))

# 可选：如果最小 cluster 太小（如 <50），可以过滤（取消注释使用）
min_size <- 50
large_clusters <- names(cluster_sizes[cluster_sizes >= min_size])
keep_cells <- visium$cell_type %in% large_clusters
visium <- visium[, keep_cells]
# 关键修复：把 cell_type 转为 factor（自动掉空水平）
visium$cell_type <- factor(visium$cell_type)

# 5. 准备 CellChat 输入（用原始 counts 层）
data.input <- GetAssayData(visium, assay = "Spatial", layer = "counts")

meta <- data.frame(
  cell_type = visium$cell_type,
  row.names = colnames(visium)
)

cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "cell_type"
)

# 6. 设置小鼠数据库（关键：用 .use 版本，并检查结构）
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

# 推荐：显式设置 ident（避免 group.by 问题）
cellchat <- setIdent(cellchat, ident.use = "cell_type")

# 打印检查 DB 是否完整（可选，但强烈建议跑一下看输出）
str(cellchat@DB)  # 应该看到 $interaction, $geneInfo 等列表，非 NULL
# 或
head(cellchat@DB$interaction)

# 7. subsetData（现在应该能过了）
cellchat <- subsetData(cellchat)  # 如果还报错，加 features = rownames(cellchat@data.signaling) 但通常不需要

# 继续后续步骤
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(
  cellchat,
  type = "triMean",
  raw.use = TRUE,
  population.size = FALSE
)

cellchat <- filterCommunication(cellchat, min.cells = 10)

# 导出
cellchat_lr <- subsetCommunication(cellchat, slot.name = "net", thresh = 0.05)
cellchat_lr_simple <- cellchat_lr[, c("source", "target", "ligand", "receptor", "prob")]
write.csv(cellchat_lr_simple, "cellchat_lr_pairs.csv", row.names = FALSE)

cat("✅ 完成！检查 cellchat@net$prob 是否有值：\n")
print(dim(cellchat@net$prob))
print(range(cellchat@net$prob, na.rm = TRUE))