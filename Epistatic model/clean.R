rm(list = ls())
library(readxl)
stem <- as.data.frame(read_excel("XM_Stem_all.xlsx"))
rownames(stem) <- stem[,1]
stem1 <-stem[,-1] 

index <- as.data.frame(read_excel("XMA-XMB.homolog.xlsx")[,1:2])
new_row_names <- c(rbind(index$A, index$B)) 
中
missing_genes <- new_row_names[!new_row_names %in% rownames(stem1)]
if (length(missing_genes) > 0) {
  stop(paste("以下基因在 stem1 中未找到:", paste(missing_genes, collapse = ", ")))
}
result <- stem1[new_row_names, , drop = FALSE]

if (nrow(result) != length(new_row_names)) {
  stop("提取的数据行数与新行名数量不匹配，请检查数据完整性！")
}
stem_dengwei <- read.csv("stem_dengwei.csv")
rownames(stem_dengwei) <- stem_dengwei[,1]
stem_dengwei1 <- stem_dengwei[,-1]
cs <- stem_dengwei1[,1:21]
ch <- stem_dengwei1[,c(1:3,22:39)]

threshold <- 0.6
get_high_zero_indices <- function(data, threshold) {
  row_zero_ratio <- rowSums(data == 0) / ncol(data) 
  high_zero_indices <- which(row_zero_ratio > threshold) 
  return(high_zero_indices)
}
cs_high_zero_indices <- get_high_zero_indices(cs, threshold)
ch_high_zero_indices <- get_high_zero_indices(ch, threshold)


all_high_zero_indices <- unique(c(cs_high_zero_indices, ch_high_zero_indices))

get_gene_pairs_indices <- function(indices) {
  paired_indices <- c()
  for (i in indices) {
    if (i %% 2 == 1) { 
      paired_indices <- c(paired_indices, i, i + 1)
    } else { 
      paired_indices <- c(paired_indices, i, i - 1)
    }
  }
  return(unique(paired_indices))  # 去重
}

final_indices_to_remove <- get_gene_pairs_indices(all_high_zero_indices)

# Step 4: 删除行
cs_filtered <- cs[-final_indices_to_remove, , drop = FALSE]
ch_filtered <- ch[-final_indices_to_remove, , drop = FALSE]

# Step 5: 输出结果
cat("原始 cs 数据行数:", nrow(cs), "\n")
cat("原始 ch 数据行数:", nrow(ch), "\n")
cat("需要删除的等位基因对数:", length(final_indices_to_remove) / 2, "\n")
cat("过滤后 cs 数据行数:", nrow(cs_filtered), "\n")
cat("过滤后 ch 数据行数:", nrow(ch_filtered), "\n")

# 保存结果
write.csv(cs_filtered, "cs.csv", row.names = TRUE)
write.csv(ch_filtered, "ch.csv", row.names = TRUE)

cat("过滤后的数据已保存为 'cs_filtered.csv' 和 'ch_filtered.csv'。\n")


