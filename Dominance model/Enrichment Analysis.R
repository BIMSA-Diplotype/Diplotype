load(file = "rdiploid_res1.RData")
load(file = "rdiploid_res2.RData")
remove_errors_recursive <- function(lst) {
  if (is.list(lst)) {
    # 对每个子元素递归调用本函数，并过滤掉 try-error 对象
    lst <- lapply(lst, remove_errors_recursive)
    lst <- lst[!sapply(lst, function(x) inherits(x, "try-error"))]
  }
  return(lst)
}
rdiploid_res1 <-remove_errors_recursive(rdiploid_res1)
gene_dat1a <- gene_dat1a[-9847,]
gene_dat1b <- gene_dat1b[-9847,]

res1 <- as.character(sapply(1:(nrow(gene_dat1a)),function(x)rdiploid_res1[[x]][[3]]))
res2 <- as.character(sapply(1:(nrow(gene_dat2a)),function(x)rdiploid_res2[[x]][[3]]))
type <- c("+/+","-/-","0/0","+/0","0/+","-/0","0/-","+/-","-/+")
get_geneid <- function(res1,gene_dat1a,gene_dat1b){
  id_res <- list()
  
  for (i in 1:9) {
    id_res[[i]] <- c( rownames( gene_dat1a)[which(res1==type[i])],rownames( gene_dat1b)[which(res1==type[i])])
    
  }
  return(id_res)
  
}

id_res1 <- get_geneid(res1,gene_dat1a,gene_dat1b)

id_res2 <- get_geneid(res2,gene_dat2a,gene_dat2b)
for (i in 1:9) {
  names(id_res1 )[i] <- type[i]
  names(id_res2)[i]  <- type[i]
  
}

saveRDS(id_res1, "id_res_ch.rds")
saveRDS(id_res2, "id_res_cs.rds")
names(id_res_ch[[1]]) <- "++"

library(org.PMVXM.eg.db)
library(clusterProfiler)
library(stringr)
library(fanyi)
library(dplyr)
library(ggplot2)
library(showtext)
library(pbapply)
library(tidyr)
library(pheatmap)
rm(list = ls())
id_res_ch=readRDS("id_res_ch.rds")
id_res_cs=readRDS("id_res_cs.rds")

############################################按照bpmfcc排序##################
library(dplyr)

#############################################################

get_heatmat2 <- function(id_gene,n_pos=20){
  # id_gene=id_res_ch
  # n_pos=15
  go_results_all <- pblapply(id_gene[1:length(id_gene)], function(gene_list) {
    result <- enrichGO(gene = gene_list,
                       OrgDb = org.PMVXM.eg.db, 
                       keyType = 'GID',
                       ont = "ALL", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 1, 
                       qvalueCutoff = 1,
                       readable = FALSE)
    
    # 只保留感兴趣的列，将 p.adjust 替换为 pvalue
    result <- result[, c("ID", "Description", "pvalue", "ONTOLOGY")]
    
    return(result)
  })
  
  go_terms_list_all <- lapply(go_results_all, function(ego) {
    ego_df <- ego[, c("Description", "pvalue", "ONTOLOGY")]  # 选择感兴趣的列
    ego_df
  })
  names(go_terms_list_all)
  type <- c("+/+","-/-","0/0","+/0","0/+","-/0","0/-","+/-","-/+")
  
  combined_df_all <- sapply(1:length(go_terms_list_all),function(c){
    if(is.null( nrow(go_terms_list_all[[c]]))){
      go_res <-NULL
      
      
    }else{
      go_res <-cbind(go_terms_list_all[[c]],genelist=rep(type[c],nrow(go_terms_list_all[[c]])))
      go_res <- go_res[order(go_res$pvalue,decreasing = F)[1:n_pos],]
      ontology_order <- c("BP" = 1, "MF" = 2, "CC" = 3)
      
      go_res<- go_res %>%
        mutate(ontology_rank = ontology_order[ONTOLOGY]) %>%  # 映射 'ONTOLOGY' 列到数字
        arrange(ontology_rank) %>%  # 按照 'ontology_rank' 排序
        select(-ontology_rank)  # 删除临时的 'on
    }
    return( go_res)
  })
  
  
  
  combined_df_all <- do.call(rbind, combined_df_all)
  
  go_matrix_all <- combined_df_all %>%
    pivot_wider(names_from = genelist, 
                values_from = pvalue, 
                values_fill = list(pvalue =NA))
  
  
  
  
  
  # go_matrix_all_filtered <- go_matrix_all %>%
  #   select(-Description, -ONTOLOGY)  # 只保留 pvalue 列
  # 
  # # 2. 计算每个 GO 术语的 p 值非 NA 的比例
  # p_value_ratios <- apply(go_matrix_all_filtered, 1, function(row) {
  #   sum(!is.na(row)) / length(row)  # 计算 p 值非 NA 的比例
  # })
  # 
  # # 将计算结果添加到原始数据中，作为一列
  # go_matrix_all$p_value_ratio <- p_value_ratios
  # 
  # # 3. 筛选出 p 值非 NA 比例较高的 GO 术语
  # # 设置一个比例阈值（例如 0.6 表示在超过 60% 的类别中有有效 p 值）
  # 
  
  high_p_value_ratio_terms <- go_matrix_all
  # ontology_order <- c("BP" = 1, "MF" = 2, "CC" = 3)
  # 
  # high_p_value_ratio_terms_sorted <- high_p_value_ratio_terms %>%
  #   mutate(ontology_rank = ontology_order[ONTOLOGY]) %>%  # 映射 'ONTOLOGY' 列到数字
  #   arrange(ontology_rank) %>%  # 按照 'ontology_rank' 排序
  #   select(-ontology_rank)  # 删除临时的 'on
  high_p_value_ratio_terms_sorted <- high_p_value_ratio_terms
  high_p_value_ratio_terms_sorted[is.na(high_p_value_ratio_terms_sorted)] <- 1
  
  mat <- as.matrix(high_p_value_ratio_terms_sorted[, c("+/+", "0/0", "+/0", "0/+", "-/0", "0/-" ,"+/-" ,"-/+")])
  rownames(mat) <- high_p_value_ratio_terms_sorted$Description
  
  mat <- cbind("+/+"=mat[,1],"-/-"=rep(1,nrow(mat)),mat[,-1])
  dim(mat)
  # 创建行注释（根据 ONTOLOGY 分类）
  annotation_row <- data.frame(
    Ontology = factor(high_p_value_ratio_terms_sorted$ONTOLOGY),
    row.names = rownames(mat)
  )
  
  
  
  
  return(list(mat=mat,annotation_row=annotation_row))
  
}
res_ch2 <- get_heatmat2(id_res_ch,n_pos=15)
res_cs2 <- get_heatmat2(id_res_cs,n_pos=15)
names(id_res_cs)
pheatmap(
  res_ch2[[1]],
  color =  colorRampPalette(c( "red","lightgray", "blue"))(1000),# 颜色映射
  scale = "none",          # 不进行行列标准化
  cluster_cols = FALSE,    # 不聚类列
  cluster_rows = FALSE,     # 按 ONTOLOGY 聚类行
  annotation_row = res_ch2[[2]],
  show_rownames = T,    # 显示行名
  show_colnames = TRUE,    # 显示列名
  fontsize_row = 8,        # 行名字体大小
  fontsize_col = 20,       # 列名字体大小
  angle_col = 0, 
  filename = "heatmap_ch.pdf",  # 指定输出路径和文件名
  width = 10,                # 自定义宽度（英寸）
  height = 12
  
)
pheatmap(
  res_cs2[[1]],
  color =  colorRampPalette(c( "red","lightgray", "blue"))(1000),# 颜色映射
  scale = "none",          # 不进行行列标准化
  cluster_cols = FALSE,    # 不聚类列
  cluster_rows = FALSE,     # 按 ONTOLOGY 聚类行
  annotation_row = res_cs2[[2]],
  show_rownames = T,    # 显示行名
  show_colnames = TRUE,    # 显示列名
  fontsize_row = 8,        # 行名字体大小
  fontsize_col = 20,       # 列名字体大小
  angle_col = 0, 
  filename = "heatmap_cs.pdf",  # 指定输出路径和文件名
  width = 10,                # 自定义宽度（英寸）
  height = 12
  
)
rownames( res_ch2[[1]])[!is.na(match(rownames( res_ch2[[1]]),rownames( res_cs2[[1]])))]
rownames( res_cs2[[1]])[!is.na(match(rownames( res_cs2[[1]]),rownames( res_ch2[[1]])))]

rownames( res_cs2[[1]])[1:10]

rownames( res_ch2[[1]])[1:10]



###########################################################################
get_heatmat <- function(id_gene,des_order=T,#是否对功能按首字母重排序
                        des_sample=T,#是否随机取点
                        n_sample=100,#随机取多少点
                        threshold =0#阈值
                        
)
{
  go_results_all <- pblapply(id_gene[1:length(id_gene)], function(gene_list) {
    result <- enrichGO(gene = gene_list,
                       OrgDb = org.PMVXM.eg.db, 
                       keyType = 'GID',
                       ont = "ALL", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05,
                       readable = FALSE)
    
    # 只保留感兴趣的列，将 p.adjust 替换为 pvalue
    result <- result[, c("GeneRatio", "Description", "pvalue", "ONTOLOGY","Count")]
    
    return(result)
  })
  
  go_terms_list_all <- lapply(go_results_all, function(ego) {
    ego_df <- ego[, c("Description", "pvalue", "ONTOLOGY")]  # 选择感兴趣的列
    ego_df
  })
  type <- c("+/+","-/-","0/0","+/0","0/+","-/0","0/-","+/-","-/+")
  
  combined_df_all <- sapply(1:length(go_terms_list_all),function(c){
    if(is.null( nrow(go_terms_list_all[[c]]))){
      go_res <-NULL
      
      
    }else{
      go_res <-cbind(go_terms_list_all[[c]],genelist=rep(type[c],nrow(go_terms_list_all[[c]])))}
    return( go_res)
  })
  
  combined_df_all <- do.call(rbind, combined_df_all)
  
  go_matrix_all <- combined_df_all %>%
    pivot_wider(names_from = genelist, 
                values_from = pvalue, 
                values_fill = list(pvalue =NA))
  
  
  
  
  if(des_order==T){
    go_matrix_all <- go_matrix_all[order(go_matrix_all$Description), ]
    
  }
  
  # 重新排序，按 GO 术语排列
  #go_matrix_all <- go_matrix_all[order(go_matrix_all$Description), ]
  go_matrix_all_filtered <- go_matrix_all %>%
    select(-Description, -ONTOLOGY)  # 只保留 pvalue 列
  
  # 2. 计算每个 GO 术语的 p 值非 NA 的比例
  p_value_ratios <- apply(go_matrix_all_filtered, 1, function(row) {
    sum(!is.na(row)) / length(row)  # 计算 p 值非 NA 的比例
  })
  
  # 将计算结果添加到原始数据中，作为一列
  go_matrix_all$p_value_ratio <- p_value_ratios
  
  # 3. 筛选出 p 值非 NA 比例较高的 GO 术语
  # 设置一个比例阈值（例如 0.6 表示在超过 60% 的类别中有有效 p 值）
  
  
  high_p_value_ratio_terms <- go_matrix_all[go_matrix_all$p_value_ratio >= threshold, ]
  ontology_order <- c("BP" = 1, "MF" = 2, "CC" = 3)
  
  high_p_value_ratio_terms_sorted <- high_p_value_ratio_terms %>%
    mutate(ontology_rank = ontology_order[ONTOLOGY]) %>%  # 映射 'ONTOLOGY' 列到数字
    arrange(ontology_rank) %>%  # 按照 'ontology_rank' 排序
    select(-ontology_rank)  # 删除临时的 'on
  high_p_value_ratio_terms_sorted[is.na(high_p_value_ratio_terms_sorted)] <- 1
  if(des_sample==T){
    high_p_value_ratio_terms_sorted <- high_p_value_ratio_terms_sorted[
      sort(sample(1:nrow(high_p_value_ratio_terms_sorted),n_sample)),]
  }
  mat <- as.matrix(high_p_value_ratio_terms_sorted[, c("+/+", "+/-", "-/+", "+/0", "-/0", "0/+", "0/-", "0/0")])
  rownames(mat) <- high_p_value_ratio_terms_sorted$Description
  
  mat <- cbind("+/+"=mat[,1],"-/-"=rep(1,nrow(mat)),mat[,-1])
  dim(mat)
  # 创建行注释（根据 ONTOLOGY 分类）
  annotation_row <- data.frame(
    Ontology = factor(high_p_value_ratio_terms_sorted$ONTOLOGY),
    row.names = rownames(mat)
  )
  return(list(mat=mat,annotation_row=annotation_row))
  
}