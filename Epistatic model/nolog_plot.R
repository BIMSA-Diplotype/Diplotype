rm(list = ls())
library(idopNetwork)
library(reshape2)
library(scales)

cs_fit_nolog <- readRDS("cs_fit_nolog.rds")
ch_fit_nolog <- readRDS("ch_fit_nolog.rds")

cs <- read.csv("cs.csv")[,-1]
ch <- read.csv("ch.csv")[,-1]


library(ggplot2)
library(grid)


stem_time_points <- c(0, 2, 4, 8, 12, 24, 48)

#gene_pairs <- list( c(1001, 1002),c(2001, 2002),c(4005, 4006),c(8005, 8006))
#gene_pairs <- list( c(1, 2),c(21, 22),c(45, 46),c(8665, 8666))
gene_pairs <- list(c(3,4),c(26669,26670),c(16873,16874),c(8099,8100))

stem_plot_data <- data.frame()
stem_average_data <- data.frame()

for (i in seq_along(gene_pairs)) {
  gene_pair <- gene_pairs[[i]]
  gene_index_1 <- gene_pair[1]
  gene_index_2 <- gene_pair[2]
  
  gene_expression_1 <- cs[gene_index_1, ]
  expression_values_1 <- lapply(seq(1, ncol(cs), by = 3), function(i) {
    as.numeric(gene_expression_1[i:(i+2)]) 
  })
  average_expression_1 <- sapply(expression_values_1, mean)

  gene_expression_2 <- cs[gene_index_2, ]
  expression_values_2 <- lapply(seq(1, ncol(cs), by = 3), function(i) {
    as.numeric(gene_expression_2[i:(i+2)]) 
  })
  average_expression_2 <- sapply(expression_values_2, mean)

  plot_data_1 <- data.frame(
    Time = rep(stem_time_points, each = 3),
    Expression = unlist(expression_values_1),
    Group = rep(stem_time_points, each = 3),
    Gene = rep(paste("Gene", gene_index_1), length(unlist(expression_values_1))),
    Pair = paste("Group", i) 
  )
  
  plot_data_2 <- data.frame(
    Time = rep(stem_time_points, each = 3),
    Expression = unlist(expression_values_2),
    Group = rep(stem_time_points, each = 3),
    Gene = rep(paste("Gene", gene_index_2), length(unlist(expression_values_2))),
    Pair = paste("Group", i) 
  )
  
  average_data_1 <- data.frame(
    Time = stem_time_points,
    AverageExpression = average_expression_1,
    Gene = paste("Gene", gene_index_1),
    Pair = paste("Group", i)  
  )
  
  average_data_2 <- data.frame(
    Time = stem_time_points,
    AverageExpression = average_expression_2,
    Gene = paste("Gene", gene_index_2),
    Pair = paste("Group", i) 
  )
  stem_plot_data <- rbind(stem_plot_data, plot_data_1, plot_data_2)
  stem_average_data <- rbind(stem_average_data, average_data_1, average_data_2)
}

y_min <- min(stem_plot_data$Expression, na.rm = TRUE)
y_max <- max(stem_plot_data$Expression, na.rm = TRUE)

library(ggplot2)
library(patchwork)

color_mapping <- unlist(lapply(seq_along(gene_pairs), function(i) {
  gene_pair <- gene_pairs[[i]]
  setNames(
    c("#ED342F", "#178642"),  # 红色和绿色
    c(paste("Gene", gene_pair[1]), paste("Gene", gene_pair[2]))  # 基因名称
  )
}))

print(color_mapping)

x_range <- range(stem_time_points)

breaks_6 <- seq(x_range[1], x_range[2], length.out = 7)

stem_plot <- ggplot(stem_plot_data, aes(x = Time, y = Expression, color = Gene)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line(data = stem_average_data, aes(x = Time, y = AverageExpression, color = Gene), size = 1.2) +
  labs(
    title = NULL,
    x = "Time (hour)",
    y = NULL
  ) +
  scale_color_manual(values = color_mapping) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_blank(), 
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "#c6ccdc"), 
    plot.background = element_rect(fill = "white"), 
    axis.line = element_line(color = "black", linewidth = 0.75),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(
    limits = c(x_range[1], x_range[2]), 
    breaks = breaks_6, 
    labels = as.character(breaks_6) 
  ) +
  facet_wrap(~ Pair, nrow = 4, scales = "free_y") +
  geom_hline(data = subset(stem_plot_data, Pair %in% unique(stem_plot_data$Pair)[1:3]), 
             aes(yintercept = 0), color = "black", linewidth = 0.75) +
  scale_y_continuous(
    limits = c(0, NA) 
  )

stem_plot


time_points <- c(0, 5, 10, 15, 20, 25, 30)


plot_data <- data.frame()
average_data <- data.frame()

for (i in seq_along(gene_pairs)) {
  gene_pair <- gene_pairs[[i]]
  gene_index_1 <- gene_pair[1]
  gene_index_2 <- gene_pair[2]
 
  gene_expression_1 <- ch[gene_index_1, ]
  expression_values_1 <- lapply(seq(1, ncol(ch), by = 3), function(i) {
    as.numeric(gene_expression_1[i:(i+2)])
  })
  average_expression_1 <- sapply(expression_values_1, mean)
  
  gene_expression_2 <- ch[gene_index_2, ]
  expression_values_2 <- lapply(seq(1, ncol(ch), by = 3), function(i) {
    as.numeric(gene_expression_2[i:(i+2)]) 
  })
  average_expression_2 <- sapply(expression_values_2, mean)
 
  plot_data_1 <- data.frame(
    Time = rep(time_points, each = 3),
    Expression = unlist(expression_values_1),
    Group = rep(time_points, each = 3),
    Gene = rep(paste("Gene", gene_index_1), length(unlist(expression_values_1))),
    Pair = paste("Group", i)
  )
  
  plot_data_2 <- data.frame(
    Time = rep(time_points, each = 3),
    Expression = unlist(expression_values_2),
    Group = rep(time_points, each = 3),
    Gene = rep(paste("Gene", gene_index_2), length(unlist(expression_values_2))),
    Pair = paste("Group", i) 
  )
  
  average_data_1 <- data.frame(
    Time = time_points,
    AverageExpression = average_expression_1,
    Gene = paste("Gene", gene_index_1),
    Pair = paste("Group", i)
  )
  
  average_data_2 <- data.frame(
    Time = time_points,
    AverageExpression = average_expression_2,
    Gene = paste("Gene", gene_index_2),
    Pair = paste("Group", i) 
  )

  plot_data <- rbind(plot_data, plot_data_1, plot_data_2)
  average_data <- rbind(average_data, average_data_1, average_data_2)
}

y_min <- min(plot_data$Expression, na.rm = TRUE)
y_max <- max(plot_data$Expression, na.rm = TRUE)

color_mapping <- unlist(lapply(seq_along(gene_pairs), function(i) {
  gene_pair <- gene_pairs[[i]]
  setNames(
    c("#ED342F", "#178642"),  # 红色和绿色
    c(paste("Gene", gene_pair[1]), paste("Gene", gene_pair[2]))  # 基因名称
  )
}))

print(color_mapping)

fruit_plot <- ggplot(plot_data, aes(x = Time, y = Expression, color = Gene)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line(data = average_data, aes(x = Time, y = AverageExpression, color = Gene), size = 1.2) +
  labs(
    title = NULL,
    x = "Temperature (℃)",
    y = NULL
  ) +
  scale_color_manual(values = color_mapping) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_blank(), 
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "#c6ccdc"), 
    plot.background = element_rect(fill = "white"), 
    axis.line = element_line(color = "black", linewidth = 0.75),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(
    limits = c(min(plot_data$Time), max(plot_data$Time)),
    breaks = time_points,
    labels = c("22", "4", "-5", "-10", "-15", "-20", "-25")
  ) +
  facet_wrap(~ Pair, nrow = 4, scales = "free_y") +
  geom_hline(data = subset(plot_data, Pair %in% unique(plot_data$Pair)[1:3]), 
             aes(yintercept = 0), color = "black", linewidth = 0.75) +
  #geom_vline(xintercept = max(time_points) + 2, color = "black", linewidth = 0.75) +
  scale_y_continuous(
    limits = c(0, NA) 
  )

fruit_plot


combined_plot1 <- stem_plot + fruit_plot

print(combined_plot1)

cs_fit_nolog <- readRDS("cs_fit_nolog.rds")
ch_fit_nolog <- readRDS("ch_fit_nolog.rds")
#gene_pairs <- list( c(1001, 1002),c(2001, 2002),c(4005, 4006),c(8005, 8006))
#gene_pairs <- list( c(1, 2),c(21, 22),c(45, 46),c(8665, 8666))
#gene_pairs <- list(c(3,4),c(26669,26670),c(16873,16874),c(8099,8100))

gene_pairs_list_1 <- list(c(1001, 1002), c(2001, 2002), c(4005, 4006), c(8005, 8006))
gene_pairs_list_2 <- list(c(1, 2), c(7, 8), c(45, 46), c(8665, 8666))
gene_pairs_list_3 <- list(c(3, 4), c(26669, 26670), c(16873, 16874), c(8099, 8100))
name_use <-  read.csv("cs.csv")
gene_names_list_1 <- lapply(gene_pairs_list_1, function(pair) {
  gene_name_1 <- name_use$X[pair[1]]
  gene_name_2 <- name_use$X[pair[2]]
  c(gene_name_1, gene_name_2) 
})

gene_names_list_2 <- lapply(gene_pairs_list_2, function(pair) {
  gene_name_1 <- name_use$X[pair[1]]
  gene_name_2 <- name_use$X[pair[2]]
  c(gene_name_1, gene_name_2) 
})

gene_names_list_3 <- lapply(gene_pairs_list_3, function(pair) {
  gene_name_1 <- name_use$X[pair[1]]
  gene_name_2 <- name_use$X[pair[2]]
  c(gene_name_1, gene_name_2) 
})
all_gene_names <- unique(c(unlist(gene_names_list_1), unlist(gene_names_list_2), unlist(gene_names_list_3)))


library(dplyr)
library(scales)
library(ggplot2)
library(reshape2)
library(patchwork)
#gene_pairs = gene_names_list_2
power_equation_plot <- function(result, gene_pairs, panel_bg_color = "lightblue") { 
  data1 = result[[2]]  # 原数据
  data2 = result[[4]]  # 拟合数据q

  df_original_all <- data.frame()
  df_fit_all <- data.frame()

  for (pair in gene_pairs) { 
    data3 <- as.matrix(data1[pair, ])
    data4 <- as.matrix(data2[pair, ]) 
    
    df_original = melt(data3)
    df_fit = melt(data4)
    
    df_original$gene_pair <- paste("Gene", pair[1], "vs", pair[2], sep = "_")
    df_fit$gene_pair <- paste("Gene", pair[1], "vs", pair[2], sep = "_")
    
    df_original$gene_source <- "Gene1"
    df_fit$gene_source <- "Gene1"
    df_original$gene_source[seq_along(df_original$value) %% 2 == 0] <- "Gene2"
    df_fit$gene_source[seq_along(df_fit$value) %% 2 == 0] <- "Gene2"

    df_original_all <- rbind(df_original_all, df_original)
    df_fit_all <- rbind(df_fit_all, df_fit)
  }
  df_original_all$gene_pair <- factor(df_original_all$gene_pair, levels = sapply(gene_pairs, function(pair) paste("Gene", pair[1], "vs", pair[2], sep = "_")))
  df_fit_all$gene_pair <- factor(df_fit_all$gene_pair, levels = sapply(gene_pairs, function(pair) paste("Gene", pair[1], "vs", pair[2], sep = "_")))

  p <- ggplot() +
    geom_point(data = df_original_all, 
               mapping = aes(x = Var2, y = value, colour = gene_source), 
               alpha = 0.6, size = 1.5) +
    geom_line(data = df_fit_all, 
              mapping = aes(x = Var2, y = value, colour = gene_source), 
              size = 1.2) +
    theme_bw() +
    facet_wrap(~gene_pair, nrow = 4, scales = "free_y") +  # 每个子图有独立 Y 轴范围
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      strip.text = element_blank(),
      strip.background = element_rect(color = "black", size = 0.75),
      panel.spacing = unit(0, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = panel_bg_color)
    ) +
    scale_color_manual(values = c("Gene1" = "#ED342F", "Gene2" = "#178642")) +
    xlab("Expression Index") +
    scale_y_continuous(
      limits = c(0, NA), 
    )
  
  p
  
  return(p)
}

#"#FFF9C4"和"#D6EAF8"
cs_fit_plot <- power_equation_plot(cs_fit_nolog, gene_names_list_3,panel_bg_color = "#f4f3ee")
ch_fit_plot <- power_equation_plot(ch_fit_nolog, gene_names_list_3,panel_bg_color = "#f4f3ee")

combined_fit_plot <- cs_fit_plot + ch_fit_plot + plot_layout(ncol = 2)

print(combined_fit_plot)

stem_plot <- stem_plot + theme(plot.margin = margin(0.2, 0, 0.2, 0.2, "cm"))

fruit_plot <- fruit_plot + theme(plot.margin = margin(1, 0, 0.2, 0.2, "cm"))

cs_fit_plot <- cs_fit_plot + theme(plot.margin = margin(0.2, 1, 0.2, 0, "cm"))

ch_fit_plot <- ch_fit_plot + theme(plot.margin = margin(0.2, 0, 0.2, 0, "cm"))
library(grid) 
combined_plot <- stem_plot | cs_fit_plot | fruit_plot| ch_fit_plot
ggsave("combined_plot.png", plot = combined_plot, width =17 , height = 10, dpi = 1200)











