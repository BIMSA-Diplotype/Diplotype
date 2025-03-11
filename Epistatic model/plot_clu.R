library(reshape2)
library(idopNetwork)
library(ggplot2)
library(parallel)
library(scales)

fun_clu_plot_byy <- function (result, best.k, label = 10, degree = 1) {
  cluster.result <- result
  times = cluster.result$Time
  times_new = seq(min(times), max(times), length = 30)
  par.mu = cluster.result$mu_par
  k = cluster.result$cluster_number
  alpha = as.numeric(table(cluster.result$cluster[, ncol(cluster.result$cluster)]))
  mu.fit = power_equation(times_new, par.mu[1:k, ])
  colnames(mu.fit) = times_new
  mu.fit = melt(as.matrix(mu.fit))
  colnames(mu.fit) = c("cluster", "x", "y")
  mu.fit$x = as.numeric(as.character(mu.fit$x))
  plot.df = cluster.result$cluster
  X = plot.df[, -ncol(plot.df)]
  plot.df$name = rownames(plot.df)
  colnames(plot.df) = c(times, "cluster", "name")
  plot.df = melt(plot.df, id.vars = c("cluster", "name"))
  colnames(plot.df) = c("cluster", "name", "x", "y")
  plot.df$x = as.numeric(as.character(plot.df$x))
  df.alpha = data.frame(cluster = 1:k, alpha = min(alpha)/alpha * 
                          degree)
  plot.df = merge(plot.df, df.alpha, by = "cluster")
  name.df = data.frame(
    label = paste0("M", 1:best.k, "(", alpha, ")"), 
    x = mean(range(times)), 
    y = sapply(1:k, function(c) max(plot.df$y[plot.df$cluster == c])) * 0.9,
    cluster = 1:best.k
  )
  p = ggplot() + 
    geom_point(plot.df, mapping = aes_string(x = "x", 
                                             y = "y", colour = factor(plot.df$cluster),
                                             group = "name"), show.legend = F, 
               shape = 1,size = 0.1,alpha = 0.2) + 
    geom_line(mu.fit, mapping = aes_string(x = "x", y = "y", colour = factor(mu.fit$cluster)), 
              size = 1.25, show.legend = F) + 
    facet_wrap(~cluster, scales = "free_y",ncol = 8) + 
    scale_alpha_manual(values = df.alpha$alpha) + 
    xlab("Habitat Index") + 
    ylab("Niche Index") + 
    theme(axis.title = element_text(size = 40)) + 
    geom_text(name.df, mapping = aes_string(x = "x", y = "y", 
                                            label = "label"), check_overlap = TRUE, size = 4) + 
    theme_bw() + 
    theme(axis.title = element_text(size = 15), 
          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 14, 
                                                                            hjust = 0), panel.spacing = unit(0, "lines"), plot.margin = unit(c(0.5, 
                                                                                                                                              0.5, 0.5, 0.5), "lines"), strip.background = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_blank(), strip.text = element_blank()) +
    scale_x_continuous(labels = scientific)  # 这里设置x轴为科学计数法
  
  if (is.null(label)) {
    p = p
  }
  else {
    xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
    ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks
    xlabel2 = parse(text = paste(label, "^", xlabel, sep = ""))
    if (ylabel[1] == 0) {
      ylabel2 = parse(text = c(0, paste(label, "^", ylabel[2:length(ylabel)], 
                                        sep = "")))
    }
    else {
      ylabel2 = parse(text = paste(label, "^", ylabel, 
                                   sep = ""))
    }
    p = p + scale_x_continuous(labels = xlabel2) + scale_y_continuous(labels = ylabel2)
  }
  
  return(p)
}


#rm(list = ls())
#s_k64 <- readRDS("s_k57.rds")
#h_k90 <- readRDS("h_k83.rds")

ss <- fun_clu_plot_byy(s_k64,best.k = 57,label = NULL, degree = 1)
hh <- fun_clu_plot_byy(h_k90,best.k = 83,label = NULL, degree = 1)
ggsave("tolerance_clu.png", plot = ss, width = 17, height = 10,dpi = 1200)
ggsave("hardiness_clu.png", plot = hh, width = 20, height = 14,dpi = 1200)








