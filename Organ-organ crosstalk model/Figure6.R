
resa <- t(sapply(1:length(rdiploid_chra),function(c)
  as.character(rdiploid_chra[[c]]$type)))
resb <- t(sapply(1:length(rdiploid_chrb),function(c)
  as.character(rdiploid_chrb[[c]]$type)))

colnames(resa)=colnames(resb)=c("fs_ss","fs_leaf","ss_leaf")
table(resa[,1])
table(resb[,1])

library(ggplot2)
library(dplyr)
library(tidyr)

# 加载必要的包
library(ggplot2)
library(dplyr)
# res_type1 <- resa[,2]
# res_type2 <- resb[,2]

get_barplot <- function(res_type1,res_type2){
data <- data.frame(
  Letter = c(names(table(res_type1))),
  Count1 = c(as.numeric(table(res_type1))),  # 第二行（朝上柱状）
  Count2 = c(as.numeric(table(res_type2)))   # 第三行（朝下柱状）
  
  
)

data <- data %>%
  mutate(
    Proportion1 = Count1 / sum(Count1) ,  # 第二行比例
    Proportion2 = -Count2 / sum(Count2)   # 第三行比例
  )

data_long <- data %>%
  pivot_longer(cols = starts_with("Count"),
               names_to = "Type",
               values_to = "Count") %>%
  mutate(
    Proportion = ifelse(Type == "Count1", Proportion1, Proportion2),
    Count = ifelse(Type == "Count2", Count, Count)  # 第三行取负数
  )

tt <- ggplot() +
  # # 绘制柱状图：使用fill映射基因型类型，并设定因子的顺序
  # geom_bar(data = data, aes(x = Letter, y = All2), fill = "#FBC02D", 
  #          stat = "identity", width = 1, color = "black") +
  geom_bar(data = data, aes(x = Letter, y = Proportion1), fill = "#E8F4C8", 
           stat = "identity", width = 1, color = "black") +
  
  # geom_bar(data = data, aes(x = Letter, y = Aa1), fill ="#7B9EBF",
  #          stat = "identity", width = 1, color = "black") +
  geom_bar(data = data, aes(x = Letter, y = Proportion2), fill ="#FFE8C2",
           stat = "identity", width = 1, color = "black") +
  
  
  geom_text(data = data_long, aes(x = Letter, y = Proportion, label = paste0(abs(Count))),
            vjust = ifelse(data_long$Type == "Count1", -0.5, 1.5), size = 10) +
  
  labs(x = NULL, y = NULL, fill = "genetype", title = NULL) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, 0.1), labels = c(0.5,0.4,0.3, 0.2,0.1,  0, 0.1, 0.2, 0.3, 0.4,0.5), limits = c(-0.5, 0.5)) +
  scale_x_discrete(limits = c("+/+","-/-","0/0","+/0","0/+","-/0","0/-","+/-","-/+")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    legend.position = "top"
  )+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()
  )+theme(axis.text.x = element_text( color="black", size=30),axis.text.y =element_text( color="black", size=30),
          panel.border = element_rect(color = "black", size = 1.5, fill = NA),
          axis.ticks.length.y = unit(0.2,"cm"),axis.ticks.length.x = unit(0.2,"cm"))

return(tt)}


fs_ssplot <- get_barplot(resa[,1],resb[,1])
fs_leafplot <- get_barplot(resa[,2],resb[,2])
ss_leafplot <- get_barplot(resa[,3],resb[,3])
pdf("barplot_fs_ss.pdf",width = 10,height = 10)
fs_ssplot
dev.off()


pdf("barplot_fs_leaf.pdf",width = 10,height = 10)
fs_leafplot
dev.off()
pdf("barplot_ss_leaf.pdf",width = 10,height = 10)
ss_leafplot
dev.off()
