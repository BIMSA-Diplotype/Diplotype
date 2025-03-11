load(file = "genetype.RData")
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
genetype1 <- genetype[-9847]
res1 <- as.character(sapply(1:(nrow(gene_dat1a)),function(c)rdiploid_res1[[c]][[3]]))
res2 <- as.character(sapply(1:(nrow(gene_dat2a)),function(c)rdiploid_res2[[c]][[3]]))


library(ggplot2)
library(dplyr)
library(tidyr)

# 加载必要的包
library(ggplot2)
library(dplyr)

data <- data.frame(
  Letter = c(names(table(res1)),"-/-"),
  Count1 = c(as.numeric(table(res1)),0),  # 第二行（朝上柱状）
  Count2 = c(as.numeric(table(res2)),0)   # 第三行（朝下柱状）
  
  
)

data <- data %>%
  mutate(
    Proportion1 = -Count1 / sum(Count1) ,  # 第二行比例
    Proportion2 = Count2 / sum(Count2)   # 第三行比例
  )
data <- data %>%
  mutate(
    All1 =data$Proportion1,
    
    AA1 = -sapply(1:9,function(c){
      length(which(genetype1[which(res1==data$Letter[c])]==0))/(length(which(genetype1==0))+length(which(genetype1==4)))+
        
        length(which(genetype1[which(res1==data$Letter[c])]==4))/(length(which(genetype1==0))+length(which(genetype1==4)))
      
    } ),
    All2 =data$Proportion2,
    
    AA2 = sapply(1:9,function(c){
      length(which(genetype[which(res2==data$Letter[c])]==0))/(length(which(genetype==0))+length(which(genetype==4)))+
      
      length(which(genetype[which(res2==data$Letter[c])]==4))/(length(which(genetype==0))+length(which(genetype==4)))
      
    } ),
  )
data_long <- data %>%
  pivot_longer(cols = starts_with("Count"),
               names_to = "Type",
               values_to = "Count") %>%
  mutate(
    Proportion = ifelse(Type == "Count1", Proportion1, Proportion2),
    Count = ifelse(Type == "Count2", Count, Count)  # 第三行取负数
  )
library(ggplot2)
library(patchwork)
pdf("bar_plot.pdf",width = 10,height = 6)
ggplot() +
  # # 绘制柱状图：使用fill映射基因型类型，并设定因子的顺序
  # geom_bar(data = data, aes(x = Letter, y = All2), fill = "#FBC02D", 
  #          stat = "identity", width = 1, color = "black") +
  geom_bar(data = data, aes(x = Letter, y = All2), fill = "#FFF9C4", 
           stat = "identity", width = 1, color = "black") +
  
  # geom_bar(data = data, aes(x = Letter, y = Aa1), fill ="#7B9EBF",
  #          stat = "identity", width = 1, color = "black") +
  geom_bar(data = data, aes(x = Letter, y = All1), fill ="#D6EAF8",
           stat = "identity", width = 1, color = "black") +
  
  
  geom_text(data = data_long, aes(x = Letter, y = Proportion, label = paste0(abs(Count))),
            vjust = ifelse(data_long$Type == "Count2", -0.5, 1.5), size = 10) +
  
  labs(x = NULL, y = NULL, fill = "genetype", title = NULL) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1), labels = c(0.4, 0.3, 0.2, 0.1, 0, 0.1, 0.2, 0.3, 0.4), limits = c(-0.4, 0.4)) +
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

dev.off()

###########################
data <- data.frame(
  Letter = c(names(table(res1)),"-/-")  # 第三行（朝下柱状）
  
  
)


data <- data %>%
  mutate(

    
    AA1 = -sapply(1:9,function(c){
      length(which(genetype1[which(res1==data$Letter[c])]==0))/(length(which(genetype1==0))+length(which(genetype1==4)))+
        
        length(which(genetype1[which(res1==data$Letter[c])]==4))/(length(which(genetype1==0))+length(which(genetype1==4)))
      
    } ),
    CountAA1 = sapply(1:9,function(c){length(which(genetype1[which(res1==data$Letter[c])]==0))+
length(which(genetype1[which(res1==data$Letter[c])]==4))}),
    

    
    AA2 = sapply(1:9,function(c){
      length(which(genetype[which(res2==data$Letter[c])]==0))/(length(which(genetype==0))+length(which(genetype==4)))+
        
        length(which(genetype[which(res2==data$Letter[c])]==4))/(length(which(genetype==0))+length(which(genetype==4)))
      
    } ),
CountAA2= sapply(1:9,function(c){length(which(genetype[which(res2==data$Letter[c])]==0))+
      length(which(genetype[which(res2==data$Letter[c])]==4))})
    
  )
data_long <- data %>%
  pivot_longer(cols = starts_with("Count"),
               names_to = "Type",
               values_to = "Count") %>%
  mutate(
    Proportion = ifelse(Type == "CountAA1", AA1, AA2),
    Count = ifelse(Type == "CountAA2", Count, Count)  # 第三行取负数
  )
pdf("bar_plot1.pdf",width = 10,height = 6)
ggplot() +

  geom_bar(data = data, aes(x = Letter, y = AA2), fill = "#FFF9C4", 
           stat = "identity", width = 1, color = "black") +
  
  geom_bar(data = data, aes(x = Letter, y = AA1), fill ="#D6EAF8",
           stat = "identity", width = 1, color = "black")  +
  geom_text(data = data_long, aes(x = Letter, y = Proportion, label = paste0(abs(Count))),
            vjust = ifelse(data_long$Type == "CountAA2", -0.5, 1.5), size = 10) +
  
  labs(x = NULL, y = NULL, fill = "genetype", title = NULL) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1), labels = c(0.4, 0.3, 0.2, 0.1, 0, 0.1, 0.2, 0.3, 0.4), limits = c(-0.4, 0.4)) +
  scale_x_discrete(limits = c("+/+","-/-","0/0","+/0","0/+","-/0","0/-","+/-","-/+"),labels =NULL,position = "top") +
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

dev.off()









