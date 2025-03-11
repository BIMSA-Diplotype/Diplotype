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
res1 <- as.character(sapply(1:(nrow(gene_dat1a)),function(x)rdiploid_res1[[x]][[3]]))
res2 <- as.character(sapply(1:(nrow(gene_dat2a)),function(x)rdiploid_res2[[x]][[3]]))


nam <- sapply(1:nrow(gene_dat1a),function(c)substr(rownames(gene_dat1a)[c], 8, 9))
table(nam)
cumsum(table(nam))

nam_w <- c(0,cumsum(table(nam)))





nam2 <- sapply(1:nrow(gene_dat2a),function(c)substr(rownames(gene_dat2a)[c], 8, 9))


nam_w2 <- c(0,cumsum(table(nam2)))



# 加载必要的包
library(dplyr)
library(ggplot2)
library(tidyr)



chromosomes <- lapply(1:8, function(i) {
  data.frame(
    chromosome = paste0("Chr", i), 
    category = c(res1[(nam_w[i]+1):nam_w[i+1]]),
    value =1:length(c(res1[(nam_w[i]+1):nam_w[i+1]]))
  )
})



all_data <- do.call(rbind, chromosomes)

# 设置间隔角度
gap_angle <- 10  # 每个染色体之间的间隔角度（单位：度）
radius1 <- 20
radius2 <- 18
line_gap <- 1
total_gap_angle <- gap_angle * 8  # 总间隔角度
total_angle <- 360
remaining_angle <- total_angle - total_gap_angle  # 剩余角度分配给染色体

chromosome_angles <- all_data  %>%            # 只计算外环（单个染色体）的点数
  group_by(chromosome) %>%
  summarise(total_points = n()) %>%      # 计算单个染色体的点数
  ungroup() %>%
  mutate(
    total_data_points = sum(total_points),          # 总点数
    proportion = total_points / total_data_points,  # 每个染色体的占比
    angle_span = proportion * remaining_angle       # 每个染色体的实际角度
  )

chromosome_angles <- chromosome_angles %>%
  arrange(chromosome) %>%  # 确保按照染色体顺序排列
  mutate(
    start_angle = cumsum(lag(angle_span + gap_angle, default = 0)),  # 每个染色体的起始角度
    end_angle = start_angle + angle_span                            # 每个染色体的结束角度
  )

# 1. 合并染色体的角度信息
all_data <- all_data %>%
  left_join(chromosome_angles, by = "chromosome")  # 将 `chromosome_angles` 中的角度信息加入到 `all_data`

# 2. 按染色体分组，计算每个点的角度和半径
all_data <- all_data %>%
  group_by(chromosome) %>%
  mutate(
    # 计算每个数据点的实际角度：
    # - value / total_points: 数据点在当前染色体中的相对位置 (0 到 1)
    # - angle_span: 当前染色体分配的总角度
    # - start_angle: 当前染色体的起始角度
    angle = value / total_points * angle_span + start_angle,
    
    # 根据环类型（Outer 或 Inner）设置不同半径
    radius = 11
  )
type <- c("+/+","-/-","+/-","-/+","+/0","-/0","0/+","0/-","0/0")
point_res=lapply(1:8,function(c)as.matrix(filter(all_data,
                                                 category==type[c])[,10]))
dat1 <- data.frame(x=rep(c(radius1,radius2),times=nrow(point_res[[1]])),
                   y=rep(as.numeric(point_res[[1]]),each=2),group=rep(1:nrow(point_res[[1]]),each=2))

dat3 <- data.frame(x=rep(c(radius1,radius2-line_gap),times=nrow(point_res[[3]])),
                   y=rep(as.numeric(point_res[[3]]),each=2),group=rep(1:nrow(point_res[[3]]),each=2))

dat4 <- data.frame(x=rep(c(radius1+line_gap,radius2),times=nrow(point_res[[4]])),
                   y=rep(as.numeric(point_res[[4]]),each=2),group=rep(1:nrow(point_res[[4]]),each=2))

dat5 <- data.frame(x=rep(c(radius2,radius2-line_gap),times=nrow(point_res[[5]])),
                   y=rep(as.numeric(point_res[[5]]),each=2),group=rep(1:nrow(point_res[[5]]),each=2))

dat6 <- data.frame(x=rep(c(radius2,radius2-line_gap),times=nrow(point_res[[6]])),
                   y=rep(as.numeric(point_res[[6]]),each=2),group=rep(1:nrow(point_res[[6]]),each=2))

dat7 <- data.frame(x=rep(c(radius1+line_gap,radius1),times=nrow(point_res[[7]])),
                   y=rep(as.numeric(point_res[[7]]),each=2),group=rep(1:nrow(point_res[[7]]),each=2))
dat8 <- data.frame(x=rep(c(radius1+line_gap,radius1),times=nrow(point_res[[8]])),
                   y=rep(as.numeric(point_res[[8]]),each=2),group=rep(1:nrow(point_res[[8]]),each=2))
siz=0.3
alp=1
gg <- ggplot()+geom_line(
  data = data.frame(dat3),  # 只绘制正常点
  aes(x , y,group=group ),
  col="green",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat1),  # 只绘制正常点
  aes(x , y,group=group ),
  col="red",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat4),  # 只绘制正常点
  aes(x , y,group=group ),
  col="green",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat5),  # 只绘制正常点
  aes(x , y,group=group ),
  col="red",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat6),  # 只绘制正常点
  aes(x , y,group=group ),
  col="blue",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat7),  # 只绘制正常点
  aes(x , y,group=group ),
  col="red",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat8),  # 只绘制正常点
  aes(x , y,group=group ),
  col="blue",alpha=alp,size = siz
)+
  scale_y_continuous(
    breaks = seq(0, 360, by = 30),  
    limits = c(0, 360)              
  )+  geom_line(
    data = data.frame(x=radius1,y=as.numeric(sapply(1:8,function(c)c(chromosome_angles[c,6:7],NA)))),  # 只绘制正常点
    aes(x , y ),
    size = 4
  )+geom_line(
    data = data.frame(x=radius2,y=as.numeric(sapply(1:8,function(c)c(chromosome_angles[c,6:7],NA)))),  # 只绘制正常点
    aes(x , y ),
    size = 4
  )+
  coord_polar(theta = "y", clip = "off") +
  scale_x_continuous(limits = c(1, 22))   +
  theme_void() 
pdf("Figure2A1.pdf",width = 50,height = 50)
gg
dev.off() 

#####################################################################
chromosomes <- lapply(1:8, function(i) {
  data.frame(
    chromosome = paste0("Chr", i), 
    category = c(res2[(nam_w2[i]+1):nam_w2[i+1]]),
    value =1:length(c(res2[(nam_w2[i]+1):nam_w2[i+1]]))
  )
})



all_data <- do.call(rbind, chromosomes)

# 设置间隔角度
gap_angle <- 10  # 每个染色体之间的间隔角度（单位：度）
radius1 <- 20
radius2 <- 18
line_gap <- 1
total_gap_angle <- gap_angle * 8  # 总间隔角度
total_angle <- 360
remaining_angle <- total_angle - total_gap_angle  # 剩余角度分配给染色体

chromosome_angles <- all_data  %>%            # 只计算外环（单个染色体）的点数
  group_by(chromosome) %>%
  summarise(total_points = n()) %>%      # 计算单个染色体的点数
  ungroup() %>%
  mutate(
    total_data_points = sum(total_points),          # 总点数
    proportion = total_points / total_data_points,  # 每个染色体的占比
    angle_span = proportion * remaining_angle       # 每个染色体的实际角度
  )

chromosome_angles <- chromosome_angles %>%
  arrange(chromosome) %>%  # 确保按照染色体顺序排列
  mutate(
    start_angle = cumsum(lag(angle_span + gap_angle, default = 0)),  # 每个染色体的起始角度
    end_angle = start_angle + angle_span                            # 每个染色体的结束角度
  )

# 1. 合并染色体的角度信息
all_data <- all_data %>%
  left_join(chromosome_angles, by = "chromosome")  # 将 `chromosome_angles` 中的角度信息加入到 `all_data`

# 2. 按染色体分组，计算每个点的角度和半径
all_data <- all_data %>%
  group_by(chromosome) %>%
  mutate(
    # 计算每个数据点的实际角度：
    # - value / total_points: 数据点在当前染色体中的相对位置 (0 到 1)
    # - angle_span: 当前染色体分配的总角度
    # - start_angle: 当前染色体的起始角度
    angle = value / total_points * angle_span + start_angle,
    
    # 根据环类型（Outer 或 Inner）设置不同半径
    radius = 11
  )
type <- c("+/+","-/-","+/-","-/+","+/0","-/0","0/+","0/-","0/0")
length(point_res[[5]])
point_res=lapply(1:8,function(c)as.matrix(filter(all_data,
                                                 category==type[c])[,10]))
dat1 <- data.frame(x=rep(c(radius1,radius2),times=nrow(point_res[[1]])),
                   y=rep(as.numeric(point_res[[1]]),each=2),group=rep(1:nrow(point_res[[1]]),each=2))

dat3 <- data.frame(x=rep(c(radius1,radius2-line_gap),times=nrow(point_res[[3]])),
                   y=rep(as.numeric(point_res[[3]]),each=2),group=rep(1:nrow(point_res[[3]]),each=2))

dat4 <- data.frame(x=rep(c(radius1+line_gap,radius2),times=nrow(point_res[[4]])),
                   y=rep(as.numeric(point_res[[4]]),each=2),group=rep(1:nrow(point_res[[4]]),each=2))

dat5 <- data.frame(x=rep(c(radius2,radius2-line_gap),times=nrow(point_res[[5]])),
                   y=rep(as.numeric(point_res[[5]]),each=2),group=rep(1:nrow(point_res[[5]]),each=2))

dat6 <- data.frame(x=rep(c(radius2,radius2-line_gap),times=nrow(point_res[[6]])),
                   y=rep(as.numeric(point_res[[6]]),each=2),group=rep(1:nrow(point_res[[6]]),each=2))

dat7 <- data.frame(x=rep(c(radius1+line_gap,radius1),times=nrow(point_res[[7]])),
                   y=rep(as.numeric(point_res[[7]]),each=2),group=rep(1:nrow(point_res[[7]]),each=2))
dat8 <- data.frame(x=rep(c(radius1+line_gap,radius1),times=nrow(point_res[[8]])),
                   y=rep(as.numeric(point_res[[8]]),each=2),group=rep(1:nrow(point_res[[8]]),each=2))
siz=0.3
alp=1
gg <- ggplot()+geom_line(
  data = data.frame(dat3),  # 只绘制正常点
  aes(x , y,group=group ),
  col="green",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat1),  # 只绘制正常点
  aes(x , y,group=group ),
  col="red",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat4),  # 只绘制正常点
  aes(x , y,group=group ),
  col="green",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat5),  # 只绘制正常点
  aes(x , y,group=group ),
  col="red",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat6),  # 只绘制正常点
  aes(x , y,group=group ),
  col="blue",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat7),  # 只绘制正常点
  aes(x , y,group=group ),
  col="red",alpha=alp,size = siz
)+geom_line(
  data = data.frame(dat8),  # 只绘制正常点
  aes(x , y,group=group ),
  col="blue",alpha=alp,size = siz
)+
  scale_y_continuous(
    breaks = seq(0, 360, by = 30),  
    limits = c(0, 360)              
  )+  geom_line(
    data = data.frame(x=radius1,y=as.numeric(sapply(1:8,function(c)c(chromosome_angles[c,6:7],NA)))),  # 只绘制正常点
    aes(x , y ),
    size = 4
  )+geom_line(
    data = data.frame(x=radius2,y=as.numeric(sapply(1:8,function(c)c(chromosome_angles[c,6:7],NA)))),  # 只绘制正常点
    aes(x , y ),
    size = 4
  )+
  coord_polar(theta = "y", clip = "off") +
  scale_x_continuous(limits = c(1, 22))   +
  theme_void() + geom_text(data = all_data %>% group_by(chromosome) %>% summarise(mid_angle = mean(angle)),
                           aes(x = radius2-3, y = mid_angle, label = chromosome), inherit.aes = FALSE, size = 40)   # 标注染色体信息

pdf("Figure2A2.pdf",width = 50,height = 50)
gg
dev.off() 





l_siz <-4
p_siz <- 3


gg <- ggplot()+  geom_line(
  data = data.frame(x=3,y=0:10),  # 只绘制正常点
  aes(x , y ),
  size = l_siz
)+geom_line(
  data = data.frame(x=5,y=0:10),  # 只绘制正常点
  aes(x , y ),
  size = l_siz
) +geom_line(
  data = data.frame(x=c(3,5),y=9),  # 只绘制正常点
  aes(x , y ),
  col="red",size = p_siz
)+geom_line(
  data = data.frame(x=c(3,5),y=8),  # 只绘制正常点
  aes(x , y ),
  col="blue",size = p_siz
)+geom_line(
  data = data.frame(x=c(3,6),y=7),  # 只绘制正常点
  aes(x , y ),
  col="green",size = p_siz
)+geom_line(
  data = data.frame(x=c(2,5),y=6),  # 只绘制正常点
  aes(x , y ),
  col="green",size = p_siz
)+geom_line(
  data = data.frame(x=c(5,6),y=5),  # 只绘制正常点
  aes(x , y ),
  col="red",size = p_siz
)+geom_line(
  data = data.frame(x=c(5,6),y=4),  # 只绘制正常点
  aes(x , y ),
  col="blue",size = p_siz
)+geom_line(
  data = data.frame(x=c(2,3),y=3),  # 只绘制正常点
  aes(x , y ),
  col="red",size = p_siz
)+geom_line(
  data = data.frame(x=c(2,3),y=2),  # 只绘制正常点
  aes(x , y ),
  col="blue",size = p_siz
) +geom_point(
  data = data.frame(x=rep(c(3,5),each=9),y=c(1:9,1:9)),  # 只绘制正常点
  aes(x , y ),
  col="black",size = 6
)+
  theme_void() 
pdf("Figure2A3.pdf",width = 4,height = 10)
gg
dev.off() 