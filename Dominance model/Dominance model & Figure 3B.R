requiredPackages = c("mvtnorm","reshape2","readxl","pbapply",'parallel','orthopolynom','glmnet','ggplot2',"Rsolnp","igraph","devtools","stringr",'patchwork',"deSolve")
for(packages in requiredPackages){
  
  require(packages,character.only = TRUE)
}

source("power fit.R")
source("ode_solving1.R")
power_equation_fit <- function(data, n=30, trans = log10) {
  data = data[,order(colSums(data))]
  if ( is.null(trans)) {
    X = colSums(data)
    trans_data = data
  } else{
    X = trans(colSums(data+1))
    trans_data = trans(data+1)
  }
  colnames(trans_data) = X
  
  
  all_model = lapply( 1:nrow(data), function(c) power_equation_all(X, trans_data[c,]))
  
  
  
  names(all_model) = rownames(data)
  no = which(sapply(all_model, length)>=1)
  all_model2 = all_model[no]
  data2 = data[no,]
  trans_data2 = trans_data[no,]
  
  new_x = seq(min(X), max(X), length = n)
  power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
  power_fit = t(vapply(all_model2, predict, newdata = data.frame(x=new_x),
                       FUN.VALUE = numeric(n), USE.NAMES = TRUE))
  
  colnames(power_fit) = new_x
  result = list(original_data = data2, trans_data = trans_data2,
                power_par = power_par, power_fit = power_fit,
                Time = X)
  return(result)
}


get_diploid <- function(i,data1,data2){
  # i=13633
  # data1=gene_dat2a
  # data2=gene_dat2b
  dat1 <- rbind(data1[i,],data2[i,])
  dat2 <- rbind(data2[i,],data1[i,])
  
  res1 <- power_equation_fit(data=dat1,trans = log)
  res2 <- power_equation_fit(data=dat2,trans = log)
  result1 <- qdODE_all(result = res1)
  result2 <- qdODE_all(result = res2)
  edep1 <- sum( (result1$predict$x[2]-result1$predict$x[1])*result1$predict[,5])
  edep2 <- sum( (result2$predict$x[2]-result2$predict$x[1])*result2$predict[,5])
  
  
  eall1 <- sum( (result1$predict$x[2]-result1$predict$x[1])*result1$predict[,3])
  eall2 <- sum( (result2$predict$x[2]-result2$predict$x[1])*result2$predict[,3])
  
  if(abs(edep1) <= eall1*0.1){
    etype1 <- "0"
  }else if((abs(edep1) > eall1*0.1)&(edep1 > 0)){
    etype1 <- "+"
  }else if((abs(edep1) > eall1*0.1)&(edep1 < 0)){
    etype1 <- "-"
  }
  
  if(abs(edep2) <= eall2*0.1){
    etype2 <- "0"
  }else if((abs(edep2) > eall2*0.1)&(edep2 > 0)){
    etype2 <- "+"
  }else if((abs(edep2) > eall2*0.1)&(edep2 < 0)){
    etype2 <- "-"
  }
  return.obj <- list(result1,result2,list(paste(etype1,etype2,sep = "/")))
  
  
  return(return.obj)
  
}

get_parllel_diploid <- function(data1,data2,thread){
  
  core.number <- 12
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix","get_diploid","data1","data2",
                      "get_legendre_par","legendre_fit","power_equation_fit","power_equation","power_equation_base","power_equation_all"), envir=environment())
  result <- parLapply(cl, 1:nrow(data1), function(c) {
    try({
      get_diploid(i = c, data1 = data1, data2 = data2)
    }, silent = TRUE)  # 捕获错误，继续执行
  })
  stopCluster(cl)
  return(result)
  
}

rdiploid_res1 <- get_parllel_diploid(data1=gene_dat1a,data2=gene_dat1b,thread = 24)


rdiploid_res2 <- get_parllel_diploid(data1=gene_dat2a,data2=gene_dat2b,thread = 12)
save(rdiploid_res1,file = "rdiploid_res1.RData")
save(rdiploid_res2,file = "rdiploid_res2.RData")
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
table(as.character(sapply(1:(nrow(gene_dat1a)-1),function(c)rdiploid_res1[[c]][[3]])))
table(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]])))



which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="+/+")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="+/+")

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="+/-")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="+/-")

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="+/0")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="+/0")


pdf("A_++.pdf",width = 11,height = 7.4)



qdODE_plot_base(rdiploid_res2[[13633]][[2]],label = 2,col="#FFF9C4",ylabel = seq(0,4,1),ylim=c(0,3.5),xlabel = seq(1,4,1),xlim = c(1,4))+
  qdODE_plot_base(rdiploid_res1[[1483]][[2]],label = 2,ylabel = seq(0,4,1),ylim=c(0,3.5),xlabel = seq(1.5,3.5,1),xlim = c(1.5,3.5))+
  
  qdODE_plot_base(rdiploid_res2[[13633]][[1]],label = 2,col="#FFF9C4",ylabel = seq(0,4,1),ylim=c(0,3.5),xlabel = seq(1,4,1),xlim = c(1,4),xl="b")+
  qdODE_plot_base(rdiploid_res1[[1483]][[1]],label = 2,ylabel = seq(0,4,1),ylim=c(0,3.5),xlabel = seq(1.5,3.5,1),xlim = c(1.5,3.5),xl="b")

dev.off() 



pdf("A_-+.pdf",width = 11,height = 7.4)

qdODE_plot_base(rdiploid_res2[[28]][[2]],label = 2,col="#FFF9C4",ylabel = seq(-0.5,1.5,0.5),ylim=c(-0.7,1.7),xlabel = seq(1,2,0.5),xlim = c(1,2))+
  qdODE_plot_base(rdiploid_res1[[39]][[2]],label = 2,ylabel = c(-4,0,4,8),ylim=c(-4,6.5),xlabel = seq(3.75,4.75,0.5),xlim = c(3.75,4.75))+
  qdODE_plot_base(rdiploid_res2[[28]][[1]],label = 2,col="#FFF9C4",ylabel = seq(-0.5,1.5,0.5),ylim=c(-0.7,1.7),xlabel = seq(1,2,0.5),xlim = c(1,2),xl="b")+
  qdODE_plot_base(rdiploid_res1[[39]][[1]],label = 2,ylabel = c(-4,0,4,8),ylim=c(-4,6.5),xlabel = seq(3.75,4.75,0.5),xlim = c(3.75,4.75),xl="b")
dev.off() 


pdf("A_0+.pdf",width = 11,height = 7.4)
qdODE_plot_base(rdiploid_res2[[29]][[2]],label = 2,col="#FFF9C4",ylabel = seq(0,4,2),ylim=c(-0.2,4.2),xlabel = seq(3.6,4.2,0.2),xlim = c(3.65,4.25))+
  qdODE_plot_base(rdiploid_res1[[45]][[2]],label = 2,ylabel = seq(0,4,2),ylim=c(-0.2,4.4),xlabel = seq(3.2,4.8,0.4),xlim = c(3.2,4.8))+
  
  qdODE_plot_base(rdiploid_res2[[29]][[1]],label = 2,col="#FFF9C4",ylabel = seq(0,4,2),ylim=c(-0.2,4.2),xlabel = seq(3.6,4.2,0.2),xlim = c(3.65,4.25),xl="b")  +
  qdODE_plot_base(rdiploid_res1[[45]][[1]],label = 2,ylabel = seq(0,4,2),ylim=c(-0.2,4.4),xlabel = seq(3.2,4.8,0.4),xlim = c(3.2,4.8),xl="b")
dev.off() 

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="-/+")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="-/+")

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="-/-")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="-/-")

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="-/0")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="-/0")

library(patchwork)

pdf("A_+-.pdf",width = 11,height = 7.4)
qdODE_plot_base(rdiploid_res2[[6]][[2]],label = 2,col="#FFF9C4",ylabel = seq(-1,2,1),ylim=c(-1.2,2.6),xlabel = seq(1.8,2.4,0.2),xlim = c(1.75,2.4))+
  qdODE_plot_base(rdiploid_res1[[13]][[2]],label = 2,ylabel = seq(0,2,1),ylim=c(-0.6,2.6),xlabel = seq(2.2,3,0.4),xlim = c(2.1,3.1))+
  
  qdODE_plot_base(rdiploid_res2[[6]][[1]],label = 2,col="#FFF9C4",ylabel = seq(-1,2,1),ylim=c(-1.2,2.6),xlabel = seq(1.8,2.4,0.2),xlim = c(1.75,2.4),xl="b")  +
  qdODE_plot_base(rdiploid_res1[[13]][[1]],label = 2,ylabel = seq(0,2,1),ylim=c(-0.6,2.6),xlabel = seq(2.2,3,0.4),xlim = c(2.1,3.1),xl="b")


dev.off() 

pdf("A_0-.pdf",width = 11,height = 7.4)


qdODE_plot_base(rdiploid_res2[[955]][[2]],label = 2,col="#FFF9C4",ylabel = seq(-1,3,1),ylim=c(-1.2,3.6),xlabel = seq(3.2,3.5,0.1),xlim = c(3.16,3.54))+
  qdODE_plot_base(rdiploid_res1[[1519]][[2]],label = 2,ylabel = seq(-1,1,1),ylim=c(-1.2,1.67),xlabel = seq(1.8,2,0.1),xlim = c(1.73,2.03))+
  qdODE_plot_base(rdiploid_res2[[955]][[1]],label = 2,col="#FFF9C4",ylabel = seq(-1,3,1),ylim=c(-1.2,3.6),xlabel = seq(3.2,3.5,0.1),xlim = c(3.16,3.54),xl="b")+
  qdODE_plot_base(rdiploid_res1[[1519]][[1]],label = 2,ylabel = seq(-1,2,1),ylim=c(-1.2,1.67),xlabel = seq(1.8,2,0.1),xlim = c(1.73,2.03),xl="b")
dev.off() 

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="0/+")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="0/+")

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="0/-")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="0/-")

which(as.character(sapply(1:nrow(gene_dat1a),function(c)rdiploid_res1[[c]][[3]]))=="0/0")
which(as.character(sapply(1:nrow(gene_dat2a),function(c)rdiploid_res2[[c]][[3]]))=="0/0")

library(patchwork)


pdf("A_+0.pdf",width = 11,height = 7.4)


qdODE_plot_base(rdiploid_res2[[58]][[2]],label = 2,col="#FFF9C4",ylabel = seq(0,3,1),ylim=c(-0.2,3.83),xlabel = seq(3.9,4.5,0.2),xlim = c(3.75,4.5))+
  qdODE_plot_base(rdiploid_res1[[48]][[2]],label = 2,ylabel = seq(0,1.6,0.4),ylim=c(-0.2,1.7),xlabel = seq(1.7,2.1,0.2),xlim = c(1.65,2.15))+
  qdODE_plot_base(rdiploid_res2[[58]][[1]],label = 2,col="#FFF9C4",ylabel = seq(0,3,1),ylim=c(-0.2,3.83),xlabel = seq(3.9,4.5,0.2),xlim = c(3.75,4.5),xl="b")+
  qdODE_plot_base(rdiploid_res1[[48]][[1]],label = 2,ylabel = seq(0,1.6,0.4),ylim=c(-0.2,1.7),xlabel = seq(1.7,2.1,0.2),xlim = c(1.65,2.15),xl="b")
dev.off() 


pdf("A_-0.pdf",width = 11,height = 7.4)


qdODE_plot_base(rdiploid_res2[[2218]][[2]],label = 2,col="#FFF9C4",ylabel = seq(-2,4,2),ylim=c(-2.2,5),xlabel = seq(3.8,4.3,0.25),xlim = c(3.8,4.3))+
  qdODE_plot_base(rdiploid_res1[[977]][[2]],label = 2,ylabel = seq(-1,3,1),ylim=c(-1,3.5),xlabel = seq(3.1,3.5,0.1),xlim = c(3.05,3.45))+
  qdODE_plot_base(rdiploid_res2[[2218]][[1]],label = 2,col="#FFF9C4",ylabel = seq(-2,4,2),ylim=c(-2.2,5),xlabel = seq(3.8,4.3,0.25),xlim = c(3.8,4.3),xl="b")+
  qdODE_plot_base(rdiploid_res1[[977]][[1]],label = 2,ylabel = seq(-1,3,1),ylim=c(-1,3.5),xlabel = seq(3.1,3.5,0.1),xlim = c(3.05,3.45),xl="b")
dev.off() 

pdf("A_00.pdf",width = 11,height = 7.4)


qdODE_plot_base(rdiploid_res2[[2]][[2]],label = 2,col="#FFF9C4",ylabel = seq(0,4,2),ylim=c(-0.3,5),xlabel = seq(4.7,5,0.1),xlim = c(4.68,5))+
  qdODE_plot_base(rdiploid_res1[[1]][[2]],label = 2,ylabel = seq(0,2,1),ylim=c(-0.1,2.6),xlabel = seq(2.9,3.2,0.1),xlim = c(2.86,3.2))+
  qdODE_plot_base(rdiploid_res2[[2]][[1]],label = 2,col="#FFF9C4",ylabel = seq(0,4,2),ylim=c(-0.3,5),xlabel = seq(4.7,5,0.1),xlim = c(4.68,5),xl="b")+
  qdODE_plot_base(rdiploid_res1[[1]][[1]],label = 2,ylabel = seq(0,2,1),ylim=c(-0.1,2.6),xlabel = seq(2.9,3.2,0.1),xlim = c(2.86,3.2),xl="b")
dev.off() 










