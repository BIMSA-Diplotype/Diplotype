requiredPackages = c("mvtnorm","reshape2","readxl","pbapply",'parallel','orthopolynom','glmnet','ggplot2',"Rsolnp","igraph","devtools","stringr",'patchwork',"deSolve")
for(packages in requiredPackages){
  
  require(packages,character.only = TRUE)
}
rm(list = ls())
source("power fit.R")
source("ode_solving1.R")
load(file = "data.RData")

datfs <- list(chra= datfs[seq(1,nrow(datfs),2),] ,chrb= datfs[seq(2,nrow(datfs),2),])
datss<- list(chra= datss[seq(1,nrow(datss),2),] ,chrb= datss[seq(2,nrow(datss),2),])
datleaf<- list(chra= datleaf[seq(1,nrow(datleaf),2),] ,chrb= datleaf[seq(2,nrow(datleaf),2),])

all(row.names(datfs)==row.names(datleaf))

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


get_diploid <- function(i,data1,data2,data3){
  # i=2870
  # data1=datfs$chra
  # data2=datss$chra
  # data3=datleaf$chra
  nam <-  paste0(c("fs_","ss_","leaf_"),substr(row.names(data1)[i],9,100))
 
 
  dat1 <- rbind(data1[i,],data2[i,],data3[i,])
  dat2 <- rbind(data2[i,],data1[i,],data3[i,])
  dat3 <- rbind(data3[i,],data1[i,],data2[i,])
  rownames(dat1) <- nam
  rownames(dat2) <- nam[c(2,1,3)]
  rownames(dat3) <- nam[c(3,1,2)]
  
  res1 <- power_equation_fit(data=dat1,trans = log)
  res2 <- power_equation_fit(data=dat2,trans = log)
  res3 <- power_equation_fit(data=dat3,trans = log)
  
  result1 <- qdODE_all(result = res1)
  result2 <- qdODE_all(result = res2)
  result3 <- qdODE_all(result = res3)
  


eall <- sapply(list(result1,result2,result3),function(result1){
  x=result1$predict$x[2]-result1$predict$x[1]
  y_all <- sum(result1$predict[,3]*x)
  y_dep1 <- sum(result1$predict[,5]*x)
  y_dep2 <- sum(result1$predict[,6]*x)
  return(c(y_all,y_dep1,y_dep2))
})
  colnames(eall) <- paste0(c("fs_","ss_","leaf_"),substr(row.names(data1)[i],9,100))
  etype <- sapply(1:3,function(c){
    if(abs(eall[2,c])<eall[1,c]*0.1){
      etype1 <- "0"
    }else if ((abs(eall[2,c]) > eall[1,c]*0.1)&(eall[2,c] > 0)){
      etype1 <- "+"
    }else if((abs(eall[2,c]) > eall[1,c]*0.1)&(eall[2,c] < 0)){
      etype1 <- "-"
    }
    if(abs(eall[3,c])<eall[1,c]*0.1){
      etype2 <- "0"
    }else if ((abs(eall[3,c]) > eall[1,c]*0.1)&(eall[3,c] > 0)){
      etype2 <- "+"
    }else if((abs(eall[3,c]) > eall[1,c]*0.1)&(eall[3,c] < 0)){
      etype2 <- "-"
    }
    return(c(etype1,etype2))
    
  })
  
  
 type=c(fs_ss=paste(etype[1,1],etype[1,2],sep = "/"),fs_leaf=paste(etype[2,1],etype[1,3],sep = "/"),ss_leaf=paste(etype[2,2],etype[2,3],sep = "/"))
  return.obj <- list(resultfs=result1,resultss=result2,resultleaf=result3,type=type)
  
  
  return(return.obj)
  
}

get_parllel_diploid <- function(data1,data2,data3,thread=12){
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix","get_diploid","data1","data2","data3",
                      "get_legendre_par","legendre_fit","power_equation_fit","power_equation","power_equation_base","power_equation_all"), envir=environment())
  result <- parLapply(cl, 1:nrow(data1), function(c) {
    try({
      get_diploid(i = c, data1 = data1, data2 = data2 ,data3 = data3)
    }, silent = TRUE)  # 捕获错误，继续执行
  })
  stopCluster(cl)
  return(result)
  
}

rdiploid_chra <- get_parllel_diploid(  data1=datfs$chra,
data2=datss$chra, data3=datleaf$chra,thread = 20)
rdiploid_chrb <- get_parllel_diploid(  data1=datfs$chrb,
                                       data2=datss$chrb, data3=datleaf$chrb,thread = 20)

save(rdiploid_chra,rdiploid_chrb,file = "rdiploid_res.RData")
