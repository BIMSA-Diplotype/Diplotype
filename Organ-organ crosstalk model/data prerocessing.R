
#######################################不去除0
rm(list = ls())

gene_dat1 <- read.csv("1.MRMfruit.FPKM.csv",row.names = 1)#cold hardiness

gene_dat2 <- read.csv("2.MRMleaf.FPKM.csv",row.names = 1)#cold stress

library(zoo)



gene_datfs <- gene_dat1[,1:24]
  gene_datss <- gene_dat1[,c(25:33)]

  gene_datleaf <-  gene_dat2

  
  n0fs <- sapply(1:nrow(gene_datfs), function(c)length(which(gene_datfs[c,]==0)))
  n0ss <- sapply(1:nrow(gene_datss), function(c)length(which(gene_datss[c,]==0)))
  n0leaf <- sapply(1:nrow(gene_datleaf), function(c)length(which(gene_datleaf[c,]==0)))
  n0 <- which((n0fs < (ncol(gene_datfs)/3))&(n0ss < (ncol(gene_datss)/3))&
          (n0leaf < (ncol(gene_datleaf)/3)))
  gene_datfs <-gene_datfs[n0,]
  gene_datss <-gene_datss[n0,]
  gene_datleaf <-gene_datleaf[n0,]

  
  gene_match <- read.csv("match.csv",header = F)#cold stress
  gene_match <- gene_match[,1:2]
  gene_match <- gene_match[intersect(   which( !is.na(match(c(gene_match[,1]),row.names(gene_datfs)))), 
               which( !is.na(match(c(gene_match[,2]),row.names(gene_datfs))))),]
  
  
  amatch <-  as.numeric( na.omit(match(c(gene_match[,1]),row.names(gene_datfs))))
  bmatch <-  as.numeric( na.omit(match(c(gene_match[,2]),row.names(gene_datfs))))



  abmatch  <- as.vector(rbind(amatch, bmatch))
  
  length(bmatch)
  

  
  datfs1 <- gene_datfs[ abmatch ,-c(1:3)]
  datss1 <- gene_datss[abmatch,]
  datleaf1 <- gene_datleaf[abmatch,]
  
 
  
  
  datss <- matrix(data=NA,nrow = nrow(datss1),ncol = 24)
  datfs <- matrix(data=NA,nrow = nrow(datss1),ncol = 24)
  datleaf <- matrix(data=NA,nrow = nrow(datss1),ncol = 24)
  row.names(datss)=  row.names(datfs)=  row.names(datleaf)=  row.names(datss1)
  colnames(datss)=  colnames(datfs)=  colnames(datleaf)= sort(rep(1:8,3))
  
  
  datss[,c(4:6,13:15,19:21)] <- as.matrix( datss1)
 
  datfs[,c(4:24)] <- as.matrix( datfs1)
  
  datleaf[,c(1:3,7:9,22:24)] <- as.matrix( datleaf1[,-(7:12)])
  
  get_data <- function(ss,fs,leaf){
   # ss<- datss[15894,]
   # fs <- datfs[15894,]
   # leaf <- datleaf[15894,]
    ss <- matrix( ss,ncol=3,byrow = TRUE)
    fs <-  matrix( fs,ncol=3,byrow = TRUE)
    leaf <- matrix( leaf,ncol=3,byrow = TRUE)
    # leaf <- sapply(1:3,function(c){
    #   y_log <- log(as.numeric(na.omit(leaf[,c])) + 1)
    # exp(spline(c(1,3,8), y_log, xout =1:8)$y)-1})
    leaf <-  sapply(1:3,function(c){
      datt <- data.frame(x=1:8,y=as.numeric((leaf[1:8,c])))
      na.approx(datt$y,x=datt$x, na.rm = FALSE)})

    # ss[2:7,] <- sapply(1:3,function(c){
    #   y_log <- log(as.numeric(na.omit(ss[,c])) + 1)
    #   exp(spline(c(2,5,7), y_log, xout =2:7)$y)-1})
    
    ss[2:7,] <-  sapply(1:3,function(c){
   datt <- data.frame(x=2:7,y=as.numeric((ss[2:7,c])))
      na.approx(datt$y,x=datt$x, na.rm = FALSE)})

    
    
    dat <- cbind(fs,ss,leaf)
    
    x=log(dat[2,]+1)
    y=log(dat[3,]+1)
    resy <- lm(x~y)
    
    
    
    # plot(x, y, pch = 16, col = "blue", main = "散点图及回归拟合线", xlab = "X", ylab = "Y")
    # 
    # # 添加拟合直线
    # abline(resx, col = "red", lwd = 2)
    x=log(dat[6,]+1)
    y=log(dat[7,]+1)
    resx <- lm(y~x)

    dat[1,1:6] <-   exp(predict(resy,newdata = data.frame(y= log(dat[2,1:6]+1))))-1
    
    
  dat[8,4:6] <-   exp(predict(resx,newdata = data.frame(x= log(dat[7,4:6]+1))))-1
  
  
  # 
  # dat[4,7:9] <-   predict(resx,newdata = data.frame(x= dat[3,7:9]))
  # dat[5,7:9] <-   predict(resx,newdata = data.frame(x= dat[4,7:9]))
  # dat[6,7:9] <-   predict(resx,newdata = data.frame(x= dat[5,7:9]))
  # dat[7,7:9] <-   predict(resx,newdata = data.frame(x= dat[6,7:9]))
  # 
  # 
  # 
  # dat[8,4:6] <-   predict(resx,newdata = data.frame(x= dat[7,4:6]))
  
  fs <- as.numeric( t(dat[,1:3]))
  ss <- as.numeric( t(dat[,4:6]))
  leaf <- as.numeric( t(dat[,7:9]))
  return(c(fs,ss,leaf))
  }
  
  dat <- t(sapply(1:nrow(datss),function(c){
    res=get_data(ss=datss[c,],fs=datfs[c,],leaf = datleaf[c,])
  }))
  dat[which(dat< 0)]=0
  datfs[,1:24] <- dat[,1:24]
  datss[,1:24] <- dat[,25:48]
  datleaf[,1:24] <- dat[,49:72]


    save(datss, datfs, datleaf,file = "data.RData")
  

    
    
    

    
    
    
  ##################################################################
  
  
  
  
  data <-  datleaf[,1:24]
  data <-as.data.frame(data)
  
  data_matrix <- log(as.matrix(data) +1 ) ###取log用这个   
  #data_matrix <- as.matrix(data) ###不取log用这个                        
  
  par(mfrow = c(4, 6) )
  
  
  for (i in 1:ncol(data_matrix)) {
    
    hist(data_matrix[, i], main = paste("Histogram of Column", i), 
         xlab = paste("Value", i), ylab = "Frequency",   
         col = "blue", border = "white", probability = TRUE)
    
    
    lines(density(data_matrix[, i]), col =rgb(1, 0, 0, 0.7), lwd = 2)  
  }
  
  
  
  
 
  
  