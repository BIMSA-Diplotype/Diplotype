rm(list = ls())
library(idopNetwork)
library(MASS)
library(mvtnorm)
library(parallel)
library(glmnet)
library(deSolve)
library(reshape2)
library(ggplot2)
library(MASS)
h_k83<- readRDS("h_k83.rds")
s_k57 <- readRDS("s_k57.rds")
power_equation <- function (x, power_par) {
  t(sapply(1:nrow(power_par), function(c) power_par[c, 1] * 
             x^power_par[c, 2]))
}
power_equation_base <- function (x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_value = min(y[y != 0])
  lmFit <- lm(log(y + runif(1, min = 0, max = min_value)) ~ 
                log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  model <- try(nls(y ~ a * x^b, start = list(a = a, b = b), 
                   control = nls.control(maxiter = 1e7, minFactor = 1e-200)))
  if ("try-error" %in% class(model)) {
    result = NULL
  }
  else {
    result = model
  }
  return(result)
}
power_equation_fit <- function (data, n = 30, trans = log10, thread = 2) {
  data = data[, order(colSums(data))]
  if (is.null(trans)) {
    X = colSums(data)
    trans_data = data
  }
  else {
    X = trans(colSums(data) + 1)
    trans_data = trans(data + 1)
  }
  colnames(trans_data) = X
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl, c("power_equation_all", "power_equation_base", 
                      "trans_data", "X"), envir = environment())
  all_model = parLapply(cl = cl, 1:nrow(data), function(c) power_equation_all(X, 
                                                                              trans_data[c, ]))
  stopCluster(cl)
  names(all_model) = rownames(data)
  no = which(sapply(all_model, length) >= 1)
  all_model2 = all_model[no]
  data2 = data[no, ]
  trans_data2 = trans_data[no, ]
  new_x = seq(min(X), max(X), length = n)
  power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), 
                       USE.NAMES = TRUE))
  power_fit = t(vapply(all_model2, predict, newdata = data.frame(x = new_x), 
                       FUN.VALUE = numeric(n), USE.NAMES = TRUE))
  colnames(power_fit) = new_x
  result = list(original_data = data2, trans_data = trans_data2, 
                power_par = power_par, power_fit = power_fit, Time = X)
  return(result)
}
power_equation_all <- function (x, y, maxit = 100) {
  result <- power_equation_base(x, y)
  iter <- 1
  while (is.null(result) && iter <= maxit) {
    iter <- iter + 1
    try(result <- power_equation_base(x, y))
  }
  return(result)
}
get_interaction <- function (data, col, reduction = FALSE) 
{
  if (nrow(data) == 2) {
    return_obj = list(ind.name = rownames(data)[col], dep.name = rownames(data)[-col], 
                      coefficient = cor(t(data))[1, 2])
  }
  else {
    data <- t(data)
    name <- colnames(data)
    y = as.matrix(data[, col])
    x = as.matrix(data[, -col])
    n <- ncol(x)
    if (reduction == TRUE) {
      vec <- abs(apply(x, 2, cor, y))
      if (all(is.na(vec))) {
        return_obj = list(ind.name = name[col], dep.name = NA, 
                          coefficient = 0)
      }
      else {
        x = x[, order(vec, decreasing = T)[1:(n/log(n))]]
      }
    }
    if (all(y == 0) | all(y == 1)) {
      return_obj = list(ind.name = name[col], dep.name = NA, 
                        coefficient = 0)
    }
    else {
      ridge_cv <- try(cv.glmnet(x = x, y = y, alpha = 0))
      if ("try-error" %in% class(ridge_cv)) {
        return_obj = list(ind.name = name[col], dep.name = NA, 
                          coefficient = 0)
      }
      else {
        ridge_cv <- cv.glmnet(x = x, y = y, type.measure = "mse", 
                              nfolds = 10, alpha = 0)
        best_ridge_coef <- abs(as.numeric(coef(ridge_cv, 
                                               s = ridge_cv$lambda.1se))[-1])
        fit <- cv.glmnet(x = x, y = y, alpha = 1, family = "gaussian", 
                         type.measure = "mse", penalty.factor = 1/best_ridge_coef, 
                         nfolds = 10, keep = TRUE, thresh = 1e-10, maxit = 1e+06)
        lasso_coef <- coef(fit, s = fit$lambda.1se)
        return_obj = list(ind.name = name[col], dep.name = lasso_coef@Dimnames[[1]][lasso_coef@i + 
                                                                                      1][-1], coefficient = lasso_coef@x[-1])
        if (length(return_obj$dep.name) == 0) {
          tmp = cor(x, y)
          return_obj$dep.name = rownames(tmp)[which.max(abs(tmp))]
          return_obj$coefficient = tmp[which.max(abs(tmp))] * 
            1/3
        }
      }
    }
  }
  return(return_obj)
}

qdODEmod <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
  })
}

qdODE_ls <- function(pars, data, Time, power_par){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  X = as.numeric(data[1,])
  fit = as.numeric(out[,2])
  ind = as.numeric(out[,(n+2)])
  sse = sum(crossprod(X-fit),sum((ind[ind<0])^2))
  return(sse)
}

qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 3, new_time = NULL, n_expand = 100){
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  out2 = data.frame(x = out[,1], y = data[1,], y.fit = out[,2],
                    ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  rownames(out2) = NULL
  
  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))
  
  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}


qdODE_all <- function(result, relationship, i, init_pars = 1, LOP_order = 3, method = "ls",
                      new_time = NULL, n_expand = 60, maxit = 1e3){
  Time = as.numeric(colnames(result$power_fit))
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result$power_fit[variable,]
  
  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    power_par = result$power_par[variable,][-1,]
    n = nrow(data)
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    if (method == "ls") {
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "L-BFGS-B",
                         lower = c(rep(-10,(length(pars_int)))),
                         upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    } else{
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                         method = "Nealder-Mead",
                         #lower = c(rep(-10,(length(pars_int)))),
                         #lower = c(0, rep(-10,(length(pars_int))-1)),
                         #upper = c(rep(10,(length(pars_int)))),
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    }
  }
  return(return.obj)
}

qdODE_parallel <- function(result, reduction = TRUE, thread = 3, maxit = 1e3){
  data = result$original_data
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, reduction = reduction))
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit","get_interaction"), envir=environment())
  result = parLapply(1:nrow(data),function(c) qdODE_all(result = result,
                                                        relationship = relationship,
                                                        i = c,
                                                        maxit = maxit
  ), cl = cl)
  stopCluster(cl)
  names(result) = rownames(data)
  names(relationship) = rownames(data)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}
fun_clu_convertbyy <- function (result, best.k) {
  cluster.result = result
  times = cluster.result$Time
  times_new = seq(min(times), max(times), length = 30)
  par.mu = cluster.result$mu_par
  colnames(par.mu) = c("a", "b")
  rownames(par.mu) = paste0("M", 1:best.k)
  k = cluster.result$cluster_number
  mu.fit = power_equation(times_new, par.mu[1:k, ])
  colnames(mu.fit) = times_new
  rownames(mu.fit) = paste0("M", 1:best.k)
  df = cluster.result$cluster
  tmp = split(df, df[, ncol(df)])
  tmp2 = lapply(tmp, function(x) {
    x["apply.omega..1..which.max."] <- NULL
    x
  })
  return_obj <- list(original_data = mu.fit, trans_data = mu.fit, 
                     power_par = par.mu, power_fit = mu.fit, Module.all = tmp2)
  return(return_obj)
}
from_matrix_to_netdata <- function(r1 = r1){
  qdode_solve_res <- qdODE_parallel(r1,reduction = TRUE,thread = 6)
  net_prepare <- function(x){
    df1 <- qdode_solve_res[["ode_result"]][[x]][[2]]
    # 假设您的数据框的名称是 df
    # 获取第五列及其后的所有列
    selected_columns <- as.matrix(df1[, 5:ncol(df1)])
    colnames(selected_columns)  = colnames(df1)[5:ncol(df1)]
    
    # 计算每一列的平均值
    column_means <- colMeans(selected_columns, na.rm = TRUE)
    
    # 获取列名
    column_names <- colnames(selected_columns)
    
    # 创建一个包含列名和平均值的数据框
    result_df1 <- data.frame(Column_Name = column_names, Mean_Value = column_means)
    indmodule <- data.frame(Indmodule = rep(colnames(df1)[4],ncol(selected_columns)))
    result_df1_1 <- as.matrix(cbind(indmodule,result_df1))
    # 打印或查看结果
    print(result_df1_1)
  }
  # 使用 sapply 将 net_prepare 函数应用于 x 从 1 到 66 的所有值
  net_data1 <- sapply(1:nrow(r1[["original_data"]]), net_prepare, simplify = FALSE)
  # 使用 do.call 和 rbind 将结果列表组合成一个大矩阵
  net_data <- do.call(rbind, net_data1)
  # 打印或查看最终的结果矩阵
  print(net_data)
  row.names(net_data) <- 1:nrow(net_data)
  net_data2 <- net_data[,c(2,1,3)]
  aaa <- as.numeric(net_data2[,3])
  zheng_or_fu <- matrix(aaa,nrow(net_data),1)
  zheng_or_fu <- cbind(zheng_or_fu, ifelse(zheng_or_fu[, 1] >= 0, 1, -1))
  zheng_or_fu[,1] <- abs(zheng_or_fu[,1])
  net_data_final <- cbind(net_data2,zheng_or_fu)[,-3]
  colnames_net_data <- c("source","target","weight","cor")
  colnames(net_data_final) <- colnames_net_data
  
  return(list(net_data_final = net_data_final, qdode_solve_res = qdode_solve_res))
}

clu_con_83 <- fun_clu_convertbyy(h_k83,83)
clu_con_83[["original_data"]] <- log10(clu_con_83[["original_data"]])
colnames(clu_con_83[["original_data"]]) <- log10(as.numeric(colnames(clu_con_83[["original_data"]])))
clu_con_83[["power_fit"]] <- clu_con_83[["original_data"]]
ode_re <- from_matrix_to_netdata(clu_con_83)







