get_legendre_matrix <- function(x,legendre_order){
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  return(legendre_matrix[,2:(legendre_order+1)])
}

#' @title use legendre polynomials to fit a given data
#' @importFrom stats lm coef
#' @param legendre_order scalar of legendre polynomials
#' @param x vector equal to the x value for legendre polynomials(in this case times)
#' @param y vector equal to the y observed data(in this case generic effect)
#' @return the polynomials coefficients
#' @examples get_legendre_par(14:1,4,1:14)
#' @export
get_legendre_par <- function(y,legendre_order,x){
  #lm_method
  legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
  return(legendre_par)
}

#' @title generate curve based on legendre polynomials
#' @importFrom orthopolynom legendre.polynomials polynomial.values scaleX
#' @param par vector of legendre polynomials coefficients
#' @param x vector equal to the x value for legendre polynomials(in this case times)
#' @return the polynomials value
#' @examples legendre_fit(rep(1,5),1:14)
#' @export
legendre_fit <- function(par,x){
  legendre_order = length(par)
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  x_interpolation <- rowSums(legendre_fit)
  return(x_interpolation)
}

get_interaction <- function(data, col, reduction = FALSE ){
  
  if (nrow(data)==2) {
    return_obj = list(ind.name = rownames(data)[col],
                      dep.name = rownames(data)[-col],
                      coefficient = cor(t(data))[1,2])
    
  } else{
    data <- t(data); name <- colnames(data)
    y = as.matrix(data[,col])
    x = as.matrix(data[,-col])
    
    n <- ncol(x)
    if (reduction == TRUE) {
      vec <- abs(apply(x, 2, cor, y))
      if (all(is.na(vec))) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)
      } else{
        x = x[,order(vec, decreasing = T)[1:(n/log(n))]]
      }
    }
    
    if ( all(y==0) |  all(y==1) ) {
      return_obj = list(ind.name = name[col],
                        dep.name = NA,
                        coefficient = 0)
    } else{
      ridge_cv <- try(cv.glmnet(x = x, y = y,alpha = 0))
      if ('try-error' %in% class(ridge_cv)) {
        return_obj = list(ind.name = name[col],
                          dep.name = NA,
                          coefficient = 0)
        
      } else{
        ridge_cv <- cv.glmnet(x = x, y = y, type.measure = "mse", nfolds = 10, alpha = 0)
        best_ridge_coef <- abs(as.numeric(coef(ridge_cv, s = ridge_cv$lambda.1se))[-1])
        
        fit <- cv.glmnet(x = x, y = y, alpha = 1, family = "gaussian", type.measure = "mse",
                         penalty.factor = 1/best_ridge_coef,
                         nfolds = 10, keep = TRUE, thresh=1e-10, maxit=1e6)
        lasso_coef <- coef(fit, s = fit$lambda.1se)
        return_obj = list(ind.name = name[col],
                          dep.name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1][-1],
                          coefficient = lasso_coef@x[-1])
        if ( length(return_obj$dep.name)==0 ) {
          tmp = cor(x,y)
          return_obj$dep.name = rownames(tmp)[which.max(abs(tmp))]
          return_obj$coefficient = tmp[which.max(abs(tmp))]*1/3
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

#' @title least-square fit for qdODE model
#' @importFrom deSolve ode
#' @param pars vector for unknown ODE parameters
#' @param data data contain independent effect as first row and dependent effect
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @return mean-square error
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

#' @title legendre polynomials fit to qdODE model
#' @importFrom deSolve ode
#' @param pars vector of qdODE parameters
#' @param data dataframe of observed data
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @param LOP_order scalar of LOP order
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @return list contain legendre polynomials parameters, qdODE values and LOP fitted values
qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 6, new_time = NULL, n_expand = 100){
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



#' @title wrapper for qdODE model
#' @importFrom stats optim
#' @param result result from power_equation_fit
#' @param relationship list contain variable selection results
#' @param i scalar for which id used for qdODE solving, must <= nrow
#' @param init_pars scalar for initial parameters
#' @param LOP_order scalar of LOP order
#' @param method scalar of qdODE solving methodm, cuurent only support least square
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_all <- function(result, init_pars = 1, LOP_order = 6, 
                      new_time = NULL, n_expand = 100, maxit = 1e3){
  Time = as.numeric(colnames(result$power_fit))
  #variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result$power_fit


    power_par = result$power_par[-1,]
    n = nrow(data)
    pars_int = c(init_pars,cor(as.numeric(data[1,]),as.numeric(data[2,])))
    
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

  
  return(return.obj)
}


#' @title wrapper for qdODE_all in parallel version
#' @param result result from power_equation_fit
#' @param reduction use n/log(n) dimension reduction
#' @param thread scales for how many threads used
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_parallel <- function(result, reduction = FALSE, thread = 2, maxit = 1e3){
  data = result$original_data
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, reduction = reduction))
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit"), envir=environment())
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

#' @title convert qdODE results to plot data
#' @importFrom reshape2 melt
#' @param result list of qdODE all
#' @importFrom stats na.omit
qdODEplot_convert <- function(result){
  data = result$predict
  n = ncol(data)
  colnames(data)[4:n] = c(paste0("ind.",colnames(data)[4]),
                          paste0("dep.",colnames(data)[5:n]))

  plot.df = melt(data, id.vars = c("x"))

  name = levels(plot.df[,2])

  ind.name = name[grep("ind", name)]
  ind.name2 = strsplit(ind.name,split = "\\.")[[1]][2]
  ind.df <- subset(plot.df, plot.df[,2] == ind.name)
  ind.df$type = "ind"
  ind.df$variable = ind.name2

  depname = levels(plot.df[,2])[grep("dep",name )]
  dep.df <- subset(plot.df, plot.df[,2] %in% depname)
  dep.df$type = "dep"
  dep.df$variable = sapply(strsplit(as.character(dep.df$variable),"\\."),"[",2)


  original.df = subset(plot.df, plot.df[,2] == "y")
  original.df$type = "original"

  fit.df = subset(plot.df, plot.df[,2] == "y.fit")
  fit.df$type = "fit"

  plot.df2 = rbind(ind.df, dep.df,fit.df)

  name.df = subset(plot.df2, plot.df[,1] == max(plot.df2[,1]))
  name.df = name.df[-nrow(name.df),]
  name.df[,2][name.df[,2] == "y.fit"] = ind.name2

  name.df = name.df[-which(name.df[,4] == "fit"),]

  name.df[,1] = name.df[,1]*1.002
  return_obj = list(plot.df2 = plot.df2,
                    name.df = name.df,
                    ind.name2 = ind.name2)
  return(return_obj)
}


#' @title plot single decompose plot
#' @import ggplot2
#' @param result list of qdODE all
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @export


qdODE_plot_base <- function(result,label = 2,ylabel, show.legend = FALSE,col= "#D6EAF8",xl="a"){
  #result=rdiploid_res1[[3]][[1]]
  result2 = qdODEplot_convert(result)
  plot.df2 = result2$plot.df2
  name.df = result2$name.df
  ind.name2 = result2$ind.name2

  p <- ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "type", colour = "type"), size = 1.1,
              show.legend = show.legend)  +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = NULL,
      labels =NULL,
      values = c("green", "blue", "red")) +
    xlab(NULL) + ylab(NULL) +
    ggtitle(NULL) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y = element_text(hjust = 0))+ theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
          axis.text.y =element_text(#face="bold", 
            color="black", size=25)  ,
          axis.text.x =element_text(#face="bold", 
            color="black", size=25) 
    )+theme(legend.position = "none",panel.background = element_rect(fill = col, color = "black"))

  if (is.null(label)) {
    p = p + scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x)*1.005))
  } else {
    xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks

    xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )
    
    # ylabel=seq(-3,1.2,0.4)
    if (0 %in% ylabel) {
      pos0 = which(ylabel==0)
      text=paste(label,"^", ylabel, sep="")
      text[pos0] = "0"
      if (any(na.omit(ylabel)<0)) {
        pos1 = which(ylabel<0)
        text[pos1] = paste0("-",label,"^",abs(ylabel[pos1]))
        ylabel2 = parse(text=text)
      } else{ylabel2 = parse(text=text)}
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    if(xl == "a"){
      p = p + theme(
                    plot.margin = unit(c(1,1,0,1),"lines"),
                    axis.ticks.length.x = unit(-0.1,"cm"))+scale_x_continuous(breaks=xlabel, labels = NULL, limits = c(min(plot.df2$x), max(plot.df2$x)*1.005)) +
        scale_y_continuous(breaks=ylabel,labels = ylabel2,limits = c(min(ylabel)-0.3,max(ylabel)+max(ylabel)*0.1))+ theme(axis.text.y = element_text(hjust = 0))
      
    }else{
    p = p +theme(
      plot.margin = unit(c(0,1,1,1),"lines"))+ scale_x_continuous(breaks=xlabel, labels = xlabel2, limits = c(min(plot.df2$x), max(plot.df2$x)*1.005)) +
      scale_y_continuous(breaks=ylabel,labels = ylabel2,limits = c(min(ylabel)-0.3,max(ylabel)+max(ylabel)*0.1))+ theme(axis.text.y = element_text(hjust = 0))
    }
    }
  return(p)
}
qdODE_plot_base <- function(result,label = 2,ylabel,ylim,xlabel,xlim, show.legend = FALSE,col= "#D6EAF8",xl="a"){
  #result=rdiploid_res1[[3]][[1]]
  result2 = qdODEplot_convert(result)
  plot.df2 = result2$plot.df2
  name.df = result2$name.df
  ind.name2 = result2$ind.name2
  
  p <- ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "type", colour = "type"), size = 1.1,
              show.legend = show.legend)  +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = NULL,
      labels =NULL,
      values = c("green", "blue", "red")) +
    xlab(NULL) + ylab(NULL) +
    ggtitle(NULL) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y = element_text(hjust = 0))+ theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
          axis.text.y =element_text(#face="bold", 
            color="black", size=25)  ,
          axis.text.x =element_text(#face="bold", 
            color="black", size=25) 
    )+theme(legend.position = "none",panel.background = element_rect(fill = col, color = "black"))
  
  if (is.null(label)) {
    p = p + scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x)*1.005))
  } else {
    #xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
    
    xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )
    
    # ylabel=seq(-3,1.2,0.4)
    if (0 %in% ylabel) {
      pos0 = which(ylabel==0)
      text=paste(label,"^", ylabel, sep="")
      text[pos0] = "0"
      if (any(na.omit(ylabel)<0)) {
        pos1 = which(ylabel<0)
        text[pos1] = paste0("-",label,"^",abs(ylabel[pos1]))
        ylabel2 = parse(text=text)
      } else{ylabel2 = parse(text=text)}
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    if(xl == "a"){
      p = p + theme(
        plot.margin = unit(c(1,1,0,1),"lines"),
        axis.ticks.length.x = unit(-0.1,"cm"))+scale_x_continuous(breaks=xlabel, labels = NULL, limits = c(xlim[1],xlim[2])) +
        scale_y_continuous(breaks=ylabel,labels = ylabel2,limits = c(ylim[1],ylim[2]))+ theme(axis.text.y = element_text(hjust = 0))
      
    }else{
      p = p +theme(
        plot.margin = unit(c(0,1,1,1),"lines"))+ scale_x_continuous(breaks=xlabel, labels = xlabel2, limits = c(xlim[1],xlim[2])) +
        scale_y_continuous(breaks=ylabel,labels = ylabel2,limits = c(ylim[1],ylim[2]))+ theme(axis.text.y = element_text(hjust = 0))
    }
  }
  return(p)
}

#' @title plot all decompose plot
#' @import ggplot2 patchwork
#' @param result list of qdODE parallel
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @param nrow scalar for subplot row number
#' @param ncol scalar for subplot column number
#' @return all effect curve decompose plot
#' @export
qdODE_plot_all <- function(result,label = 10, show.legend = TRUE, nrow = NULL, ncol = NULL){
  p = lapply(result$ode_result, qdODE_plot_base, label = label, show.legend = show.legend)
  #p = lapply(1:length(result$ode_result), function(c) qdODE_plot_base(result$ode_result[[c]], label = label))
  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))

  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Niche Index", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(p, nrow = nrow, ncol = ncol)) +
    plot_annotation(caption = "Habitat Index",
                    theme = theme(plot.caption = element_text(size = 20,hjust=.5))) +
    plot_layout(widths = c(.05, 1),guides = 'collect')

  return(pp)
}

#' @title plot single decompose plot for two data
#' @import ggplot2
#' @param result1 list of qdODE all for first data
#' @param result2 list of qdODE all for second data
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @param remove.label to remove x and y label
#' @export
biqdODE_plot_base <- function(result1, result2, label = 10, show.legend = FALSE, remove.label = FALSE){
  resulta = qdODEplot_convert(result1)
  resultb = qdODEplot_convert(result2)
  plot.df1 = resulta$plot.df2
  name.df1 = resulta$name.df
  ind.name2 = resulta$ind.name2
  name.df1$x = name.df1$x*0.99

  plot.df2 = resultb$plot.df2
  name.df2 = resultb$name.df
  name.df2$x = name.df2$x*0.99

  lower1 = min(plot.df1[,3])
  upper1 = max(plot.df1[,3])
  lower2 = min(plot.df2[,3])
  upper2 = max(plot.df2[,3])
  y_min = round(min(lower1, lower2),1)-0.05
  y_max = round(max(upper1, upper2),1)+0.05

  p1 = ggplot(plot.df1, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    geom_text(name.df1, mapping = aes_string(label = "variable", colour = "type",
                                             x = "x", y = "value"), size = 3,show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) + theme_bw() +
    xlab("Habitat Index") + ylab("Niche Index") +
    theme(axis.text.y = element_text(hjust = 0)) + xlab(NULL)+
    ggtitle(ind.name2) + theme(plot.title = element_text(hjust = 1))


  s1 = p1 + scale_y_continuous(limits = c(y_min, y_max))

  p2 = ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1,
              show.legend = show.legend) +
    geom_text(name.df2, mapping = aes_string(label = "variable", colour = "type",
                                             x = "x", y = "value"), size = 3,show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) + theme_bw() +
    theme(axis.text.y = element_text(hjust = 0)) + xlab(NULL)

  s2 = p2 + scale_y_continuous(limits = c(y_min, y_max))

  if (is.null(label)) {
    p1 = p1 + scale_x_continuous(limits = c(min(plot.df1$x), max(plot.df1$x)*1.005))+
      theme(axis.title.x = element_text(vjust=1))
    p2 = p2 + scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x)*1.005))
  } else {
    xlabel1 = ggplot_build(s1)$layout$panel_params[[1]]$x.sec$breaks
    xlabel2 = ggplot_build(s2)$layout$panel_params[[1]]$x.sec$breaks

    ylabel = ggplot_build(s1)$layout$panel_params[[1]]$y.sec$breaks

    xlabel1.2 = parse(text= paste(label,"^", xlabel1, sep="") )
    xlabel2.2 = parse(text= paste(label,"^", xlabel2, sep="") )
    if (0 %in% ylabel) {
      pos0 = which(ylabel==0)
      text=paste(label,"^", ylabel, sep="")
      text[pos0] = "0"
      if (any(na.omit(ylabel)<0)) {
        pos1 = which(ylabel<0)
        text[pos1] = paste0("-",label,"^",abs(ylabel[pos1]))
        ylabel2 = parse(text=text)
      } else{ylabel2 = parse(text=text)}
    } else{
      ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
    }
    p1 = p1 + scale_x_continuous(labels = xlabel1.2, limits = c(min(xlabel1)-0.05, max(xlabel1)+0.05)) +
      scale_y_continuous(limits = c(y_min,y_max),labels = ylabel2) +
      theme(axis.text.y = element_text(hjust = 0))+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))+
      theme(axis.title.x = element_text(vjust=1))

    p2 = p2 + scale_x_continuous(labels = xlabel2.2, limits = c(min(xlabel2)-0.05, max(xlabel2)+0.05)) +
      scale_y_continuous(limits = c(y_min,y_max),labels = ylabel2) +
      theme(axis.text.y = element_text(hjust = 0)) +
      theme(axis.text.y = element_blank(), axis.ticks.length.y = unit(-0.1,"cm")) +
      ylab(NULL)+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))

  }

  if (remove.label == TRUE) {
    p1 = p1 + theme(axis.title=element_blank())
    p2 = p2 + theme(axis.title=element_blank())
  }

  pp = p1+p2
  return(pp)
}


#' @title plot all decompose plot for two data
#' @import ggplot2
#' @param result1 list of qdODE all for first data
#' @param result2 list of qdODE all for second data
#' @param label relabel x and y label due to log-transform, set 10 as default
#' @param show.legend to show legend
#' @param remove.label to remove x and y label
#' @param nrow scalar for subplot row number
#' @param ncol scalar for subplot column number
#' @export
biqdODE_plot_all <- function(result1, result2, label = 10, show.legend = FALSE,
                             remove.label = TRUE, nrow = NULL, ncol = NULL){
  n = length(result1$ode_result)
  p = lapply(1:n, function(c)
    biqdODE_plot_base(result1$ode_result[[c]], result2$ode_result[[c]],
                      label = label, show.legend = show.legend, remove.label = remove.label))

  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))

  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Niche Index", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(p, nrow = nrow, ncol = ncol)) +
    plot_annotation(caption = "Habitat Index",
                    theme = theme(plot.caption = element_text(size = 20,hjust=.5))) +
    plot_layout(widths = c(.03, 1),guides = 'collect')

  return(pp)
}
