requiredPackages = c("mvtnorm","reshape2","readxl","pbapply",'parallel','orthopolynom','glmnet','ggplot2',"Rsolnp","igraph","devtools","stringr",'patchwork',"deSolve")
for(packages in requiredPackages){
  
  require(packages,character.only = TRUE)
}
rm(list = ls())

load(file = "rdiploid_res.RData")
source("power fit.R")
source("ode_solving1.R")
source("network_reconstruction.R")

res_err <-  c(which(sapply(rdiploid_chra, function(x) inherits(x, "try-error"))),
              which(sapply(rdiploid_chrb, function(x) inherits(x, "try-error"))))

rdiploid_chra <- rdiploid_chra[-res_err]
rdiploid_chrb <- rdiploid_chrb[-res_err]






which(as.character(sapply(1:length(rdiploid_chra),function(c)rdiploid_chra[[c]][[4]][1]))=="+/+")

which(as.character(sapply(1:length(rdiploid_chra),function(c)rdiploid_chra[[c]][[4]][1]))=="+/+")
t(sapply(1:length(rdiploid_chra),function(c)rdiploid_chra[[c]][[4]]))

which(sapply(1:length(rdiploid_chra),function(c){
  if(all(rdiploid_chra[[c]][[4]]==c("+/-","+/-","+/-"))){
  return(1)}else{  return(0)}}
  )==1)
which(sapply(1:length(rdiploid_chra),function(c){
  if(all(rdiploid_chra[[c]][[4]]==c("+/+","+/+","+/-"))){
    return(1)}else{  return(0)}}
)==1)


         #     531, 5,21,25
library(patchwork)

pdf("curve1.pdf",width = 11,height = 7)
qdODE_plot_base(rdiploid_chra[[531]][[3]],label = 2,ylabel = seq(-2,4,2),ylim=c(-2.1,4.5),xlabel = seq(1,4,1),xlim = c(1,4))+qdODE_plot_base(rdiploid_chra[[531]][[1]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[531]][[2]],label = NULL)+plot_layout(ncol = 1)

qdODE_plot_base(rdiploid_chra[[531]][[3]],label = NULL)+qdODE_plot_base(rdiploid_chra[[531]][[1]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[531]][[2]],label = NULL)+plot_layout(ncol = 1)
pdf("06aG01695.pdf",width = 10,height = 11)
qdODE_plot_base(rdiploid_chra[[10418]][[3]],label = 2,ylabel = seq(-2,3,2),ylim=c(-2.2,3.6),xlabel = seq(2,3.5,0.5),xlim = c(1.6,3.8),col="#E8F4C8")+
  qdODE_plot_base(rdiploid_chra[[10418]][[1]],label = 2,ylabel = seq(-2,3,2),ylim=c(-2.2,3.6),xlabel = seq(2,3.5,0.5),xlim = c(1.6,3.8),col="#E8F4C8")+
  qdODE_plot_base(rdiploid_chra[[10418]][[2]],label = 2,ylabel = seq(-2,3,2),ylim=c(-2.2,3.6),xlabel = seq(2,3.5,0.5),xlim = c(1.6,3.8),col="#E8F4C8",xl = "b")+

qdODE_plot_base(rdiploid_chrb[[10418]][[3]],label = 2,ylabel = seq(-0.5,1.5,0.5),ylim=c(-0.8,1.8),xlabel = seq(1.4,2.4,0.2),xlim = c(1.33,2.37))+
  qdODE_plot_base(rdiploid_chrb[[10418]][[1]],label = 2,ylabel = seq(-0.5,1.5,0.5),ylim=c(-0.8,1.8),xlabel = seq(1.4,2.4,0.2),xlim = c(1.33,2.37))+
  qdODE_plot_base(rdiploid_chrb[[10418]][[2]],label = 2,ylabel = seq(-0.5,1.5,0.5),ylim=c(-0.8,1.8),xlabel = seq(1.4,2.4,0.2),xlim = c(1.33,2.37),xl = "b")+plot_layout(ncol = 2,byrow = F)

dev.off() 


qdODE_plot_base(rdiploid_chra[[5]][[3]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[5]][[1]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[5]][[2]],label = NULL)+plot_layout(ncol = 1)
  
qdODE_plot_base(rdiploid_chrb[[5]][[3]],label = NULL)+
  qdODE_plot_base(rdiploid_chrb[[5]][[1]],label = NULL)+
  qdODE_plot_base(rdiploid_chrb[[5]][[2]],label = NULL)+plot_layout(ncol = 1)


qdODE_plot_base(rdiploid_chra[[316]][[3]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[316]][[1]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[316]][[2]],label = NULL)+plot_layout(ncol = 1)


pdf("01aG00010.pdf",width = 10,height = 11)
qdODE_plot_base(rdiploid_chra[[5]][[3]],label = 2,ylabel = seq(0,2,1),ylim=c(-0.4,2.3),xlabel = seq(1.8,2.8,0.2),xlim = c(1.75,2.75),col="#E8F4C8")+
  qdODE_plot_base(rdiploid_chra[[5]][[1]],label = 2,ylabel = seq(0,2,1),ylim=c(-0.4,2.3),xlabel = seq(1.8,2.8,0.2),xlim = c(1.75,2.75),col="#E8F4C8")+
  qdODE_plot_base(rdiploid_chra[[5]][[2]],label = 2,ylabel = seq(0,2,1),ylim=c(-0.4,2.3),xlabel = seq(1.8,2.8,0.2),xlim = c(1.75,2.75),col="#E8F4C8",xl = "b")+
  
  qdODE_plot_base(rdiploid_chrb[[5]][[3]],label = 2,ylabel = seq(-0.4,1.2,0.4),ylim=c(-0.3,1.3),xlabel = seq(1.2,1.8,0.2),xlim = c(1.15,1.85))+
  qdODE_plot_base(rdiploid_chrb[[5]][[1]],label = 2,ylabel = seq(-0.4,1.2,0.4),ylim=c(-0.3,1.3),xlabel = seq(1.2,1.8,0.2),xlim = c(1.15,1.85))+
  qdODE_plot_base(rdiploid_chrb[[5]][[2]],label = 2,ylabel = seq(-0.4,1.2,0.4),ylim=c(-0.3,1.3),xlabel = seq(1.2,1.8,0.2),xlim = c(1.15,1.85),xl = "b")+plot_layout(ncol = 2,byrow = F)

dev.off() 



qdODE_plot_base(rdiploid_chra[[316]][[3]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[316]][[1]],label = NULL)+
  qdODE_plot_base(rdiploid_chra[[316]][[2]],label = NULL)+plot_layout(ncol = 1)


qdODE_plot_base(rdiploid_chrb[[316]][[3]],label = NULL)+
  qdODE_plot_base(rdiploid_chrb[[316]][[1]],label = NULL)+
  qdODE_plot_base(rdiploid_chrb[[316]][[2]],label = NULL)+plot_layout(ncol = 1)


pdf("01aG00590.pdf",width = 10,height = 11)
qdODE_plot_base(rdiploid_chra[[316]][[3]],label = 2,ylabel = seq(0,3,1),ylim=c(-0.3,3.8),xlabel = seq(3.4,4.2,0.2),xlim = c(3.35,4.25),col="#E8F4C8")+
  qdODE_plot_base(rdiploid_chra[[316]][[1]],label = 2,ylabel = seq(0,3,1),ylim=c(-0.3,3.8),xlabel = seq(3.4,4.2,0.2),xlim = c(3.35,4.25),col="#E8F4C8")+
  qdODE_plot_base(rdiploid_chra[[316]][[2]],label = 2,ylabel = seq(0,3,1),ylim=c(-0.3,3.8),xlabel = seq(3.4,4.2,0.2),xlim = c(3.35,4.25),col="#E8F4C8",xl = "b")+
  
  qdODE_plot_base(rdiploid_chrb[[316]][[3]],label = 2,ylabel = seq(-1,3,1),ylim=c(-1.3,3.9),xlabel = seq(3.4,4,0.3),xlim = c(3.15,4.25))+
  qdODE_plot_base(rdiploid_chrb[[316]][[1]],label = 2,ylabel = seq(-1,3,1),ylim=c(-1.3,3.9),xlabel = seq(3.4,4,0.3),xlim = c(3.15,4.25))+
  qdODE_plot_base(rdiploid_chrb[[316]][[2]],label = 2,ylabel = seq(-1,3,1),ylim=c(-1.3,3.9),xlabel = seq(3.4,4,0.3),xlim = c(3.15,4.25),xl = "b")+plot_layout(ncol = 2,byrow = F)

dev.off() 

