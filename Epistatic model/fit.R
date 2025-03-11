rm(list = ls())
library(idopNetwork)
cs <- read.csv("D:/一元抗寒/clean_data/cs.csv")
ch <- read.csv("D:/一元抗寒/clean_data/ch.csv")
rownames(cs) <- cs[,1]
rownames(ch) <- ch[,1]
cs <- cs[,-1]
ch <- ch[,-1]
cs_fit <- power_equation_fit(cs,n = 30,trans = NULL,thread = 12)
saveRDS(cs_fit,"D:/一元抗寒/cs_fit_nolog.rds")
ch_fit <- power_equation_fit(ch,n = 30,trans = NULL,thread = 12)
saveRDS(ch_fit,"D:/一元抗寒/ch_fit_nolog.rds")










