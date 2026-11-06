setwd("C:/Users/2016391/OneDrive - UAB/Documentos/GitHub/paper_missings/Review 2/code")
       
library(data.table)
library(magrittr)

dt <- read.csv("antidepressant_data2.csv")

source("estimators2.r")

setDT(dt)

mnar_estimator (dt[,],col_t=c("t"),cols_w = c("x","g"), apply_bounds=FALSE, 
                estimator=NULL, n_sample=100, delta=TRUE)
