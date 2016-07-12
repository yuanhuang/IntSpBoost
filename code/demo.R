#=====================================================================================================
# program: demo.R                                                                                     
# note: code submission for paper "Promoting similarity of model sparsity structures                  
#                                  in integrative analysis of cancer genetic data".                   
#
# purpose: illustrate how to implement the proposed method with provided datasets
#          covariates datasets x1,x2,x3 
#          logobstime datasets y1,y2,y3
#          censor indicator status1 status2 status3 (1 if events, 0 if censor)
# steps: 
#    (1) source R files to load pre-defined functions
#        source("functions.R")  
#    (2) import datasets.  Demo datasets are provided for illurstration with dim=1000,n1=n2=n3=100.
#    (3) pre-process datasets:
#        get weighted data for AFT model (d.w)
#    (4) specify tunings.
#    (5) implement the algorithms.
#=====================================================================================================

rm(list=ls()) 

# load functions that need to be called by the proposed method in main.R
library(Matrix)   # install.packages("Matrix")
source("function.R")

# import datasets.
x1_read = as.matrix(read.table("../demoData/x1.txt",header = TRUE,sep=","))
x2_read = as.matrix(read.table("../demoData/x2.txt",header = TRUE,sep=","))
x3_read = as.matrix(read.table("../demoData/x3.txt",header = TRUE,sep=","))

y1_read = unlist(read.table("../demoData/y1.txt",header = TRUE,sep=","))
y2_read = unlist(read.table("../demoData/y2.txt",header = TRUE,sep=","))
y3_read = unlist(read.table("../demoData/y3.txt",header = TRUE,sep=","))

status1_read = unlist(read.table("../demoData/status1.txt",header = TRUE,sep=","))
status2_read = unlist(read.table("../demoData/status2.txt",header = TRUE,sep=","))
status3_read = unlist(read.table("../demoData/status3.txt",header = TRUE,sep=","))

# process datasets.

d1.w         = AFT.getWeightedData(x=x1_read,obstime=y1_read,delta.ifobserve=status1_read)            
d2.w         = AFT.getWeightedData(x=x2_read,obstime=y2_read,delta.ifobserve=status2_read)            
d3.w         = AFT.getWeightedData(x=x3_read,obstime=y3_read,delta.ifobserve=status3_read)

xList.use    = list(d1.w[["x.use"]],d2.w[["x.use"]],d3.w[["x.use"]])
yList.use    = list(d1.w[["y.use"]],d2.w[["y.use"]],d3.w[["y.use"]])

# specify the tunings.

lambda_percent.use = 0.1                               

# implement the integrative sparse boosting.

fit = int.spBoost(xList=xList.use,
                  yList=yList.use,
                  lambda_percent=lambda_percent.use)
  

############################################# END of file  #################################################################

