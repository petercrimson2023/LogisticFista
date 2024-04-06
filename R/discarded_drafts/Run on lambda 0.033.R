# Loading the Data
load("x_matrix.RData")
load("y.RData")

#Testing on the small dataset with PCA dimension reduction
#x_matrix = as.matrix(read.csv("adata_pca.csv")[,1:50])

# including functions and packages
library(Rcpp)
library(caret)
library(ggplot2)
library(patchwork)
library(nnet)
sourceCpp("fista.cpp")
source("evaluation.R")



n = dim(x_matrix)[1]
p = dim(x_matrix)[2]
k = length(table(y))


beta0 = matrix(0.1,k*p,1)
lambda = 0.033
L_init = 1


y_one_hot=nnet::class.ind(y)

class_list = c("acinar","alpha","beta","delta","ductal","gamma","Class Combined")

result = fista(lambda,L_init,beta0,x_matrix,y_one_hot,n=n,p=p,k=k)

plot_function(result,dir_name)

y_predict = softmax_predict(result$theta,x_matrix,y,p,k)

accuracy_function(result,x_matrix,y,p,k)

plot_function(result)


library(mclust)
library(aricode)

# NMI

nmi_result = NMI(y,y_predict)
nmi_result

# ANI

ani_result = adjustedRandIndex(y,y_predict)
ani_result

