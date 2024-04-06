
library(Rcpp)
library(RcppArmadillo)

sourceCpp("test.cpp")


adata_nn = read.csv("adata_nn.csv",header=T)

y_raw = adata_nn[,"celltype"]
y_num = as.numeric(factor(y_raw))

p = ncol(adata_nn)-2
X_raw = adata_nn[,1:p]
X = matrix(unlist(X_raw),nrow=nrow(X_raw),ncol=p)
X = cbind(rep(1,nrow(X_raw)),X)
p = p+1
y = ifelse(y_num %in% c(2,7, 8, 10, 11, 12, 13),15,y_num)

### recoding y into 0-6

for (i in 1:length(y))
{
  if(y[i] == 1){
    y[i] = 0
  }else if(y[i] %in% c(3, 4, 5, 6)){
    y[i] = y[i]-2
  }else if(y[i] == 9){
    y[i] = 5
  }else
  {
    y[i] = 6
  }
}

k = length(table(y))
n = nrow(X)

beta0 = matrix(0.1,k*p,1)
lambda = 0.045
L_init = 1

y_one_hot=nnet::class.ind(y)


result = fistaCpp(lambda,L_init,beta0,X,y_one_hot,n=n,p=p,k=k)

print(result$theta)

print(paste0("iteration times: ",result$iter_times))

source("fista_self_try.R")

plot_loss(result$loss)

y_predict = softmax_predict(result$theta,X,y,p,k)

print(paste("general accuracy:",sum(y_predict==y)/length(y)))

result_matrix = as.matrix(table(y_predict,y))

print(diag(result_matrix)/colSums(result_matrix))

