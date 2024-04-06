adata_pca = read.csv("adata_pca.csv",header=T)
source("fista_self_try.R")


y_raw = adata_pca[,"celltype"]
y_num = as.numeric(factor(y_raw))


p = ncol(adata_pca)-2
X_raw = adata_pca[,1:p]
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
#lambda = 0.551096
lambda = 0.1
L_init = 1

result = fista(f,g,X,y,k,p,n,p_y,Q,lambda,beta0,L_init,max_iter=10000,eita=1.5,eps=1e-6,loss_compute=TRUE)

print(matrix(result$theta,nrow=p,ncol=k))

plot_loss(result$loss)

y_predict = softmax_predict(result$theta,X,y,p,k)

print(sum(y_predict==y)/n)









