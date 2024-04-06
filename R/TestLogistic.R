set.seed(12345)
n <- 30
p <- 2
k = 2
X <- matrix(rnorm(n*p),n,p)
y <- sample(c(0,1),nrow(X),replace=TRUE)
beta <- matrix(c(1,2,3,4),p*k,1)

arrays = cbind(X%*%beta[(1+0*p):(p+0*p)],X%*%beta[(1+1*p):(p+1*p)])
py = t(apply(arrays,1,softmax))
y = apply((py==apply(py,1,max)),1,which)-1



beta0 <- matrix(1,k*p,1) # initial starting vector
lambda <- 5

softmax <- function(x) {
  exp_x <- exp(x - max(x))
  return(exp_x / sum(exp_x))
}



library(nnet)

f <- function(beta) {
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_one_hot <- class.ind(y)  # 需要nnet包
  loss <- -sum(y_one_hot * log(softmax_scores)) / n
  return(loss)
}

gradf <- function(beta) {
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_one_hot <- class.ind(y)  # 需要nnet包
  grad_matrix <- t(X) %*% (softmax_scores - y_one_hot) / n
  return(matrix(grad_matrix, nrow = p * k, ncol = 1))
}


# f1<-function(beta) {
#   arrays = cbind(X%*%beta[(1+0*p):(p+0*p)],X%*%beta[(1+1*p):(p+1*p)])
#   return(-sum(log(t(apply(arrays,1,softmax)))*cbind(y==0,y==1))/n)
# }
# 
# gradf1<-function(beta){
#   X_T = t(X)
#   arrays = cbind(X%*%beta[(1+0*p):(p+0*p)],X%*%beta[(1+1*p):(p+1*p)])
#   p_arrays = (1-t(apply(arrays,1,softmax)))*cbind(y==0,y==1)
#   t= rbind(X_T%*%p_arrays[,1],X_T%*%p_arrays[,2])/n
#   return(-t)
# }



beta_old = 0.5
beta_new = beta0
i=1

while (sum(abs(beta_new-beta_old))>0.0001) {
  print(i)
  beta_old = beta_new
  beta_new = beta_old - 0.1*gradf(beta_old)
  i=i+1
}

c(beta_new[3]-beta_new[1],beta_new[4]-beta_new[2])




