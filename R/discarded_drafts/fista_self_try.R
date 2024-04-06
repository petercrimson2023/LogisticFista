# set.seed(12345)
# n <- 500
# p <- 30
# k = 3
# X <- matrix(rnorm(n*p),n,p)
# y <- sample(c(0,1,2),nrow(X),replace=TRUE)
# 
# beta <- matrix(c(1,2,3,4,3,2),p*k,1)
# beta0 <- matrix(1,k*p,1) # initial starting vector
# lambda <- 1e-5

### simple test case

# n <- 4
# p <- 2
# k = 3
# X <- matrix(c(1,1,1,1,3,2,5,4),n,p)
# y <- matrix(c(0,1,2,1),nrow(X))
# 
# beta <- matrix(c(1,2,-2,3,0,1),p*k,1)
# beta0 <- matrix(0.1,k*p,1) # initial starting vector
# lambda <- 0.2
# L = 1

softmax <- function(x) 
{
  exp_x <- exp(x - max(x))
  return(exp_x / sum(exp_x))
}

#library(nnet)

f <- function(beta,X,y,n,p,k) 
{
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_one_hot <- class.ind(y)  # 需要nnet包
  loss <- -sum(y_one_hot * log(softmax_scores+1e-22)) / n
  return(loss)
}

gradf <- function(beta,X,y,n,p,k)
{
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_one_hot <- class.ind(y)  # 需要nnet包
  grad_matrix <- t(X) %*% (softmax_scores - y_one_hot) / n
  return(matrix(grad_matrix, nrow = p * k, ncol = 1))
}

g <- function(beta,lambda) { lambda*norm(as.matrix(beta),'1') }

gradg <- function(beta,tau,lambda) { sign(beta)*(sapply(abs(beta) - tau*lambda,
                                                        FUN=function(x) {max(x,0)})) }

p_y = function(lambda,L,theta,X,y,n,p,k)
{
  u = theta - 1/L * gradf(theta,X,y,n,p,k)
  return(gradg(u,1/L,lambda))
}

Q = function(theta1,theta2,X,y,lambda,L,n,p,k)
{
  arg = f(theta2,X,y,n,p,k)+t(theta1-theta2) %*% gradf(theta2,X,y,n,p,k) +L/2 * t(theta1-theta2) %*% (theta1-theta2)+g(theta1,lambda)
  return(arg)
}


softmax_predict=function(beta,X,y,p,k)
{
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_predict <- apply(softmax_scores,1,which.max)
  return(y_predict-1)
}


fista = function(Loss,g,X,y,k,p,n,p_y,Q,lambda,theta0,L_init,max_iter=10000,eps=1e-10,eita=1.2,loss_compute=FALSE)
{
  # Initialization
  L_old = L_init
  gama = theta0
  theta_old = theta0
  t_old = 1
  loss_list = c(100)
  
  # Iteration
  
  condition = TRUE
  times=1
  
  while(condition)
  {
    ik = 1
    smallest_ik_condition = T
    
    while (smallest_ik_condition)
    {
      L_bar = eita^ik * L_old
      
      p_l_gama = p_y(lambda,L_bar,gama,X,y,n,p,k)
      fvalue = Loss(p_l_gama,X,y,n,p,k)+g(p_l_gama,lambda)
      qvalue = Q(p_l_gama,gama,X,y,lambda,L_bar,n,p,k)
      
      
      if (fvalue <= qvalue)
      {
        smallest_ik_condition = F
      }
      else
      {
        ik = ik+1
      }
    }
    
    loss_list = c(loss_list,fvalue)
    
    L_new = L_bar
    theta_new = p_y(lambda,L_new,gama,X,y,n,p,k)
    tk = (1+sqrt(1+4*t_old^2))/2
    gama = theta_new + (t_old-1)/tk * (theta_new-theta_old)
    t_old = tk
    L_old = L_new
    
    
    if(!loss_compute){
      end_loop_condition = times>max_iter || (max(abs(theta_new-theta_old))<eps)
    }else{
      end_loop_condition = times>max_iter || max(abs(theta_new-theta_old))<eps || (abs(loss_list[times+1]-loss_list[times])<1e-3)
    }
    
    
    
    
    if (end_loop_condition)
    {
      condition = F
    }
    else
    {
      # print(c(times,max(abs(theta_new-theta_old)),fvalue))
      times = times+1
      theta_old = theta_new
    }
    
  }
  
  result_list = list(theta = theta_new,loss = loss_list,iter_times = times)
  
  return(result_list)
  
}

create_file = function(lambda)
{
  dir_name = paste("result_",lambda,sep="")
  dir.create(dir_name)
}

plot_function=function(result,dir_name)
{
  g1 = ggplot(data.frame(x=1:21215,y=result$theta[,1]),aes(x=x,y=y))+
    geom_point()+
    ggtitle("Class 1")
  
  g2 = ggplot(data.frame(x=1:21215,y=result$theta[,2]),aes(x=x,y=y))+
    geom_point()+
    ggtitle("Class 3")
  
  g3 = ggplot(data.frame(x=1:21215,y=result$theta[,3]),aes(x=x,y=y))+
    geom_point()+
    ggtitle("Class 4")
  
  g4 = ggplot(data.frame(x=1:21215,y=result$theta[,4]),aes(x=x,y=y))+
    geom_point()+
    ggtitle("Class 5")
  
  g5 = ggplot(data.frame(x=1:21215,y=result$theta[,5]),aes(x=x,y=y))+
    geom_point()+
    ggtitle("Class 6")
  
  g6 = ggplot(data.frame(x=1:21215,y=result$theta[,6]),aes(x=x,y=y))+
    geom_point()+
    ggtitle("Class 9")
  
  g7 = ggplot(data.frame(x=1:21215,y=result$theta[,7]),aes(x=x,y=y))+
    geom_point()+
    ggtitle("Class combined")
  
  g8 = ggplot(data.frame(x=2:length(result$loss),loss=result$loss[2:length(result$loss)]),
              aes(x=x,y=loss))+
    geom_line()+
    ggtitle("Loss")
  
  combined_plot=g1+g2+g3+g4+g5+g6+g7+g8+plot_layout(ncol=3)
  
  ggsave(file.path(dir_name,"combined_plot.png"), plot = combined_plot, width = 10, height = 8, units = "in")
  
}

write_coefficients=function(result,dir_name)
{
  
  class_list = c("Class 1","Class 3","Class 4","Class 5","Class 6","Class 9","Class Combined")
  
  for (i in 1:7){
    cat("Non Zero coefficients for predicting",class_list[i],"\n")
    
    non_zero_indecis = c(1:21215)[(abs(result$theta[,i])>1e-7)]
    
    file_name = file.path(dir_name,paste("non_zero_coefficients_for_predicting_",class_list[i],".txt",sep=""))
    
    write.table(non_zero_indecis,file=file_name,row.names = F,col.names = F)
  }
}

write_confusion = function(result,dir_name)
{
  class_list = c("Class 1","Class 3","Class 4","Class 5","Class 6","Class 9","Class Combined")
  
  y_predict = softmax_predict(result$theta,x_matrix,y,p,k)
  
  confusion_matrix = confusionMatrix(factor(y_predict),factor(y))$byClass[,c("Sensitivity","Specificity","F1","Balanced Accuracy")]
  
  rownames(confusion_matrix) = class_list
  
  write.csv(confusion_matrix*100,file=file.path(dir_name,"confusion_matrix.csv"))
}





# Test on iris data
# data("iris")
# full_data = iris
# full_data$Species = as.numeric(full_data$Species)
# X = matrix(unlist(full_data[,1:4]),nrow=nrow(full_data),ncol=4)
# X = cbind(rep(1,nrow(full_data)),X)
# y = full_data[,5]-1
# 
# k = length(table(y))
# p = ncol(X)
# beta0 = matrix(0.1,k*p,1) # initial starting vector
# lambda = 0.1
# L_init = 1
# n = nrow(X)
# 
# result = fista(f,g,X,y,k,p,n,p_y,Q,lambda,beta0,L_init,max_iter=10000)
# 
# y_predict = softmax_predict(result$theta,X,y)
# 
# 
# 
# plot(1:length(result$loss),result$loss,type='l',xlab='iteration',ylab='loss')
# print(sum(y_predict==y)/n)
# print(result$theta)
# 
# precision_list = c()
# 
# for (i in 1:100){
#   lambda_i =i/1000
#   result = fista(f,g,X,y,k,p,n,p_y,Q,lambda_i,beta0,L_init,max_iter=10000)
#   y_predict = softmax_predict(result$theta,X,y)
#   precsion = sum(y_predict==y)/n
#   precision_list = c(precision_list,precsion)
# }
# 
# plot(1:100,precision_list,type='l',xlab='lambda',ylab='precision')




