# load("adata_assays_originalexp_x.RData")
#
# load("adata_assays.Rdata")

#load("adata_assays_originalexp_counts.RData")

# x_matrix = as.matrix(count_data,nrow=21215,ncol=6321)
#
# x_matrix = t(x_matrix)

#save(x_matrix,file="x_matrix.RData")

#X = x_matrix;rm(x_matrix);rm(count_data);gc()

load("x_matrix.RData")

#load("adata_Seurat.Rdata")

library(Rcpp)
library(caret)
library(ggplot2)
library(patchwork)
sourceCpp("test.cpp")
source("fista_self_try.R")



test_n =n




temp_small_data = read.csv("adata_nn.csv")
x_small = as.matrix(temp_small_data[,1:32])
n = dim(x_small)[1]
p = dim(x_small)[2]
y_raw = temp_small_data$celltype
y_num = as.numeric(factor(y_raw))
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


beta0 = matrix(0.1,k*p,1)
lambda = 0.05
L_init = 1

y_one_hot=nnet::class.ind(y)

# sample_index = sample(1:n,test_n)
#
# X = x_matrix[sample_index,]
# y_one_hot = y_one_hot[sample_index,]
# y = y[sample_index]


result = fistaCpp(lambda,L_init,beta0,x_small,y_one_hot,n=test_n,p=p,k=k)

save(result,file="result_fistaCpp_lamba_0.05.RData")



y_predict = softmax_predict(result$theta,x_small,y,p,k)

print(paste("general accuracy:",sum(y_predict==y)/length(y)))

#result_matrix = as.matrix(table(y_predict,y))



print(confusionMatrix(factor(y_predict),factor(y)))

cofusion_matrix = confusionMatrix(factor(y_predict),factor(y))

for (i in 1:7){
  cat("Non Zero coefficients for predicting class ",i,"\n")
  print(c(1:21215)[(abs(result$theta[,i])>1e-7)])
}

png("loss_plot.png", width = 800, height = 600)
print(plot_loss(result$loss))
dev.off()


## Ploting the coefficients for each class

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

g1+g2+g3+g4+g5+g6+g7+g8+plot_layout(ncol=3)



