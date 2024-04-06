# evaluation and drawing functions

library(ggplot2)
library(patchwork)
library(caret)
library(mclust)
library(aricode)
library(RColorBrewer)
library(umap)
library(dplyr)
library(nnet)

softmax <- function(x)
{
  exp_x <- exp(x - max(x))
  return(exp_x / sum(exp_x))
}

softmax_predict=function(beta,X,y,p,k)
{
  beta_matrix = matrix(beta, nrow = p, ncol = k)
  scores = X %*% beta_matrix
  softmax_scores = t(apply(scores, 1, softmax))
  y_predict = apply(softmax_scores,1,which.max)
  return(y_predict-1)
}

tokenize_y = function(y_raw)
{
  y_num = as.numeric(factor(y_raw))
  y = ifelse(y_num %in% c(2,7, 8, 10, 11, 12, 13),15,y_num)
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
  y_one_hot=class.ind(y)
  return(y_one_hot)
}

plot_thetas=function(result,dir_name = NULL, save_file = FALSE)
{
  p = length(result$theta[,1])

  g1 = ggplot(data.frame(x = 1:p, y = result$theta[, 1]), aes(x = x, y =
                                                                y)) +
    geom_point() +
    ggtitle("acinar")

  g2 = ggplot(data.frame(x = 1:p, y = result$theta[, 2]), aes(x = x, y =
                                                                y)) +
    geom_point() +
    ggtitle("alpha")

  g3 = ggplot(data.frame(x = 1:p, y = result$theta[, 3]), aes(x = x, y =
                                                                y)) +
    geom_point() +
    ggtitle("beta")

  g4 = ggplot(data.frame(x = 1:p, y = result$theta[, 4]), aes(x = x, y =
                                                                y)) +
    geom_point() +
    ggtitle("delta")

  g5 = ggplot(data.frame(x = 1:p, y = result$theta[, 5]), aes(x = x, y =
                                                                y)) +
    geom_point() +
    ggtitle("ductal")

  g6 = ggplot(data.frame(x = 1:p, y = result$theta[, 6]), aes(x = x, y =
                                                                y)) +
    geom_point() +
    ggtitle("gamma")

  g7 = ggplot(data.frame(x = 1:p, y = result$theta[, 7]), aes(x = x, y =
                                                                y)) +
    geom_point() +
    ggtitle("Class combined")

  g8 = ggplot(data.frame(x = 2:length(result$loss), loss = result$loss[2:length(result$loss)]),
              aes(x = x, y = loss)) +
    geom_line() +
    ggtitle("Loss")

  combined_plot = g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8+
    plot_layout(ncol = 2, nrow = 4)

  if(save_file == TRUE)
  {
    ggsave(file.path(dir_name,"combined_plot.png"), plot = combined_plot, width = 10, height = 8, units = "in")

  }

  print(combined_plot)

}

accuracy_function = function(result,x_matrix,y,p,k)
{
  class_list = c("acinar","alpha","beta","delta","ductal","gamma","Class Combined")

  y_predict = softmax_predict(result$theta,x_matrix,y,p,k)

  confusion_matrix = confusionMatrix(factor(y_predict),factor(y))$byClass[,c("Sensitivity","Specificity","F1","Balanced Accuracy")]

  rownames(confusion_matrix) = class_list

  #write.csv(confusion_matrix*100,file=file.path(dir_name,"confusion_matrix.csv"))

  return(list(confustion = confusion_matrix*100,
              nmi = NMI(y,y_predict),
              ani = adjustedRandIndex(y,y_predict)))
}

plot_umap = function(result,x_matrix,y,p,k)
{
  color_palette = brewer.pal(n = k, name = "Dark2")

  class_list = class_list = c("acinar","alpha","beta","delta","ductal","gamma","Class Combined")

  y_name = sapply(as.numeric(y),function(i){
    return(class_list[i+1])
  })

  y_pred = softmax_predict(result$theta,x_matrix,y,p,k)%>% as.numeric()

  y_pred_name = sapply(y_pred,function(i){
    return(class_list[i+1])
  })

  umap_pca = umap(x_matrix)

  df_pca_model = data.frame(
    X1=umap_pca$layout[,1],
    X2=umap_pca$layout[,2],
    label=factor(y_pred_name)
  )

  df_pca_original = data.frame(
    X1=umap_pca$layout[,1],
    X2=umap_pca$layout[,2],
    label=factor(y_name)
  )

  ggplot()+
    geom_point(data=df_pca_model,aes(x=X1,y=X2,color=label),alpha=0.2)+
    geom_point(data=df_pca_original,aes(x=X1,y=X2,color=label),alpha=0.2)+
    scale_color_manual(values = color_palette) +
    labs(title="UMAP for L1 softmax")
}

