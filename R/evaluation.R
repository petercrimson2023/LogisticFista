#' @import ggplot2
#' @import patchwork
#' @import caret
#'
#' Softmax Transformation
#'
#' This function applies the softmax transformation to a numeric vector,
#' converting it into a vector of probabilities that sum to 1.
#'
#' @param x Numeric vector to be transformed.
#' @return A numeric vector of probabilities.
#' @examples
#' x <- c(1, 2, 3)
#' softmax(x)
#'
#' Softmax-based Prediction
#'
#' Predicts class labels for a dataset using softmax regression.
#' Requires regression coefficients, feature matrix, labels, and class details.
#'
#' @param beta Regression coefficients.
#' @param X Feature matrix.
#' @param y True labels.
#' @param p Number of features.
#' @param k Number of classes.
#' @return Predicted class labels.
#' @examples
#' # Assuming beta, X, y, p, and k are defined
#' softmax_predict(beta, X, y, p, k)
#'
#' Plotting Softmax Regression Results
#'
#' Visualizes the results of softmax regression by generating a series of plots
#' for each class showing the regression coefficients and a plot for the loss function.
#'
#' @param result Results of the softmax regression.
#' @examples
#' # Assuming result is defined
#' plot_function(result)
#'
#' Accuracy Metrics for Softmax Regression
#'
#' Calculates and returns accuracy metrics like Sensitivity, Specificity, F1 Score,
#' and Balanced Accuracy for a softmax regression model.
#'
#' @param result Results of the softmax regression.
#' @param x_matrix Input feature matrix.
#' @param y True labels.
#' @param p Number of features.
#' @param k Number of classes.
#' @return Confusion matrix and accuracy metrics.
#' @examples
#' # Assuming result, x_matrix, y, p, and k are defined
#' accuracy_function(result, x_matrix, y, p, k)
#'
#'
#' @export
#'






library(ggplot2)
library(patchwork)
library(caret)


softmax <- function(x)
{
  exp_x <- exp(x - max(x))
  return(exp_x / sum(exp_x))
}



softmax_predict=function(beta,X,y,p,k)
{
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_predict <- apply(softmax_scores,1,which.max)
  return(y_predict-1)
}

plot_function=function(result)
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

  #ggsave(file.path(dir_name,"combined_plot.png"), plot = combined_plot, width = 10, height = 8, units = "in")

  print(combined_plot)

}

accuracy_function = function(result,x_matrix,y,p,k)
{
  class_list = c("acinar","alpha","beta","delta","ductal","gamma","Class Combined")

  y_predict = softmax_predict(result$theta,x_matrix,y,p,k)

  confusion_matrix = confusionMatrix(factor(y_predict),factor(y))$byClass[,c("Sensitivity","Specificity","F1","Balanced Accuracy")]

  rownames(confusion_matrix) = class_list

  #write.csv(confusion_matrix*100,file=file.path(dir_name,"confusion_matrix.csv"))

  return(confusion_matrix*100)
}
