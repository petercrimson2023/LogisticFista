current_path <- getwd()
dir_list = list.files(current_path, full.names = TRUE)


class_1_df <- data.frame(
  lambda = numeric(),
  F1 = numeric(),
  balanced_accuracy = numeric(),
  n_list = numeric(),
  Iteration_times = numeric(),
  final_loss = numeric()
)

class_3_df <- data.frame(
  lambda = numeric(),
  F1 = numeric(),
  balanced_accuracy = numeric(),
  n_list = numeric(),
  Iteration_times = numeric(),
  final_loss = numeric()
)

class_4_df <- data.frame(
  lambda = numeric(),
  F1 = numeric(),
  balanced_accuracy = numeric(),
  n_list = numeric(),
  Iteration_times = numeric(),
  final_loss = numeric()
)

class_5_df <- data.frame(
  lambda = numeric(),
  F1 = numeric(),
  balanced_accuracy = numeric(),
  n_list = numeric(),
  Iteration_times = numeric(),
  final_loss = numeric()
)

class_6_df <- data.frame(
  lambda = numeric(),
  F1 = numeric(),
  balanced_accuracy = numeric(),
  n_list = numeric(),
  Iteration_times = numeric(),
  final_loss = numeric()
)

class_9_df <- data.frame(
  lambda = numeric(),
  F1 = numeric(),
  balanced_accuracy = numeric(),
  n_list = numeric(),
  Iteration_times = numeric(),
  final_loss = numeric()
)

class_combined_df <- data.frame(
  lambda = numeric(),
  F1 = numeric(),
  balanced_accuracy = numeric(),
  n_list = numeric(),
  Iteration_times = numeric(),
  final_loss = numeric()
)




for (dir_name in dir_list[1:length(dir_list)-1]) {
  path <- dir_name
  print(path)

  param_value <- numeric(length(path))
  

  

  match <- regexpr("[0]+\\.[0-9]+", path)
  if (match != -1) {
    param_value <-
      as.numeric(substr(path, match, match + attr(match, "match.length") - 1))
  }
  
  print(param_value)

load(file=file.path(path,"result.RData"))
betas = result$theta

final_loss = result$loss[length(result$loss)]
Iteration_times=length(result$loss)-1
n_list = list()

confusion_matrix = read.csv(file=file.path(path,"confusion_matrix.csv"))
F1 = confusion_matrix$F1
balanced_accuracy = confusion_matrix$Balanced.Accuracy
for (i in 1:7)
{
  n_list[length(n_list)+1] = sum(abs(betas[,i])>1e-7)
}
n_list=unlist(n_list)

class_1_df = rbind(class_1_df, 
                   data.frame(param_value,
                              F1[1],
                              balanced_accuracy[1],
                              n_list[1],
                              Iteration_times,
                              final_loss))


class_3_df = rbind(class_3_df,
                   data.frame(param_value,
                              F1[2],
                              balanced_accuracy[2],
                              n_list[2],
                              Iteration_times,
                              final_loss))


class_4_df = rbind(class_4_df,
                   data.frame(param_value,
                              F1[3],
                              balanced_accuracy[3],
                              n_list[3],
                              Iteration_times,
                              final_loss))


class_5_df = rbind(class_5_df,
                   data.frame(param_value,
                              F1[4],
                              balanced_accuracy[4],
                              n_list[4],
                              Iteration_times,
                              final_loss))

class_6_df = rbind(class_6_df,
                   data.frame(param_value,
                              F1[5],
                              balanced_accuracy[5],
                              n_list[5],
                              Iteration_times,
                              final_loss))


class_9_df = rbind(class_9_df,
                   data.frame(param_value,
                              F1[6],
                              balanced_accuracy[6],
                              n_list[6],
                              Iteration_times,
                              final_loss))


class_combined_df = rbind(class_combined_df,
                          data.frame(param_value,
                                     F1[7],
                                     balanced_accuracy[7],
                                     n_list[7],
                                     Iteration_times,
                                     final_loss))
cat("finish",param_value,"\n")
rm(result)
}


colnames(class_1_df) <- c("param_value", "F1", "balanced_accuracy", "n_list", "Iteration_times", "final_loss")
colnames(class_3_df) <- c("param_value", "F1", "balanced_accuracy", "n_list", "Iteration_times", "final_loss")
colnames(class_4_df) <- c("param_value", "F1", "balanced_accuracy", "n_list", "Iteration_times", "final_loss")
colnames(class_5_df) <- c("param_value", "F1", "balanced_accuracy", "n_list", "Iteration_times", "final_loss")
colnames(class_6_df) <- c("param_value", "F1", "balanced_accuracy", "n_list", "Iteration_times", "final_loss")
colnames(class_9_df) <- c("param_value", "F1", "balanced_accuracy", "n_list", "Iteration_times", "final_loss")
colnames(class_combined_df) <- c("param_value", "F1", "balanced_accuracy", "n_list", "Iteration_times", "final_loss")


library(ggplot2)
library(patchwork)

class_plot = function(class_df, file_name, point_size = 0.5, save_condition = TRUE) {
  g1 = ggplot(class_df, aes(x = param_value, y = F1)) +
    geom_line(color = "blue") +  
    geom_point(size = point_size, color = "red") + 
    labs(title = "F1 vs lambda",
         x = "lambda",
         y = "F1")
  
  g2 = ggplot(class_df, aes(x = param_value, y = balanced_accuracy)) +
    geom_line(color = "green") +  
    geom_point(size = point_size, color = "orange") +  
    labs(title = "balanced_accuracy vs lambda",
         x = "lambda",
         y = "balanced_accuracy")
  
  g3 = ggplot(class_df, aes(x = param_value, y = n_list)) +
    geom_line(color = "purple") +  
    geom_point(size = point_size, color = "pink") +  
    labs(title = "n_list vs lambda",
         x = "lambda",
         y = "n_list")
  
  g4 = ggplot(class_df, aes(x = param_value, y = Iteration_times)) +
    geom_line(color = "brown") +  
    geom_point(size = point_size, color = "gray") + 
    labs(title = "Iteration_times vs lambda",
         x = "lambda",
         y = "Iteration_times")
  
  g5 = ggplot(class_df, aes(x = param_value, y = final_loss)) +
    geom_line(color = "blue") +  
    geom_point(size = point_size, color = "red") +  
    labs(title = "final_loss vs lambda",
         x = "lambda",
         y = "final_loss")
  
  g_finale = g1 + g2 + g3 + g4 + g5 + plot_layout(ncol = 3) +
    ggtitle(paste("Class Plot: ", file_name))
  
  if (save_condition) {
    ggsave(paste0(file_name, ".png"), plot = g_finale, width = 10, height = 8, dpi = 400)
  } else {
    print(g_finale)
  }
}



# 
class_plot(class_1_df,"acinar")
class_plot(class_3_df,"alpha")
class_plot(class_4_df,"beta")
class_plot(class_5_df,"delta")
class_plot(class_6_df,"ductal")
class_plot(class_9_df,"gamma")
class_plot(class_combined_df,"class_combined")
# 


# extracting results with lambda between 0.026 to 0.040

# class_1_df_026_040 = class_1_df[class_1_df$param_value >= 0.026 & class_1_df$param_value <= 0.040, ]
# class_3_df_026_040 = class_3_df[class_3_df$param_value >= 0.026 & class_3_df$param_value <= 0.040, ]
# class_4_df_026_040 = class_4_df[class_4_df$param_value >= 0.026 & class_4_df$param_value <= 0.040, ]
# class_5_df_026_040 = class_5_df[class_5_df$param_value >= 0.026 & class_5_df$param_value <= 0.040, ]
# class_6_df_026_040 = class_6_df[class_6_df$param_value >= 0.026 & class_6_df$param_value <= 0.040, ]
# class_9_df_026_040 = class_9_df[class_9_df$param_value >= 0.026 & class_9_df$param_value <= 0.040, ]
# class_combined_df_026_040 = class_combined_df[class_combined_df$param_value >= 0.026 & class_combined_df$param_value <= 0.040, ]
# 
# 
# class_plot(class_1_df_026_040,"class_1_026_040",save_condition = TRUE)
# class_plot(class_3_df_026_040,"class_3_026_040",save_condition = TRUE)
# class_plot(class_4_df_026_040,"class_4_026_040",save_condition = TRUE)
# class_plot(class_5_df_026_040,"class_5_026_040",save_condition = TRUE)
# class_plot(class_6_df_026_040,"class_6_026_040",save_condition = TRUE)
# class_plot(class_9_df_026_040,"class_9_026_040",save_condition = TRUE)
# class_plot(class_combined_df_026_040,"class_combined_026_040",save_condition = TRUE)



