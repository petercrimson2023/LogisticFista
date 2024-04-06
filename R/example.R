
Rcpp::sourceCpp("cpp/fista.cpp")
load("example_data/y.RData")
source("R/eval.R")

# Load data
x_matrix = read.csv("example_data/adata_pca.csv")[,1:50] %>% as.matrix()
y_true = as.factor(y)
y_one_hot=nnet::class.ind(y)

# Parameters
n = nrow(x_matrix)
p = ncol(x_matrix)
k = length(unique(y_true))
beta0 = matrix(0.1,k*p,1)
lambda = 0.045
L_init = 1

# Fitting L1 multilevel logistic regression
result = fista(lambda,L_init,beta0,x_matrix,y_one_hot,n=n,p=p,k=k)

# Evaluation and drawing
plot_thetas(result)
accuracy_function(result,x_matrix,y_true,p,k)
plot_umap(result,x_matrix,y_true,p,k)



