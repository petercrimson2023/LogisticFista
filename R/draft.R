#load("adata_Seurat.Rdata")

length(adata_Seurat$celltype)

y_name = adata_Seurat$celltype

cell_type_name = rownames(table(y_name))

#load("x_matrix.RData")

gene_name = colnames(x_matrix)


# Generating ANI NMI UMAP

#load("y.RData")

#load("x_matrix.RData")

y_true = as.factor(y)

source("fista_self_try.R")

#load("./eita_1_1/eita_1.1/result_list/result_0.033/result.RData")

p = ncol(x_matrix)
k = length(unique(y_true))

beta = result$theta

y_pred = softmax_predict(beta,x_matrix,y,p,k)

library(mclust)
library(aricode)

# NMI

nmi_result = NMI(y_true,y_pred)
nmi_result

# ANI

ani_result = adjustedRandIndex(y_true,y_pred)
ani_result

# UMAP

#library(umap)

#umap_result = umap(x_matrix)


save(umap_result,file="umap_result.RData")

save(y_pred,file="y_pred_0.033.RData")

# load("y_pred_0.033.RData")
# load("umap_result.RData")
# load("y.RData")

umap_model = umap_result$layout
head(umap_model)
rm(umap_result)

umap_original = umap_model
head(umap_original)

library(ggplot2)

library(RColorBrewer)


num_labels = length(unique(y))
color_palette <- brewer.pal(n =num_labels, name = "Dark2")




df_model = data.frame(X1=umap_model[,1],X2=umap_model[,1],label=factor(y_pred))

df_original = data.frame(X1=umap_original[,1],X2=umap_original[,2],label=factor(y))

ggplot()+
  geom_point(data=df_model,aes(x=X1,y=X2,color=label),alpha=0.2)+
  geom_point(data=df_original,aes(x=X1,y=X2,color=label),alpha=0.2)+
  scale_color_manual(values = color_palette) +
  theme_minimal()+
  labs(title="UMAP for L1 logistic")


x_pca = read.csv("adata_pca.csv")

cell_type_name = names(table(x_pca$celltype))

x_pca = x_pca[,1:50]

class_list = c("Class 1","Class 3","Class 4","Class 5","Class 6","Class 9","Class Combined")

class_list = c(cell_type_name[1],cell_type_name[3],cell_type_name[4],cell_type_name[5],cell_type_name[6],cell_type_name[9],"combined")

y_name = sapply(y,function(x){
  return(class_list[x+1])
})

y_pred_name = sapply(y_pred,function(x){
  return(class_list[x+1])
})

umap_pca = umap(x_pca)

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
  theme_minimal()+
  labs(title="UMAP for L1 softmax")






