library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(patchwork)

# Failed Data reading part
# preprocessed_data = Read10X_h5(file = "preprocessed_adata.h5ad")
# 
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# 
# Convert("preprocessed_adata.h5ad", "preprocessed_adata.h5seurat")
# 
# seuratObject <- LoadH5Seurat("preprocessed_adata.h5seurat",meta.data = FALSE,misc = FALSE)


# Read with zellkonverter
# library(SeuratDisk)
# library(patchwork)

# Installing zellkonverter
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("zellkonverter")


ad = zellkonverter::readH5AD("preprocessed_adata.h5ad")
adata_Seurat <- as.Seurat(ad, counts = "X", data = NULL)
rm(ad);gc();


head(adata_Seurat)

check_basic = function(X,condition=F,table_condition=T){
  print("head")
  print(head(X))
  print("tail")
  print(tail(X))
  print("dim")
  print(dim(X))

  print("length")
  print(length(X))
  print("class")
  print(class(X))
  print("str")
  print(str(X))
  print("ncol")
  print(ncol(X))
  print("nrow")
  print(nrow(X))
  print("levels")
  print(levels(X))
  print("summary")
  print(summary(X))
  
  # print("table")
  # print(table(X))
  if(table_condition)
  {
    print("table")
    print(table(X))
  }
  if(condition)
  {
    print("X[1:5,1:5]")
    print(X[1:5,1:5])
  }

}
slots <- slotNames(adata_Seurat)

# 遍历槽并应用 summary() 函数
for (slot in slots) {
  # 获取槽的内容
  slot_content <- slot(adata_Seurat, slot)
  cat("size of the object ",slot," ",lobstr::obj_size(slot_content),"\n")
  # 检查是否为列表
  if (is.list(slot_content)) {
    # 对列表应用 summary()
    
    cat("Summary of slot", slot, ":\n")
    print(summary(slot_content))
  }
}


save(adata_Seurat,file="adata_Seurat.Rdata")

adata_assays = slot(adata_Seurat,"assays")

save(adata_assays,file="adata_assays.Rdata")


save(adata_assays$originalexp$data@x,file="adata_assays_originalexp_x.Rdata")

x_data = adata_assays$originalexp$data@x

save(x_data,file="adata_assays_originalexp_x.RData")

count_data = adata_assays$originalexp@counts

save(count_data,file="adata_assays_originalexp_counts.RData")

# check_basic(adata_assays)


# adata_Seurat[["percent.mt"]] <- PercentageFeatureSet(adata_Seurat, pattern = "^MT-")
# 
# VlnPlot(adata_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# FeatureScatter(adata_Seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt")
