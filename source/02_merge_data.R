setwd("E:/shiyiLab/WeiboNiu_data/singleCellAnalysis20231204/04.scRNA.data.merge/")
# merge data from same group
library(Seurat)
library(tidyverse)

# load data
sample_meta <- readRDS("../03.scRNA.data.QC/03_sample_meta.rds")
sample_seurat_list <- readRDS("../03.scRNA.data.QC/03_sample_seurat_obj_list.rds")

#remove C3-B
#sample_meta %>% filter(name != "C3-B") -> sample_meta
#sample_seurat_list[sample_meta$name] -> sample_seurat_list

# merged <- merge(x = sample_seurat_list[[1]],
#                 y = sample_seurat_list[2:11],
#                 add.cell.ids = sample_meta$name,
#                 project = "TSC_CON")



# normalize data, find variable genes, scale data, run PCA
sample_seurat_list %>%
  lapply(., function(obj) {
    obj %>%
      NormalizeData() %>%
      FindVariableFeatures(nfeatures = 3000) %>%
      ScaleData() %>%
      RunPCA()
  }) -> sample_seurat_list

# select integration features (genes)
integration_features <- SelectIntegrationFeatures(sample_seurat_list,
                                                  nfeatures = 2000
)

# use integration features to do the PCA again
# regress out the batch effect and prepare for merging
# sample_seurat_list %>%
#   lapply(., function(obj) {
#     obj %>%
#       RunPCA(features = integration_features)
#   }) -> sample_seurat_list

# integrate data
integration_anchors <- 
  FindIntegrationAnchors(sample_seurat_list,
                         anchor.features = integration_features,
                        # reduction = "rpca"
                        )

merged <- IntegrateData(integration_anchors)
DefaultAssay(merged) <- "integrated"

merged %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(reduction = "pca", dims = 1:50) -> merged

DimPlot(merged,
        group.by = "orig.ident", label = T)
ggsave("./04_umap_all_sample.png",
       width = 9, height = 7)

DimPlot(merged,
        group.by = "orig.ident", label = T, split.by = "group")
ggsave("./04_umap_all_sample_split_by_group.png",
       width = 12, height = 7)


DimPlot(merged,
        group.by = "group", label = T)
ggsave("./04_umap_group.png",
       width = 9, height = 7)

saveRDS(merged, "./04_merged.rds",
        compress = F)

save.image(file = '04_merge_data.RData')
