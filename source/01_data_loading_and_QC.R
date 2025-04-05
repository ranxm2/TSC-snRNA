setwd("E:/shiyiLab/WeiboNiu_data/singleCellAnalysis20231204/03.scRNA.data.QC/")
# load data from cellranger and perform basic QC
library(Seurat)
library(tidyverse)
library(data.table)
library(scDblFinder)

### read sample meta information
sample.info <- openxlsx::read.xlsx(
  "E:/shiyiLab/WeiboNiu_data/singleCellAnalysis20231204/data Table 2.1 Sample Sequencing Statistics_snRNAseq_snATACseq_WN.xlsx",
  sheet = 1)
sample_meta <- sample.info[c(2,4,5,7), 2:3]
sample_meta$group <- c('PSZ-6', 'PSZ-7-edge', 'PSZ-7-tuber', 'PSZ-7-outside')
#colnames(sample_meta) <- c("name", "group")
#sample_meta$name <- sub("_", "-", sample_meta$name)

### read cellranger out into Seurat object
lapply(1:nrow(sample_meta), function(i) {
  print(i)
  h5 <- sprintf(
    "../01.cellrangerOuts/%s/filtered_feature_bc_matrix.h5",
    sample_meta$Sample.ID[i])
  print(h5)
  obj <- CreateSeuratObject(
    counts = Read10X_h5(h5),
    project = sample_meta$Sample.ID[i]
  )
  obj$group <- sample_meta$group[i]
  # # add sample name to cell id
  obj <- RenameCells(obj, add.cell.id = sample_meta$Sample.ID[i])
  obj
}) -> sample_seurat_list
names(sample_seurat_list) <- sample_meta$Sample.ID


### Data QC
lapply(1:nrow(sample_meta), function(i) {
  print(i)
  obj <- sample_seurat_list[[i]]
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # before QC plot
  VlnPlot(obj,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3
  )
  ggsave(sprintf("./%s_before_QC.png", sample_meta$Sample.ID[i]),
         height = 5, width = 8
  )
  
  # doublets detection
  sce <- scDblFinder(as.SingleCellExperiment(obj))
  assertthat::are_equal(colnames(sce), colnames(obj))
  obj$scDblFinder.class <- sce$scDblFinder.class
  
  # QC
  obj <- subset(obj,
                subset =
                  nFeature_RNA > 200 &
                  nFeature_RNA < 10000 &
                  percent.mt < 20 &
                  scDblFinder.class == "singlet"
  )
  
  # after QC plot
  VlnPlot(obj,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3
  )
  ggsave(sprintf("./%s_after_QC.png", sample_meta$Sample.ID[i]),
         height = 5, width = 8
  )
  obj
}) -> sample_seurat_list
names(sample_seurat_list) <- sample_meta$Sample.ID



# save clean data
saveRDS(sample_seurat_list,
        file = "./03_sample_seurat_obj_list.rds",
        compress = F,
)
saveRDS(sample_meta,
        file = "./03_sample_meta.rds"
)
