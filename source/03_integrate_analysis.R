setwd('E:/shiyiLab/WeiboNiu_data/singleCellAnalysis20231204/05.integrate.analysis/')
library(tidyverse)
library(Seurat)

merged <- readRDS("../04.scRNA.data.merge/04_merged.rds")
sample_meta <- readRDS("../03.scRNA.data.QC/03_sample_meta.rds")

#remove C3-B
#sample_meta %>% filter(name != "C3-B") -> sample_meta
#merged <- subset(merged, subset = orig.ident != "C3-B")
merged@meta.data$group <- factor(merged@meta.data$group, ordered = TRUE, 
                                 levels = c('PSZ-6', 'PSZ-7-tuber', 'PSZ-7-edge', 'PSZ-7-outside'))
DefaultAssay(merged)

FeaturePlot(merged, "TMSB10", split.by = "group")
ggsave("05_TMSB10.pdf",
       width = 10, height = 5)


merged %>%
  FindNeighbors() %>%
  FindClusters(resolution = seq(0.1, 1, 0.1)) -> merged
saveRDS(merged, "05_merged_with_clustering.rds",
        compress = F)

DimPlot(merged, 
        label = T,
        split.by = "group")
options(mc.cores = 5)
pbmcapply::pbmclapply(seq(0.1, 1, 0.1), function(res){
  #res <- 0.1
  resolution <- sprintf("integrated_snn_res.%s", res)
  merged$seurat_clusters <- merged[[resolution]]
  Idents(merged) <- "seurat_clusters"
  DimPlot(merged, 
          label = T,
          split.by = "group")
  ggsave(
    sprintf("05_umap_resolution_%s.pdf", res),
    width = 8, height = 5)
  
  merged@meta.data %>%
    group_by(seurat_clusters, group) %>%
    count() %>%
    print(n = 50) %>%
    rename(n_cells = n) %>%
    write.csv(sprintf("05_cluster_size_resolution_%s.csv", res),
              row.names = F, quote = F)
  
  ## Markers
  FindAllMarkers(merged, assay = "RNA",
                 only.pos = T, densify = T) -> cluster.markers
  write.csv(cluster.markers,
            file = sprintf("05_cluster_gene_markers_resolution_%s.csv", res))
})

# convert the cluster gene markers csv to xlsx
lapply(seq(0.1, 1, 0.1), function(res){
  csvfile <- sprintf("05_cluster_gene_markers_resolution_%s.csv", res)
  csvdata <- read.table(csvfile, sep = ',', header = T)
  #print(csvdata[1:2, ])
  lapply(unique(csvdata$cluster), function(cl){
    clusterData <- csvdata[which(csvdata$cluster==cl), ]
    clusterData
  }) -> xlsxdata
  names(xlsxdata) <- paste0('cluster', unique(csvdata$cluster))
  xlsxdata %>% 
    openxlsx::write.xlsx(file = sprintf('05_cluster_gene_markers_resolution_%s.xlsx', res))
})