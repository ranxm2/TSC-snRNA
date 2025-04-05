setwd("E:/shiyiLab/WeiboNiu_data/singleCellAnalysis20231204/06.cell.type.annotation/")
library(Seurat)
library(tidyverse)
library(openxlsx)

# use cell marker db to generate cell annotation
# this is for reference only, may need manual updates

# load cell marker database
# http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html
# http://117.50.127.228/CellMarker/CellMarker_download.html
markers_db <- read.xlsx("http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Seq.xlsx",
                        check.names = F, sheet = 1)
unique(markers_db$tissue_type)
markers_db %>% filter(
  # tissue_type == "Brain",
  species == "Human",
                      cancer_type == "Normal") -> markers_db
markers_db %>% group_by(cell_name) %>%
  summarise(total = n()) -> cell_type_marker_count

### load data
merged05 <- readRDS("../05.integrate.analysis/05_merged_with_clustering.rds")
sample_meta <- readRDS("../03.scRNA.data.QC/03_sample_meta.rds")
#sample_meta %>% filter(name != "C3-B") -> sample_meta

lapply(seq(0.1, 1, 0.1), function(res){
  merged <- merged05
  merged$seurat_clusters <- merged[[sprintf("integrated_snn_res.%s", res)]]
  Idents(merged) <- "seurat_clusters"
  DimPlot(merged, label = T, split.by = "group")
  
  cluster_markers <- read.csv(
    sprintf("../05.integrate.analysis/05_cluster_gene_markers_resolution_%s.csv", res))
  table(table(cluster_markers$gene))
  
  # significant markers
  cluster_markers %>%
    filter(p_val_adj < 0.05, 
           avg_log2FC > 0) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 500) -> sig_markers
  
  # expression data
  DefaultAssay(merged) <- "RNA"
  FetchData(merged, c(unique(sig_markers$gene), "seurat_clusters")) -> expr_data
  expr_data %>%
    pivot_longer(!seurat_clusters, names_to = "gene", values_to = "expr") %>%
    group_by(seurat_clusters, gene) %>%
    summarise(avg_expr = round(mean(expr), 4)) -> avg_expr_data
  
  lapply(unique(cluster_markers$cluster), function(cl){
    # cl <- 0
    print(cl)
    sig_markers %>%
      filter(cluster == cl) %>%
      arrange(desc(avg_log2FC)) -> cl_marker_genes
    
    markers_db %>%
      filter(Symbol %in% cl_marker_genes$gene) %>%
      left_join(cl_marker_genes, by = c("Symbol" = "gene")) %>%
      group_by(cell_name) %>%
      summarise(n = n(),
                markers = paste(unique(Symbol), collapse = ",")) %>%
      left_join(cell_type_marker_count) %>%
      mutate(proportion = n / total) %>%
      arrange(desc(proportion)) -> res_anno
    res_anno$cluster <- cl
    
    # add gene logfc to the annotation res
    print("adding marker gene logfc....")
    sapply(res_anno$markers, function(genes){
      genes <- genes %>% str_split(",", simplify = T) %>% .[1,]
      cluster_markers %>%
        mutate(avg_log2FC = round(avg_log2FC, 4)) %>%
        filter(cluster == cl,
               gene %in% genes) %>%
        arrange(match(gene, genes)) %>%
        select(avg_log2FC) %>%
        unlist() %>% paste0(., collapse = ",")
    }) -> res_anno$log2FC
    
    res_anno$log2FC %>% sapply(function(x){
      str_split(x ,",", simplify = T) %>% .[1,] %>%
        as.numeric() %>% sum()
    }) -> res_anno$log2FC_sum
    
    # add gene avg expr to the annotation res
    print("adding marker gene expr....")
    sapply(res_anno$markers, function(genes){
      genes <- genes %>% str_split(",", simplify = T) %>% .[1,]
      avg_expr_data %>%
        filter(seurat_clusters == cl, gene %in% genes) %>%
        arrange(match(gene, genes)) %>%
        ungroup() %>%
        select(avg_expr) %>%
        unlist() %>% paste0(., collapse = ",")
    }) -> res_anno$avg_expr
    
    res_anno$avg_expr %>% sapply(function(x){
      str_split(x ,",", simplify = T) %>% .[1,] %>%
        as.numeric() %>% sum()
    }) -> res_anno$avg_expr_sum
    
    res_anno %>%
      arrange(desc(avg_expr_sum))
  }) -> cell.type.annotation
  
  names(cell.type.annotation) <- paste0("cluster ", unique(cluster_markers$cluster))
  openxlsx::write.xlsx(
    cell.type.annotation,
    file = sprintf("./cell_type_annotation_resolution_%s.xlsx", res)
  )
  
  return(0)
})