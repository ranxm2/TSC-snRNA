setwd('E:/shiyiLab/WeiboNiu_data/singleCellAnalysis20231204/07.cell.cluster.DEGs/')
require(tidyverse)
require(Seurat)
require(ggpubr)
######## load data: merged and cell annotation ###########
merged <- readRDS('../05.integrate.analysis/05_merged_with_clustering.rds')
merged$seurat_clusters <- merged[["integrated_snn_res.0.1"]]
Idents(merged) <- 'seurat_clusters'
merged$celltype_group <- paste(Idents(merged), merged$group, sep = '_')
group.combn <- combn(c('PSZ-6', 'PSZ-7-outside', 'PSZ-7-edge', 
                       'PSZ-7-tuber'), 2)
#group.combn
cell_annotation <- 
  openxlsx::read.xlsx("../todolist20240108/cell_type_annotation_resolution_0.1 update20240109.xlsx", 
                      sheet = 1, colNames = F)
cell_annotation$X1 <- as.factor(cell_annotation$X1)

DimPlot(merged, label = T, split.by = 'group')
#require(ggplot2)
merged@meta.data %>%
  group_by(seurat_clusters, group) %>%
  count() %>%
  group_by(group) %>%
  mutate(percent = n / sum(n) * 100) %>%
  left_join(cell_annotation, by = c("seurat_clusters" = "X1")) %>%
  ggplot(aes(x = group, y = percent, fill = X2)) +
  geom_bar(stat = "identity", width = 0.8) + 
  labs(fill = "Cell") +
  theme_bw(base_size = 14,
           base_family = "Arial") +
  theme(axis.text = element_text(color = "black"))

########## find Marker genes for every cell cluster ##########
#merged$celltype_group
Idents(merged) <- 'celltype_group'
DefaultAssay(merged) <- 'RNA'
# test code
if (FALSE){
FindMarkers(merged, ident.1 = '0_PSZ-6',ident.2 = '0_PSZ-7-tuber', 
            verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
FindMarkers(merged, ident.2 = '0_PSZ-6',ident.1 = '0_PSZ-7-tuber', 
            verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
FindMarkers(merged, ident.1 = '0_PSZ-6', ident.2 = '0_PSZ-7-outside', 
            verbose = FALSE, test.use = 'wilcox', min.pct = 0.1)
}
#### find the degs for 6 comparison and save the result to xlsx
for(i in 1:6){
  g2 <- group.combn[1,i]
  g1 <- group.combn[2,i]
  lapply(seq(0,10), function(j){
    i2 <- paste(as.character(j), g2, sep = '_')
    i1 <- paste(as.character(j), g1, sep = '_')
    print(paste(i1, i2, sep = ' to '))
    deg <- FindMarkers(merged, ident.1 = i1, ident.2 = i2, 
                       verbose = FALSE, test.use = 'wilcox',
                       min.pct = 0.1, logfc.threshold = 0)
    deg
  }) -> deg.list
  names(deg.list) <- paste(cell_annotation$X1, 
                           cell_annotation$X2, 
                           sep = '_')
  if(FALSE){
  lapply(seq(0, 10), function(j){
    n <- paste(cell_annotation$X1[j+1], 
               cell_annotation$X2[j+1], 
               sep = '_')
    print(n)
    deg <- deg.list[[n]]
    deg[which(deg$p_val_adj< 0.05), ]
  }) -> deg.list_f0
  names(deg.list_f0) <- paste(cell_annotation$X1, 
                              cell_annotation$X2, 
                              sep = '_')
  deg.list_f0 %>% 
    openxlsx::write.xlsx(file = sprintf('DEG.20240125/07_DEGs_%s_%s_11cluster_lfc0.xlsx',
                                        g1, g2),
                         rowNames = T)
  lapply(seq(0, 10), function(j){
    n <- paste(cell_annotation$X1[j+1], 
               cell_annotation$X2[j+1], 
               sep = '_')
    print(n)
    deg <- deg.list[[n]]
    deg[which(deg$p_val_adj< 0.05 & 
                abs(deg$avg_log2FC)>0.01), ]
  }) -> deg.list_f001
  names(deg.list_f001) <- paste(cell_annotation$X1, 
                                cell_annotation$X2, 
                                sep = '_')
  deg.list_f001 %>% 
    openxlsx::write.xlsx(file = sprintf('DEG.20240125/07_DEGs_%s_%s_11cluster_lfc001.xlsx',
                                        g1, g2),
                         rowNames = T)
  }
  
  lapply(seq(0, 10), function(j){
    n <- paste(cell_annotation$X1[j+1], 
               cell_annotation$X2[j+1], 
               sep = '_')
    print(n)
    deg <- deg.list[[n]]
    deg[which(deg$p_val_adj< 0.05 &
                abs(deg$avg_log2FC)>0.25), ]
  }) -> deg.list_f025
  names(deg.list_f025) <- paste(cell_annotation$X1, 
                                cell_annotation$X2, 
                                sep = '_')
  deg.list_f025 %>% 
    openxlsx::write.xlsx(file = sprintf('DEG.20240125/07_DEGs_%s_%s_11cluster_lfc025.xlsx',
                                        g1, g2), 
                         rowNames = T)
  deg.list %>%
    saveRDS(file = sprintf("DEG.20240125/07_DEGs_%s_%s_11cluster.rds", 
                           g1, g2), compress = F)
}

##### find DEGs for psz-7-tuber and psz-6
lapply(seq(0,10), function(i){
  #i
  i2 <- paste(as.character(i), 'PSZ-6', sep = "_")
  i1 <- paste(as.character(i), 'PSZ-7-tuber', sep = '_')
  deg <- FindMarkers(merged, ident.1 = i1, ident.2 = i2,
                     verbose = FALSE, test.use = 'wilcox', 
                     min.pct = 0.1, logfc.threshold = 0)
  deg
}) -> deg.psz7tuber.psz6.11cluster
names(deg.psz7tuber.psz6.11cluster) <- paste(cell_annotation$X1, 
                                             cell_annotation$X2, 
                                             sep = '_')
deg.psz7tuber.psz6.11cluster %>%
  saveRDS(file = '07_DEGs_PSZ-7-tuber_PSZ-6_11cluster.rds', 
          compress = F)
deg.psz7tuber.psz6.11cluster %>% 
  openxlsx::write.xlsx(file = 'DEG.PSZ-7-tuber.PSZ-6.11clusters.xlsx', 
                       rowNames = TRUE)
deg.psz7tuber.psz6.11cluster %>% 
  openxlsx::write.xlsx(file = 'DEG.PSZ-7-tuber.PSZ-6.11clusters_RNA.xlsx', 
                       rowNames = TRUE)

# psz-7-tuber psz-6 up DEG
lapply(seq(0,10), function(i){
  #i
  i2 <- paste(as.character(i), 'PSZ-6', sep = "_")
  i1 <- paste(as.character(i), 'PSZ-7-tuber', sep = '_')
  deg <- FindMarkers(merged, ident.1 = i1, ident.2 = i2,
                     verbose = FALSE, test.use = 'wilcox', 
                     min.pct = 0.1, logfc.threshold = 0)
  deg[which(deg$avg_log2FC>0.25 & deg$p_val_adj<0.05), ]
}) -> deg.psz7tuber.psz6.11cluster.up
names(deg.psz7tuber.psz6.11cluster.up) <- paste(cell_annotation$X1, 
                                             cell_annotation$X2, 
                                             sep = '_')
deg.psz7tuber.psz6.11cluster.up %>%
  saveRDS(file = '07_DEGs_PSZ-7-tuber_PSZ-6_11cluster.up.rds', 
          compress = F)
deg.psz7tuber.psz6.11cluster.up %>% 
  openxlsx::write.xlsx(file = 'DEG.PSZ-7-tuber.PSZ-6.11clusters.up.xlsx', 
                       rowNames = TRUE)

# psz-7-tuber psz-6 down
lapply(seq(0,10), function(i){
  #i
  i2 <- paste(as.character(i), 'PSZ-6', sep = "_")
  i1 <- paste(as.character(i), 'PSZ-7-tuber', sep = '_')
  deg <- FindMarkers(merged, ident.1 = i1, ident.2 = i2,
                     verbose = FALSE, test.use = 'wilcox', 
                     min.pct = 0.1, logfc.threshold = 0)
  deg[which(deg$avg_log2FC< -0.25 & deg$p_val_adj<0.05), ]
}) -> deg.psz7tuber.psz6.11cluster.down
names(deg.psz7tuber.psz6.11cluster.down) <- paste(cell_annotation$X1, 
                                             cell_annotation$X2, 
                                             sep = '_')
deg.psz7tuber.psz6.11cluster.down %>%
  saveRDS(file = '07_DEGs_PSZ-7-tuber_PSZ-6_11cluster.down.rds', 
          compress = F)
deg.psz7tuber.psz6.11cluster.down %>% 
  openxlsx::write.xlsx(file = 'DEG.PSZ-7-tuber.PSZ-6.11clusters.down.xlsx', 
                       rowNames = TRUE)

# psz-7-outside psz-6
lapply(seq(0,10), function(i){
  #i
  i2 <- paste(as.character(i), 'PSZ-6', sep = "_")
  i1 <- paste(as.character(i), 'PSZ-7-outside', sep = '_')
  deg <- FindMarkers(merged, ident.1 = i1, ident.2 = i2,
                     verbose = FALSE, test.use = 'wilcox', 
                     min.pct = 0.1, logfc.threshold = 0)
  deg
}) -> deg.psz7outside.psz6.11cluster
names(deg.psz7outside.psz6.11cluster) <- paste(cell_annotation$X1, 
                                               cell_annotation$X2, 
                                               sep = '_')
deg.psz7outside.psz6.11cluster %>%
  saveRDS(file = '07_DEGs_PSZ-7-outside_PSZ-6_11cluster.rds', 
          compress = F)
deg.psz7outside.psz6.11cluster %>% 
  openxlsx::write.xlsx(file = 'DEG.PSZ-7-outside.PSZ-6.11clusters.xlsx', 
                       rowNames = TRUE)
# psz-7-outside psz-6 up
lapply(seq(0,10), function(i){
  #i
  i2 <- paste(as.character(i), 'PSZ-6', sep = "_")
  i1 <- paste(as.character(i), 'PSZ-7-outside', sep = '_')
  deg <- FindMarkers(merged, ident.1 = i1, ident.2 = i2,
                     verbose = FALSE, test.use = 'wilcox', 
                     min.pct = 0.1, logfc.threshold = 0)
  deg[which(deg$avg_log2FC>0.25 & deg$p_val_adj<0.05), ]
}) -> deg.psz7outside.psz6.11cluster.up
names(deg.psz7outside.psz6.11cluster.up) <- paste(cell_annotation$X1, 
                                               cell_annotation$X2, 
                                               sep = '_')
deg.psz7outside.psz6.11cluster.up %>%
  saveRDS(file = '07_DEGs_PSZ-7-outside_PSZ-6_11cluster.up.rds', 
          compress = F)
deg.psz7outside.psz6.11cluster.up %>% 
  openxlsx::write.xlsx(file = 'DEG.PSZ-7-outside.PSZ-6.11clusters.up.xlsx', 
                       rowNames = TRUE)
# psz-7-outside psz-6 down
lapply(seq(0,10), function(i){
  #i
  i2 <- paste(as.character(i), 'PSZ-6', sep = "_")
  i1 <- paste(as.character(i), 'PSZ-7-outside', sep = '_')
  deg <- FindMarkers(merged, ident.1 = i1, ident.2 = i2,
                     verbose = FALSE, test.use = 'wilcox', 
                     min.pct = 0.1, logfc.threshold = 0)
  deg[which(deg$avg_log2FC< -0.25 & deg$p_val_adj<0.05), ]
}) -> deg.psz7outside.psz6.11cluster.down
names(deg.psz7outside.psz6.11cluster.down) <- paste(cell_annotation$X1, 
                                                  cell_annotation$X2, 
                                                  sep = '_')
deg.psz7outside.psz6.11cluster.down %>%
  saveRDS(file = '07_DEGs_PSZ-7-outside_PSZ-6_11cluster.down.rds', 
          compress = F)
deg.psz7outside.psz6.11cluster.down %>% 
  openxlsx::write.xlsx(file = 'DEG.PSZ-7-outside.PSZ-6.11clusters.down.xlsx', 
                       rowNames = TRUE)

######### enrichment for every cell cluster ##############
require(gprofiler2)
require(ggsci)
i <- gost(rownames(deg.psz6.psz7tuber.11cluster$`0_Oligodendrocyte`) %>% unique(),
     organism = "hsapiens", 
     correction_method = "fdr", evcodes = T)
# enrichment for every cell clusters, comparison6, up, down, abs(lfc)>0.25
for (i in 1:6){
  g2 <- group.combn[1,i]
  g1 <- group.combn[2,i]
  print(paste(g1, g2, sep = ' to '))
  deg.list <- readRDS(file = sprintf('DEG.20240125/07_DEGs_%s_%s_11cluster.rds',
                                    g1, g2))
  # up enrichment
  lapply(seq(0, 10), function(j){
    s <- paste(as.character(j), cell_annotation$X2[j+1], sep = "_")
    print(s)
    degtable <- deg.list[[s]]
    degtable <- degtable[which(degtable$avg_log2FC>0.25 &
                                 degtable$p_val_adj<0.05), ]
    query <- gost(query = rownames(degtable) %>%
                    unique(), organism = 'hsapiens',
                  correction_method = 'fdr', evcodes = T)
    query.res <- query[['result']]
    query.res
  }) -> query.res.list
  names(query.res.list) <- paste(cell_annotation$X1, 
                                 cell_annotation$X2, 
                                 sep = '_')
  query.res.list %>%
    openxlsx::write.xlsx(file = sprintf('DEG.20240125/07_enrichment_%s_%s_11cluster_up_lfc025.xlsx', 
                                        g1, g2))
  # down enrichment 
  lapply(seq(0, 10), function(j){
    s <- paste(as.character(j), cell_annotation$X2[j+1], sep = "_")
    print(s)
    degtable <- deg.list[[s]]
    degtable <- degtable[which(degtable$avg_log2FC < -0.25 &
                                 degtable$p_val_adj<0.05), ]
    query <- gost(query = rownames(degtable) %>%
                    unique(), organism = 'hsapiens',
                  correction_method = 'fdr', evcodes = T)
    query.res <- query[['result']]
    query.res
  }) -> query.res.list
  names(query.res.list) <- paste(cell_annotation$X1, 
                                 cell_annotation$X2, 
                                 sep = '_')
  query.res.list %>%
    openxlsx::write.xlsx(file = sprintf('DEG.20240125/07_enrichment_%s_%s_11cluster_down_lfc025.xlsx', 
                                        g1, g2))
}

# plot the enrichment result for every comparison and cluster, down and up abs(LFC)>0.25
for (i in 1:6){
  g2 <- group.combn[1, i]
  g1 <- group.combn[2, i]
  print(paste(g1, g2, sep = ' to '))
  for (ud in c('up', 'down')){
    enrichment_file <- 
      sprintf('DEG.20240125/07_enrichment_%s_%s_11cluster_%s_lfc025.xlsx',
              g1, g2, ud)
    readxl::excel_sheets(enrichment_file) %>%
      lapply(., function(s){
        openxlsx::read.xlsx(enrichment_file, 
                            sheet = s) 
      }) -> enrichment_list
    readxl::excel_sheets(enrichment_file) -> names(enrichment_list)
    if (ud=='up'){
      hc <- "#a50f15"
      lc <- "#fc9272"
    }else if (ud=='down'){
      hc <- "#132b43"
      lc <- "#56b1f7"
    }
    lapply(readxl::excel_sheets(enrichment_file), function(c){
      en <- enrichment_list[[c]]
      if (is.null(en)){
        print('empty data!')
        return(NULL)
      }
      en %>%
        filter(source == 'GO:BP') %>%
        slice_min(order_by = p_value, n = 15) %>%
        head(15) %>%
        #filter(term_name %in% TSC_ASD_pathway$V1) %>%
        arrange(desc(p_value)) %>%
        mutate(term_name = factor(term_name, levels = term_name)) %>%
        ggplot(aes(term_name, -log10(p_value))) +
        geom_bar(aes(fill = -log10(p_value)), stat = "identity", show.legend = F) +
        scale_fill_gradient(high = hc, low = lc) +
        ylab("-log10(FDR)") +
        xlab("") +
        coord_flip() +
        theme_classic(base_size = 12, base_family = "Arial") +
        theme(axis.title = element_text(color = "black"),
              axis.text = element_text(color = "black", size = 12))
      
      export::graph2pdf(last_plot(), font = "Arial",
                        sprintf("DEG.20240125/07_enrichment_%s_%s_%s_%s_lfc025.pdf",
                                g1, g2, c, ud),
                        width = 7, height = 5)
    })
  }
}

# the go enrichment for DEG up psz-6 tsc-tuber
deg.psz7tuber.psz6.11cluster <- readRDS(
  file = 'DEG.20240125/07_DEGs_PSZ-7-tuber_PSZ-6_11cluster.rds')
# lfc >0 
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  print(s)
  degtable <- deg.psz7tuber.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC>0 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  if (FALSE){
  query.res %>%
    filter(source == "GO:BP") %>%
    slice_min(order_by = p_value, n = 20) %>%
    arrange(desc(intersection_size / term_size)) %>%
    mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
    ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
    geom_point(aes(size = intersection_size)) +
    scale_color_gsea() +
    ggtitle(label = "Pathway Enrichment Analysis",
            subtitle = sprintf("Top 20 Pathways")) +
    xlab("gene ratio") +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = "black"))
  #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-tuber_up_go_enrichment.pdf'), 
  #       width = 10, height = 5)
  }
  query.res
}) -> query.res.cellclusters.psz7tuber.psz6.up.lfc0
names(query.res.cellclusters.psz7tuber.psz6.up.lfc0) <- paste(cell_annotation$X1, 
                                                         cell_annotation$X2, 
                                                         sep = '_')
query.res.cellclusters.psz7tuber.psz6.up.lfc0 %>%
  openxlsx::write.xlsx(file = 'DEG.20240125/07_enrichment_PSZ-7-tuber_PSZ-6_11cluster_up_lfc0.xlsx',
                       rowNames = T)

# lfc > 0.01
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  print(s)
  degtable <- deg.psz7tuber.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC>0.01 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  if (FALSE){
    query.res %>%
      filter(source == "GO:BP") %>%
      slice_min(order_by = p_value, n = 20) %>%
      arrange(desc(intersection_size / term_size)) %>%
      mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
      ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
      geom_point(aes(size = intersection_size)) +
      scale_color_gsea() +
      ggtitle(label = "Pathway Enrichment Analysis",
              subtitle = sprintf("Top 20 Pathways")) +
      xlab("gene ratio") +
      theme_bw() +
      theme(axis.text = element_text(size = 12, color = "black"))
    #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-tuber_up_go_enrichment.pdf'), 
    #       width = 10, height = 5)
  }
  query.res
}) -> query.res.cellclusters.psz7tuber.psz6.up.lfc001
names(query.res.cellclusters.psz7tuber.psz6.up.lfc001) <- paste(cell_annotation$X1, 
                                                              cell_annotation$X2, 
                                                              sep = '_')
query.res.cellclusters.psz7tuber.psz6.up.lfc001 %>%
  openxlsx::write.xlsx(file = 'DEG.20240125/07_enrichment_PSZ-7-tuber_PSZ-6_11cluster_up_lfc001.xlsx',
                       rowNames = T)

# lfc > 0.25
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  print(s)
  degtable <- deg.psz7tuber.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC>0.25 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  if (FALSE){
  query.res %>%
    filter(source == "GO:BP") %>%
    slice_min(order_by = p_value, n = 20) %>%
    arrange(desc(intersection_size / term_size)) %>%
    mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
    ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
    geom_point(aes(size = intersection_size)) +
    scale_color_gsea() +
    ggtitle(label = "Pathway Enrichment Analysis",
            subtitle = sprintf("Top 20 Pathways")) +
    xlab("gene ratio") +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = "black"))
  #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-tuber_up_go_enrichment.pdf'), 
  #       width = 10, height = 5)
  }
  query.res
}) -> query.res.cellclusters.psz7tuber.psz6.up.lfc025
names(query.res.cellclusters.psz7tuber.psz6.up.lfc025) <- paste(cell_annotation$X1, 
                                                         cell_annotation$X2, 
                                                         sep = '_')
#query.res.cellclusters.psz7tuber.psz6.up.lfc025 %>% 
#  saveRDS(file = '07_cellclustersDEG_PSZ-7-tuber_PSZ-6_up_goquery.rds',
#                  compress = F)
query.res.cellclusters.psz7tuber.psz6.up.lfc025 %>%
  openxlsx::write.xlsx(file = 'DEG.20240125/07_enrichment_PSZ-7-tuber_PSZ-6_11cluster_up_lfc025.xlsx', 
                       rowNames = T)

# the go enrichment for DEG down PSZ-6 TSC-tuber
# lfc < 0
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  print(s)
  degtable <- deg.psz7tuber.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC < 0 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  if (FALSE){
    query.res %>%
      filter(source == "GO:BP") %>%
      slice_min(order_by = p_value, n = 20) %>%
      arrange(desc(intersection_size / term_size)) %>%
      mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
      ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
      geom_point(aes(size = intersection_size)) +
      scale_color_gsea() +
      ggtitle(label = "Pathway Enrichment Analysis",
              subtitle = sprintf("Top 20 Pathways")) +
      xlab("gene ratio") +
      theme_bw() +
      theme(axis.text = element_text(size = 12, color = "black"))
    #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-tuber_up_go_enrichment.pdf'), 
    #       width = 10, height = 5)
  }
  query.res
}) -> query.res.cellclusters.psz7tuber.psz6.down.lfc0
names(query.res.cellclusters.psz7tuber.psz6.down.lfc0) <- paste(cell_annotation$X1, 
                                                                cell_annotation$X2, 
                                                                sep = '_')
query.res.cellclusters.psz7tuber.psz6.down.lfc0 %>%
  openxlsx::write.xlsx(file = 'DEG.20240125/07_enrichment_PSZ-7-tuber_PSZ-6_11cluster_down_lfc0.xlsx',
                       rowNames = T)
# lfc < -0.01
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  print(s)
  degtable <- deg.psz7tuber.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC < -0.01 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  if (FALSE){
    query.res %>%
      filter(source == "GO:BP") %>%
      slice_min(order_by = p_value, n = 20) %>%
      arrange(desc(intersection_size / term_size)) %>%
      mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
      ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
      geom_point(aes(size = intersection_size)) +
      scale_color_gsea() +
      ggtitle(label = "Pathway Enrichment Analysis",
              subtitle = sprintf("Top 20 Pathways")) +
      xlab("gene ratio") +
      theme_bw() +
      theme(axis.text = element_text(size = 12, color = "black"))
    #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-tuber_up_go_enrichment.pdf'), 
    #       width = 10, height = 5)
  }
  query.res
}) -> query.res.cellclusters.psz7tuber.psz6.down.lfc001
names(query.res.cellclusters.psz7tuber.psz6.down.lfc001) <- paste(cell_annotation$X1, 
                                                                cell_annotation$X2, 
                                                                sep = '_')
query.res.cellclusters.psz7tuber.psz6.down.lfc001 %>%
  openxlsx::write.xlsx(file = 'DEG.20240125/07_enrichment_PSZ-7-tuber_PSZ-6_11cluster_down_lfc001.xlsx',
                       rowNames = T)
# lfc < -0.25
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  print(s)
  degtable <- deg.psz7tuber.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC< -0.25 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  if (FALSE){
  query.res %>%
    filter(source == "GO:BP") %>%
    slice_min(order_by = p_value, n = 20) %>%
    arrange(desc(intersection_size / term_size)) %>%
    mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
    ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
    geom_point(aes(size = intersection_size)) +
    scale_color_gsea() +
    ggtitle(label = "Pathway Enrichment Analysis",
            subtitle = sprintf("Top 20 Pathways")) +
    xlab("gene ratio") +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = "black"))
  #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-tuber_down_go_enrichment.pdf'), 
  #       width = 10, height = 5)
  }
  query.res
}) -> query.res.cellclusters.psz7tuber.psz6.down.lfc025
names(query.res.cellclusters.psz7tuber.psz6.down.lfc025) <- paste(cell_annotation$X1, 
                                                         cell_annotation$X2, 
                                                         sep = '_')
#query.res.cellclusters.psz7tuber.psz6.down.lfc025 %>% 
#  saveRDS(file = '07_cellclustersDEG_PSZ-7-tuber_PSZ-6_down_goquery.rds',
#          compress = F)
query.res.cellclusters.psz7tuber.psz6.down.lfc025 %>%
  openxlsx::write.xlsx(file = 'DEG.20240125/07_enrichment_PSZ-7-tuber_PSZ-6_11cluster_down_lfc025.xlsx',
                       rowNames = T)

# the go enrichment for DEG up PSZ-6 PSZ-7-outside
#deg.psz6.psz7outside.11cluster
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  degtable <- deg.psz7outside.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC>0.25 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  query.res %>%
    filter(source == "GO:BP") %>%
    slice_min(order_by = p_value, n = 20) %>%
    arrange(desc(intersection_size / term_size)) %>%
    mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
    ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
    geom_point(aes(size = intersection_size)) +
    scale_color_gsea() +
    ggtitle(label = "Pathway Enrichment Analysis",
            subtitle = sprintf("Top 20 Pathways")) +
    xlab("gene ratio") +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = "black"))
  #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-outside_up_go_enrichment.pdf'), 
  #       width = 10, height = 5)
  query.res
}) -> query.res.cellclusters.psz7outside.psz6.up
names(query.res.cellclusters.psz7outside.psz6.up) <- paste(cell_annotation$X1, 
                                                         cell_annotation$X2, 
                                                         sep = '_')
query.res.cellclusters.psz7outside.psz6.up %>% 
  saveRDS(file = '07_cellclustersDEG_PSZ-7-outside_PSZ-6_up_goquery.rds',
          compress = F)
query.res.cellclusters.psz7outside.psz6.up %>%
  openxlsx::write.xlsx(file = '07_cellclustersDEG_PSZ-7-outside_PSZ-6_up_goquery.xlsx')

# the go enrichment for DEG down PSZ-6 PSZ-7-outside
#deg.psz6.psz7outside.11cluster
lapply(seq(0,10), function(i){
  s <- paste(as.character(i), cell_annotation$X2[i+1], sep = '_')
  degtable <- deg.psz7outside.psz6.11cluster[[s]]
  degtable <- degtable[which(degtable$avg_log2FC< -0.25 & degtable$p_val_adj<0.05), ]
  query <- gost(query = rownames(degtable) 
                %>% unique(), organism = "hsapiens", 
                correction_method = "fdr", evcodes = T)
  query.res <- query[["result"]]
  query.res %>%
    filter(source == "GO:BP") %>%
    slice_min(order_by = p_value, n = 20) %>%
    arrange(desc(intersection_size / term_size)) %>%
    mutate(term_name = factor(term_name, levels = rev(term_name))) %>%
    ggplot(aes(intersection_size / term_size, term_name, color = -log10(p_value))) +
    geom_point(aes(size = intersection_size)) +
    scale_color_gsea() +
    ggtitle(label = "Pathway Enrichment Analysis",
            subtitle = sprintf("Top 20 Pathways")) +
    xlab("gene ratio") +
    theme_bw() +
    theme(axis.text = element_text(size = 12, color = "black"))
  #ggsave(paste0('07_', as.character(i), '_cellClusterDEG_PSZ-6_TSC-outside_down_go_enrichment.pdf'), 
  #       width = 10, height = 5)
  query.res
}) -> query.res.cellclusters.psz7outside.psz6.down
names(query.res.cellclusters.psz7outside.psz6.down) <- paste(cell_annotation$X1, 
                                                           cell_annotation$X2, 
                                                           sep = '_')
query.res.cellclusters.psz7outside.psz6.down %>% 
  saveRDS(file = '07_cellclustersDEG_PSZ-7-outside_PSZ-6_down_goquery.rds',
          compress = F)
query.res.cellclusters.psz7outside.psz6.down %>%
  openxlsx::write.xlsx(file = '07_cellclustersDEG_PSZ-7-outside_PSZ-6_down_goquery.xlsx')

