#set the work path
setwd('E:/shiyiLab/WeiboNiu_data/singleCellAnalysis20231204/08.post.annotation.plot/')
#
library(Seurat)
library(tidyverse)
library(ggpubr)
library(export)

# using resolution 0.1
cluster_annotation <-
  openxlsx::read.xlsx("group.rename/cell_type_annotation_resolution_0.1 update20240109.xlsx",
                      sheet = 1, colNames = F)
cluster_annotation$X1 <- as.factor(cluster_annotation$X1)

merged <- readRDS("../05.integrate.analysis/05_merged_with_clustering.rds")
merged$group.rename <- ifelse(merged$group=='PSZ-6', 'CTRL', 
                              ifelse(merged$group=='PSZ-7-tuber', 'TSC-tuber',
                                     ifelse(merged$group=='PSZ-7-edge', 'TSC-edge', 
                                            'TSC-outside')))
merged$group.rename <- factor(merged$group.rename, ordered = TRUE, 
                              levels = c('CTRL', 'TSC-tuber', 'TSC-edge', 'TSC-outside'))

merged$seurat_clusters <- merged$integrated_snn_res.0.1
Idents(merged) <- "seurat_clusters"
merged$seurat_clusters %>% table()

##### UMAP
DimPlot(merged, split.by = "group.rename", label = T)
  export::graph2pdf(last_plot(),
                  "group.rename/umap.resolution.0.1.pdf",
                  width = 9, height = 5)

# add figures on 20241031
  DimPlot(merged, label = T)
  export::graph2pdf(last_plot(), 
                    'add_figures_20241031/allcell.umap.resolution.0.1.pdf', 
                    width = 9, height = 5)
  export::graph2pdf(last_plot(), 
                    'add_figures_20241031/allcell.umap.resolution.0.1.w4.5.pdf', 
                    width = 4.5, height = 5)
  export::graph2pdf(last_plot(), 
                    'add_figures_20241031/allcell.umap.resolution.0.1.w3.pdf', 
                    width = 3, height = 5)
  
  merged$group.rename %>% table()
  merged %>% subset(group.rename == 'CTRL') %>%
    DimPlot(label = T)
  export::graph2pdf(last_plot(), 
                    'add_figures_20241031/ctrl.umap.resolution.0.1.pdf', 
                    width = 9, height = 5)

######################## marker gene scatter plot
cluster_annotation$X3 %>% paste(collapse = ",") %>%
  str_remove_all(" ") %>%
  str_split(",") %>% .[[1]] -> marker_genes

marker_genes <- marker_genes[!marker_genes %in% c("TBR1", "SATB2")]

marker_genes <- "NEFM"
DefaultAssay(merged) <- "RNA"
expr <- FetchData(merged, c(marker_genes, "seurat_clusters", "group.rename"))
# TBR1,  CUX1,  POU3F2,  SATB2,  LHX1
expr %>%
  #mutate(group = ifelse(group == "CON", "CTRL", "TSC")) %>%
  left_join(cluster_annotation, by = c(seurat_clusters = "X1")) %>%
  group_by(X2, group.rename) %>%
  mutate(expr = NEFM) %>%
  ggbarplot(x = "X2", y = "expr", add = "mean_se", fill = "group.rename",
            position = position_dodge(), width = 0.6,
            xlab = "", ylab = "Average Expression") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#00BFFF", "#FFD700")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        axis.text = element_text(color = "black")) +
  stat_compare_means(aes(group = group.rename), label = "p.signif",
                     label.y = 1.5)
export::graph2pdf(last_plot(),
      sprintf("%s_average_expr.pdf",marker_genes),
      width = 8, height = 6)

DotPlot(object = merged, features = unique(marker_genes)) + 
  scale_y_discrete(labels = cluster_annotation$X2) +
  theme_bw(base_size = 16, base_family = "Arial") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"))
export::graph2pdf(last_plot(),
                  "update20240127/maker_gene_dotplot.pdf",
                  width = 15, height = 7)

###################### feature plot of marker genes
# Creating a character vector
my_vector <- c("RELN","PTGDS","S100B","NNAT","CHL1","CNTN4","RBFOX1","PTN","CNTN6","NEFL","CHCHD2",
               "ANXA2","SLC1A3","GRIN2B","NEUROD2","GAD2","GABBR2")
my_vector <- c("TGFB2", "NLGN1", "TSLP", "DKK1","BMP4",
                 "NFIA", "SOX9", "RUNX2", "STAT3", "RORB", "NEFL", "NEFM")
openxlsx::read.xlsx('snRNAseq_marker_genes.xlsx', colNames = FALSE)$X1 -> my_vector
openxlsx::read.xlsx('group.rename/snRNAseq_marker_genes.xlsx', 
                    sheet = 1, colNames = FALSE)$X1 -> my_vector_1
openxlsx::read.xlsx('group.rename/snRNAseq_marker_genes.xlsx', 
                    sheet = 2, colNames = FALSE)$X1 -> my_vector_2

DefaultAssay(merged) <- "RNA"
lapply(my_vector_1, function(gene){
  print(gene)
  FeaturePlot(merged, gene, split.by = "group.rename")
  export::graph2pdf(last_plot(),
                    sprintf("group.rename/snRNAseq_marker_genes_sheet1/feature_plot_%s.pdf",gene),
                    width = 12, height = 5)
})

lapply(my_vector_2, function(gene){
  print(gene)
  FeaturePlot(merged, gene, split.by = "group.rename")
  export::graph2pdf(last_plot(),
                    sprintf("group.rename/snRNAseq_marker_genes_sheet2/feature_plot_%s.pdf",gene),
                    width = 12, height = 5)
})

my_vector_3 <- c("CCL3", 'SYT1', 'SLC1A2')
lapply(my_vector_3, function(gene){
  print(gene)
  FeaturePlot(merged, gene, split.by = "group.rename")
  export::graph2pdf(last_plot(), 
                    sprintf("group.rename/snRNAseq_marker_genes_sheet2/feature_plot_%s.pdf", gene),
                    width = 12, height = 5)
})

# feature plot of genes in new table;make the point on the topside
openxlsx::read.xlsx('markers20240928/Microglia_markers.xlsx', colNames = F)$X1 -> my_vector
DefaultAssay(merged) <- 'RNA'
lapply(my_vector, function(gene){
  print(gene)
  FeaturePlot(merged, gene, split.by = 'group.rename', order = T)
  export::graph2pdf(last_plot(), 
                    sprintf('markers20240928/feature_plot_%s_newcols.pdf', gene), 
                    width = 12, height = 5)
})
#FeaturePlot(merged, gene, split.by = 'group.rename', cols = c('grey', 'red'))

# add the feature plot in vlnplot for LPL gene in microgila
merged %>% subset(subset= seurat_clusters=='2') -> merged_microgila
merged_microgila@meta.data %>% dim()
DefaultAssay(merged_microgila) <- "RNA"
merged_microgila$LPL <- FetchData(merged_microgila, 'LPL')
#merged_microgila$LPL %>% log10() -> merged_microgila$lpllog
merged_microgila@meta.data %>% 
  ggplot(aes(x = group.rename, y = LPL, fill = group.rename)) +
  geom_bar(stat = 'summary', fun = 'mean', 
           fill = c('#00bfc1', '#f9766e', '#f9766e', '#f9766e')) +
  geom_jitter(size = 0.8) +
  geom_signif(comparisons = list(c('CTRL', 'TSC-tuber'), 
                                 c('CTRL', 'TSC-edge'), 
                                 c('CTRL', 'TSC-outside')), 
              test = 't.test',  
              #tip_length = c(0,0,0,0,0,0),
              y_position = c(4.4, 4.8, 5.2), 
              annotations = c('****')) +
  theme_classic2() +
  ylab('LPL Normalized Expression') + 
  xlab('') + NoLegend() -> p.lpl.jitter.withbar
ggsave('add_figures_20241031/LPL.jitter.withbar_bluevsred.pdf', plot = p.lpl.jitter.withbar, 
       bg = 'transparent', width = 9, height = 5)

merged_microgila@meta.data %>% 
  ggplot(aes(x = group.rename, y = LPL, fill = group.rename)) +
  geom_bar(stat = 'summary', fun = 'mean', )+
           #fill = c('#00bfc1', '#f9766e', '#f9766e', '#f9766e')) +
  geom_jitter(size = 0.8) +
  geom_signif(comparisons = list(c('CTRL', 'TSC-tuber'), 
                                 c('CTRL', 'TSC-edge'), 
                                 c('CTRL', 'TSC-outside')), 
              test = 't.test',  
              #tip_length = c(0,0,0,0,0,0),
              y_position = c(4.4, 4.8, 5.2), 
              annotations = c('****')) +
  theme_classic2() +
  ylab('LPL Normalized Expression') + 
  xlab('') + NoLegend() -> p.lpl.jitter.withbar
ggsave('add_figures_20241031/LPL.jitter.withbar_autocolor.pdf', plot = p.lpl.jitter.withbar, 
       bg = 'transparent', width = 9, height = 5)

merged_microgila@meta.data %>%
  ggplot(aes(x = group.rename, y = LPL, fill = group.rename, 
             group = group.rename)) +
  geom_jitter(aes(color = group.rename), size = 0.8) +
  scale_color_manual(values = c('#00bfc1', '#f9766e', '#f9766e', '#f9766e')) +
  geom_signif(comparisons = list(c('CTRL', 'TSC-tuber'), 
                                 c('CTRL', 'TSC-edge'), 
                                 c('CTRL', 'TSC-outside')), 
              test = 't.test',  
              #tip_length = c(0,0,0,0,0,0),
              y_position = c(4.4, 4.8, 5.2), 
              annotations = c('****')) +
  theme_classic2() +
  ylab('LPL Normalized Expression') + 
  xlab('') + NoLegend() -> p.lpl.jitter.nobar
ggsave('add_figures_20241031/LPL.jitter.nobar_bluevsred.pdf', plot = p.lpl.jitter.nobar, 
       bg = 'transparent', width = 9, height = 5)

merged_microgila@meta.data %>%
  ggplot(aes(x = group.rename, y = LPL, fill = group.rename, 
             group = group.rename)) +
  geom_jitter(aes(color = group.rename), size = 0.8) +
  #scale_color_manual(values = c('#00bfc1', '#f9766e', '#f9766e', '#f9766e')) +
  geom_signif(comparisons = list(c('CTRL', 'TSC-tuber'), 
                                 c('CTRL', 'TSC-edge'), 
                                 c('CTRL', 'TSC-outside')), 
              test = 't.test',  
              #tip_length = c(0,0,0,0,0,0),
              y_position = c(4.4, 4.8, 5.2), 
              annotations = c('****')) +
  theme_classic2() +
  ylab('LPL Normalized Expression') + 
  xlab('') + NoLegend() -> p.lpl.jitter.nobar
ggsave('add_figures_20241031/LPL.jitter.nobar_autocolor.pdf', plot = p.lpl.jitter.nobar, 
       bg = 'transparent', width = 9, height = 5)

VlnPlot(merged, features = c('LPL'), group.by = 'group.rename', 
        pt.size = 0) + NoLegend() +
  geom_signif(comparisons = list(c('CTRL', 'TSC-tuber'), 
                                 c('CTRL', 'TSC-edge'), 
                                 c('CTRL', 'TSC-outside')), 
              y_position = c(4.4, 5.0, 5.6)) +
  ylim(-2, 6)
export::graph2pdf(last_plot(), 
                  'add_figures_20241031/LPL.vlnplot.pvalue.pdf', 
                  width = 9, height = 5)

####################### cluster DEG
merged$deg_group <- paste0(merged$group, "_", merged$seurat_clusters)
DefaultAssay(merged) <- "RNA"
Idents(merged) <- "deg_group"
lapply(unique(merged$seurat_clusters), function(cl){
  print(cl)
  g1 <- paste0("CON_", cl)
  g2 <- paste0("TSC_", cl)
  deg <- FindMarkers(merged, ident.1 = g1, ident.2 = g2)
  deg$cluster <- cl
  deg$gene <- rownames(deg)
  deg
}) -> deg_list
names(deg_list) <- paste0("cluster ", unique(merged$seurat_clusters))
openxlsx::write.xlsx(deg_list, 
                    "cluster_DEGs.xlsx")

lapply(deg_list, function(deg){
  deg %>%
    filter(p_val_adj < 0.05)
}) -> deg_list

############### enrichment barplot
list.files("TSC-organoid/results/SingleCell/post_annotation/cluster_DEG_enrichment/") %>%
  grep("WN.xlsx$", ., value = T) %>%
  grep("up", ., value = T) %>%
  lapply(., function(x){
    sprintf("TSC-organoid/results/SingleCell/post_annotation/cluster_DEG_enrichment/%s", x) %>%
      openxlsx::read.xlsx(., sheet = 2) -> dat
    dat %>%
      arrange(desc(p_value)) %>%
      mutate(term_name = factor(term_name, levels = term_name)) %>%
      ggplot(aes(term_name, -log10(p_value))) +
      geom_bar(aes(fill = -log10(p_value)), stat = "identity", show.legend = F) +
      # scale_fill_gradient(high = "#a50f15", low = "#fc9272") +
      scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
      ylab("-log10(FDR)") +
      xlab("") +
      coord_flip() +
      theme_classic(base_size = 14) +
      theme(axis.title = element_text(color = "black"),
            axis.text = element_text(color = "black"))
    
    export::graph2pdf(last_plot(),
                      sprintf("TSC-organoid/results/SingleCell/post_annotation/cluster_DEG_enrichment/%s.pdf", x),
                      width = 8, height = 5, font = "Arial")
  })

############# number of DEG bar plot
# old code
readxl::excel_sheets(
  "../07.cell.cluster.DEGs/DEG.PSZ-7-tuber.PSZ-6.11clusters.xlsx"
) %>% lapply(., function(x){
  openxlsx::read.xlsx("../07.cell.cluster.DEGs/DEG.PSZ-7-tuber.PSZ-6.11clusters.xlsx",
                      sheet = x)
}) -> cluster_degs
readxl::excel_sheets(
  "../07.cell.cluster.DEGs/DEG.PSZ-7-tuber.PSZ-6.11clusters.xlsx"
) -> names(cluster_degs)
# update 20240127
## plot the barplot for 6 comparison
g1 <- 'PSZ-7-outside'
g2 <- 'PSZ-6'
readxl::excel_sheets(
  sprintf(
    '../07.cell.cluster.DEGs/DEG.20240125/07_DEGs_%s_%s_11cluster_lfc025.xlsx',
    g1, g2)
) %>% lapply(., function(x){
  openxlsx::read.xlsx(
    sprintf(
      '../07.cell.cluster.DEGs/DEG.20240125/07_DEGs_%s_%s_11cluster_lfc025.xlsx',
      g1, g2), sheet = x)
}) -> cluster_degs
readxl::excel_sheets(
  sprintf(
    '../07.cell.cluster.DEGs/DEG.20240125/07_DEGs_%s_%s_11cluster_lfc025.xlsx', 
    g1, g2)
) -> names(cluster_degs)

lapply(cluster_degs, function(x){
  sum(x$p_val_adj < 0.05 & x$avg_log2FC > 0.25) -> up
  sum(x$p_val_adj < 0.05 & x$avg_log2FC < -0.25) -> down
  c(up, down)
}) %>% Reduce(rbind, .) -> deg_count
as.data.frame(deg_count) -> deg_count
colnames(deg_count) <- c("up", "down")
deg_count$cluster <- names(cluster_degs)
deg_count$deg_count <- deg_count$up + deg_count$down

sapply(deg_count$cluster, function(x){
  str_extract(x, "\\d+") %>% as.numeric()
}) -> deg_count$cluster_numeric
deg_count$C <- paste0("C", deg_count$cluster_numeric)
deg_count %>%
  mutate(label = paste0(deg_count, " (up:", up, ",down:", down, ")")) -> deg_count
deg_count
deg_count %>%
  arrange(desc(cluster_numeric)) %>%
  mutate(C = factor(C, levels = C)) %>%
  ggplot(aes(C, deg_count)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#99d8c9") +
  geom_text(aes(label = label), position = position_nudge(y = 240)) +
  coord_flip() + xlab("") + ylab("Cluster Specific DEGs") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_family = "Arial", base_size = 13) +
  theme(axis.text = element_text(color = "black", size = 13)) +
  ylim(0,2900)

export::graph2pdf(last_plot(),
  sprintf(
    "update20240127/08_cluster_DEGs_count_%s_%s_11cluster_lfc025.pdf",
    g1, g2), width = 14, height = 4
)

########### microglia enrich path way bar plot
# fig 2C, enrich bar plot
TSC_up_enrichment <- openxlsx::read.xlsx('group.rename/07_cellclustersDEG_PSZ-7-tuber_PSZ-6_up_goquery_WN.xlsx', 
                                         sheet = 1)
TSC_down_enrichment <- openxlsx::read.xlsx('group.rename/07_cellclustersDEG_PSZ-7-tuber_PSZ-6_down_goquery_WN.xlsx', 
                                           sheet = 1)
# up enrichment 
TSC_up_enrichment %>%
  #filter(term_name %in% TSC_ASD_pathway$V1) %>%
  arrange(desc(p_value)) %>%
  mutate(term_name = factor(term_name, levels = term_name)) %>%
  ggplot(aes(term_name, -log10(p_value))) +
  geom_bar(aes(fill = -log10(p_value)), stat = "identity", show.legend = F) +
  scale_fill_gradient(high = "#a50f15", low = "#fc9272") +
  ylab("-log10(FDR)") +
  xlab("") +
  coord_flip() +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 12))

# ggsave("RNA-Seq/results/1_TSC_microglia_PBS/plots/Fig2C_TSC_ASD_enrichment.png", width = 7, height = 5)
export::graph2pdf(last_plot(), font = "Arial",
                  "group.rename/Fig1C_TSC_up_enrichment.pdf",
                  width = 7, height = 5)

# down enrichment 
TSC_down_enrichment %>%
  #filter(term_name %in% TSC_ASD_pathway$V1) %>%
  arrange(desc(p_value)) %>%
  mutate(term_name = factor(term_name, levels = term_name)) %>%
  ggplot(aes(term_name, -log10(p_value))) +
  geom_bar(aes(fill = -log10(p_value)), stat = "identity", show.legend = F) +
  scale_fill_gradient(high = "#132b43", low = "#56b1f7") +
  ylab("-log10(FDR)") +
  xlab("") +
  coord_flip() +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 12))

# ggsave("RNA-Seq/results/1_TSC_microglia_PBS/plots/Fig2C_TSC_ASD_enrichment.png", width = 7, height = 5)
export::graph2pdf(last_plot(), font = "Arial",
                  "group.rename/Fig1D_TSC_down_enrichment.pdf",
                  width = 7, height = 5)