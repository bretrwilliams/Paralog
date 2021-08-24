require(dplyr)
require(magrittr)
require(tidyr)
require(tibble)
require(purrr)
require(taigr)
require(celllinemapr)
require(ggrepel)
require(plotly)
require(RColorBrewer)
options(stringsAsFactors = F)
#===== Heatmap of MAPK pathway =====
# Visualizations for % of hits in each cell line, consistent with Manhattan
fdr_Scores <- readRDS("Data/LatestScores.rds")$fdr_sensitive_lethality
fdr_Scores%<>%replace(is.na(.), 1)%>%.[grep("AAVS1|TRIM",rownames(fdr_Scores),invert = T,value=T),]
require(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene.data <- getBM(attributes=c('hgnc_symbol','go_id'),filters = 'go', values = 'GO:0000165', mart = ensembl)
MAPK_genes <- gene.data%>%
  filter(go_id=="GO:0000165")%>%
  distinct()%>%
  pull(hgnc_symbol)
gsea_MAPK <- read.delim("Data/MAPK_GSEA_geneset.txt")[-1,1]
MAPK_allG <- union(MAPK_genes,gsea_MAPK)
MAPK_pair <- data.frame(pair =rownames(fdr_Scores))%>%
  separate(col = "pair",into = c("gene1","gene2"),sep = ";",remove = F)%>%
  filter(gene1 %in% MAPK_allG & gene2 %in% MAPK_allG)%>%
  pull(pair)
columnOrders <- c("PK1_PANCREAS",
                  "PATU8988S_PANCREAS",
                  "A549_LUNG",
                  "MELJUSO_SKIN",
                  "HS944T_SKIN",
                  "IPC298_SKIN",
                  "HS936T_SKIN",
                  "MEL202_UVEA",
                  "MEWO_SKIN",
                  "HSC5_SKIN",
                  "GI1_CENTRAL_NERVOUS_SYSTEM")
heatmapD <- fdr_Scores[MAPK_pair,rev(columnOrders)]
colnames(heatmapD)%<>%mgsub::mgsub(.,c("\\_.*","PATU8988S"),c("","PATU"))
save_pheatmap_pdf <- function(x, filename, width=4.5, height=7.5) {
  pdf(filename, width = width, height = height)
  grid::grid.draw(x$gtable)
  dev.off()
}
pheatmap::pheatmap(t(heatmapD),
                   color = colorRampPalette(
                     brewer.pal(n = 9, name ="GnBu"))(100)%>%rev(),
                   cluster_rows = F,
                   cellwidth = 5.5,
                   cellheight = 18,
                   fontsize_row = 5,
                   fontsize_col = 4,
                   clustering_method = "single",
                   treeheight_col = 0,
                   show_colnames= TRUE,
                   angle_col = "90",
                   main = "MAPK-pathway related paralog pairs")%>%
  save_pheatmap_pdf(., paste0("Fig S4a Heatmap for MAPK pathway/FigS4a_heatmap_MAPK_fdr_SingleClustered_",Sys.Date(),
                              ".pdf"),height = 6,width = 20)