require(dplyr)
require(pheatmap)
Model <- readRDS("Data/LatestModel.rds")
Score <- readRDS("Data/LatestScores.rds")
senscore <- Score$sensitive_lethality
lfc <- Model$Input$LFC
lfcMean <- lfc%>%
  as.data.frame()%>%
  tibble::rownames_to_column("rowname")%>%
  inner_join(
    Model$Input$guide.pair.annot%>%
      filter(grepl("PRKC",Aureus_gene)&
               grepl("PRKC",Pyogenes_gene)),by="rowname")%>%
  mutate(pair = purrr::pmap_chr(list(Aureus_gene,Pyogenes_gene),
                                ~paste(sort(c(...)),collapse = ";")))%>%
  select(-rowname,-Aureus_gene,-Pyogenes_gene)%>%
  group_by(pair)%>%
  summarise_all(mean)%>%
  as.data.frame()
lfcMeanMat <- lfcMean%>%select(-pair)%>%as.matrix
rownames(lfcMeanMat) <- lfcMean$pair
### Heatmap ####
fdrD <- -log10(Score$fdr_sensitive_lethality[grep("PRKC[A-Z];PRKC[A-Z]",rownames(senscore)),])
fdrD[which(fdrD > 2, arr.ind=TRUE)] <- 2
heatmapFDR<- pheatmap(fdrD,
                       cellwidth = 15,
                       cellheight = 10,
                       fontsize = 6,
                      angle_col = "45")
fdrrow <- rownames(fdrD[heatmapFDR$tree_row[["order"]],])
fdrcol <- colnames(fdrD[,heatmapFDR$tree_col[["order"]]])
heatmapLFC <- pheatmap(lfcMeanMat[fdrrow,fdrcol],
                       cluster_cols=F,
                       cluster_rows = F,
                       cellwidth = 15,
                       cellheight =10,
                       fontsize = 6,angle_col = "45",
                       color=colorRampPalette(c("#D73027", "#FC8D59", "#FEE090",
                                                "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4"))(100))
save_pheatmap_pdf <- function(x, filename, width=5,height=7.5) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmapFDR, "Fig S3f-g PRKC Heatmap/FigS3f_Heatmap_FDR.pdf")
save_pheatmap_pdf(heatmapLFC, "Fig S3f-g PRKC Heatmap/FigS3g_heatmap_LFC_nocluster.pdf")
