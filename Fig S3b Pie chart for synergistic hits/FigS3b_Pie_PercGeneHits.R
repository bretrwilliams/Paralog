require(dplyr)
require(magrittr)
require(tidyr)
require(tibble)
require(purrr)
require(taigr)
require(celllinemapr)
require(ggrepel)
require(plotly)
options(stringsAsFactors = F)
#===== Pie chart: # of genes detected in DKO ======
Model <- readRDS("Data/LatestModel.rds")
fdr_Scores <- readRDS("Data/LatestScores.rds")$fdr_sensitive_lethality
fdr_Scores%<>%replace(is.na(.), 1)%>%.[grep("AAVS1|TRIM",rownames(fdr_Scores),invert = T,value=T),]
hits_pool <- fdr_Scores%>%
  as.data.frame()%>%
  rownames_to_column("pair")%>%
  filter(!grepl("AAVS1|TRIM",pair))%>%
  gather(Cell.line,fdr,-pair)%>%
  filter(fdr < 0.05)%>%
  mutate(Cell.line = gsub("\\_.*","",Cell.line))
genes_hits <- lapply(unique(hits_pool$pair),strsplit,split=";")%>%unlist()%>%unique()
all_genes <- Model$Input$guide.pair.annot%>%
  mutate(pair = pmap_chr(list(Aureus_gene,Pyogenes_gene),~paste0(sort(c(...)),collapse = ";")))%>%
  filter(grepl("AAVS1",pair))%>%
  mutate(gene = gsub("AAVS1;|;AAVS1","",pair))%>%
  pull(gene)%>%unique()
p_hit <- round(length(genes_hits)*100/length(all_genes),2)
require(scales)
data.frame(
  group = c("Hits", "Non-hits"),
  value = c(p_hit,100-p_hit))%>%
  ggplot(., aes(x="", y=value, fill=group))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    scale_fill_manual(values = c("#faaf3f","#1b75bb"))+
    geom_text(aes(y = c(90,45), 
                label = c("Synergy \n 22% \n (n = 731)",
                          "No synergy \n 78% \n (n = 2553)")),
              size=5,family="Helvetica")+
  theme(
    plot.title =  element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")+
  ggsave(filename = "Fig S3b Pie chart for synergistic hits/FigS3b_PieChart_Hits_genes.pdf",width=6,height=6)
