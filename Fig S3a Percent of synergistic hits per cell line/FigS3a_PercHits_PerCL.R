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
#===== Percentage of synergistic hits per cl =====
# Visualizations for % of hits in each cell line, consistent with Manhattan
fdr_Scores <- readRDS("Data/LatestScores.rds")$fdr_sensitive_lethality
fdr_Scores%<>%replace(is.na(.), 1)%>%.[grep("AAVS1|TRIM",rownames(fdr_Scores),invert = T,value=T),]
hit_percent <- apply(fdr_Scores,2,function(s) sum(s<0.05)*100/nrow(fdr_Scores))%>%
  as.data.frame()%>%
  set_colnames("Percent")%>%
  rownames_to_column("Cell.line")%>%
  mutate(Cell.line = gsub("\\_.*","",Cell.line))
hit_percent$Cell.line%<>%factor(.,levels = hit_percent%>%arrange(Percent)%>%pull(Cell.line),ordered = T)
ggplot(hit_percent, aes(x=Cell.line, y=Percent)) +
  geom_bar(stat="identity",fill="black")+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text.y = element_text(color='black'),
        axis.text.x = element_text(color='black',angle=90,hjust = 1,vjust=0.5),
        axis.title.x = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  ylim(c(0,7))+
  labs(y="Percent of synergistic hits")+
  ggsave(filename = "Fig S3a Percent of synergistic hits per cell line/FigS3a_Percent_synergistichits_fdr0.05_trimRemoved.pdf",width=8,height=8)
