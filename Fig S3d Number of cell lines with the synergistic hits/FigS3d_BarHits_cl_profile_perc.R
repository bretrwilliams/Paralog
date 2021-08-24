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
#===== Synergistic pairs in each cell line, binned by shared/unique =====
fdr_Scores <- readRDS("Data/LatestScores.rds")$fdr_sensitive_lethality
fdr_Scores%<>%replace(is.na(.), 1)%>%.[grep("AAVS1|TRIM",rownames(fdr_Scores),invert = T,value=T),]
hits_pool <- fdr_Scores%>%
  as.data.frame()%>%
  rownames_to_column("pair")%>%
  filter(!grepl("AAVS1|TRIM",pair))%>%
  gather(Cell.line,fdr,-pair)%>%
  filter(fdr < 0.05)%>%
  mutate(Cell.line = gsub("\\_.*","",Cell.line))
hits_profile <- hits_pool%>%
  group_by(pair)%>%
  summarize(n=n())%>%
  mutate(hits_in = case_when(n==1~"1 cell line",
                             n<=3~"2-3 cell lines",
                             n>3 & n<=6~"4-6 cell lines",
                             n>6~"More than 7 cell lines")%>%
           factor(.,levels=rev(c("More than 7 cell lines","4-6 cell lines","2-3 cell lines","1 cell line"))))%>%
  select(-n)%>%
  distinct()%>%
  group_by(hits_in)%>%
  summarize(n=n())%>%
  mutate(perc = round(n*100/sum(n),2))
ggplot(hits_profile,aes(x=perc, y=hits_in)) + 
  geom_bar(position="stack", stat="identity",fill="#faaf3f")+
  scale_x_continuous(limits = c(0,100),breaks = seq(0,100,25))+
  theme(
    plot.title = element_text(color='black', hjust = 0.5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(color='black',size=rel(1.3)),
    legend.position = "right"
  )+labs(y=NULL,x="Percent of synergistic pairs")
ggsave(filename = "Fig S3d Number of cell lines with the synergistic hits/FigS3d_Hits_cl_profile_perc.pdf",width=10,height=4)

