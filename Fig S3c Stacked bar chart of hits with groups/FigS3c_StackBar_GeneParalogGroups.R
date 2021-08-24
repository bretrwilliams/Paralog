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
#===== Synergistic pairs binned by protein class and number of paralogs in family =====
paralog_library <- readRDS("Data/annotated_redundome.rds")
fdr_Scores <- readRDS("Data/LatestScores.rds")$fdr_sensitive_lethality
fdr_Scores%<>%replace(is.na(.), 1)%>%.[grep("AAVS1|TRIM",rownames(fdr_Scores),invert = T,value=T),]
hits_pool <- fdr_Scores%>%
  as.data.frame()%>%
  rownames_to_column("pair")%>%
  filter(!grepl("AAVS1|TRIM",pair))%>%
  gather(Cell.line,fdr,-pair)%>%
  filter(fdr < 0.05)%>%
  mutate(Cell.line = gsub("\\_.*","",Cell.line))
familyPerc <- paralog_library%>%
  mutate(pair = pmap_chr(list(Gene_1_HUGO,Gene_2_HUGO),~paste0(sort(c(...)),collapse = ";")))%>%
  mutate(pair = gsub("PRSS46","PRSS46P",pair))%>%
  filter(pair %in% rownames(fdr_Scores))%>%
  select(Gene_Family,pair)%>%
  distinct()%>%
  filter(!grepl("AAVS1|TRIM",pair))%>%
  mutate(Gene_Family =  gsub("propellor","propeller",Gene_Family)%>%
           gsub("Hyrdolase","Hydrolase",.)%>%
           gsub("apoptotic","-apoptotic",.))%>%
  mutate(Gene_Family = case_when(grepl("apoptotic",Gene_Family)~"Apoptotic",
                                 Gene_Family %in% c("B-propeller","Achilles Paralog Dependent")~"Others",
                                 TRUE~Gene_Family))%>%
  filter(Gene_Family!="Others")
familyPerc%<>%
  mutate(hits = case_when(pair %in% hits_pool$pair~"Yes",TRUE~"No"))%>%
  group_by(Gene_Family,hits)%>%
  summarize(n=n())
familyPerc%>%
  inner_join(familyPerc%>%
               group_by(Gene_Family)%>%
               summarize(n_total = sum(n)),by="Gene_Family")%>%
  mutate(perc = round(n*100/n_total,2))%>%
  mutate(Gene_Family = factor(Gene_Family,levels = rev(c("Transferase","Ligase","Hydrolase","Apoptotic"))))%>%
  ggplot(.,aes(fill=hits, y=perc, x=Gene_Family)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("#1b75bb","#faaf3f"))+
  coord_flip()+
  theme(
    plot.title = element_text(color='black', hjust = 0.5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color='black',size=8,family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.position = "none"
  )+labs(x=NULL,y=NULL)+
  ggsave(filename = "Fig S3c Stacked bar chart of hits with groups/FigS3c_Hits_genefamily_perc_annot.pdf",width=8,height=5)
  
require(igraph)
Groups_of_pairs <- data.frame(pair = rownames(fdr_Scores))%>%
  mutate(symbol1 = gsub(";.*","",pair),
         symbol2 = gsub(".*;","",pair))%>%
  select(-pair)%>%
  graph_from_data_frame()%>%
  components()
Paralog_lst <- split(names(Groups_of_pairs$membership),Groups_of_pairs$membership)
Paralog_groups <- lapply(names(Paralog_lst),function(x) dat = 
                           data.frame(memeber = as.character(Paralog_lst[[x]]),stringsAsFactors = F)%>%
                           mutate(group = x))%>%bind_rows()
group_size <- Paralog_groups%>%
  distinct()%>%
  group_by(group)%>%summarise(grp_size=n())
group_info <- data.frame(pair = rownames(fdr_Scores))%>%
  mutate(symbol1 = gsub(";.*","",pair),
         symbol2 = gsub(".*;","",pair))%>%
  inner_join(Paralog_groups,by=c("symbol1" = "memeber"))%>%
  inner_join(group_size,by="group")

group_info%<>%
  mutate(hits = case_when(pair %in% hits_pool$pair~"Yes",TRUE~"No"))%>%
  mutate(grp_size = case_when(grp_size<6~"Less than 6 genes",
                              grp_size>10~"6-10 genes",
                              TRUE~"More than 10 genes")%>%
           factor(.,levels=c("More than 10 genes",
                             "6-10 genes",
                             "Less than 6 genes")))%>%
  group_by(grp_size,hits)%>%
  summarize(n=n())

group_info%>%
  inner_join(group_info%>%
               group_by(grp_size)%>%
               summarize(n_total = sum(n)),by="grp_size")%>%
  mutate(perc = round(n*100/n_total,2))%>%
  ggplot(.,aes(fill=hits, y=perc, x=grp_size)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("#1b75bb","#faaf3f"))+
  coord_flip()+
  theme(
    plot.title = element_text(color='black', hjust = 0.5,family = "Helvetica"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color='black',size=8,family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.position = "none"
  )+labs(x=NULL,y=NULL)
ggsave(filename = "Fig S3c Stacked bar chart of hits with groups/FigS3c_Hits_grpsize_perc_annot.pdf",width=8,height=5)
