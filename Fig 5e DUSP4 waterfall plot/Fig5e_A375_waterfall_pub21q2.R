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
dname <- 'public-21q2-110d'
dver <- 12
ceres <- load.from.taiga(data.name=dname, data.version=dver, data.file='Achilles_gene_effect')
sample.info <- load.from.taiga(data.name=dname, data.version=dver,data.file='sample_info')
mut <- load.from.taiga(data.name=dname, data.version=dver, data.file='CCLE_mutations')
mut_annot <- mut%>%
  filter(Hugo_Symbol %in% c("BRAF") &isCOSMIChotspot %in% c("True","TRUE",TRUE))%>%
  select(Hugo_Symbol,DepMap_ID)%>%
  distinct()
Mut_ceres <- ceres%>%
  as.data.frame()%>%
  tibble::rownames_to_column("DepMap_ID")%>%
  select(DepMap_ID,`DUSP4 (1846)`)%>%
  magrittr::set_colnames(c("DepMap_ID","DUSP4"))%>%
  left_join(mut_annot,by="DepMap_ID")%>%
  arrange(-DUSP4)%>%
  mutate(rank=1:n())%>%
  inner_join(sample.info%>%select(DepMap_ID,clName = CCLE_Name,strpName = stripped_cell_line_name),by="DepMap_ID")%>%
  mutate(Mutations = case_when(is.na(Hugo_Symbol)~"WT",TRUE~Hugo_Symbol))
Mut_ceres$Mutations%<>%factor(.,levels=c("WT","BRAF"))
cols <- c("WT" = "#273891","BRAF" = "#f7941d") 
w <-ggplot(Mut_ceres, 
       aes(x=rank, y=DUSP4, fill=Mutations, color=Mutations)) +
  geom_rect(aes(xmin = -Inf,xmax=Inf,ymin=-Inf,ymax = -1),fill = "lightpink", alpha = 0.01,color="transparent")+
  geom_rect(aes(xmin = -Inf,xmax=Inf,ymin=-1,ymax = Inf),fill = "gray91", alpha = 0.03,color="transparent")+
  scale_colour_manual(values = cols,name="",labels = c("BRAF wild-type", "BRAF mutant"))+
  scale_fill_manual(values = cols,name="",labels = c("BRAF wild-type", "BRAF mutant"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        text = element_text(color='black'),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        strip.background = element_rect(colour="black", fill="gray81", size=1, linetype="solid"),
        axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.05, 0.05),
        legend.justification = c("left", "bottom"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6))+
  labs(x=NULL,y="DUSP4 CERES")+
  geom_bar(stat="identity", width=0.5, position = position_dodge(width=0.4))+
  geom_text_repel(
    data = Mut_ceres%>%filter(strpName %in% c("A375","A375SKINCJ1","A375SKINCJ3", "HT144", "HT144SKINFV1")),
    aes(x = rank,y=c(0.1,0,-0.8,0,0),label = c("A375 \nParental",
                                          "HT144 \nParental",
                                          "HT144 \nDr",
                                          "A375 \nTDr",
                                          "A375 \nTDSr"),fill=NULL),
    force = 1,
    nudge_x=c(-5,-5,-20,-10,5),
    nudge_y=c(0.2,-0.3,-0.3,0.3,0.3),
    color=c("black","black","darkolivegreen4","dodgerblue4","firebrick"),
    size = 3,
    arrow = arrow(length = unit(0.01, 'npc'),type = "closed"),
    fontface = 'bold', 
    family = 'Times',inherit.aes = F)+
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)))

ggsave(w,width = 5,height=5,filename = paste0("Fig 5e DUSP4 waterfall plot/","Fig5e_A375_ceres_waterfall_pub21q2_nonsilent_",Sys.Date(),".pdf"),dpi = 2000)


