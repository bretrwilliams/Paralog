require(taigr)
require(tidyr)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(magrittr)
import::from(tibble,rownames_to_column)
ceres <- load.from.taiga(data.name='public-21q2-110d', data.version=12, data.file='Achilles_gene_effect')
sample.info <- load.from.taiga(data.name='public-21q2-110d', data.version=12, data.file='sample_info')
ht144_cls <- sample.info%>%
  select(CCLE_Name,DepMap_ID)%>%
  filter(CCLE_Name %in% c("HT144_SKIN","HT144_SKIN_FV1_RESISTANT"))
##### For HT144 #####
ht144 <- ceres%>%
  .[ht144_cls$DepMap_ID,]%>%
  as.data.frame()%>%
  rownames_to_column("DepMap_ID")%>%
  gather(gene,ceres21,-DepMap_ID)%>%
  mutate(gene = gsub(" \\(.*","",gene))%>%
  inner_join(ht144_cls,by="DepMap_ID")%>%
  select(-DepMap_ID)%>%
  spread(CCLE_Name,ceres21)%>%
  mutate(diff = HT144_SKIN_FV1_RESISTANT-HT144_SKIN)%>%
  arrange(diff)%>%
  mutate(index=1:n())
ht144%<>%left_join(ht144%>%arrange(diff)%>%slice(1:10)%>%mutate(note="Neg")%>%
            rbind(ht144%>%arrange(-diff)%>%slice(1:10)%>%mutate(note="Pos"))%>%
            select(gene,note),
          by=c("gene"))%>%
  mutate(note=case_when(is.na(note)~"Others",TRUE~note))

ggplot(ht144, aes(x=index, y=diff,color=note))+
  scale_colour_manual(values = c(Neg="darkgreen",Pos="black",Others="gray30"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        text = element_text(color='black'),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        strip.background = element_rect(colour="black", fill="gray81", size=1, linetype="solid"),
        legend.position = "none")+
  labs(x="Gene Index",y="Delta CERES (Resistant - Parental)",
       title = "HT144: Dr versus parental",hjust=0.5)+
  geom_point(alpha=0.6)+
  geom_text(
    data = ht144%>%filter(note!="Others"),
            aes(label = gene,fill=NULL),
            size = rel(2.8),
            colour=c(rep("darkgreen",10),rep("black",10)),
            fontface = 'bold', 
            family = 'Times',
    nudge_x = c(rep(1500,10),rep(-1300,10)),
    nudge_y = c(c(0,-0.03,-0.015,0.02,seq(0.03,0.025+0.04*5,length.out = 6)),
                c(-0.08,-0.04,-0.025,0,0.03,
                  0.04,0.06,0.08,0.1,0.05))
    )+
  geom_hline(yintercept = 0, linetype="dashed", size=0.5)
ggsave(width = 6,height=6,filename = paste0("Fig S6c HT144 ranked depedency/FigS6c_HT144_scatter_",Sys.Date(),".pdf"),dpi = 2000)
