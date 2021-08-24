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
ceres <- load.from.taiga(data.name='public-21q2-110d', data.version=12, data.file='Achilles_gene_effect')
sample.info <- load.from.taiga(data.name='public-21q2-110d', data.version=12, data.file='sample_info')
A375_cls <- sample.info%>%
  select(CCLE_Name,DepMap_ID)%>%
  filter(CCLE_Name %in% c("A375_SKIN","A375_SKIN_CJ1_RESISTANT","A375_SKIN_CJ3_RESISTANT"))
cl3c <- load.from.taiga(data.name='public-21q2-110d', data.version=12, data.file='Achilles_gene_effect')%>%
  .[A375_cls$DepMap_ID,]%>%
  as.data.frame()%>%
  rownames_to_column("DepMap_ID")%>%
  gather(gene,ceres21,-DepMap_ID)%>%
  mutate(gene = gsub(" \\(.*","",gene))%>%
  inner_join(A375_cls,by="DepMap_ID")%>%
  select(-DepMap_ID)%>%
  spread(CCLE_Name,ceres21)%>%
  mutate(m1p = A375_SKIN_CJ1_RESISTANT-A375_SKIN,
         m3p = A375_SKIN_CJ3_RESISTANT-A375_SKIN)
D <- cl3c%>%
  select(m1p,m3p,gene)%>%
  gather(situ,diff,-gene)%>%
  mutate(situ=case_when(situ=="m1p"~"A375: TDr versus parental",
                        situ=="m3p"~"A375: TDSr versus parental"))%>%
  group_by(situ)%>%
  arrange(diff)%>%
  mutate(index=1:n())
D%<>%
  left_join(D%>%group_by(situ)%>%arrange(diff)%>%slice(1:5)%>%
              as.data.frame()%>%
              mutate(note=c(rep("CJ1_neg5",5),rep("CJ3_neg5",5)))%>%
              select(gene,situ,note),
            by=c("gene","situ"))%>%
  mutate(note=case_when(is.na(note)~"Others",TRUE~note))

ggplot(D, aes(x=index, y=diff,color=note)) + 
  facet_wrap(~situ,nrow=1,scales = "free")+
  xlim(c(0,18353))+ylim(-1.5,1.5)+
  scale_colour_manual(values = c(CJ3_neg5="firebrick",CJ1_neg5="dodgerblue4",Others="gray30"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        text = element_text(color='black'),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Helvetica"),
        axis.title = element_text(color='black',family="Helvetica"),
        strip.background = element_blank(),
        strip.text = element_text(family="Helvetica",face="bold",size=15),
        legend.position = "none")+
  labs(x="Gene Index",y="Delta CERES \n (Resistent - Parental)")+
  geom_point(alpha=0.6,size=0.8)+
  geom_hline(yintercept=0,linetype="dashed", size=0.5)+
  geom_text_repel(data = D%>%group_by(situ)%>%arrange(diff)%>%slice(1:5)%>%
              rbind(D%>%group_by(situ)%>%arrange(-diff)%>%slice(1:5)),
            aes(label = gene,fill=NULL),
            box.padding = 0.1,
            point.padding = 0.5,
            size = 3,
            colour=c(rep("dodgerblue4",5),rep("firebrick",5),rep("black",10)),
            fontface = 'bold', 
            family = 'Helvetica')
ggsave(width = 10,height=5,filename = "Fig 5c-d A375 ranked depedency/Fig5cd_A375_Delta_Ceres_scatter_public21q2.pdf",dpi = 2000)


