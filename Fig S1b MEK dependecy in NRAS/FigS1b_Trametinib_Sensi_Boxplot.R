require(taigr)
require(ggplot2)
require(ggpubr)
library(dplyr)
library(tidyr)
#===== Data Loading and processing =====
mut <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='CCLE_mutations')
mut_annot <- mut%>%
  select(Hugo_Symbol,isCOSMIChotspot,DepMap_ID)%>%
  filter(Hugo_Symbol =="NRAS" & isCOSMIChotspot %in% c("True","TRUE",TRUE))%>%
  group_by(DepMap_ID)%>%
  summarise(Mutations = paste(sort(unique(Hugo_Symbol)),collapse=";"))
### Drug Screens ###
AUC <- load.from.taiga(data.name='ctrp-v2-9f98', data.version=1, data.file='v20.data.curves_post_qc')
maps <- load.from.taiga(data.name='ctrp-v2-9f98', data.version=1, data.file='v20.meta.per_compound')%>%
  select(master_cpd_id,cpd_name)
cl <- load.from.taiga(data.name='ctrp-v2-9f98', data.version=1, data.file='v20.meta.per_cell_line')%>%select(ccl_name,master_ccl_id)
expr <- load.from.taiga(data.name='ctrp-v2-9f98', data.version=1, data.file='v20.meta.per_experiment')%>%select(experiment_id,master_ccl_id)
cl.annot <- load.from.taiga(data.name='other-ccle2-c93e', data.version=1, data.file='Cell_lines_annotations_20181226')
cl.annot%<>%
  select(CCLE_ID,DepMap_ID)%>%
  mutate(ccl_name = gsub("\\_.*","",CCLE_ID))
clmap <- cl%>%
  inner_join(expr,by="master_ccl_id")%>%
  left_join(cl.annot,by="ccl_name")
AUC_Drug <- AUC%>%
  select(experiment_id,auc = area_under_curve,master_cpd_id)%>%
  inner_join(maps%>%filter(cpd_name %in% c("trametinib","selumetinib")),by="master_cpd_id")%>%
  select(-master_cpd_id)%>%
  distinct()%>%
  left_join(clmap,by="experiment_id")%>%
  select(interest = cpd_name,Measure=auc,DepMap_ID)%>%
  left_join(mut_annot,by="DepMap_ID")
### CERES ###
ceres <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='Achilles_gene_effect')
G <- c("NRAS (4893)","MAP2K1 (5604)","MAP2K2 (5605)")
Mut_ceres <- ceres[,G]%>%
  magrittr::set_colnames(gsub(" \\(.*","",G))%>%
  as.data.frame()%>%
  tibble::rownames_to_column("DepMap_ID")%>%
  left_join(mut_annot,by="DepMap_ID")
Mut_ceres%<>%
  mutate(Mutations = case_when(is.na(Mutations)~"WT",
                               TRUE~Mutations))
gatheredD <- Mut_ceres%>%
  gather(interest,Measure,-DepMap_ID,-Mutations)%>%
  group_by(interest)%>%
  arrange(-Measure)%>%
  mutate(rank=1:nrow(Mut_ceres))%>%
  mutate(source="CERES")%>%
  select(Mutations,Measure,interest,source,rank)
gatheredD2 <- AUC_Drug%>%
  group_by(interest)%>%
  arrange(-Measure)%>%
  mutate(rank=1:n())%>%
  mutate(source="AUC")%>%
  select(Mutations,Measure,interest,source,rank)
gatheredD2%<>%
  mutate(Mutations = case_when(is.na(Mutations)~"WT",
                               TRUE~Mutations))
#===== Plotting =====
b1 <- ggplot(gatheredD2%>%
               mutate(interest  = gsub("tra","Tra",interest))%>%
               filter(interest=="Trametinib")%>%
               mutate(Mut = factor(case_when(Mutations=="NRAS"~"MUT",TRUE~"WT"),
                                   levels=c("WT","MUT"))), 
             aes(x=Mut, y=Measure,fill=Mut)) + 
  geom_boxplot()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),name="Mutations")+
  facet_wrap(~interest,nrow=1)+
  theme(plot.title = element_text(color='black', hjust = 0.25),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = .5, fill = NA),
        text = element_text(color='black'),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',size=rel(1.5)),
        axis.title= element_text(color='black',size=rel(1.5)),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=rel(1.5)))+
  labs(x=NULL,y="Dose AUC")+
  stat_compare_means(method = "t.test",label =  "p.format")
b2 <- ggplot(gatheredD%>%
               filter(interest %in% c("NRAS","MAP2K1","MAP2K2"))%>%
               mutate(interest =factor(interest,levels = c("NRAS","MAP2K1","MAP2K2")))%>%
               mutate(Mut = factor(case_when(Mutations=="NRAS"~"MUT",TRUE~"WT"),
                                   levels=c("WT","MUT"))), 
             aes(x=Mut, y=Measure,fill=Mut)) + 
  geom_boxplot()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),name="Mutations")+
  facet_wrap(~interest,nrow=1)+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        text = element_text(color='black'),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',size=rel(1.5)),
        axis.title= element_text(color='black',size=rel(1.5)),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=rel(1.5)),
        legend.position = "none")+
  labs(x=NULL,y="CERES")+
  geom_hline(yintercept = -1, linetype="dashed", size=0.2)+
  stat_compare_means(method = "t.test",label =  "p.format")
ggsave(egg::ggarrange(b1, b2,
                   ncol=2,widths =  c(2,6),
                   labels = c("Compound:","sgRNA:")),
         width = 12,height = 7.5,dpi = 2000,
       filename = paste0("Fig S1b MEK dependecy in NRAS/FigS1b_MEK_dep_boxplot_",Sys.Date(),".pdf"))
