require(taigr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)
##### 1) Data loading and processing #####
CCLE.mutations <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='CCLE_mutations')
Achilles.gene.effect <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='Achilles_gene_effect')
CCLE.expression <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='CCLE_expression')
sample.info <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='sample_info')
proteomics <- load.from.taiga(data.name='total-proteome--5c50',data.version=1, data.file='protein_quant_current_normalized')
mut_annot <- CCLE.mutations%>%
  select(Hugo_Symbol,isCOSMIChotspot,DepMap_ID)%>%
  filter(Hugo_Symbol %in% c("BRAF","NRAS") & isCOSMIChotspot %in% c("True","TRUE",TRUE))%>%
  group_by(DepMap_ID)%>%
  summarise(Genes = paste(sort(unique(Hugo_Symbol)),collapse=";"))
DepD <- Achilles.gene.effect[,c("DUSP4 (1846)","DUSP6 (1848)")]%>%
  as.data.frame()%>%
  magrittr::set_colnames(c("DUSP4","DUSP6"))%>%
  tibble::rownames_to_column("DepMap_ID")%>%
  gather(Dgene,Dependency,-DepMap_ID)%>%
  mutate(source="CERES")
CCLE <- CCLE.expression[,c("DUSP4 (1846)","DUSP6 (1848)")]%>%
  as.data.frame()%>%
  magrittr::set_colnames(c("DUSP4","DUSP6"))%>%
  tibble::rownames_to_column("DepMap_ID")%>%
  gather(Egene,expr,-DepMap_ID)%>%
  inner_join(DepD,by="DepMap_ID")%>%
  left_join(mut_annot,by="DepMap_ID")%>%
  mutate(Situ = purrr::pmap_chr(list(Egene,Dgene),~paste0("x:",c(...)[2],
                                                          ", ","y:",c(...)[1])))%>%
  mutate(ExprSource = "CCLE")
prote <- proteomics%>%
  filter(Gene_Symbol %in% c("DUSP4","DUSP6"))%>%
  t()%>%
  as.data.frame(stringsAsFactors = FALSE)%>%
  magrittr::set_colnames(c("DUSP4","DUSP6"))%>%
  tibble::rownames_to_column("cl")%>%
  filter(grepl("_TenP",cl))%>%
  mutate(cclename = gsub("_TenP.*","",cl))%>%
  inner_join(sample.info%>%select(DepMap_ID,cclename=CCLE_Name),by="cclename")%>%
  select(DepMap_ID,DUSP4,DUSP6)
prote%<>%
  gather(Egene,expr,-DepMap_ID)%>%
  inner_join(DepD,by="DepMap_ID")%>%
  left_join(mut_annot,by="DepMap_ID")%>%
  mutate(Situ = purrr::pmap_chr(list(Egene,Dgene),~paste0("x:",c(...)[2],
                                                          ", ","y:",c(...)[1])))%>%
  mutate(ExprSource = "Proteomics")%>%
  filter(source!="DEMETER")%>%
  mutate_at(vars(Dependency,expr),as.numeric)
##### 2) Plotting #####
pA<- CCLE%>%
  rbind(prote)%>%
  filter(!Situ %in% c("x:DUSP6, y:DUSP4","x:DUSP6, y:DUSP6"))%>%
  mutate(Mutation = factor(case_when(is.na(Genes)~"WT",
                              TRUE~Genes),
                           levels=c("BRAF","NRAS","WT")))%>%
  group_by(ExprSource,source,Situ)%>%
  do(plot = if(unique(.$ExprSource)=="CCLE"){
    ggplot(data=.,aes(Dependency,expr,colour = Mutation)) +
      scale_color_manual(values=c("firebrick4", "turquoise4", "grey50"),name = "Mutations")+
      geom_point(data=.%>%filter(Mutation=="WT"),alpha=0.7)+
      geom_point(data=.%>%filter(Mutation!="WT"),alpha=0.7)+
      theme(plot.title = element_text(color='black', hjust = 0.5),
             plot.background = element_blank(),
             panel.background = element_blank(),
             panel.border = element_rect(color = "black", size = 0.25, fill = NA),
             text = element_text(color='black'),
             panel.grid = element_blank(),
             axis.text = element_text(color='black'),
             strip.background = element_rect(colour="black", fill="gray81", 
                                             size=1, linetype="solid"),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.position = c(0.98, 0.02),
            legend.justification = c("right", "bottom"),
            legend.box.just = "right")+
      guides(color=guide_legend(override.aes = list(size=rel(5),alpha=1)))+
       geom_vline(xintercept = 0, linetype="dashed")+
      annotate(geom="text", x= -1.3, y=9.3, 
               label=paste0("p = ",formatC(lm(data = .,Dependency~expr)%>%summary()%>%coefficients()%>%.[2,4],
                                           digits = 3,format = "G")),color="black")+
      labs(x=paste0(unique(.$Dgene)," ",unique(.$source)),y=paste0(unique(.$Egene)," expression"))+
      xlim(c(-2,1.5))+
      ylim(c(0,10))}else{
        ggplot(data=.,aes(Dependency,expr,colour = Mutation)) +
          scale_color_manual(values=c("firebrick4", "turquoise4", "grey50"),name = "Mutations")+
          geom_point(data=.%>%filter(Mutation=="WT"),alpha=0.7)+
          geom_point(data=.%>%filter(Mutation!="WT"),alpha=0.7)+
          theme(plot.title = element_text(color='black', hjust = 0.5),
                plot.background = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color = "black", size = 0.25, fill = NA),
                text = element_text(color='black'),
                panel.grid = element_blank(),
                axis.text = element_text(color='black'),
                strip.background = element_rect(colour="black", fill="gray81", 
                                                size=1, linetype="solid"),
                legend.background = element_blank(),
                legend.key = element_blank(),
                legend.position = c(0.98, 0.02),
                legend.justification = c("right", "bottom"),
                legend.box.just = "right")+
          annotate(geom="text", x= -1.2, y=3.8, 
                   label=paste0("p = ",formatC(lm(data = .,Dependency~expr)%>%summary()%>%coefficients()%>%.[2,4],
                                               digits = 3,format = "G")),color="black")+
          guides(color=guide_legend(override.aes = list(size=rel(5),alpha=1)))+
          geom_vline(xintercept = 0, linetype="dashed")+
          labs(x=paste0(unique(.$Dgene)," ",unique(.$source)),
               y=paste0(unique(.$Egene)," proteomics"))+
          xlim(c(-1.5,1.5))+ylim(c(-4,4))})

pBD <- Achilles.gene.effect[,c("DUSP4 (1846)","BRAF (673)")]%>%
  as.data.frame()%>%
  magrittr::set_colnames(c("DUSP4","BRAF"))%>%
  tibble::rownames_to_column("DepMap_ID")%>%
  left_join(mut_annot,by="DepMap_ID")%>%
  mutate(Mutation = factor(case_when(is.na(Genes)~"WT",TRUE~Genes),levels=c("BRAF","NRAS","WT")))
pB <- ggplot(data=pBD,aes(DUSP4,BRAF,colour = Mutation)) +
  annotate(geom="text", x= -1.5, y=0.8,
           label=paste0("p = ",formatC(lm(data=pBD,DUSP4~BRAF)%>%summary()%>%coefficients()%>%.[2,4],
                                       digits = 3,format = "G")),color="black")+
  scale_color_manual(values=c("firebrick4", "turquoise4", "grey50"),name = "Mutations")+
  geom_point(data=.%>%filter(Mutation=="WT"),alpha=0.7)+
  geom_point(data=.%>%filter(Mutation!="WT"),alpha=0.7)+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        text = element_text(color='black'),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        strip.background = element_rect(colour="black", fill="gray81", 
                                        size=1, linetype="solid"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right")+
  guides(colour = guide_legend(override.aes = list(size=rel(5))))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(x="DUSP4 CERES",y="BRAF CERES")+
  xlim(c(-2,1.5))+
  ylim(c(-3,1))

interestedG <- sapply(c("MAPK1 ","SOX10 ","DUSP4 ","MITF ","MAP2K1 ","BRAF "),function(s) grep(s,colnames(Achilles.gene.effect),value=T))%>%as.character()
tabD <- cor(Achilles.gene.effect[,interestedG],method = "pearson")[,"BRAF (673)"]%>%
  round(.,3)%>%
  as.data.frame()%>%
  magrittr::set_colnames("Correlation")%>%
  tibble::rownames_to_column("Gene")%>%
  filter(!grepl("BRAF",Gene))%>%
  mutate(Gene=gsub(" \\(.*","",Gene))
interestedG2 <- sapply(c("BRAF ","SOX10 ","MAPK1 ","DUSP4 ","MAP2K1 "),function(s) grep(s,colnames(Achilles.gene.effect),value=T))%>%as.character()
indi_BRAF <- ifelse(rownames(Achilles.gene.effect) %in% (mut_annot%>%filter(grepl("BRAF",Genes))%>%pull(DepMap_ID)),1,0)
tabD2 <- apply(Achilles.gene.effect[,interestedG2],2,function(s) ltm::biserial.cor(s, indi_BRAF, use = c("all.obs"), level = 2))%>%
  round(.,3)%>%
  as.data.frame()%>%
  magrittr::set_colnames("Correlation")%>%
  tibble::rownames_to_column("Gene")%>%
  mutate(Gene=gsub(" \\(.*","",Gene))
tab1 <- ggtexttable(tabD, rows = NULL, theme = ttheme("mBlue"))%>%
  tab_add_title(text = "Top correlated CRISPR \n gene to BRAF knockout", face = "plain", size = 10,
                just="left") 
tab2 <- ggtexttable(tabD2, rows = NULL, theme = ttheme("mBlue"))%>%
  tab_add_title(text = "Top correlated CRISPR \n gene to BRAF mutation", face = "plain", size = 10,
                just="left") 
tabs <- ggpubr::ggarrange(tab1, tab2,  ncol = 2, nrow = 1,align="h")
##### 3) Final combination #####
ggpubr::ggarrange(ggpubr::ggarrange(tabs, pB,  ncol = 2, nrow = 1,labels = c("a","b"),align="h"),
          ggpubr::ggarrange(plotlist=pA$plot, nrow = 2,ncol=2,align="hv",labels=c("c","","d","")),
          nrow = 2,heights =c(1.3,3))
ggsave(width = 8,height = 8,dpi = 2000,filename = paste0("Fig S5 DUSP4 AVANA correlations/FigS5_DUSP46_BRAF_NRAS_Scatter_",Sys.Date(),".pdf"))



