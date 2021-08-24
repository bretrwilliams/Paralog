require(dplyr)
require(ggplot2)
require(ggrepel)
require(ggpubr)
require(taigr)
require(egg)
import::from(magrittr,"%<>%","set_colnames","set_rownames","set_names")
import::from(tibble,"column_to_rownames","rownames_to_column")
import::from(tidyr,"spread","unite","gather","separate")
#===== ggplot2 extended utils =====
discrete_gradient_pal <- function(colours, bins = 5) {
  ramp <- scales::colour_ramp(colours)
  
  function(x) {
    if (length(x) == 0) return(character())
    
    i <- floor(x * bins)
    i <- ifelse(i > bins-1, bins-1, i)
    ramp(i/(bins-1))
  }
}
scale_fill_discrete_gradient <- 
  function(..., colours, bins = 5, 
           na.value = "grey50", 
           guide = "colourbar", 
           aesthetics = "fill", colors)  {
    colours <- if (missing(colours)) 
      colors
    else colours
    continuous_scale(
      aesthetics,
      "discrete_gradient",
      discrete_gradient_pal(colours, bins),
      na.value = na.value,
      guide = guide,
      ...
    )
  } 
#===== Fig Data loading =====
Model<- readRDS("Data/LatestModel.rds")
LFC_prism <- load.from.taiga(data.name='primary-screen-e5c7', data.version=8, data.file='primary_replicate_collapsed_logfold_change')
trt_prism <- load.from.taiga(data.name='primary-screen-e5c7', data.version=8, data.file='primary_replicate_collapsed_treatment_info')
sampleInfo <- sample.info <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='sample_info')
CCLE.mutations <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='CCLE_mutations')
#===== Fig 2b =====
cellline_order <- c( "GI1","MEWO","PATU8988S","MELJUSO","A549","PK1","HS944T","IPC298")
h1 <- data.frame(CDK4 = Model$y["CDK4",],
                 CDK6 = Model$y["CDK6",],
                 CDK4_CDK6 = Model$s["CDK4;CDK6",])%>%
  rownames_to_column("cl")%>%
  mutate(cl = gsub("\\_.*","",cl))%>%
  filter(cl %in% cellline_order)%>%
  mutate(CDK4_CDK6 = CDK4_CDK6+CDK4+CDK6)%>%
  set_colnames(c("cl","CDK4","CDK6","CDK4:CDK6"))%>%
  gather(pair,LFC,-cl)%>%
  mutate(pair = factor(pair,levels = c("CDK4","CDK6","CDK4:CDK6")),
         cl = factor(cl,levels = rev(cellline_order)))%>%
  ggplot(data = ., aes(x=pair, y=cl, fill=LFC)) + 
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "top")+
  scale_y_discrete(expand=c(0,0),position = "right")+
  labs(x="",y="")+
  theme(
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,hjust = -0.05,face="bold",size=rel(1.05)),
    axis.text.y = element_blank(),
    legend.position = "bottom")+
  coord_fixed()+
  scale_fill_gradient2(low="deepskyblue4",
                       high="slategray",
                       mid = "white",
                       limits = c(-3,3),
                       breaks = c(-3,0,3),name="          LFC")+
  guides(fill = guide_colorbar(title.position = "bottom",frame.colour = "black"))
# PRISM screen data loading and wrangling
rowmap <- sampleInfo%>%
  select(DepMap_ID,ccle = stripped_cell_line_name)%>%
  filter(ccle %in% cellline_order)
colmap <- trt_prism%>%
  select(column_name,name)%>%
  distinct()%>%
  filter(name %in% c("palbociclib","ribociclib"))
h2 <- LFC_prism[,colmap$column_name]%>%
  as.data.frame()%>%
  rownames_to_column("DepMap_ID")%>%
  inner_join(rowmap,by="DepMap_ID")%>%
  select(-DepMap_ID)%>%
  gather(column_name,prism,-ccle)%>%
  inner_join(colmap,by="column_name")%>%
  select(-column_name)%>%
  mutate(name = stringr::str_to_title(name))%>%
  mutate(name = factor(name,levels = c("Palbociclib","Ribociclib")),
         ccle = factor(ccle,levels = rev(cellline_order)))%>%
  group_by(name)%>%
  arrange(prism)%>%
  mutate(r_c = 1:n())%>%
  ggplot(data = ., aes(x=name, y=ccle, fill=r_c)) + 
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "top")+
  scale_y_discrete(expand=c(0,0),position = "right")+
  labs(x="",y="")+
  theme(
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,hjust = -0.05,face="bold",size=rel(1)),
    axis.text.y = element_blank(),
    legend.position = "bottom")+
  coord_fixed()+
  scale_fill_discrete_gradient(limits = c(1,8),
                               breaks = c(1,8), 
                               labels = c("Lo","Hi"),
                               colors = rev(c("#eaeaea","#bababa","#a3a4a6","#9aa0a4","#5d819d","#4c7a9c","#085b94","#015393")),
                               bins = 8,
                               guide = guide_colourbar(frame.colour = "black", 
                                                       ticks.colour = "black",
                                                       barwidth=rel(5)),
                               name="          PRISM\n            LFC")+
  guides(fill = guide_colorbar(title.position = "bottom",frame.colour = "black"))

h3 <- CCLE.mutations[,-1]%>%
  filter(Hugo_Symbol %in% c("RB1"))%>%
  select(Hugo_Symbol,DepMap_ID,Variant_annotation)%>%
  right_join(rowmap,by="DepMap_ID")%>%
  mutate(annot = case_when(Variant_annotation=="damaging"~"dmgMut",
                           TRUE~"others")%>%factor())%>%
  mutate(gene = "RB1 status",
         ccle = factor(ccle,levels = rev(cellline_order)))%>%
  select(ccle,gene,annot)%>%
  ggplot(data = ., aes(x=gene, y=ccle, fill=annot)) + 
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "top")+
  scale_y_discrete(expand=c(0,0),position = "right")+
  labs(x="",y="")+
  theme(
    axis.text = element_text(face="bold",size=rel(0.8)),
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,hjust = -0.05),
    legend.position = "none")+
  coord_fixed()+
  scale_fill_manual(values=c("#f05a29","#bbbdbf"))
ggarrange(h1,h2,h3,nrow = 1,widths= c(3,2,1))%>%
  ggsave(width = 5,height=5,
         filename = paste0("Fig 2b-c CRISPR_VS_PRISM//",
                           "Fig2b_Pheatmap_Screen_PRISM_CDK46_",Sys.Date(),".pdf"))
#===== Fig 2c =====
cellline_order <- c( "GI1","MEWO","IPC298","MELJUSO","HS944T","PATU8988S","A549","PK1")
pairs <- c("MAP2K1;MAP2K2","MAPK1;MAPK3","ARAF;BRAF","ARAF;RAF1","BRAF;RAF1")
screend <- data.frame(pairs= pairs)%>%
  separate(pairs,into=c("g1","g2"),remove=F)%>%
  gather(side,gene,-pairs)%>%
  inner_join(Model$y%>%as.data.frame()%>%rownames_to_column("gene")%>%gather(cl,lfc,-gene),by="gene")%>%
  select(-gene)%>%
  spread(side,lfc)%>%
  inner_join(Model$s%>%as.data.frame()%>%rownames_to_column("pairs")%>%gather(cl,lfc,-pairs),by=c("pairs","cl"))%>%
  mutate(pairLFC = g1+g2+lfc)%>%
  mutate(cl = gsub("\\_.*","",cl))%>%
  filter(cl %in% cellline_order)%>%
  mutate(cl = factor(cl,levels = rev(cellline_order)))%>%
  select(-lfc)%>%
  gather(label,LFC,-pairs,-cl)%>%
  separate(pairs,into=c("g1","g2"),remove=F)%>%
  mutate(label = case_when(label=="g1" ~ g1,
                           label=="g2"~g2,
                           label=="pairLFC"~pairs))%>%
  mutate(og = case_when(grepl(";",label)~"pair",TRUE~"gene")%>%factor(.,levels=c("gene","pair")))%>%
  mutate(part = case_when(pairs == "MAP2K1;MAP2K2"~"p1",
                          pairs == "MAPK1;MAPK3"~"p2",
                          TRUE~"p3")%>%factor(.,levels = c("p1","p2","p3")))%>%
  select(cl,label,LFC,part,og)%>%
  mutate(label = gsub(";"," - ",label))

fixed_order <- screend%>%
  select(label,part,og)%>%
  distinct()%>%
  arrange(part,og)%>%
  pull(label)
h1 <- ggplot(screend%>%
               mutate(label = factor(label,levels = fixed_order)), aes(x=label,y=cl,fill=LFC))+
  geom_tile(color = "black",lwd = 0.3,linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "top")+
  scale_y_discrete(expand=c(0,0),position = "right")+
  labs(x="",y="")+
  facet_grid(~part,drop=T,scales="free_x",space = "free_x")+
  theme(
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.spacing = unit(2, "lines"),
    axis.text.x = element_text(angle = 90,hjust = -0.05,face="bold",size=rel(1.05)),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom")+
  scale_fill_gradient2(low="#015393",
                       high="#eaeaea",
                       mid = "white",
                       limits = c(-2.4,2.4),
                       breaks = c(-2,0.5,2),
                       name=paste0(paste0(rep("           ",5),collapse = ""),"LFC"))+
  guides(fill = guide_colorbar(title.position = "bottom",frame.colour = "black",barwidth=rel(25)))

# PRISM screen data loading and wrangling
rowmap <- sampleInfo%>%
  select(DepMap_ID,ccle = stripped_cell_line_name)%>%
  filter(ccle %in% cellline_order)
colmap <- trt_prism%>%
  select(column_name,name)%>%
  distinct()%>%
  filter(name %in% c("trametinib","PD-0325901"))
h2 <- LFC_prism[,colmap$column_name]%>%
  as.data.frame()%>%
  rownames_to_column("DepMap_ID")%>%
  inner_join(rowmap,by="DepMap_ID")%>%
  select(-DepMap_ID)%>%
  gather(column_name,prism,-ccle)%>%
  inner_join(colmap,by="column_name")%>%
  select(-column_name)%>%
  mutate(name = gsub("tra","Tra",name))%>%
  mutate(name = factor(name,levels = c("Trametinib","PD-0325901")),
         ccle = factor(ccle,levels = rev(cellline_order)))%>%
  group_by(name)%>%
  arrange(prism)%>%
  mutate(r_c = 1:n())%>%
  ggplot(data = ., aes(x=name, y=ccle, fill=r_c)) + 
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "top")+
  scale_y_discrete(expand=c(0,0),position = "right")+
  labs(x="",y="")+
  theme(
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,hjust = -0.05,face="bold",size=rel(1)),
    axis.text.y = element_blank(),
    legend.position = "bottom")+
  coord_fixed()+
  scale_fill_discrete_gradient(limits = c(1,8),
                               breaks = c(1,8), 
                               labels = c("Lo","Hi"),
                               colors = rev(c("#eaeaea","#bababa","#a3a4a6","#9aa0a4","#5d819d","#4c7a9c","#085b94","#015393")),
                               bins = 8,
                               guide = guide_colourbar(frame.colour = "black", 
                                                       ticks.colour = "black",
                                                       barwidth=rel(5)),
                               name="          PRISM\n            LFC")+
  guides(fill = guide_colorbar(title.position = "bottom",frame.colour = "black"))


h3 <- CCLE.mutations[,-1]%>%
  filter(Hugo_Symbol %in% c("KRAS","NRAS"))%>%
  select(gene = Hugo_Symbol,DepMap_ID,isCOSMIChotspot)%>%
  right_join(expand.grid(rowmap$DepMap_ID,c("KRAS","NRAS"))%>%
               set_colnames(c("DepMap_ID","gene"))%>%
               inner_join(rowmap,by="DepMap_ID"),by=c("DepMap_ID","gene"))%>%
  mutate(annot = case_when(isCOSMIChotspot %in% c("True","TRUE",TRUE)~"cosmic",
                           TRUE~"others")%>%factor())%>%
  mutate(ccle = factor(ccle,levels = rev(cellline_order)),
         gene = paste0(gene," status"))%>%
  select(ccle,gene,annot)%>%
  ggplot(data = ., aes(x=gene, y=ccle, fill=annot)) + 
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "top")+
  scale_y_discrete(expand=c(0,0),position = "right")+
  labs(x="",y="")+
  theme(
    axis.text = element_text(face="bold",size=rel(0.8)),
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,hjust = -0.05),
    legend.position = "none")+
  coord_fixed()+
  scale_fill_manual(values=c("#f05a29","#bbbdbf"))
ggarrange(h1,h2,h3,nrow = 1,widths= c(6.2,0.8,0.8))%>%
  ggsave(width = 7.8,height=5,
         filename = paste0("Fig 2b-c CRISPR_VS_PRISM//",
                           "Fig2c_Pheatmap_Screen_PRISM_MEK_",Sys.Date(),".pdf"))



