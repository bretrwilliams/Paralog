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
source("Utility_Functions/gene_pair.R")
source("Utility_Functions/gemini_avg_lfc.R")
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
avgLFC <- gemini_avg_lfc(Model, LFC_center = "mean")
sampleInfo <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='sample_info')
CCLE.mutations <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='CCLE_mutations')
#===== Fig 3b =====
cellline_order <- c( "GI1","HSC5","MEWO","HS944T","IPC298","MELJUSO","HS936T","MEL202","PK1","PATU8988S","A549")
dusps <- grep("DUSP[1-9];DUSP[1-9]$|AAVS1;DUSP[1-9]$",rownames(avgLFC),value=T)%>%sort()
h1 <- avgLFC[dusps,]%>%
  as.data.frame()%>%
  rownames_to_column("pair")%>%
  gather(cl,LFC,-pair)%>%
  mutate(cl = gsub("\\_.*","",cl)%>%factor(.,levels = rev(cellline_order)),
         pair = factor(pair,levels=dusps))%>%
  ggplot(data = ., aes(x=pair, y=cl, fill=LFC)) + 
  geom_tile(color = "black",lwd = 0.3,linetype = 1)+
  geom_hline(yintercept = c(1.5*2.33,1.5*3,1.5*5.66),color="black",size=rel(2))+
  scale_x_discrete(expand=c(0,0),position = "bottom")+
  scale_y_discrete(expand=c(0,0),position = "left")+
  labs(x="",y="")+
  theme(
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,hjust = -0.05,face="bold",size=rel(1.05)),
    axis.text.y = element_text(face="bold",size=rel(1.05)),
    legend.position = "right")+
  coord_fixed()+
  scale_fill_gradient2(low="#021893",
                       high="#c0c0c0",
                       mid = "#d4d4d4",
                       limits = c(-2.7,3),
                       breaks = c(-2,0,2),
                       name="LFC")+
  guides(fill = guide_colorbar(title.position = "top",
                               barheight =rel(10),
                               frame.colour = "black",ticks = F))+
  geom_rect(xmin=1.5*15.6, xmax=1.5*16.4, ymin=1.5*3, ymax=1.5*5.66, color="red",fill=NA)
rowmap <- sampleInfo%>%
  select(DepMap_ID,ccle = stripped_cell_line_name)%>%
  filter(ccle %in% cellline_order)
mut_gene <- c("NRAS","KRAS","GNAQ","NF1")
h2 <- CCLE.mutations[,-1]%>%
  filter(Hugo_Symbol %in% mut_gene)%>%
  select(gene = Hugo_Symbol,DepMap_ID,Variant_annotation)%>%
  right_join(expand.grid(rowmap$DepMap_ID,mut_gene)%>%
               set_colnames(c("DepMap_ID","gene"))%>%
               inner_join(rowmap,by="DepMap_ID"),by=c("DepMap_ID","gene"))%>%
  mutate(annot = case_when(Variant_annotation!="silent"~"nonSilent",
                           TRUE~"others"))%>%
  mutate(ccle = factor(ccle,levels = rev(cellline_order)),
         gene = factor(paste0(gene," status"),levels =paste0(mut_gene," status")))%>%
  select(ccle,gene,annot)%>%
  group_by(ccle,gene)%>%
  summarise(annot = ifelse(n()==2,"nonSilent",annot)%>%
              factor(.,levels=c("nonSilent","others")))%>%
  ggplot(data = ., aes(x=gene, y=ccle, fill=annot)) + 
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "bottom")+
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
  scale_fill_manual(values=c("orange","#bbbdbf"))+
  geom_hline(yintercept = c(1.5*2.33,1.5*3,1.5*5.66),color="black",size=rel(1.5))
ggarrange(h1,h2,nrow = 1,widths= c(8,1))%>%
  ggsave(width = 9,height=5,
         filename = paste0("Fig 3b DUSP heatmap/",
                           "Fig3b_Pheatmap_Sreen_Mutation_DUSP_",Sys.Date(),".pdf"))
#===== Fig 3b (Inferred LFC) =====
cellline_order <- c( "GI1","HSC5","MEWO","HS944T","IPC298","MELJUSO","HS936T","MEL202","PK1","PATU8988S","A549")
dusps <- grep("DUSP[1-9];DUSP[1-9]$",rownames(avgLFC),value=T)%>%sort()
dusaps_map <- data.frame(pairs= dusps)%>%separate(pairs,into=c("g1","g2"),remove=F)
pairLFC <- (Model$s[dusaps_map$pairs,]+Model$y[dusaps_map$g1,]+Model$y[dusaps_map$g2,])%>%
  as.data.frame()%>%
  set_rownames(dusaps_map$pairs)%>%
  rownames_to_column("label")%>%
  gather(cl,lfc,-label)
geneLFC <- Model$y[paste0("DUSP",1:9),]%>%
  as.data.frame()%>%
  rownames_to_column("label")%>%
  gather(cl,lfc,-label)%>%
  mutate(label = paste0("AAVS1;",label))
Infscreend <- rbind(pairLFC,geneLFC)%>%
  mutate(cl = gsub("\\_.*","",cl))%>%
  filter(cl %in% cellline_order)%>%
  mutate(cl = factor(cl,levels = rev(cellline_order)))%>%
  mutate(label = factor(label,levels=c(paste0("AAVS1;DUSP",1:9),dusps)))
h1 <- ggplot(data = Infscreend, aes(x=label, y=cl, fill=lfc)) + 
  geom_tile(color = "black",lwd = 0.3,linetype = 1)+
  geom_hline(yintercept = c(1.5*2.33,1.5*3,1.5*5.66),color="black",size=rel(2))+
  scale_x_discrete(expand=c(0,0),position = "bottom")+
  scale_y_discrete(expand=c(0,0),position = "left")+
  labs(x="",y="")+
  theme(
    axis.ticks=element_blank(),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,hjust = -0.05,face="bold",size=rel(1.05)),
    axis.text.y = element_text(face="bold",size=rel(1.05)),
    legend.position = "right")+
  coord_fixed()+
  scale_fill_gradient2(low="#021893",
                       high="#c0c0c0",
                       mid = "#d4d4d4",
                       limits = c(-2.7,3),
                       breaks = c(-2,0,2),
                       name="LFC")+
  guides(fill = guide_colorbar(title.position = "top",
                               barheight =rel(10),
                               frame.colour = "black",ticks = F))+
  geom_rect(xmin=1.5*15.6, xmax=1.5*16.4, ymin=1.5*3, ymax=1.5*5.66, color="red",fill=NA)

rowmap <- sampleInfo%>%
  select(DepMap_ID,ccle = stripped_cell_line_name)%>%
  filter(ccle %in% cellline_order)
mut_gene <- c("NRAS","KRAS","GNAQ","NF1")
h2 <- CCLE.mutations[,-1]%>%
  filter(Hugo_Symbol %in% mut_gene)%>%
  select(gene = Hugo_Symbol,DepMap_ID,Variant_annotation,isCOSMIChotspot)%>%
  right_join(expand.grid(rowmap$DepMap_ID,mut_gene)%>%
               set_colnames(c("DepMap_ID","gene"))%>%
               inner_join(rowmap,by="DepMap_ID"),by=c("DepMap_ID","gene"))%>%
  mutate(annot = case_when((gene=="NF1" & Variant_annotation=="damaging")|isCOSMIChotspot %in% c("True","TRUE",TRUE)~"mut",
                           TRUE~"others"))%>%
  mutate(ccle = factor(ccle,levels = rev(cellline_order)),
         gene = factor(paste0(gene," status"),levels =paste0(mut_gene," status")))%>%
  select(ccle,gene,annot)%>%
  group_by(ccle,gene)%>%
  summarise(annot = ifelse(n()==2,"mut",annot)%>%
              factor(.,levels=c("mut","others")))%>%
  ggplot(data = ., aes(x=gene, y=ccle, fill=annot)) + 
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1)+
  scale_x_discrete(expand=c(0,0),position = "bottom")+
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
  scale_fill_manual(values=c("orange","#bbbdbf"))+
  geom_hline(yintercept = c(1.5*2.33,1.5*3,1.5*5.66),color="black",size=rel(1.5))
ggarrange(h1,h2,nrow = 1,widths= c(8,1))%>%
  ggsave(width = 9,height=5,
         filename = paste0("Fig 3b DUSP heatmap/",
                           "Fig3b_Inferred_Pheatmap_Sreen_Mutation_DUSP_",Sys.Date(),".pdf"))


