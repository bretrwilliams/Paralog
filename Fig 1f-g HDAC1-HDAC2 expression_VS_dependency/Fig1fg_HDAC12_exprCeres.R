require(dplyr)
require(ggplot2)
require(ggrepel)
require(ggpubr)
import::from(magrittr,"%<>%","set_colnames","set_rownames","set_names")
import::from(tibble,"column_to_rownames","rownames_to_column")
import::from(tidyr,"spread","unite","gather","separate")
HDAC2expr_HDAC1ceres <- read.csv("Data/HDAC2 log2 Expression Internal 20Q1 vs HDAC1 CERES Achilles Avana 20Q1 Public CERES.csv")
HDAC1expr_HDAC2ceres <- read.csv("Data/HDAC1 log Expression Internal 20Q1 vs HDAC2 CERESAchilles Avana 20Q1 Public CERES.csv")
#===== Fig 1f =====
Model<- readRDS("Data/LatestModel.rds")
cellline_order <- c( "MELJUSO","GI1","MEL202", "PK1",  "MEWO","HS944T","IPC298","A549","HSC5","HS936T", "PATU8988S")
h1 <- data.frame(HDAC1 = Model$y["HDAC1",],
                 HDAC2 = Model$y["HDAC2",],
                 HDAC1_HDAC2 = Model$s["HDAC1;HDAC2",])%>%
  rownames_to_column("cl")%>%
  mutate(cl = gsub("\\_.*","",cl))%>%
  mutate(HDAC1_HDAC2 = HDAC1_HDAC2+HDAC1+HDAC2)%>%
  set_colnames(c("cl","HDAC1","HDAC2","HDAC1:HDAC2"))%>%
  gather(pair,LFC,-cl)%>%
  mutate(pair = factor(pair,levels = c("HDAC1","HDAC2","HDAC1:HDAC2")),
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
    axis.text.x = element_text(angle = 90,hjust = -0.05),
    axis.text.y = element_blank(),
    legend.position = "bottom")+
  coord_fixed()+labs(x="",y="")+
  scale_fill_gradient2(low="deepskyblue4",
                       high="slategray",
                       mid = "white",
                       limits = c(-1.6,1.6),
                       breaks = c(-1.5,0,1.5),name="          LFC")+
  guides(fill = guide_colorbar(title.position = "bottom"))

h2 <- HDAC1expr_HDAC2ceres[,3:4]%>%
  set_colnames(c("expr","cl"))%>%
  filter(cl %in% cellline_order)%>%
  mutate(gene="HDAC1              ")%>%
  rbind(HDAC2expr_HDAC1ceres[,3:4]%>%
          set_colnames(c("expr","cl"))%>%
          filter(cl %in% cellline_order)%>%
          mutate(gene="HDAC2              "))%>%
  mutate(gene = factor(gene,levels = c("HDAC1              ","HDAC2              ")),
         cl = factor(cl,levels = rev(cellline_order)))%>%
  ggplot(data = ., aes(x=gene, y=cl, fill=expr)) + 
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
    axis.text.x = element_text(angle = 90,hjust = -0.05),
    legend.position = "bottom")+
  coord_fixed()+labs(x="",y="")+
  scale_fill_gradient2(low = "grey", mid = "floralwhite", high = "firebrick3", midpoint =5.5,
                       limits=c(5,8),
                       breaks=c(5,5.5,8),labels=c("Lo","","Hi"),name="        mRNA")+
  guides(fill = guide_colorbar(title.position = "bottom"))

ggarrange(h1,h2,nrow = 1,widths= c(1,1.2))%>%
  ggsave(width = 3.4,height=5,
         filename = paste0("Fig 1f-g HDAC1-HDAC2 expression_VS_dependency/",
                           "Fig1f_Pheatmap_Expr_Ceres_HDAC12_",Sys.Date(),".pdf"))

#===== Fig 1g =====
HDAC1expr_HDAC2ceresP <- ggplot(data = HDAC1expr_HDAC2ceres[,c(2:4,7)]%>%
                                     set_colnames(c("ceres","expr","cl","Paralog"))%>%
                                     arrange(Paralog),
                                   aes(x = ceres,y = expr, color = Paralog )) +
  geom_text_repel(data=.%>%filter(cl=="MEL202"),aes(label=cl),box.padding = 0.8)+
  annotate(geom="text", x= -1.2, y=9.8, label=paste0("p = ",formatC(cor.test(HDAC1expr_HDAC2ceres[,2],HDAC1expr_HDAC2ceres[,3])$p.val,digits = 3,format = "G")),color="black")+
  geom_point(size=0.75, alpha=0.4) +
  labs(title = "Avana DepMap Public 20Q1",
       x = "HDAC2 CERES", y = "HDAC1 expression")+
  scale_color_manual(values = c("Yes" = "red", "No" =  "grey"))+
  xlim(-1.5, 1)+
  ylim(4,10)+
  geom_vline(xintercept = 0, linetype="dashed", size = 0.25) +
  theme(plot.title = element_text(color='black'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position = "none")

HDAC2expr_HDAC1ceresP= ggplot(data = HDAC2expr_HDAC1ceres[,c(2:4,7)]%>%
                               set_colnames(c("ceres","expr","cl","Paralog"))%>%
                               arrange(Paralog), aes(x = ceres, y = expr, color = Paralog )) +
  geom_point(size=0.75, alpha=0.4) +
  labs(title = "",
       x = "HDAC1 CERES", y = "HDAC2 expression")+
  annotate(geom="text", x= -1.2, y=9.8, label=paste0("p = ",formatC(cor.test(HDAC2expr_HDAC1ceres[,2],HDAC2expr_HDAC1ceres[,3])$p.val,digits = 3,format = "G")),color="black")+
  scale_color_manual(values = c("Yes" = "red", "No" =  "grey"))+
  xlim(-1.5, 1)+
  ylim(4,10)+
  geom_vline(xintercept = 0, linetype="dashed", size = 0.25) +
  theme(plot.title = element_text(color='black'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position = "none")
ggarrange(HDAC1expr_HDAC2ceresP,HDAC2expr_HDAC1ceresP,nrow = 2,heights = c(1.1,1))%>%
  ggsave(width = 5,height=5,
         filename = paste0("Fig 1f-g HDAC1-HDAC2 expression_VS_dependency/","Fig1g_Scatter_Expr_Ceres_HDAC12_",Sys.Date(),".pdf"))
