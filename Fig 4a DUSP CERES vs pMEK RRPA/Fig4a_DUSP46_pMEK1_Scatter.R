library(plyr)
library(dplyr)
library(ggplot2)
library(gemini)
library(ggrepel)
library(ggpubr)
DUSP46_MEK1_RRPA <- read.csv("Data/pMEK_ DUSP46_CERSE_20Q1Public.csv")[,c(2:4,14)]%>%
  set_colnames(c("DUSP6_ceres","DUSP4_ceres","MEK1_RPPA","COSMIC mutation"))%>%
  mutate(`COSMIC mutation` = gsub(" Hotspot","",`COSMIC mutation`))
###### Scatter plot for DUSP4 and phosphoMEK RRPA ###### 
DUSP4_dependency = ggplot(data = DUSP46_MEK1_RRPA, aes(x =DUSP4_ceres,y = MEK1_RPPA, color = `COSMIC mutation`)) +
  geom_rect(aes(xmin = -Inf,xmax=-1,ymin=-Inf,ymax = Inf),fill = "lightpink", alpha = 0.03,color="transparent")+
  geom_rect(aes(xmin = -1,xmax=Inf,ymin=-Inf,ymax = Inf),fill = "gray91", alpha = 0.03,color="transparent")+
  geom_point(size=1, alpha=0.8) +
  labs(x = "DUSP4 CERES", y = "pMEK1 (S217/S221)",title="Avana Public 20Q1 CRISPR CCLE")+
  scale_color_manual(values = c("BRAF" = "red",
                                "KRAS" =  "#529EFF",
                                "NRAS" =  "dark green",
                                "Other" = "dark gray"),name=" COSMIC mutation")+
  xlim(-2, 1)+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")+
  theme(plot.title = element_text(color='black'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.02, 0.02),
        legend.justification = c("left", "bottom"),
        legend.box.just = "left",
        legend.spacing.x = unit(0.01, 'cm'),
        legend.spacing.y = unit(0.01, 'cm'))+
  guides(color=guide_legend(nrow=2, byrow=TRUE,override.aes = list(size=3,alpha=1)))+
  annotate(geom="text", x= -1.88, y=2.3, 
           label=paste0("p = ",formatC(lm(data = DUSP46_MEK1_RRPA,MEK1_RPPA~DUSP4_ceres)%>%summary()%>%coefficients()%>%.[2,4],
                                       digits = 3,format = "G")),color="black")

###### Scatter plot for DUSP6 and phosphoMEK RRPA ###### 
DUSP6_dependency = ggplot(data = DUSP46_MEK1_RRPA, aes(x =DUSP6_ceres,y = MEK1_RPPA, color = `COSMIC mutation`)) +
  geom_rect(aes(xmin = -Inf,xmax=-1,ymin=-Inf,ymax = Inf),fill = "lightpink", alpha = 0.03,color="transparent")+
  geom_rect(aes(xmin = -1,xmax=Inf,ymin=-Inf,ymax = Inf),fill = "gray91", alpha = 0.03,color="transparent")+
  geom_point(size=1, alpha=0.8) +
  labs(x = "DUSP6 CERES", y = "pMEK1 (S217/S221)",title="          ")+
  scale_color_manual(values = c("BRAF" = "red",
                                "KRAS" =  "#529EFF",
                                "NRAS" =  "dark green",
                                "Other" = "dark gray"),name=" COSMIC mutation")+
  xlim(-2, 1)+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")+
  theme(plot.title = element_text(color='black'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.02, 0.02),
        legend.justification = c("left", "bottom"),
        legend.box.just = "left",
        legend.spacing.x = unit(0.01, 'cm'),
        legend.spacing.y = unit(0.01, 'cm'))+
  guides(color=guide_legend(nrow=2, byrow=TRUE,override.aes = list(size=3,alpha=1)))+
  annotate(geom="text", x= -1.88, y=2.3, 
           label=paste0("p = ",formatC(lm(data = DUSP46_MEK1_RRPA,MEK1_RPPA~DUSP6_ceres)%>%summary()%>%coefficients()%>%.[2,4],
                                       digits = 3,format = "G")),color="black")
ggarrange(DUSP4_dependency,DUSP6_dependency,nrow = 1)%>%
  ggsave(width = 15,height=3.5,
         filename = paste0("Fig 4a DUSP CERES vs pMEK RRPA/Fig4a_Scatter_pMEK1_DUSP46_",Sys.Date(),".pdf"))
