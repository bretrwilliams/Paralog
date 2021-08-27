##### libraries #####
library(tidyverse)
library(taigr)
library(celllinemapr)
library(mixtools)
library(ggrepel)
source("Utility_Functions/gene_pair.R")
source("Utility_Functions/gemini_avg_lfc.R")
##### Model and avg_LFC #####
Model <- readRDS("Data/LatestModel_ERKi.rds")
avg_LFC <- gemini_avg_lfc(Model, LFC_center = "median")
##### visulize LFCs #####
avg_LFC_row2col <- as.data.frame(avg_LFC)%>%
  rownames_to_column(var = "pairs") %>%
  mutate(pairs = gsub(";","-",pairs))%>%
  mutate(random_index = sample(n()))%>%
  mutate(dir_r = case_when(abs(MELJUSO_SKIN_ERKi - MELJUSO_SKIN_DMSO)>0.5~"Responsive",
                           TRUE~"Normal"))
g <- ggplot(data = avg_LFC_row2col, aes(x = random_index, y = MELJUSO_SKIN_ERKi - MELJUSO_SKIN_DMSO)) +
  geom_point(aes(fill=dir_r),alpha = 0.8,color="grey2",shape=21)+
  scale_fill_manual(values=c("black","firebrick2"))+
  geom_text_repel(data = avg_LFC_row2col%>%
                    filter(pairs %in% c("AAVS1-DUSP4",
                                        "DUSP4-DUSP6",
                                        "AAVS1-DUSP6",
                                        "DUSP1-DUSP4",
                                        "DUSP2-DUSP4",
                                        "DUSP4-DUSP5",
                                        "DUSP6-DUSP7")),
                  aes(label = pairs),
                  size =c(rep(3.5,5),5,3.5),
                  color = c(rep("black",5),"firebrick","black"),
                  fontface = 'bold',
                  show.legend=F,
                  point.padding = unit(1, "lines"),
                  box.padding = unit(0.8, "lines"),
                  segment.color = 'grey50',colour="black") +
  geom_abline(slope = 0,intercept = 0, color="cornflowerblue", linetype="dashed", size=1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ylim(c(-2,4))+
  labs(x="Gene Index",y = "Sensitive <- Sensitizer index -> Resistance",title="DMSO vs SCH984 (ERK inhibitor)")
ggsave(g, filename = paste0("Fig 3h ERKi resistance screen/Fig3h_DMSO_ERKi_SensitivityIndex_ReplicatesMerged_",Sys.Date(),".pdf"), 
       width = 6, height = 6)
