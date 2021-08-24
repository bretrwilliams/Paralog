library(gemini)
library(tidyverse)
library(taigr)
library(celllinemapr)
library(ggpmisc)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
source("Utility_Functions/ceres_utils.R")
source("Utility_Functions/gene_pair.R")
source("Utility_Functions/gemini_avg_lfc.R")
#===== Data loading and preprocessing =====
Hart_DepMap_essential_list<- readRDS("Data/essential_genes.rds")
Hart_nonessentials <- read.csv("Data/NEGv1.txt",sep = "\t")$Gene
Model <- readRDS("Data/LatestModel.rds")
avg_LFC <- gemini_avg_lfc(Model, LFC_center = "mean")
lfc_long <- avg_LFC %>% 
  as.data.frame() %>%
  rownames_to_column("gene_pair") %>%
  filter(grepl("AAVS1", gene_pair)) %>%
  mutate(gene = gsub("(;|AAVS1)", "", gene_pair)) %>%
  select(-gene_pair) %>%
  gather(ccle_name, lfc, -gene)
Achilles.replicate.map <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='Achilles_replicate_map')
lfc.cls <- Achilles.replicate.map%>%
  filter(grepl("ACH",DepMap_ID))%>%
  mutate(ccle = purrr::pmap_chr(list(DepMap_ID),celllinemapr::arxspan.to.ccle))%>%
  filter(ccle %in% unique(lfc_long$ccle_name))
avana.lfc <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='Achilles_logfold_change')
avana.guide.map <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='Achilles_guide_map')
avana.sub <- avana.guide.map%>%
  mutate(gene =gsub(" \\(.*","",gene))%>%
  filter(gene %in% lfc_long$gene)%>%select(sgrna,gene)
lfc.all <- avana.lfc%>%
  as.data.frame()%>%
  tibble::rownames_to_column("sgrna")%>%
  gather(replicate_ID,lfc,-sgrna)%>%
  inner_join(lfc.cls%>%select(replicate_ID,ccle),by="replicate_ID")%>%
  inner_join(avana.sub,by="sgrna")%>%
  group_by(gene,ccle)%>%
  summarise(lfcAvana = mean(lfc))
lfc.data <- lfc.all%>%
  inner_join(lfc_long, by = c("gene", "ccle" = "ccle_name")) %>%
  separate(ccle, c("cellline", "organ"), sep = "_", extra = "merge", remove = FALSE) %>%
  mutate(Hart = ifelse(gene %in% Hart_DepMap_essential_list$Hart_essentials, "essential", 
                       ifelse(gene %in% Hart_nonessentials, "non-essential", "other")))
#===== Correlation between Avana LFC and AAVS1:GeneX highlighting HartEssentialGene===========
g <- ggplot(lfc.data%>%arrange(order(Hart)), aes(x =  lfcAvana, y = lfc, color = Hart)) + 
  geom_point(size=0.5,data = lfc.data%>%filter(Hart=="other"), alpha=0.75) +
  geom_point(size = 0.5,data = lfc.data%>%filter(Hart!="other"), alpha=0.75)+
  geom_text(data=lfc.data%>%filter(Hart!="other"),
            aes(x = -4.7, y = 2.2, color = NULL, 
                label = paste("r =", signif(cor(lfc, lfcAvana, use = "pairwise.complete"), 3))),
            show.legend = FALSE) +
  geom_smooth(data = lfc.data%>%filter(Hart!="other"),
              color = "#999999", method = 'lm') + theme(text = element_text(size = 14)) +
  scale_color_manual(values=c("#FF9900", "#0066CC", "grey50"),
                     name = "   Gene essentiality", 
                     labels = c("Pan-essential genes", "Non-essential genes","Others"))+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 0.25, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position = c(0.96, 0.1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  labs(y="Log-fold change (sgGeneA : sgAAVS1)",
       x = "AVANA Log-fold change")+
  ylim(c(-5,2.3))+
  geom_hline(aes(yintercept = 0),linetype="dashed")+
  geom_vline(aes(xintercept= 0),linetype="dashed")+
  xlim(c(-5,2.3))
ggsave(g,filename = paste0("Fig S2b Avana_correlation/FigS2b_Corr_Scatter_Avana_20Q1public_LFC_", Sys.Date(),".pdf"),
         width = 7, height = 7)

