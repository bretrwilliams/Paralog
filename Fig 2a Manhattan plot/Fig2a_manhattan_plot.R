library(gemini)
library(tidyverse)
library(taigr)
library(celllinemapr)
library(mixtools)
require(ggrepel)
require(ggnewscale)

source("Utility_Functions/gene_pair.R")
source("Utility_Functions/gemini_avg_lfc.R")
Model <- readRDS("Data/LatestModel.rds")
Scores <- readRDS("Data/LatestScores.rds")
avg_LFC <- gemini_avg_lfc(Model, LFC_center = "mean")
read_synth_lethal <- function(filename) {
  read_csv(filename) %>%
    mutate(dep_pair = ifelse(Dep_Gene < Loss_Gene, 
                             paste0(gsub(" \\(.*", "", Dep_Gene), ";", gsub(" \\(.*", "", Loss_Gene)),
                             paste0(gsub(" \\(.*", "", Loss_Gene), ";", gsub(" \\(.*", "", Dep_Gene)))) %>%
    mutate(dep_gene1 = ifelse(Dep_Gene < "AAVS1", 
                              paste0(gsub(" \\(.*", "", Dep_Gene), ";", "AAVS1"),
                              paste0("AAVS1", ";", gsub(" \\(.*", "", Dep_Gene)))) %>%
    mutate(dep_gene2 = ifelse(Loss_Gene < "AAVS1", 
                              paste0(gsub(" \\(.*", "", Loss_Gene), ";", "AAVS1"),
                              paste0("AAVS1", ";", gsub(" \\(.*", "", Loss_Gene))))
}
avana_synth_lethal_genetic_hit <- read_synth_lethal("Data/Synthetic lethal paralog hit.csv")
##### ggplots #####
score_df = -log10(Scores$fdr_sensitive_lethality) %>%
	as.data.frame() %>%
	rownames_to_column('genepair') %>%
	tidyr::gather(cell.line, GEMINI.score, -genepair)

lfc_df = avg_LFC %>%
  as.data.frame() %>%
	rownames_to_column('genepair') %>%
	tidyr::gather(cell.line, LFC, -genepair)

t_score = -log10(0.05) # thresholds for FDR

plot_df <- score_df %>% 
  inner_join(lfc_df, by = c("genepair", "cell.line")) %>%
  filter(complete.cases(.)) %>%
  arrange(genepair, cell.line)

PRKCE_members = row.names(avg_LFC)[grepl("PRKCE", row.names(avg_LFC)) & !(grepl("PKN", row.names(avg_LFC))) &
                                    !(grepl("AAVS1", row.names(avg_LFC)))]

cellline_order <- c("HSC5", "MEWO", "MELJUSO", "IPC298", "HS936T", "HS944T", "PATU8988S", "PK1", "A549", "GI1","MEL202")
# top 3 per line
plot_df <- plot_df %>%
  mutate(Top_3 = 0)
for (i in unique(plot_df$cell.line)){
  plot_df$Top_3[plot_df$cell.line == i &
                  plot_df$genepair %in% rownames(Scores$fdr_sensitive_lethality)[order(Scores$fdr_sensitive_lethality[,i])[1:3]]] = 1
}

plot_df %<>%
  mutate(cell.line2 = gsub("_.*", "", cell.line)) %>%
  mutate(highlight = case_when(genepair %in% PRKCE_members~"PRKCE_pairs",
                               genepair=="CDK4;CDK6"~"CDK4_CDK6",
                               Top_3==1~"Top_pairs",
                               TRUE~"others"))

plot_df$cell.line2 <- factor(plot_df$cell.line2, levels = cellline_order)

gg <- ggplot(data = plot_df, aes(x = genepair, y = GEMINI.score, color = LFC, label = genepair)) + 
	facet_grid(~ cell.line2, scales = 'free_x', space = 'free_x', switch = 'x') + 
  scale_y_sqrt(expand = expansion(mult = c(0, 0.01)),   
               sec.axis = dup_axis(breaks = 0),
               labels = function(x){paste(x, "-")})+
  geom_point() + 
  scale_colour_gradientn(colours=c("red", "#CCCCCC", "#333333" )) +
  geom_hline(yintercept = t_score, color = 'black', linetype = "dashed") +
  geom_hline(yintercept = c(-0.1,26), color = 'black') +
  new_scale_color()+
  geom_text_repel(data=.%>%
                    filter(highlight!="others" & GEMINI.score>t_score),aes(color=highlight),
                  size = 4, force = 6,fontface="bold") +
  scale_color_manual(values = c("dodgerblue4", "firebrick1", "grey"),guide = 'none')+
	xlab("") +
	ylab("-log10(GEMINI synergy FDR)") +
  theme_classic()+
	theme(
	  legend.position = c(0.05, 0.98),
	  legend.justification = c("left", "top"),
	  legend.box.just = "left",
	  axis.ticks = element_blank(), 
	  axis.title.y = element_text(size = rel(1.5)),
	  axis.title.y.right = element_blank(),
	  axis.text.x= element_blank(), 
	  axis.text.y.right = element_blank(),                      
	  axis.text.y = element_text(margin = margin(r = 0)), 
	  panel.spacing = unit(0, "mm"),                      
	  strip.background = element_blank(),
		strip.text.x = element_text(size = rel(1.5), colour = "black", angle = 0))
ggsave(gg, filename = paste0("Fig 2a Manhattan plot/Fig2a_ManhattanOfFDR_",Sys.Date(),".pdf"), width = 15, height = 11)
# Background segmentation layer difficult to add