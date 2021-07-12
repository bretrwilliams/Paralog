library(gemini)
library(tidyverse)
library(taigr)
library(celllinemapr)
library(mixtools)
library(dplyr)

##### input data #####
# run this in a bash shell (creates a symbolic link): ln -s /your/path/in/Dropbox/to/the/Paralog ~/Paralog
setwd("~/Paralog/")
source("Rfunctions/gene_pair.R")
source("Rfunctions/gemini_avg_lfc.R")

Model <- readRDS("Gemini_Model/LatestModel.rds")
avg_LFC <- gemini_avg_lfc(Model, LFC_center = "mean")
#InferredLFC <- read.csv("~/Sellers Lab Dropbox/Takahiro Ito/NRAS - Paralog screens/Paralog dependency screen V1.0/Paralog screen V1.0/Data csv/InferredLFC.csv", row.names=1)
  
##### Rugplot (paralog pairs vs LFC) for genetic interactions identified from Avana, Hart non-essential and essential genes #####

# read avana data
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

avana_synth_lethal_genetic_hit <- read_synth_lethal("data/Synthetic lethal paralog hit.csv")
Hart_essential <- read.delim("Rfunctions/CEGv2.txt", stringsAsFactors = F) %>% 
  mutate(hart_pair = ifelse(GENE < "AAVS1", 
                            paste0(GENE, ";", "AAVS1"),
                            paste0("AAVS1", ";", GENE))) 
Hart_nonessential <- read.delim("Rfunctions/NEGv1.txt", stringsAsFactors = F) %>% 
  mutate(hart_pair = ifelse(Gene< "AAVS1", 
                            paste0(GENE, ";", "AAVS1"),
                            paste0("AAVS1", ";", Gene))) 
#avana_synth_lethal_genetic_hit <- read_synth_lethal("data/Synthetic lethal paralog all parameter hit.csv")

# dataframe for rugplot
df_rug <- avg_LFC %>%
  as.data.frame() %>%
  rownames_to_column('genepair') %>%
  tidyr::gather(cell.line, LFC, -genepair) %>%
  mutate(avana_lethal_pair = genepair %in% avana_synth_lethal_genetic_hit$dep_pair) %>%
  mutate(avana_lethal_gene = genepair %in% avana_synth_lethal_genetic_hit$dep_gene1 | genepair %in% avana_synth_lethal_genetic_hit$dep_gene2) %>%
  mutate(hart_nonessential = genepair %in% Hart_nonessential$hart_pair) %>%
  mutate(hart_essential = genepair %in% Hart_essential$hart_pair) %>%
  mutate(Index = 1:n()) %>% 
  mutate(avana_lethal_gene_2 = ifelse(avana_lethal_gene, 1,0)) %>%
  mutate(avana_lethal_pair_2 = ifelse(avana_lethal_pair, 2,0)) %>%
  mutate(hart_nonessential_2 = ifelse(hart_nonessential, 3,0)) %>%
  mutate(hart_essential_2 = ifelse(hart_essential, 4,0))


# df_rug <- InferredLFC %>%
#   as.data.frame() %>%
#   rownames_to_column('genepair') %>%
#   tidyr::gather(cell.line, LFC, -genepair) %>%
#   mutate(avana_lethal_pair = genepair %in% avana_synth_lethal_genetic_hit$dep_pair) %>%
#   mutate(avana_lethal_gene = genepair %in% avana_synth_lethal_genetic_hit$dep_gene1 | genepair %in% avana_synth_lethal_genetic_hit$dep_gene2) %>%
#   mutate(hart_nonessential = genepair %in% Hart_nonessential$hart_pair) %>%
#   mutate(hart_essential = genepair %in% Hart_essential$hart_pair) %>%
#   mutate(Index = 1:n()) %>% 
#   mutate(avana_lethal_gene_2 = ifelse(avana_lethal_gene, 1,0)) %>%
#   mutate(avana_lethal_pair_2 = ifelse(avana_lethal_pair, 2,0)) %>%
#   mutate(hart_nonessential_2 = ifelse(hart_nonessential, 3,0)) %>%
#   mutate(hart_essential_2 = ifelse(hart_essential, 4,0))



# rugplot function  
rugplot <- function(dat, feature_x, feature_y, feature_rug1, feature_rug2, feature_rug3, feature_rug4){
  cols <- c("1" = "#003366","2" = "#003366", "3" = "#006600", "4" = "#006600","0" = "white")
  min_y <- min(dat[[feature_y]], na.rm = T)
  dat[[feature_x]] <- factor(dat[[feature_x]], levels = dat[[feature_x]][order(dat[[feature_y]], decreasing = F)])
  
  gg = ggplot(data = dat, aes_string(x = feature_x, y = feature_y)) +
    geom_point(color = "black", size = 0.5, alpha = .25) +
    geom_segment(data = dat, mapping = aes(x = dat[[feature_x]], xend = dat[[feature_x]],
                                           y = min_y - 0.1*(as.numeric(dat[[feature_rug1]])),
                                           yend = min_y - 0.3*(as.numeric(dat[[feature_rug1]])),
                                           color = factor(dat[[feature_rug1]])), alpha = 0.25,inherit.aes = F) +
    geom_segment(data = dat, mapping = aes(x = dat[[feature_x]], xend = dat[[feature_x]],
                                           y = min_y - 0.2*(as.numeric(dat[[feature_rug2]])),
                                           yend = min_y - 0.3*(as.numeric(dat[[feature_rug2]])),
                                           color = factor(dat[[feature_rug2]])), alpha = 0.25, inherit.aes = F) +
    geom_segment(data = dat, mapping = aes(x = dat[[feature_x]], xend = dat[[feature_x]],
                                           y = min_y - 0.7/3*(as.numeric(dat[[feature_rug3]])),
                                           yend = min_y - 0.3*(as.numeric(dat[[feature_rug3]])),
                                           color = factor(dat[[feature_rug3]])), alpha = 0.25, inherit.aes = F) +
    geom_segment(data = dat, mapping = aes(x = dat[[feature_x]], xend = dat[[feature_x]],
                                           y = min_y - 1/4*(as.numeric(dat[[feature_rug4]])),
                                           yend = min_y - 0.3*(as.numeric(dat[[feature_rug4]])),
                                           color = factor(dat[[feature_rug4]])), alpha = 0.25, inherit.aes = F) +
   scale_color_manual(name="Type",values = cols, breaks=c("0","1","2", "3", "4"),
                      labels = c("others", "Synthetic lethal pair with sgAAVS1", "Synthetic lethal pair dKO", "Non-essential", "Essential"),limits=c("1","2","3","4")) +
    scale_x_discrete(expand = c(0.05,0.05)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          # panel.border = element_rect(color = "black", size = 1, fill = NA),
          text = element_text(color='black', size = 18),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'))
  
  return(gg)
}

# plot and save
gg <- rugplot(dat = df_rug, feature_x = "Index", feature_y = "LFC",
              feature_rug1 = "avana_lethal_gene_2",
              feature_rug2 = "avana_lethal_pair_2",
              feature_rug3 = "hart_nonessential_2",
              feature_rug4 = "hart_essential_2")
ggsave(gg, filename = paste0("~/Paralog/RugPlot/AggregatedAcrossAllLines_inferredLFC_",Sys.Date(),".pdf"), width = 14, height = 8)

##### define statistical diff between the four rows and our null  #####
# # non cell line specific pairs using CCLE expression
# CCLE_RNAseq <- load.from.taiga(data.name='internal-19q3-4c37', data.version=8, data.file='CCLE_expression_full') %>%
#   as.data.frame(stringsAsFactors = FALSE) %>%
#   rownames_to_column("arxspan.id") %>%
#   gather(gene, log2tpm, -arxspan.id) %>%
#   mutate(cellline = celllinemapr::arxspan.to.ccle(arxspan.id, ignore.problems = TRUE)) %>%
#   #filter(cellline %in% Model$Input$sample.annot[, Model$sample.column.name]) %>%
#   dplyr::select(-arxspan.id) %>%
#   mutate(symbol = gsub(" \\(.*", "", gene)) %>%
#   dplyr::select(-gene)
# 
# LOG2_TPM_ABS_THRESHOLD = 0.1
# percentNotExpressed_df <- CCLE_RNAseq %>%
#   group_by(symbol) %>%
#   summarise(percentage_of_cell = sum(log2tpm < LOG2_TPM_ABS_THRESHOLD)/n())
# 
# Tr = 0.8
# percentNotExpressed_df_Tr <- percentNotExpressed_df %>% filter(percentage_of_cell > Tr)
# NotExpressed_Paralog <- sapply(rownames(Model$s), function(x){
#   g = strsplit(x,split = ";")[[1]][1]
#   h = strsplit(x,split = ";")[[1]][2]
#   return(ifelse(g %in% percentNotExpressed_df_Tr$symbol & h %in% percentNotExpressed_df_Tr$symbol &
#                   !grepl("TRIM",g) & !grepl("TRIM",h), TRUE, FALSE))
# })
# NotExpressed_Paralog <- names(NotExpressed_Paralog)[NotExpressed_Paralog==T]

# calculating the p-values
t.test(df_rug$LFC[df_rug$avana_lethal_pair], df_rug$LFC[df_rug$avana_lethal_gene], alternative = "less")
wilcox.test(df_rug$LFC[df_rug$avana_lethal_pair], df_rug$LFC[df_rug$avana_lethal_gene])

t.test(df_rug$LFC[df_rug$hart_essential], df_rug$LFC[df_rug$hart_nonessential], alternative = "less")
wilcox.test(df_rug$LFC[df_rug$hart_essential], df_rug$LFC[df_rug$hart_nonessential])





  
##### Rugplot (LFC density) for genetic interactions identified from Avan #####
# rugplot_density <- function(dat, feature_x, feature_y, feature_rug1, feature_rug2){
#   cols <- c("1" = "tomato4","2" = "lightseagreen","0" = "white")
#   min_y <- 0
#   dat[[feature_x]] <- factor(dat[[feature_x]], levels = dat[[feature_x]][order(dat[[feature_y]], decreasing = F)])
#   
#   gg = ggplot(data = dat, aes_string(x = feature_y)) +
#     geom_density(color = "black", size = 2.5) +
#     geom_segment(data = dat, mapping = aes(x = dat[[feature_y]], xend = dat[[feature_y]],
#                                            y = min_y - 0.1*(as.numeric(dat[[feature_rug1]])),
#                                            yend = min_y - 0.3*(as.numeric(dat[[feature_rug1]])),
#                                            color = factor(dat[[feature_rug1]])), inherit.aes = F) +
#     geom_segment(data = dat, mapping = aes(x = dat[[feature_y]], xend = dat[[feature_y]],
#                                            y = min_y - 0.2*(as.numeric(dat[[feature_rug2]])),
#                                            yend = min_y - 0.3*(as.numeric(dat[[feature_rug2]])),
#                                            color = factor(dat[[feature_rug2]])), inherit.aes = F) +
#     scale_color_manual(name="Type",values = cols, breaks=c("0","1","2"),
#                        labels = c("others", "avana_lethal_gene", "avana_lethal_pair"),limits=c("1","2")) +
#     scale_x_discrete(expand = c(0.05,0.05)) +
#     theme(plot.title = element_text(color='black', hjust = 0.5),
#           plot.background = element_blank(),
#           panel.background = element_blank(),
#           panel.border = element_rect(color = "black", size = 1, fill = NA),
#           text = element_text(color='black', size = 18),
#           panel.grid = element_blank(),
#           axis.text = element_text(color='black'))
#   
#   return(gg)
# }
# 
# # plot and save
# gg <- rugplot_density(dat = df_rug, feature_x = "Index", feature_y = "LFC",
#               feature_rug1 = "avana_lethal_gene_2",
#               feature_rug2 = "avana_lethal_pair_2")
