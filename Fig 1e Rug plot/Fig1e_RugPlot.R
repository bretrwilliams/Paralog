library(gemini)
library(tidyverse)
library(taigr)
library(celllinemapr)
library(mixtools)
library(dplyr)

##### input data #####
source("Utility_Functions/gene_pair.R")
source("Utility_Functions/gemini_avg_lfc.R")
Model <- readRDS("Data/LatestModel.rds")
avg_LFC <- gemini_avg_lfc(Model, LFC_center = "mean")
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
avana_synth_lethal_genetic_hit <- read_synth_lethal("Data/Synthetic lethal paralog hit.csv")
Hart_essential <- read.delim("Data/CEGv2.txt", stringsAsFactors = F) %>% 
  mutate(hart_pair = ifelse(GENE < "AAVS1", 
                            paste0(GENE, ";", "AAVS1"),
                            paste0("AAVS1", ";", GENE))) 
Hart_nonessential <- read.delim("Data/NEGv1.txt", stringsAsFactors = F) %>% 
  mutate(hart_pair = ifelse(Gene< "AAVS1", 
                            paste0(GENE, ";", "AAVS1"),
                            paste0("AAVS1", ";", Gene))) 
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
  mutate(avana_lethal_gene_2 = ifelse(avana_lethal_gene, 3,0)) %>%
  mutate(avana_lethal_pair_2 = ifelse(avana_lethal_pair, 4,0)) %>%
  mutate(hart_nonessential_2 = ifelse(hart_nonessential, 1,0)) %>%
  mutate(hart_essential_2 = ifelse(hart_essential, 2,0))
# rugplot function  
rugplot <- function(dat, feature_x, feature_y, feature_rug1, feature_rug2, feature_rug3, feature_rug4){
  cols <- c("0" = "white", "1" = "#006600", "2" = "#006600","3" = "#003366","4" = "#003366")
  min_y <- min(dat[[feature_y]], na.rm = T)
  dat[[feature_x]] <- factor(dat[[feature_x]], levels = dat[[feature_x]][order(dat[[feature_y]], decreasing = F)])
  g1 = ggplot(data = dat, aes_string(x = feature_x, y = feature_y)) +
    geom_point(color = "black", size = 0.5, alpha = .25) +
    geom_hline(yintercept = 0,linetype="dashed",color="grey50")+
    theme(axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = rel(0.5), fill = NA),
          panel.grid = element_blank(),
          legend.position = "none")+
    labs(title="11 cell lines")+
    ylim(c(-4.1,2.1))
  rugfs <- c(feature_rug1,feature_rug2,feature_rug3,feature_rug4)
  shifts <- c(0.1,0.2,0.7/3,0.25)
  rugsD <- lapply(c(1:4),function(s){
      dat%>%
        mutate(xloc = dat[[feature_x]])%>%
        mutate(yloc1 = min_y - (shifts[s])*(as.numeric(dat[[rugfs[s]]])))%>%
        mutate(yloc2 = min_y - 0.3*(as.numeric(dat[[rugfs[s]]])))%>%
        mutate(col = factor(dat[[rugfs[s]]]))%>%
        select(xloc,yloc1,yloc2,col)
    })%>%bind_rows()
  g2 <- ggplot() +
    geom_segment(data = rugsD, mapping = aes(x = xloc, xend = xloc,y = yloc1,yend = yloc2,color =col), alpha = 0.4,inherit.aes = F) +
    scale_color_manual(values = cols, breaks=c("0","1","2", "3", "4"),limits=c("1","2","3","4"))+
    scale_x_discrete(expand = expansion(mult = c(0, 0)))+
    scale_y_continuous(expand = expansion(mult = c(0, 0)),limits = c(min(rugsD$yloc2),min_y - 0.1),
                       breaks = c(-4.75,-4.45,-4.15,-3.86),
                       labels = c("Paralog lethal dKO",
                                  "Paralog lethal-sgAAVS1",
                                  "Pan-essential-sgAAVS1",
                                  "Non-essential-sgAAVS1"))+
    geom_hline(yintercept = sort(union(rugsD$yloc1,rugsD$yloc2))[-c(1,10)],color="black",size=rel(0.3))+
    theme(axis.text.x=element_blank(), 
          axis.ticks=element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = rel(0.5), fill = NA),
          text = element_text(color='black', size = rel(4)),
          panel.grid = element_blank(),
          legend.position = "none")+
    labs(x="Depleted          <------    ------>             Enriched",y="")
  return(list(g1,g2))
}
# plot and save
gg <- rugplot(dat = df_rug, feature_x = "Index", feature_y = "LFC",
              feature_rug1 = "hart_nonessential_2",
              feature_rug2 = "hart_essential_2",
              feature_rug3 = "avana_lethal_gene_2",
              feature_rug4 = "avana_lethal_pair_2")
ggsave(gg[[1]], filename = paste0("Fig 1e Rug plot/Fig1e_Waterfall_LFC_11cls",Sys.Date(),".pdf"), width = 6, height = 5)
ggsave(gg[[2]], filename = paste0("Fig 1e Rug plot/Fig1e_RugMatch_Waterfall_11cls",Sys.Date(),".pdf"), width = 9, height = 5)
# calculating the p-values from one-sided t-test
t.test(df_rug$LFC[df_rug$avana_lethal_pair], df_rug$LFC[df_rug$avana_lethal_gene], alternative = "less")
t.test(df_rug$LFC[df_rug$hart_essential], df_rug$LFC[df_rug$hart_nonessential], alternative = "less")
