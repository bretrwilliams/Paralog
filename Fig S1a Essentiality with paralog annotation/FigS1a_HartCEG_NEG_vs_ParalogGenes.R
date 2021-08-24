library(tidyverse)
library(taigr)
library(UpSetR)

Achilles.gene.effect <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='Achilles_gene_effect')
Avana_genes <- data.frame(gene = colnames(Achilles.gene.effect)) %>%
  separate(gene, c("symbol", "geneid"), sep = " ") %>%
  mutate(geneid = as.numeric(gsub("(\\(|\\))", "", geneid))) %>%
  unique()
hugo <- read_tsv("Data/all_HUGO.txt") %>%
  filter(Status == "Approved") %>%
  select(geneid = `NCBI Gene ID(supplied by NCBI)`, HGNC_ID = `HGNC ID`)
CEG <- read_tsv("Data/CEGv2.txt") %>% inner_join(hugo, by = "HGNC_ID") %>% filter(geneid %in% Avana_genes$geneid)
NEG <- read_tsv("Data/NEGv1.txt") %>% rename(geneid = ENTREZ_ID) %>% filter(geneid %in% Avana_genes$geneid)
paralog_raw <- read_csv("Data/ensembl-paralogs_v2-ensembl-paralogs.csv")
paralog_genes_split <-  paralog_raw %>%
  filter(Percent_Sequence_Identity > 00) %>%
  separate(Query_ID, c("symbol_q", "geneid_q"), sep = " ") %>%
  mutate(geneid_q = as.numeric(gsub("(\\(|\\))", "", geneid_q))) %>%
  separate(Paralog_ID, c("symbol_p", "geneid_p"), sep = " ") %>%
  mutate(geneid_p = as.numeric(gsub("(\\(|\\))", "", geneid_p)))
paralog_genes <- data.frame(geneid = c(paralog_genes_split$geneid_q, paralog_genes_split$geneid_p),
                            symbol = c(paralog_genes_split$symbol_q, paralog_genes_split$symbol_p)) %>%
  unique() %>%
  filter(geneid %in% Avana_genes$geneid)
non_paralog_genes <- Avana_genes %>%
  filter(!geneid %in% paralog_genes$geneid)
UpSetR::upset(
  fromList(
    list(CEG = CEG$geneid, 
         NEG = NEG$geneid,
         Avana_genes = Avana_genes$geneid,
         paralog_genes = paralog_genes$geneid,
         nonparalog_genes = non_paralog_genes$geneid)
  ),
  text.scale = c(1.3, 1.0, 1, 1, 1.5, 1.5)
)
numerator1_CEG <- length(intersect(CEG$geneid, paralog_genes$geneid)) / length(paralog_genes$geneid)
numerator2_CEG <- length(intersect(CEG$geneid, non_paralog_genes$geneid)) / length(non_paralog_genes$geneid)
denominator_CEG <- numerator1_CEG + numerator2_CEG
numerator1_NEG <- length(intersect(NEG$geneid, paralog_genes$geneid)) / length(paralog_genes$geneid)
numerator2_NEG <- length(intersect(NEG$geneid, non_paralog_genes$geneid)) / length(non_paralog_genes$geneid)
denominator_NEG <- numerator1_NEG + numerator2_NEG
fraction_data <- bind_rows(
  data.frame(gene_set = "essential genes", 
             basis = c("paralog genes", "non paralog genes"),
             percent = 100 * c(numerator1_CEG/denominator_CEG, numerator2_CEG/denominator_CEG)),
  data.frame(gene_set = "non essential genes", 
             basis = c("paralog genes", "non paralog genes"),
             percent = 100 * c(numerator1_NEG/denominator_NEG, numerator2_NEG/denominator_NEG))
)

p <- ggplot(fraction_data%>%mutate(gene_set = case_when(gene_set=="non essential genes"~"Non",TRUE~"Pan")
                                     %>%factor(.,levels=c("Non","Pan"))), aes(x = gene_set, y = percent, fill = basis)) + 
  geom_bar(stat = "identity") + 
  ylab("Percent of genes") + 
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)),breaks = c(0,25,50,75,100),limits = c(0,140))+
  scale_fill_manual(values = c("dodgerblue4","firebrick"),labels = c("No paralog","Paralog"),name="") +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = c(0.02, 0.95),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6),
  )+
  guides(fill= guide_legend(reverse = TRUE))

ggsave(p, filename = paste0(
  "Fig S1a Essentiality with paralog annotation/",
  "FigS1a_CEG_NEG_paralog_non_paralog_",Sys.Date(),".pdf"), width = 3, height = 6)
