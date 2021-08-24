# Paralog analysis
library(tidyverse)
library(magrittr)
library(dplyr)
library(taigr)
# Install taigr with https://github.com/broadinstitute/taigr
### EXPRESSION
# Define non-expressed genes in CCLE
CCLE_RNAseq <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='CCLE_expression_full', transpose = TRUE)
rownames(CCLE_RNAseq) <- gsub(" \\(.*", "", rownames(CCLE_RNAseq))
TPM_THRESHOLD = 2
percentExpressed <- apply(CCLE_RNAseq, 1, function(x) sum(x > TPM_THRESHOLD, na.rm = T)/length(x[!is.na(x)]))
percentExpressed_df <-
	data.frame(
		percentage_of_cells = percentExpressed,
		gene = rownames(CCLE_RNAseq),
		stringsAsFactors = F
	)
# Read in paralog groups
paralog_df <- readRDS("Data/paralog_list_df.rds")
paralog_lst <- lapply(unique(paralog_df$group), function(x){
	return(unique(paralog_df$gene_symbol[paralog_df$group==x]))
})
paralog_expression_df <- merge(percentExpressed_df, paralog_df, by.y = "gene_symbol", by.x = "gene", all.y  = T)
THRESHOLD <- 0.5
group_n_expressed = sapply(unique(paralog_expression_df$group), function(x){
  ccle_exp = paralog_expression_df$percentage_of_cells[paralog_expression_df$group==x]
  percent = sum(ccle_exp > THRESHOLD, na.rm = T)
  return(percent)
})
group_summary <- data.frame(
  members_expressed = group_n_expressed,
  group_size = sapply(paralog_lst[unique(paralog_expression_df$group)], length),
  group = unique(paralog_expression_df$group),
  stringsAsFactors = F
)
group_summary$n_expressed <- ""
group_summary$n_expressed[group_summary$members_expressed==0] <- "0 gene"
group_summary$n_expressed[group_summary$members_expressed==1] <- "1 gene"
group_summary$n_expressed[group_summary$members_expressed==2] <- "2 genes"
group_summary$n_expressed[group_summary$members_expressed>2] <- ">2 genes"
df = group_summary %>%
  count(group_size, n_expressed)
df$group_size <- as.factor(df$group_size)
df$n_expressed <- factor(df$n_expressed, levels = c("0 gene",
                                                    "1 gene",
                                                    "2 genes",
                                                    ">2 genes"))
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}
reduceLabel <- function(x,l){
  x <- sort(unique(x))
  x <- as.numeric(as.character(x))
  seqFlags <- x >l & x< max(x)
  x <- as.character(x)
  x[seqFlags] <- "."
  x
}
op <- ggplot(data = df, aes(x = group_size, y = n, fill = n_expressed)) +
  geom_col(colour = "black",key_glyph = "polygon3",size=rel(0.2))
sp <- op+
  scale_fill_manual(values = c("#535355","#808083","#e5e6e6","#db3b33"),name = "mRNA expression in CCLE")+
  scale_x_discrete(labels = reduceLabel(df$group_size,l=10))+
  scale_y_continuous(limits = c(0,1600), expand = expansion(mult = c(0, .1))) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = rel(0.5), fill = NA),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=rel(1.5),family = "Helvetica"),
        legend.position = c(0.15, 0.98),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text=element_text(size=rel(1.5),family = "Helvetica"),
        legend.direction = "vertical",
        legend.title = element_text(size=rel(1.5),family = "Helvetica"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.8, "cm"),
        legend.spacing.x = unit(0.05, "cm")
        )+
  labs(x="Number of genes per paralog family",y = " Number of \n paralog families")
ggsave(sp,filename = paste0("Fig 1b Number of Paralog members with expression/Fig1b_n_expressed_in_",
                         THRESHOLD*100,"_percent_of_ccle.pdf"), dpi = 1200, width = 8, height = 3)
