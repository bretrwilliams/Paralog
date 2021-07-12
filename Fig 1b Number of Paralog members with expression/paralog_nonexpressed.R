# Paralog analysis
library(tidyverse)
library(magrittr)
library(dplyr)
library(taigr)

### EXPRESSION
# Define non-expressed genes in CCLE
#CCLE_RNAseq <- load.from.taiga(data.name='internal-19q2-9504', data.version=24, data.file='CCLE_expression', transpose = TRUE)
CCLE_RNAseq <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='CCLE_expression_full', transpose = TRUE)
rownames(CCLE_RNAseq) <- gsub(" \\(.*", "", rownames(CCLE_RNAseq))
colnames(CCLE_RNAseq) <- celllinemapr::arxspan.to.ccle(colnames(CCLE_RNAseq), ignore.problems = TRUE)
#View(data.frame(rownames(CCLE_RNAseq ),CCLE_RNAseq[,grep("MELJUSO_SKIN",colnames(CCLE_RNAseq))]))

# CCLE_RNAseq <- read.delim("~/Documents/Data/CCLE_RNAseq_rsem_genes_tpm_20180929.txt", sep = '\t', stringsAsFactors = F)
# 
# genes <- CCLE_RNAseq[,1]
# genes <- sapply(genes, function(s) strsplit(s, split = '.', fixed = T)[[1]][1])
# 
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# bm <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
# 
# genes.hugo <- bm$hgnc_symbol[match(bm$ensembl_gene_id, genes)]
# 
# new_ids <- bm$hgnc_symbol[match(genes, bm$ensembl_gene_id, nomatch = NA)]
# CCLE_RNAseq[,1] <- new_ids
# CCLE_RNAseq <- CCLE_RNAseq[,-2]

TPM_THRESHOLD = 2
percentExpressed <- apply(CCLE_RNAseq, 1, function(x) sum(x > TPM_THRESHOLD, na.rm = T)/length(x[!is.na(x)]))

percentExpressed_df <-
	data.frame(
		percentage_of_cells = percentExpressed,
		gene = rownames(CCLE_RNAseq),
		stringsAsFactors = F
	)

# Read in paralog groups
paralog_df <- readRDS("~/Dropbox/paralog_AndreasMahdi/Paralog_Paper/Fig1b Number of Paralog members with expression/paralog_list_df.rds")

paralog_lst <- lapply(unique(paralog_df$group), function(x){
	return(unique(paralog_df$gene_symbol[paralog_df$group==x]))
})

paralog_expression_df <- merge(percentExpressed_df, paralog_df, by.y = "gene_symbol", by.x = "gene", all.y  = T)
for(THRESHOLD in seq(0.1, 1, by = 0.1)){
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
	group_summary$n_expressed[group_summary$members_expressed==0] <- "None expressed"
	group_summary$n_expressed[group_summary$members_expressed==1] <- "1 expressed"
	group_summary$n_expressed[group_summary$members_expressed==2] <- "2 expressed"
	group_summary$n_expressed[group_summary$members_expressed>2] <- "3+ expressed"
	
	
	df = group_summary %>%
		count(group_size, n_expressed)
	
	df$group_size <- as.factor(df$group_size)
	df$n_expressed <- factor(df$n_expressed, levels = c("None expressed",
	                                                    "1 expressed",
	                                                    "2 expressed",
	                                                    "3+ expressed"))
	ggplot(data = df, aes(x = group_size, y = n, fill = n_expressed)) +
		geom_bar(stat = 'identity') +
		theme(plot.background = element_blank(), panel.background = element_blank()) +
		ggsave(filename = paste0("~/Dropbox/paralog_AndreasMahdi/Paralog_Paper/Fig1b Number of Paralog members with expression/n_expressed_in_",THRESHOLD*100,"_percent_of_ccle.pdf"), dpi = 320, width = 8, height = 6)
}

#### Load in paralog screen data, identify hits, characterize expression #####
# library(gemini)
# Model <- readRDS("~/Dropbox/GEMINI/Paralog/Gemini_Model/2019-04-12_AllLines_merged_model_sdx_1_priorshape_0.5.rds")
# Scores <- gemini_score(Model, pc_threshold = -1)
# Scores$strong <- subset(Scores$strong, !grepl("TRIM", rownames(Scores$strong)))
# Scores$sensitive_lethality <- subset(Scores$sensitive_lethality, !grepl("TRIM", rownames(Scores$sensitive_lethality)))
# Scores$sensitive_recovery <- subset(Scores$sensitive_recovery, !grepl("TRIM", rownames(Scores$sensitive_recovery)))
# 
# 
# sample_hit_annot <- list()
# for(sample in colnames(Scores$sensitive_lethality)){
# 	x = Scores$sensitive_lethality[,sample]
# 	top50 = rownames(Scores$sensitive_lethality)[order(x, decreasing = T)[1:50]]
# 	hit.annot <- data.frame(
# 		hits = top50,
# 		cell_line = sample,
# 		g = sapply(top50, function(s) strsplit(s, split = ';')[[1]][1]),
# 		h = sapply(top50, function(s) strsplit(s, split = ';')[[1]][2])
# 	)
# 	hit.annot$g_exp_pct <- paralog_expression_df$percentage_of_cells[match(hit.annot$g, paralog_expression_df$gene)]
# 	hit.annot$h_exp_pct <- paralog_expression_df$percentage_of_cells[match(hit.annot$h, paralog_expression_df$gene)]
# 	sample_hit_annot[[sample]] <- hit.annot
# }
# 
# sample_hit_annot_merged <- bind_rows(sample_hit_annot)
# sample_hit_annot_merged$contains_trim <- grepl("TRIM", sample_hit_annot_merged$g) | grepl("TRIM", sample_hit_annot_merged$h)
# 
# ggplot(sample_hit_annot_merged, aes(x = g_exp_pct, y = h_exp_pct, shape = cell_line, color = contains_trim)) +
# 	geom_point(size = 3) +
# 	ggpmisc::stat_quadrant_counts(xintercept = 0.5, yintercept = 0.5, geom = 'label_npc') +
# 	scale_shape_manual(values = c(15, 16, 17, 18, 19,20, 4, 8, 10, 13)) +
# 	ggsave("~/Dropbox/GEMINI/Paralog/Expression/ExpressionPercentageCCLE_withTRIM.pdf", width = 8, height = 6, dpi = 320)
# 
# ggplot(subset(sample_hit_annot_merged, !contains_trim), aes(x = g_exp_pct, y = h_exp_pct, shape = cell_line, color = contains_trim)) +
# 	geom_point(size = 3) +
# 	ggpmisc::stat_quadrant_counts(xintercept = 0.5, yintercept = 0.5, geom = 'label_npc') +
# 	scale_shape_manual(values = c(15, 16, 17, 18, 19,20, 4, 8, 10, 13)) +
# 	ggsave("~/Dropbox/GEMINI/Paralog/Expression/ExpressionPercentageCCLE_withoutTRIM.pdf", width = 8, height = 6, dpi = 320)
# 
