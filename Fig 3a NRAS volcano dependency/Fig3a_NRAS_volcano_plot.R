##### required libraries #####
library(limma)
library(magrittr)
library(tibble)
library(plyr)
library(dplyr)
library(plotly)
library(ggplot2)
library(gemini)
library(ggrepel)
library(ggpubr)
library(tidyverse)
library(taigr)
library(celllinemapr)
library(mixtools)

##### source codes #####
source("Utility_Functions/gene_pair.R")
source("Utility_Functions/gemini_avg_lfc.R")

##### LFC and GEMINI model and GEMINI scores #####
Model <- readRDS("Data/LatestModel.rds")
Score <- readRDS("Data/LatestScores.rds")
LFC <- gemini_avg_lfc(Model, LFC_center = "mean")
##### removing TRIM family #####
Score <- lapply(Score, function(x){
  if (is.matrix(x)){
    return(x[!grepl("TRIM", rownames(x)),])
  } else {
    return(x)
  }
})
Sensitive_filter = Score$sensitive_lethality
Strong_filter = Score$strong
##### cell line selection #####
Volcanocells <- c("HS944T_SKIN", "IPC298_SKIN", "MELJUSO_SKIN","HS936T_SKIN", "GI1_CENTRAL_NERVOUS_SYSTEM","HSC5_SKIN") 
Sensitive_filter <- Sensitive_filter[,colnames(Sensitive_filter) %in% Volcanocells]
Sensitive_filter[is.na(Sensitive_filter)] <- 0
LFC_filter <- LFC[!grepl("AAVS1", rownames(LFC)),]
LFC_filter <- LFC_filter[,colnames(LFC_filter) %in% Volcanocells]
LFC_filter <- LFC_filter[match(rownames(Sensitive_filter),rownames(LFC_filter)),]
##### limma #####
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL, target_type = 'Gene') {
	require(limma)
	require(magrittr)
	require(tibble)
	require(plyr)
	require(dplyr)
	
	udata <- which(!is.na(vec))
	if (!is.numeric(vec)) {
		pred <- factor(vec[udata])
		stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
		n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
		n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
		min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
	} else {
		pred <- vec[udata]
		min_samples <- colSums(!is.na(mat[udata,]))
	}
#there must be more than one unique value of the independent variable
	if (length(unique(pred)) <= 1) {
		return(NULL)
	}
#if using covariates add them as additional predictors to the model
	if (!is.null(covars)) {
		if (!is.data.frame(covars)) {
			covars <- data.frame(covars)
		}
		combined <- covars[udata,, drop = FALSE]
		combined[['pred']] <- pred
		form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
		design <- model.matrix(form, combined)
		design <- design[, colSums(design) != 0, drop = FALSE]
	} else {
		design <- model.matrix(~pred)
	}
	if (!is.null(weights)) {
		if (is.matrix(weights)) {
			weights <- t(weights[udata,])
		} else{
			weights <- weights[udata]
		}
	}
	fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
	fit <- limma::eBayes(fit)
	targ_coef <- grep('pred', colnames(design), value = TRUE)
	results <- limma::topTable(fit, coef = targ_coef, number = Inf)
	
	if (colnames(results)[1] == 'ID') {
		colnames(results)[1] <- target_type
	} else {
		results %<>% rownames_to_column(var = target_type)
	}
	results$min_samples <- min_samples[results[[target_type]]]
	
	two_to_one_sided <- function(two_sided_p, stat, test_dir) {
	
#helper function for converting two-sided p-values to one-sided p-values
		one_sided_p <- two_sided_p / 2
		if (test_dir == 'right') {
			one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
		} else {
			one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
		}
		return(one_sided_p)
	}
	results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
													 'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
	results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
							   p.right = two_to_one_sided(p.value, EffectSize, 'right'),
							   q.left = two_to_one_sided(q.value, EffectSize, 'left'),
							   q.right = two_to_one_sided(q.value, EffectSize,'right'))
	return(results)
}
##### Identification of dependency for NRAS #####
NRASmut<- c("HS944T_SKIN", "IPC298_SKIN", "MELJUSO_SKIN", "HS936T_SKIN")
NRASwt <- c("GI1_CENTRAL_NERVOUS_SYSTEM","HSC5_SKIN")
indicator_gemini <- as.numeric(colnames(Sensitive_filter) %in% NRASmut)
indicator_LFC <- as.numeric(colnames(LFC_filter) %in% NRASwt)
res.gemini = run_lm_stats_limma(mat = t(Sensitive_filter), vec = indicator_gemini)
res.lfc = run_lm_stats_limma(mat = t(LFC_filter), vec = indicator_LFC )
res.simple.gemini <- res.gemini[,c(1,2,5)] %>% set_colnames(c("Gene","Gemini_EffectSize","Gemini_p.value"))
res.simple.lfc <- res.lfc[,c(1,2,5)] %>% set_colnames(c("Gene","lfc_EffectSize","lfc_p.value"))
res.simple.merged <- inner_join(res.simple.gemini, res.simple.lfc, by = "Gene")%>%
  mutate(textcol = case_when(Gene %in% c("YWHAE;YWHAZ","DUSP4;DUSP6")~1,TRUE~0))
##### ggplot #####
gg = ggplot(data = res.simple.merged, 
            aes(x = lfc_EffectSize, y = -log10(lfc_p.value), label = Gene)) +
  geom_point(aes(size = Gemini_EffectSize, color = Gemini_EffectSize),alpha=0.75) +
  scale_radius(limits = c(-1,1),range = c(1, 5),breaks=c(-1,0,0.5,1),name="GEMINI\neffect size")+
  scale_y_sqrt() +
  scale_color_gradient2(low="gray23", mid = "gray50", high="red",
                        limits = c(-1,1),breaks=c(-1,0,0.5,1),name="GEMINI\neffect size")+
  guides(color= guide_legend(reverse = TRUE), size=guide_legend(reverse = TRUE))+
  geom_text_repel(
    data = subset(res.simple.merged, Gene %in% c("YWHAE;YWHAZ", "BRAF;RAF1","CDK4;CDK6","DUSP4;DUSP6")),
    aes(label = Gene),
    size = c(rel(5),rel(5),rel(3),rel(3)),
    color = c("tomato","tomato","black","black"),
    box.padding = unit(0.8, "lines"),
    point.padding = unit(0.8, "lines"),
    fontface="bold")+
    labs(title = "NRAS-mutant dependency",
           x = "Enriched NRAS mut <- LFC effect size -> Depleted NRAS mut", y = "-log10(P value)")+
    xlim(-2.5, 2.5)+
    theme(aspect.ratio = 1, 
          panel.border = element_rect(color = "black", size = 0.25, fill = NA),
          panel.grid = element_blank(),
          panel.background = element_blank(), 
          axis.text.x=element_text(colour="black"),
          legend.position = c(0.98, 0.02),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.key = element_rect(fill = NA),
          legend.box.background = element_rect(colour = "black"))
ggsave(gg, filename = paste0("Fig 3a NRAS volcano dependency/Fig3a_NRAS_volcano_Plot.pdf"), width = 6, height = 6)

