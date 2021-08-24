require(dplyr)
require(magrittr)
require(tidyr)
require(purrr)
require(taigr)
require(ggrepel)
require(plotly)
import::from(tibble,"column_to_rownames","rownames_to_column")
options(stringsAsFactors = F)
#===== In-Context Synergy: Volcano plots =====
gemini_avg_lfc <- function(Model) {
  stopifnot("gemini.model" %in% class(Model))
  anno <- Model$Input$guide.pair.annot
  colnames(anno) <- c("key", "g1", "g2")
  Model$Input$LFC %>% 
    as.data.frame() %>%
    rownames_to_column("key") %>%
    gather(cellline, lfc, -key) %>%
    left_join(anno, by = "key") %>%
    mutate(dep_pair = ifelse(g1 < g2, 
                             paste0(gsub(" \\(.*", "", g1), Model$pattern_join, gsub(" \\(.*", "", g2)),
                             paste0(gsub(" \\(.*", "", g2), Model$pattern_join, gsub(" \\(.*", "", g1)))) %>%
    group_by(dep_pair, cellline) %>%
    dplyr::summarize(lfc = mean(lfc)) %>%
    ungroup()%>%
    spread(cellline, lfc) %>%
    column_to_rownames("dep_pair") %>%
    as.matrix()
}
Score <- readRDS("Data/LatestScores.rds")$sensitive_lethality
Score%<>%replace(is.na(.), 0)%>%.[grep("AAVS1|TRIM",rownames(Score),invert = T,value=T),]
Model <- readRDS("Data/LatestModel.rds")
LFC <- gemini_avg_lfc(Model)
LFC%<>%.[rownames(Score),]
CCLE.gene.cn <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='CCLE_gene_cn')
CCLE.mutations <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='CCLE_mutations')
sample.info <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='sample_info')
screenCL <- sample.info%>%
  select(DepMap_ID,CCLE_cl=CCLE_Name)%>%
  filter(CCLE_cl %in% colnames(Score))
cn_data <- screenCL%>%
  inner_join(CCLE.gene.cn[screenCL$DepMap_ID,c("CDKN2A (1029)","MYC (4609)")]%>%
               as.data.frame()%>%
               set_colnames(c("CDKN2A_cn","MYC_cn"))%>%
               rownames_to_column("DepMap_ID"),by="DepMap_ID")
mutation_data <- CCLE.mutations[,-1]%>%
  filter(Hugo_Symbol %in% c("TP53","KRAS","NRAS"))%>%
  select(Hugo_Symbol,isTCGAhotspot,isCOSMIChotspot,DepMap_ID)%>%
  inner_join(cn_data%>%select(DepMap_ID,CCLE_cl),by="DepMap_ID")
mut_knras <- mutation_data%>%filter(Hugo_Symbol %in% c("KRAS","NRAS"))%>%pull(CCLE_cl)%>%unique()
mut_kras <- mutation_data%>%filter(Hugo_Symbol %in% c("KRAS"))%>%pull(CCLE_cl)%>%unique()
mut_tp53 <- mutation_data%>%filter(Hugo_Symbol %in% c("TP53"))%>%pull(CCLE_cl)%>%unique()
deplet_CDKN2A <- cn_data%>%filter(CDKN2A_cn < quantile(CCLE.gene.cn[,"CDKN2A (1029)"],0.3))%>%pull(CCLE_cl)
amplify_MYC <- cn_data%>%filter(MYC_cn > quantile(CCLE.gene.cn[,"MYC (4609)"],0.7))%>%pull(CCLE_cl)
indicators <- list(mut_knras,mut_tp53,deplet_CDKN2A,amplify_MYC)
labels <- c("KRAS+NRAS_mutation","TP53_mutation","CDKN2A_depletion","MYC_amplification")
highSets <- list(c("SEPHS1;SEPHS2","YWHAE;YWHAZ","CDK4;CDK6"),
                 c("USP32;USP4","MDM2;MDM4","USP25;USP28","CTDSP2;CTDSPL2","TBL1X;TBL1XR1"),
                 c("TPTE;TPTE2","HTRA2;HTRA4"),
                 c("ACSL3;ACSL4","MAPK14;MAPK8","MAPK13;MAPK14"))
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
volcano_plot_paralog <- function(Score,LFC,mut,allcl,g_subset=NULL,t=NULL){
  Score%<>%.[,allcl]
  LFC%<>%.[,allcl]
  indicator_gemini <- as.numeric(colnames(Score) %in% mut)
  indicator_LFC <- as.numeric(colnames(LFC) %in% setdiff(allcl,mut))
  res.gemini = run_lm_stats_limma(mat = t(Score), vec = indicator_gemini)
  res.lfc = run_lm_stats_limma(mat = t(LFC), vec = indicator_LFC)
  res.simple.gemini <- res.gemini[,c(1,2,5)] %>% set_colnames(c("Gene","Gemini_EffectSize","Gemini_p.value"))
  res.simple.lfc <- res.lfc[,c(1,2,5)] %>% set_colnames(c("Gene","lfc_EffectSize","lfc_p.value"))
  res.simple.merged <- inner_join(res.simple.gemini, res.simple.lfc, by = "Gene")
  gg = ggplot(data = res.simple.merged, aes(x = lfc_EffectSize, y = -log10(lfc_p.value), label = Gene)) +
    geom_point(aes(size = Gemini_EffectSize, color = Gemini_EffectSize),alpha=0.75) +
    scale_radius(limits = c(-1,1),range = c(1, 5),breaks=c(-1,0,0.5,1),name="GEMINI\neffect size")+
    scale_y_sqrt() +
    scale_color_gradient2(low="gray23", mid = "gray50", high="red",
                          limits = c(-1,1),breaks=c(-1,0,0.5,1),name="GEMINI\neffect size")+
    guides(color= guide_legend(reverse = TRUE), size=guide_legend(reverse = TRUE))+
    geom_text_repel(
      data = subset(res.simple.merged, Gene %in% g_subset),
      aes(label = Gene),
      size = rel(4),
      color = "black",
      box.padding = unit(0.8, "lines"),
      point.padding = unit(0.8, "lines"),
      fontface="bold")+
    labs(title = t,
         x = "LFC effect size", y = "-log10(P value)")+
    xlim(-2.5, 2.5)+
    geom_vline(xintercept = c(-1,1),linetype="dashed",color="gray40")+
    theme(aspect.ratio = 1, 
          panel.border = element_rect(color = "black", size = 0.25, fill = NA),
          panel.grid = element_blank(),
          panel.background = element_rect(color = "transparent", fill = "gray91"), 
          axis.text.x=element_text(colour="black"),
          legend.position = c(0.98, 0.02),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.key = element_rect(fill = NA),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black",fill="transparent"))
  return(list(sigD = res.simple.merged,volcanop = gg))
}
volp_list <- lapply(1:length(labels),function(s){
         p <- volcano_plot_paralog(Score=Score,LFC=LFC,mut=indicators[[s]],allcl=colnames(Score),
                                   g_subset =highSets[[s]], 
                                   t=paste0(gsub("\\_.*","",labels[s])," dependencies"))$volcanop
         ggsave(p,filename = paste0("Fig S4b Volcano plots for in-context synergy/FigS4b_Volcano_",labels[s],"_dependency.pdf"),
                width=6.5,height=6.5)
         return(p)
       })
ggpubr::ggarrange(plotlist = volp_list[c(2,1,4,3)],ncol = 2,nrow=2)
ggsave(filename = paste0("Fig S4b Volcano plots for in-context synergy/FigS4b_CombVolcano_InContext_dependency",Sys.Date(),".pdf"),
       width=13,height=13)