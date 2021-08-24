require(dplyr)
require(magrittr)
require(tidyr)
require(tibble)
require(purrr)
require(taigr)
require(celllinemapr)
require(ggrepel)
require(plotly)
options(stringsAsFactors = F)
#===== Function loading =====
left_right_comparison <- function(D,annot,posctrl,type="single",FUN="mean"){
  colnames(annot) <- c("rowname","symbol1","symbol2")
  if(type=="single"){
    mapping_12 <-annot%>%
      filter(symbol1=="AAVS1"|symbol2=="AAVS1")%>%
      mutate(Side = case_when(symbol1=="AAVS1"~"side2",TRUE~"side1"),
             symbol = case_when(symbol1=="AAVS1"~symbol2,TRUE~symbol1))%>%
      select(rowname,Side,symbol)%>%distinct()
    desireDesign <- mapping_12%>%group_by(symbol)%>%summarise(n=n())%>%filter(n==6)%>%pull(symbol)
    mapping_12%<>%filter(symbol %in% desireDesign)
  }else{
    mapping_12 <-annot%>%
      filter(!(symbol1=="AAVS1"|symbol2=="AAVS1"))%>%
      mutate(symbol = purrr::pmap_chr(list(symbol1,symbol2),~paste0(sort(c(...)),collapse = ";")))%>%
      mutate(Lead_g = gsub(";.*","",symbol))%>%
      mutate(Side = case_when(Lead_g==symbol1~"side1",TRUE~"side2"))%>%
      select(Side,symbol,rowname)
    desireDesign <- mapping_12%>%group_by(symbol)%>%summarise(n=n())%>%filter(n==18)%>%pull(symbol)
    mapping_12%<>%filter(symbol %in% desireDesign)
  }
  compare_data <- lapply(colnames(D),function(cl){
    D[,cl]%>%
      as.data.frame()%>%
      set_colnames(cl)%>%
      mutate(rowname = rownames(D))%>%
      inner_join(mapping_12,by="rowname")%>%
      group_by(symbol,Side)%>%
      summarise_at(vars(cl),FUN,na.rm=T)%>%
      spread(Side,cl)%>%
      na.omit()%>%
      mutate(Essentiality = case_when(symbol %in% posctrl~"yes",TRUE~"no"))%>%
      mutate(source = cl)})%>%bind_rows()
  return(compare_data)
}
left_right_scatter <- function(inputlist,type,FUN,fillColor="skyblue3",rowno=2,comb=F){
  if(type=="single"){
    lr_compare<- lapply(inputlist,function(s) left_right_comparison(D =s$LFC,annot = s$annot,posctrl = s$posctrl_gene,type=type,FUN=FUN)) 
  }else{
    lr_compare<- lapply(inputlist,function(s) left_right_comparison(D =s$LFC,annot = s$annot,posctrl = s$posctrl_pair,type=type,FUN=FUN))  
  }
  lr_compare%<>%bind_rows()%>%as.data.frame()
  min_lim <- floor(apply(lr_compare%>%select(side1,side2),2,min,na.rm=T)%>%min())
  max_lim <- ceiling(apply(lr_compare%>%select(side1,side2),2,max,na.rm=T)%>%max())
  lims <- c(min_lim,max_lim)
  tailedX <-ifelse(type=="single"," of Left guides paired with AAVS1",
                   " of geneX: geneY")
  tailedY <-ifelse(type=="single"," of Right guides paired with AAVS1",
                   " of geneY: geneX")
  corsD <- lr_compare%>%
    group_by(source)%>%
    summarize(values = cor(side1, side2, use = "pairwise.complete",method="pearson"))%>%
    as.data.frame()%>%
    mutate(cors =  paste("r =", signif(values, 2)))%>%
    arrange(values)
  corsD%<>%mutate(
    source = factor(source,levels = 
                      (corsD%>%arrange(values)%>%pull(source)),
                    ordered = TRUE))
  
  lr_compare%<>%mutate(source=factor(source,levels = 
                                       (corsD%>%arrange(values)%>%pull(source)),
                                     ordered = TRUE))
  p <- ggplot(lr_compare, aes(x=side1, y=side2,color=Essentiality)) +
    geom_abline(intercept = 0)+ 
    geom_point(size=1.3,data = .%>%filter(Essentiality=="no"), alpha=0.8) +
    geom_point(size = 1.3,data = .%>%filter(Essentiality=="yes"), alpha=0.8)+
    scale_color_manual(values=c("grey50","red3"),labels = c( "Others","Pan-essential genes"))+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          legend.key = element_rect(fill = NA),
          legend.position ="none")+
    xlim(lims)+ylim(lims)+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    labs(x=paste0(stringi::stri_trans_totitle(FUN),tailedX),
         y = paste0(stringi::stri_trans_totitle(FUN),tailedY))
  if(comb){
    p <- p+
      geom_text(aes(x = min_lim+1.5, y = max_lim-0.5, color = NULL,label=
                      paste("r =", signif( cor(side1, side2, use = "pairwise.complete",method="pearson"), 2))),
                show.legend = FALSE,size=4.5,inherit.aes = F,family="Times")
  }else{
    p <- p+
      facet_wrap(~source,nrow=rowno)+
      geom_text(data = corsD,aes(x = min_lim+1.5, y = max_lim-0.5, color = NULL,label=cors),
                show.legend = FALSE,size=4.5,inherit.aes = F,family="Times")+
      theme(strip.text.x = element_text(
        size = 12, color = "white", family = "Times"),
        strip.background = element_rect(
          color="black", fill=fillColor, size=1, linetype="solid"))
  }
  return(p)}
#===== Pyo-Aur balance =====
Model <- readRDS("Data/LatestModel.rds")
V1list <- list(V1 = list(LFC = Model$Input$LFC%>%
                           set_colnames(gsub("\\_.*","",colnames(Model$Input$LFC))),
                         annot = Model$Input$guide.pair.annot,posctrl_gene=c("")))
g <- left_right_scatter(inputlist=V1list,type="double",FUN="mean",fillColor="skyblue3",comb = T)+
  labs(x="Log-fold change (sgGeneA:sgGeneB)",
       y="Log-fold change (sgGeneB:sgGeneA)")
ggsave(g,filename = "Fig S2c Aureus-Pyogenes average pair LFC/FigS2c_Aur_pyo_compatison_pair_mean_combined.pdf",width=8,height=8)
