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
#===== Expression versus synergy =====
CCLE.expression <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='CCLE_expression')
sample.info <- load.from.taiga(data.name='public-20q1-c3b6', data.version=10, data.file='sample_info')
expr_cl <- CCLE.expression[sample.info%>%
  filter(CCLE_Name %in% colnames(fdr_Scores))%>%
  pull(DepMap_ID),]%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("gene")%>%
  gather(DepMap_ID,log2tpm,-gene)%>%
  mutate(gene = gsub(" \\(.*","",gene))%>%
  inner_join(sample.info%>%select(stripped_cell_line_name,DepMap_ID),by="DepMap_ID")%>%
  select(gene,log2tpm,Cell.line = stripped_cell_line_name)
fdr_Scores <- readRDS("Data/LatestScores.rds")$fdr_sensitive_lethality
fdr_Scores%<>%replace(is.na(.), 1)%>%.[grep("AAVS1|TRIM",rownames(fdr_Scores),invert = T,value=T),]
expr_profile <- fdr_Scores%>%
  as.data.frame()%>%
  rownames_to_column("pair")%>%
  gather(Cell.line,fdr,-pair)%>%
  mutate(Cell.line = gsub("\\_.*","",Cell.line))%>%
  filter(!grepl("AAVS1|TRIM",pair))%>%
  separate(pair,into = c("symbol1","symbol2"),sep = ";",remove = F)%>%
  mutate(signif = case_when(fdr < 0.05~"yes",
                            fdr>0.2~"no",
                            TRUE~"others"))%>%
  select(-fdr)%>%
  gather(side,gene,-pair,-Cell.line,-signif)%>%
  inner_join(expr_cl,by=c("gene","Cell.line"))%>%
  mutate(exprLevel = case_when(log2tpm < log2(5)~0,TRUE~1))%>%
  group_by(pair,Cell.line,signif)%>%
  summarize(exprno=sum(exprLevel))%>%
  mutate(exprno = case_when(exprno<=1~paste0(exprno," gene"),
                            TRUE~paste0(exprno," genes")))%>%
  group_by(signif,exprno)%>%
  summarise(n=n())%>%
  filter(signif !="others")

expr_profile%>%
  inner_join(expr_profile%>%
               group_by(signif)%>%
               summarize(n_total = sum(n)),by="signif")%>%
  mutate(perc = round(n*100/n_total,2))%>%
  mutate(exprno = factor(exprno))%>%
  ggplot(.,aes(fill=exprno, y=perc, x=signif)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("black","gray70","#db3b33"),name = "mRNA expression in CCLE")+
  scale_y_continuous(limits = c(0,130),breaks = seq(0,100,25))+
  scale_x_discrete(labels = c("Non \n synergistic","Synergistic"))+
  theme(
    plot.title = element_text(color='black', hjust = 0.5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color='black',size=rel(1.5),family = "Helvetica"),
    axis.title = element_text(color='black',size=rel(1.5),family = "Helvetica"),
    legend.position = c(0.05, 0.98),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.text=element_text(size=rel(1.5),family = "Helvetica"),
    legend.direction = "vertical",
    legend.title = element_text(size=rel(1.5),family = "Helvetica"),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.x = unit(0.05, "cm")
  )+labs(x=NULL,y="Percentage")
ggsave(filename = "Fig S3e CCLE expression versus synergy/FigS3e_Hits_ExprCategory_stack.pdf",width=5,height=8)