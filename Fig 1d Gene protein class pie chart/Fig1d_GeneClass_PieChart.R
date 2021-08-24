library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
geneV1 <- readRDS("Data/Genes_V1.rds")
annot <- readRDS("Data/GeneralAnnot.rds")
#PieDonut(acs,aes(Dx,smoking),r1=0.7,r0=0.3,showPieName=F)
#===== Updated Analysis =====
spitAnnot <- annot%>%
  filter(grepl("\\|",protein_class))%>%
  separate(protein_class,sep="\\|",into=c("v1","v2","v3"),remove=F)%>%
  select(symbol,"v1","v2","v3")%>%
  gather(vs,General_Term,-symbol)%>%
  na.omit()%>%
  select(-vs)%>%
  mutate(General_Term = stringr::str_to_title(General_Term))
enriched_annot<- annot%>%
  filter(!grepl("\\|",protein_class))%>%
  select(symbol,General_Term)%>%
  distinct()%>%
  mutate(General_Term = stringr::str_to_title(General_Term))%>%
  rbind(spitAnnot)%>%
  unique()%>%
  na.omit()
supp_annot <- geneV1%>%
  filter(!Genes %in% enriched_annot$symbol)%>%
  filter(Family == "Achilles Paralog Dependent")%>%
  select(symbol = Genes,General_Term = Family)%>%
  mutate(General_Term = stringr::str_to_title(General_Term))
reannot <- enriched_annot%>%
  rbind(.,supp_annot)%>%
  mutate(Class= case_when(General_Term %in% c("Kinase","Methyltransferase","Acetyltransferase","Transferase")~"Transferase",
                          General_Term %in% c("Phosphatase","Protease","Hydrolase","Lipase")~"Hydrolase",
                          General_Term %in% c("Ubiquitin","Ligase")~"Ligase",
                          General_Term=="Achilles Paralog Dependent"~"Other",
                          TRUE~"NC"))%>%
  filter(Class!="NC")%>%
  mutate(g1in = case_when(symbol %in% geneV1$Genes~'Included',TRUE~"n"))
require(webr)
pdf("Fig 1d Gene protein class pie chart/Fig1d_GeneProteinClass_DocutPie.pdf",width = 5,height = 5)
reannot%>%
  group_by(Class,g1in)%>%
  summarize(n=n())%>%
  mutate(Class =factor(Class,levels = c("Ligase","Hydrolase","Transferase","Other")))%>%
  PieDonut(.,aes(pies=Class,donuts=g1in,count=n),
           ratioByGroup=T,showPieName=FALSE,
           r0=0.4,r1=0.8,r2=1.2,
           showRatioPie = F,
           pieLabelSize = 2.5,
           donutLabelSize = 3,
           color = "transparent")
dev.off()