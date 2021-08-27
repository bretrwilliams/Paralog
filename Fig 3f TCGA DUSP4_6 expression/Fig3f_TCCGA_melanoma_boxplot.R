library(readr)
library(tidyverse)
library(data.table)

#===== Data Loading =====
#clinical
SKCM_clinical_raw<- read_delim("Data/data_clinical_patient.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(SKCM_clinical_raw) <- SKCM_clinical_raw[4,]#set colname to a easier to operate format
SKCM_clinical <-SKCM_clinical_raw[5:(nrow(SKCM_clinical_raw)),]#remove extra header

#mutation
SKCM_mutation_raw <- read_delim("Data/data_mutations_extended.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
SKCM_mutation_raw_cleaned <- SKCM_mutation_raw[,c("Hugo_Symbol","Variant_Classification","HGVSp_Short","IMPACT","Tumor_Sample_Barcode","Protein_position")]

#create lists of TCGA patients with different mutations
BRAF_mut_tumor <- SKCM_mutation_raw_cleaned %>% filter(Hugo_Symbol == "BRAF", Variant_Classification =="Missense_Mutation", IMPACT %in% c("MODERATE","HIGH"), Protein_position %in% c(600,601,599)) #select tumors with BRAF mutations
NRAS_mut_tumor <- SKCM_mutation_raw_cleaned %>% filter(Hugo_Symbol == "NRAS", Variant_Classification =="Missense_Mutation", IMPACT %in% c("MODERATE","HIGH"),Protein_position %in% c(12,13,61)) #select tumors with NRAS mutations
NF1_mut_tumor <- SKCM_mutation_raw_cleaned %>% filter(Hugo_Symbol == "NF1", Variant_Classification =="Missense_Mutation", IMPACT %in% c("MODERATE","HIGH")) #select tumors with NRAS mutations
mutation_yes <- unique(SKCM_mutation_raw$Tumor_Sample_Barcode) # filter samples with mutation data
normal_barcode <- unique(SKCM_mutation_raw$Matched_Norm_Sample_Barcode)

#RNA seq
SKCM_RNA_seq_raw <- read_delim("Data/data_RNA_Seq_v2_expression_median.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

SKCM_RNA_seq_raw_noNA <- SKCM_RNA_seq_raw %>% drop_na(Hugo_Symbol)
SKCM_RNA_seq_raw_noNA_traw<- t(SKCM_RNA_seq_raw_noNA)
colnames(SKCM_RNA_seq_raw_noNA_traw) <-  SKCM_RNA_seq_raw_noNA_traw[1,]
SKCM_RNA_seq_0<-SKCM_RNA_seq_raw_noNA_traw[3:(nrow(SKCM_RNA_seq_raw_noNA_traw)),]#remove extra header
SKCM_RNA_seq <- as.data.frame(SKCM_RNA_seq_0) %>%rownames_to_column(var = "ID")

#RPPA
SKCM_RPPA_raw <- read_table2("Data/data_rppa.csv", na = "0")#note: data are manually corrected to remove unexpected whitespaces causing parsing error
SKCM_RPPA_raw <-cbind(SKCM_RPPA_raw[,1],apply(SKCM_RPPA_raw[,-1],2,as.numeric))
SKCM_RPPA_raw[is.na(SKCM_RPPA_raw)] <- 0
SKCM_RPPA_raw_noNA_traw<- t(SKCM_RPPA_raw)
colnames(SKCM_RPPA_raw_noNA_traw) <-  SKCM_RPPA_raw_noNA_traw[1,]
SKCM_RPPA_raw_noNA_traw <- as.data.frame(SKCM_RPPA_raw_noNA_traw[-1,])
SKCM_RPPA <- SKCM_RPPA_raw_noNA_traw %>% rownames_to_column(var = "ID") %>% mutate_if(is.factor,list(~as.character(as.numeric)))
SKCM_RPPA[is.na(SKCM_RPPA)] <- 0
SKCM_RPPA.1 <- as.data.frame(cbind(SKCM_RPPA[,1],apply(SKCM_RPPA[,-1],2,as.numeric)))

#merge expression and RPPA
paralog_merged <- inner_join((SKCM_RNA_seq[,c("ID","DUSP4","DUSP6","YWHAE","YWHAZ","RAF1")]),SKCM_RPPA[,c("ID","\"MAP2K1|MEK1_pS217_S221","\"MAPK1_MAPK3|MAPK_pT202_Y204")],by="ID") %>% filter(ID %in% mutation_yes) 

# add BRAF & NRAS & NF1 mutation status
BRAF_mut <- paralog_merged$ID %in% BRAF_mut_tumor$Tumor_Sample_Barcode
NRAS_mut <- paralog_merged$ID %in% NRAS_mut_tumor$Tumor_Sample_Barcode
NF1_mut <- paralog_merged$ID %in% NF1_mut_tumor$Tumor_Sample_Barcode
mutation_status <- c()
mutation_status_nf1 <- c()
#===== Data Wrangling =====
#merge expression+RPPA+ mutation
paralog_merged.all <- paralog_merged %>% add_column(BRAF_mut) %>% add_column(NRAS_mut) %>% add_column(NF1_mut)
paralog_merged.all[,2:8] <- sapply(paralog_merged.all[,2:8], function(x){as.numeric(as.character(x))})
colnames(paralog_merged.all)[7:8] <- c("MEK1_pS117_S221_RPPA","MAPK_pT202_Y204")

mutation_status <- paralog_merged.all$BRAF_mut*5+paralog_merged.all$NRAS_mut*1
mutation_status <- plyr::mapvalues(mutation_status, from=c(0, 1, 5, 6), to=c("None","NRAS hotspot","BRAF hotspot","BRAF+NRAS hotspot"))

mutation_status_MAPK <-  paralog_merged.all$BRAF_mut*5+paralog_merged.all$NRAS_mut*1+paralog_merged.all$NF1_mut*2
mutation_status_MAPK <- plyr::mapvalues(mutation_status_MAPK, from=c(0, 1, 5, 6,2,7,3), to=c("No MAPK Hotspot","NRAS hotspot","BRAF hotspot","multiple hotspots","NF1 hotspot","multiple hotspots","multiple hotspots"))

#merged dataset with expression and mutation status
paralog_merged.all <- paralog_merged.all %>% add_column(mutation_status)%>%add_column(mutation_status_MAPK)
paralog_merged.all$mutation_status_MAPK <- 
  factor(paralog_merged.all$mutation_status_MAPK, 
         levels = c("No MAPK Hotspot","NRAS hotspot",
                    "BRAF hotspot","NF1 hotspot","multiple hotspots"))
paralog_merged.all%<>%
  select(DUSP4,DUSP6,mutation_status_MAPK)%>%
  filter(mutation_status_MAPK %in% c("No MAPK Hotspot","NRAS hotspot","BRAF hotspot"))%>%
  mutate(MutStatus = case_when(mutation_status_MAPK=="No MAPK Hotspot"~"MAPK WT",
                                          mutation_status_MAPK=="BRAF hotspot"~"BRAF MUT",
                                          mutation_status_MAPK=="NRAS hotspot"~"NRAS MUT")%>%
           factor(.,levels = c("MAPK WT","NRAS MUT","BRAF MUT")))%>%
  select(-mutation_status_MAPK)%>%
  gather(dusp,expr,-MutStatus)
#===== Plotting =====
paralog_merged.all%<>%unite("mutDUSP", c("dusp","MutStatus"),remove=F)
levelComb <- paralog_merged.all%>%select(dusp,MutStatus,mutDUSP)%>%
  distinct()%>%arrange(dusp,MutStatus)%>%pull(mutDUSP)
paralog_merged.all%<>%mutate(mutDUSP = factor(mutDUSP,levels = levelComb))
my_comparisons <- list(c("DUSP4_MAPK WT","DUSP4_NRAS MUT"),
                       c("DUSP4_MAPK WT","DUSP4_BRAF MUT"),
                       c("DUSP6_MAPK WT","DUSP6_NRAS MUT"),
                       c("DUSP6_MAPK WT","DUSP6_BRAF MUT"))
wil_box<- ggplot(paralog_merged.all, aes(y=log2(expr), x=mutDUSP,color=MutStatus)) +
  geom_rect(aes(xmin = 3.5,xmax=Inf,ymin=-Inf,ymax = 16),
            fill = "gray80", alpha = 0.02,color="transparent")+
  geom_boxplot(show.legend = FALSE)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth=0.05,aes(fill=MutStatus))+
  scale_color_manual(values=c("black","deepskyblue","firebrick"),labels = c("MAPK WT (n=77)",
                                                                            "NRAS MUT (n=78)",
                                                                            "BRAF MUT (n=137)"))+
  scale_fill_manual(values=c("black","deepskyblue","firebrick"),labels = c("MAPK WT (n=77)",
                                                                           "NRAS MUT (n=78)",
                                                                           "BRAF MUT (n=137)"))+
  stat_compare_means(label = "p.format",comparisons = my_comparisons,
                     label.y = rep(c(14.5,15.2),2))+ 
  labs(x="",y="mRNA expression log2(RSEM)")+
  scale_x_discrete(breaks = c("DUSP4_NRAS MUT",
                              "DUSP6_NRAS MUT"),
                   labels = c("DUSP4","DUSP6"))+
  scale_y_continuous(limits = c(7,18),breaks = c(7.5,10,12.5,15),
                     expand = expansion(mult = c(0.01,0)))+
  theme(
    plot.background=element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = .5, fill = NA),
    axis.ticks.x = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.position = c(0.001, 0.999),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.text=element_text(size=rel(0.8)),
    legend.direction = "vertical",
    legend.title = element_blank())+
  labs(title  = "TCGA melanoma")
ggsave(wil_box,filename = "Fig 3f TCGA DUSP4_6 expression/Fig3f_TCGA_mel_DUSP46_wilcox.pdf",
       dpi = 1000,width = 4, height = 6)


