# Library distribution - ETP diagnostic plot
Paralog <- readRDS("Data/LatestModel.rds")
Paralog_pDNA <- as.vector(Paralog$Input$counts[,"pDNA"])
CDKO <- readRDS("Data/CDKO.rds")
CDKO_pDNA <- as.vector(CDKO$Input$counts[,"pDNA"])
ShenMali <- readRDS("Data/shen-mali.rds")
ShenMali_pDNA <- as.vector(ShenMali$Input$counts[,grep("T3", colnames(ShenMali$Input$counts))])
ZhaoMali <- readRDS("Data/zhao-mali.rds")
ZhaoMali_pDNA <- as.vector(ZhaoMali$Input$counts[,grep("d3", colnames(ZhaoMali$Input$counts))])
BigPapi <- readRDS("Data/BigPapi.rds")
BigPapi_pDNA <- as.vector(BigPapi$Input$counts[,"pDNA"])
ETP_counts <- list(
  `CDKO` = CDKO_pDNA,
  `Shen-Mali` = ShenMali_pDNA,
  `Zhao-Mali` = ZhaoMali_pDNA,
  `Big Papi` = BigPapi_pDNA,
  `Paralog - pDNA` = Paralog_pDNA
)

source("Utility_Functions/ETP_diagnostic_plot.R")
g = ETP_diagnostic_plot(ETP_counts,fix="Digenic Paralog")$g
g = g+
  scale_color_manual(values=c("black",
                              "red",
                              "blue",
                              "purple",
                              "darkgreen",
                              "darkgoldenrod1"),name="Combinatorial Screens")+
  theme(plot.title = element_text(color='black'),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(color='black', size = 14),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.key = element_rect(fill = "white")) + 
  guides(color = guide_legend(order = 1,
                              family="Helvetica",
                              override.aes = list(size=3.5,fill=NA)))+
  labs(x = "Fraction of library, ranked by abundance", y = "Cumulative fraction of reads")
ggsave(g,filename = "Fig S1d Library distribution/FigS1d_pDNA_distribution.pdf", width = 7, height = 7)
