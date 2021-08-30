Model <- readRDS("Data/LatestModel.rds")
inputRep <- Model$Input
inputRep$replicate.map$samplename <- inputRep$replicate.map$colname
inputRep%<>%gemini::gemini_calculate_lfc(Input = ., sample.column.name = 'samplename')
LFCRepCors <- cor(inputRep$LFC, method = 'pearson')
ori_Name <- c("REP","Meljuso","TIGI1_004","TIPK1","TIMEL202_003","HS944")
sub_Name <- c("Rep","MELJUSO","GI1","PK1","MEL202","HS944T")
colnames(LFCRepCors)%<>%mgsub::mgsub(.,ori_Name,sub_Name)
rownames(LFCRepCors)%<>%mgsub::mgsub(.,ori_Name,sub_Name)
pdf(paste0("Fig S2a Paralog screen replicate correlation/FigS2a_Replicate_LFC_correlation_",Sys.Date(),".pdf"), width = 10, height = 9)
pheatmap::pheatmap(mat = LFCRepCors,
                   color = colorRampPalette(c("white","gray", "#FFAAAA", "#FF6666"))(100))
dev.off()