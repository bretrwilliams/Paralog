require(dplyr)
require(ggplot2)
import::from(magrittr,"%<>%","set_colnames","set_rownames","set_names")
recombRates <- openxlsx::read.xlsx("Fig S1e Recombination Rate/Recombination rate.xlsx")[1:12,c(1,4)]%>%
  set_colnames(c("cl","frac"))%>%
  mutate(cl = gsub("\\_.*","",cl))%>%
  mutate(cl = factor(cl,levels = .$cl))
ggplot(recombRates,aes(x = cl, y = frac)) + 
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = expansion(mult=c(0,0)),limits = c(0,100))+
  labs(x="",y="Mismatch reads (%)")+
  theme(
    axis.text = element_text(face="bold",size=rel(0.8)),
    axis.ticks=element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 0.25, fill = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = paste0("Fig S1e Recombination Rate/","FigS1e_RecombinationRate_",Sys.Date(),".pdf"),
       width = 8,height = 4)