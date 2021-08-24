#' Diagnostic Plot for pDNA counts
#'@author Sidharth Jain
#'@param ETP_counts either a numeric vector or list of numeric vectors containing the total pDNA read counts.
#'If the list is named, each named element will be plotted as a separate line.  See description for usage.
#'@description
# ETP <- read.delim("~/Downloads/scores-JD-AACK-redo_BT9TR-full.txt/scores-JD-AACK-redo_BT9TR-full.txt", sep = '\t', header = T, stringsAsFactors = F)
# Redundome_ETP <- rowSums(ETP[,c(-1,-2)])
# 
# Papi_ETP <- read.delim("Papi/Supplementary Tables/Supplementary Table 3 SynLet screening data.txt")
# Papi_ETP <- Papi_ETP[Papi_ETP$Cell.Line==Papi_ETP$Cell.Line[1] & Papi_ETP$Time.Point=="Day 21", "pDNA.Reads"]
# 
# ETP_diagnostic_plot(list(Papi = Papi_ETP, Redundome = Redundome_ETP))

ETP_diagnostic_plot <- function(ETP_counts,fix = NULL){
	library(ggplot2)
	library(tidyverse)
	if(is.list(ETP_counts)){
		df <- pbmcapply::pbmclapply(1:length(ETP_counts), function(i){
		   	c <- ETP_counts[[i]]
		   	name <- names(ETP_counts)[i]
		   	etp <- sort(c, decreasing = T)
		   	cumfrac <- cumsum(etp) / sum(etp)
		   	ranketp <- (1:length(etp)) / length(etp)
		   	dat = data.frame(cumulative_fraction = cumfrac,
		   					relative_rank = ranketp,
		   					screen = name)
	   		return(dat)
	   },mc.cores = 8) %>%
	   	bind_rows()
	}else if(is.vector(ETP_counts)){
		etp <- sort(ETP_counts, decreasing = T)
		cumfrac <- cumsum(etp) / sum(etp)
		ranketp <- (1:length(etp)) / length(etp)
		df = data.frame(cumulative_fraction = cumfrac,
						relative_rank = ranketp,
						screen = "")
	}
	
	else{
		stop("Unable to process ETP_counts in current format.")
	}
	
	df.control = data.frame(cumulative_fraction = seq(0,1,0.001),
							relative_rank = seq(0,1,0.001),
							screen = "Optimal Screen")
	
	df = rbind(df.control, df)
	ordered_screen <- df%>%
	  filter(cumulative_fraction>=0.9 & cumulative_fraction<0.9+1e-4)%>%
	  group_by(screen)%>%
	  arrange(relative_rank)%>%
	  slice(1)%>%
	  arrange(-relative_rank)%>%
	  pull(screen)
	fixed_screen <- c(c("Optimal Screen",fix),ordered_screen[!ordered_screen %in% c("Optimal Screen",fix)])
	g = ggplot(data = df%>%
	             mutate(screen = factor(screen,levels = fixed_screen)), 
	           aes(x = relative_rank, y = cumulative_fraction, color = screen)) +
		geom_smooth(size=1)
	return(list(g=g,screenDistri = ordered_screen))
}
