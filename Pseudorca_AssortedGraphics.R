## Network and Histogram graphics for Pseudorca social cluster work ##
## written by AEH, September 2022 for Cascadia Research Collective ##

# Data for the histograms/networks is generated in the Pseudorca_RcomCalculation_RCode script, also included in this repository #
# In this script, that data is plotted using ggplot2 #

# load up required packages #
library(ggplot2)
library(tidyverse)
library(lubridate)
library(forcats)
library(ggpattern)

## Histograms of Ncom from bootstrap replicates ##
# Load up data #
bootstrapdata <- read.csv("Histogram CSVs/Bootstrapped_ncoms_allmethods.csv")

# review data
summary(bootstrapdata)

# Basic code for each histogram #
# Just adjust parameters for each one as needed #
ggplot() + 
  geom_bar(data=bootstrapdata, aes(x=WT_ncom), stat="count", fill="darkgrey") +
  scale_x_continuous(breaks=c(1:9), expand=c(0,0), limits=c(0,10)) +
  scale_y_continuous(breaks=c(100, 200, 300, 400, 500, 600), limits=c(0,700), expand=c(0,0)) +
  theme_classic() + theme(axis.text=element_text(color="black", size=14), text=element_text(size=14)) +
  xlab("Number of communities in bootstrap replicates") + ylab("Count") +
  ggtitle("C") + theme(plot.title=element_text(hjust=-.1))
ggsave("Pc_Dist2PQ2_19992021_Seen5_WT_BootstrapNcomHistogram_v2.pdf", dpi=600)

# with standardized scales across all methods #
ggplot() + 
  geom_bar(data=bootstrapdata, aes(x=CL_ncom), stat="count", fill="darkgrey") +
  scale_x_continuous(breaks=c(5,10,15,20,25), expand=c(0,0), limits=c(0,26)) +
  scale_y_continuous(breaks=c(100, 200, 300, 400, 500,600, 700, 800, 900, 1000), limits=c(0,1025), expand=c(0,0)) +
  theme_classic() + theme(axis.text=element_text(color="black", size=14), text=element_text(size=14)) +
  xlab("Number of communities in bootstrap replicates") + ylab("Count")

ggsave("Pc_Dist2PQ2_19992021_Seen5_CL_BootstrapNcomHistogram_standardized.jpg", dpi=600)  

# stacked barplot #
# reformat data #
combinedhistogram <- data.frame(ncoms=c(1,3,4,5,6,7,8,9), 
                                count=c(1,8,499,290,151,33,15,3), 
                                method="Leading.eigenvector.community()")
EB <- data.frame(ncoms=c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,23,25), 
                 count=c(1,304,204,142,102,79,50,44,24,15,11,6,9,1,2,1,2,1,1,1), 
                 method="Edge.betweenness.community()")
combinedhistogram <- bind_rows(combinedhistogram, EB)
FG <- data.frame(ncoms=c(3,4,5,6,7,8), 
                count=c(10,792,146,39,11,2), 
                method="Fastgreedy.community()")
combinedhistogram <- bind_rows(combinedhistogram, FG)
WT <- data.frame(ncoms=c(4,5,6,7,8,9), 
                 count=c(639,245,87,15,7,7), 
                 method=c("Walktrap.community()"))
combinedhistogram <- bind_rows(combinedhistogram, WT)
LPC <- data.frame(ncoms=c(4,5,6,7,8,9,10), 
                  count=c(441,374,120,50,11,2,2), 
                  method="Label.propagation.community()")
combinedhistogram <- bind_rows(combinedhistogram, LPC)
CL <- data.frame(ncoms=c(4,5,6,7,8,9), 
                 count=c(660,269,54,14,2,1), 
                 method="Cluster_louvain()")
combinedhistogram <- bind_rows(combinedhistogram, CL)
# make the histogram #
ggplot() + geom_col(data=combinedhistogram, 
                    aes(x=ncoms,y=count, fill=method), 
                    position=position_dodge2(preserve = "single"), width=1) +
  scale_fill_viridis_d() +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0), limits=c(0,800)) +
  theme_classic() + theme(axis.text=element_text(color="black", size=14), legend.position = "NONE", text=element_text(size=14)) +
  xlab("Number of communities in bootstrap replicates") + ylab("Count")
ggsave("Pc_Dist2PQ2_19992021_Seen5_BootstrapNcomHistogram_Combined_withoutlegend.jpg", dpi=600)






## Formula and script to customize edge weights for weighted networks ##
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
edges$weight <- round_df(edges$weight, 1)

edgeslight_10 <- filter(edges, weight < 0.11)
edgeslight_20 <- filter(edges, weight < 0.21) %>%
  filter(weight > 0.10)
edgeslight_30 <- filter(edges, weight < 0.31) %>%
  filter(weight > 0.20)
edgeslight_40 <- filter(edges, weight < 0.41) %>%
  filter(weight > 0.30)
edgeslight_50 <- filter(edges, weight < 0.51) %>%
  filter(weight > 0.40)
edgeslight_60 <- filter(edges, weight < 0.61) %>%
  filter(weight > 0.50)
edgeslight_70 <- filter(edges, weight < 0.71) %>%
  filter(weight > 0.60)
edgeslight_80 <- filter(edges, weight < 0.81) %>%
  filter(weight > 0.70)
edgeslight_90 <- filter(edges, weight < 0.91) %>%
  filter(weight > 0.80)
edgeslight_100 <- filter(edges, weight < 1.01) %>%
  filter(weight > 0.90)


# Code for weighted edges for weighted networks with thin lines #
geom_segment(data=edgeslight_10, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +
  geom_segment(data=edgeslight_20, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.5, color="grey") +
  geom_segment(data=edgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="grey") +
  geom_segment(data=edgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="grey") +
  geom_segment(data=edgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="grey") +
  geom_segment(data=edgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="grey") +
  geom_segment(data=edgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="grey") +
  geom_segment(data=edgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="grey") +
  geom_segment(data=edgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="grey") +
  geom_segment(data=edgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="grey")
  
  
  
  
  
  

## Full Networks ##
# Each network is reproduced below, for ease of future reproduction #
# C1= Yellow circles, shape 21 ,#fde725 
# C2=Green triangles, shape 24, #35b779 
# C3/C5=Blue squares, shape 22, #31688e
# C4=Purple diamonds, shape 23, #440154
# Additional cluster 1 = white circles, shape 21, #FFFFFF
# Additional cluster 2 = black circles, shape 21, #000000
# for weighted networks, put size=weight in aes() call, for unweighted, put size=0.25 outside aes() call #
# upload nodes and edges #
nodes <- read.csv("Network CSVs/Pc_Dist2PQ2_19992021_Seen5_fullnetwork_nodecoords_clustermembership.csv")
edges <- read.csv("Network CSVs/Pc_Dist2PQ2_19992021_Seen5_fullnetwork_edgelist.csv")

# Break down cluster members into individual items #
nodes_LEC1 <- filter(nodes, leccluster==1)
nodes_LEC2 <- filter(nodes, leccluster==2)
nodes_LEC3 <- filter(nodes, leccluster==3)
nodes_LEC4 <- filter(nodes, leccluster==4)
nodes_LEC5 <- filter(nodes, leccluster==5)
nodes_EB1 <- filter(nodes, ebcluster==1)
nodes_EB2 <- filter(nodes, ebcluster==2)
nodes_EB3 <- filter(nodes, ebcluster==3)
nodes_EB4 <- filter(nodes, ebcluster==4)
nodes_EB5 <- filter(nodes, ebcluster==5)
nodes_EB6 <- filter(nodes, ebcluster==6)
nodes_FG1 <- filter(nodes, fgcluster==1)
nodes_FG2 <- filter(nodes, fgcluster==2)
nodes_FG3 <- filter(nodes, fgcluster==3)
nodes_FG4 <- filter(nodes, fgcluster==4)
nodes_WC1 <- filter(nodes, wccluster==1)
nodes_WC2 <- filter(nodes, wccluster==2)
nodes_WC3 <- filter(nodes, wccluster==3)
nodes_WC4 <- filter(nodes, wccluster==4)
nodes_LPC1 <- filter(nodes, lpccluster==1)
nodes_LPC2 <- filter(nodes, lpccluster==2)
nodes_LPC3 <- filter(nodes, lpccluster==3)
nodes_LPC4 <- filter(nodes, lpccluster==4)
nodes_CL1 <- filter(nodes, clcluster==1)
nodes_CL2 <- filter(nodes, clcluster==2)
nodes_CL3 <- filter(nodes, clcluster==3)
nodes_CL4 <- filter(nodes, clcluster==4)



# LEC full network plot #
ggplot() + 
  geom_segment(data=edges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y),size=0.25, color="grey") +
  geom_point(data=nodes_LEC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=nodes_LEC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_LEC3, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=nodes_LEC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  geom_point(data=nodes_LEC5, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("PC_Dist2PQ2_19992021_Seen5_LEC_FullNetwork_Unweighted.jpg", dpi=600)
# LEC full network plot - weighted #
ggplot() +
  geom_segment(data=edgeslight_10, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="#808B96") +
  geom_segment(data=edgeslight_20, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.5, color="#808B96") +
  geom_segment(data=edgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=edgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=edgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=edgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=edgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=edgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=edgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=edgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=nodes_LEC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=nodes_LEC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_LEC3, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=nodes_LEC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  geom_point(data=nodes_LEC5, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_LEC_FullNetwork_Weighted_v4.jpg", dpi=600)
# EB full network plot #
ggplot() + 
  geom_segment(data=edges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +
  geom_point(data=nodes_EB1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=nodes_EB2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_EB3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_EB4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=nodes_EB5, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=nodes_EB6, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_EB_FullNetwork_Unweighted.jpg", dpi=600)
# EB full network plot - weighted #
ggplot() +
  geom_segment(data=edgeslight_10, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="#808B96") +
  geom_segment(data=edgeslight_20, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.5, color="#808B96") +
  geom_segment(data=edgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=edgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=edgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=edgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=edgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=edgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=edgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=edgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=nodes_EB1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=nodes_EB2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_EB3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_EB4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=nodes_EB5, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=nodes_EB6, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_EB_FullNetwork_Weighted_v2.jpg", dpi=600)  
  
  # FG full network plot #
ggplot() +
  geom_segment(data=edges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y, size=weight), color="grey") +  
  geom_point(data=nodes_FG1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=nodes_FG2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_FG3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_FG4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_FG_FullNetwork_Unweighted.jpg", dpi=600)
# FG full network plot - weighted #
ggplot() +
  geom_segment(data=edgeslight_10, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="#808B96") +
  geom_segment(data=edgeslight_20, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.5, color="#808B96") +
  geom_segment(data=edgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=edgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=edgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=edgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=edgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=edgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=edgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=edgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=nodes_FG1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=nodes_FG2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_FG3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_FG4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_FG_FullNetwork_Weighted_v2.jpg", dpi=600)  
  
  # WC full network plot #
ggplot() +
  geom_segment(data=edges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +
  geom_point(data=nodes_WC1, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=nodes_WC2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_WC3, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_WC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_WC_FullNetwork_Unweighted.jpg", dpi=600)

# WC full network plot - weighted #
ggplot() +
  geom_segment(data=edgeslight_10, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="#808B96") +
  geom_segment(data=edgeslight_20, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.5, color="#808B96") +
  geom_segment(data=edgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=edgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=edgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=edgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=edgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=edgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=edgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=edgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=nodes_WC1, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=nodes_WC2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_WC3, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_WC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_WC_FullNetwork_Weighted_v2.jpg", dpi=600)
  
# LPC full network plot #
ggplot() +
  geom_segment(data=edges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +  
  geom_point(data=nodes_LPC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=nodes_LPC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_LPC3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_LPC4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")  
ggsave("Pc_Dist2PQ2_19992021_Seen5_LPC_FullNetwork_Unweighted.jpg", dpi=600)

# LPC full network plot - weighted #
ggplot() +
  geom_segment(data=edgeslight_10, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="#808B96") +
  geom_segment(data=edgeslight_20, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.5, color="#808B96") +
  geom_segment(data=edgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=edgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=edgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=edgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=edgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=edgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=edgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=edgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=nodes_LPC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=nodes_LPC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=nodes_LPC3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_LPC4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "") 
ggsave("Pc_Dist2PQ2_19992021_Seen5_LPC_FullNetwork_Weighted.jpg", dpi=600)

# CL full network plot #
ggplot() +
  geom_segment(data=edges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +  
  geom_point(data=nodes_CL1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=nodes_CL2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_CL3, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=nodes_CL4, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  theme_void() +
  theme(legend.position = "") 
ggsave("Pc_Dist2PQ2_19992021_Seen5_CL_FullNetwork_Unweighted.jpg", dpi=600)
# CL full network plot - weighted #
ggplot() +
  geom_segment(data=edgeslight_10, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="#808B96") +
  geom_segment(data=edgeslight_20, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.5, color="#808B96") +
  geom_segment(data=edgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=edgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=edgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=edgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=edgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=edgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=edgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=edgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=nodes_CL1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=nodes_CL2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=nodes_CL3, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=nodes_CL4, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  theme_void() +
  theme(legend.position = "") +
  ggtitle("A") + theme(plot.title=element_text(hjust=0.01))
ggsave("Pc_Dist2PQ2_19992021_Seen5_CL_FullNetwork_Weighted_A.pdf", dpi=600)






## Filtered Networks ##
# Each network is reproduced below, for ease of future reproduction #
# C1= Yellow circles, shape 21 ,#fde725 
# C2=Green triangles, shape 24, #35b779 
# C3/C5=Blue squares, shape 22, #31688e
# C4=Purple diamonds, shape 23, #440154
# Additional cluster 1 = white circles, shape 21, #FFFFFF
# Additional cluster 2 = black circles, shape 21, #000000
# for weighted networks, put size=weight in aes() call, for unweighted, put size=0.25 outside aes() call #
filterededges$weight <- round_df(filterededges$weight, 1)

filterededgeslight_30 <- filter(filterededges, weight < 0.31) %>%
  filter(weight > 0.20)
filterededgeslight_40 <- filter(filterededges, weight < 0.41) %>%
  filter(weight > 0.30)
filterededgeslight_50 <- filter(filterededges, weight < 0.51) %>%
  filter(weight > 0.40)
filterededgeslight_60 <- filter(filterededges, weight < 0.61) %>%
  filter(weight > 0.50)
filterededgeslight_70 <- filter(filterededges, weight < 0.71) %>%
  filter(weight > 0.60)
filterededgeslight_80 <- filter(filterededges, weight < 0.81) %>%
  filter(weight > 0.70)
filterededgeslight_90 <- filter(filterededges, weight < 0.91) %>%
  filter(weight > 0.80)
filterededgeslight_100 <- filter(filterededges, weight < 1.01) %>%
  filter(weight > 0.90)




# upload nodes and edges #
filterednodes <- read.csv("Network CSVs/Pc_Dist2PQ2_19992021_Seen5_filterednetwork_nodecoords_clustermembership.csv")
filterededges <- read.csv("Network CSVs/Pc_Dist2PQ2_19992021_Seen5_filterednetwork_edgelist.csv")

# Break down cluster members into individual items #
filterednodes_LEC1 <- filter(filterednodes, leccluster==1)
filterednodes_LEC2 <- filter(filterednodes, leccluster==2)
filterednodes_LEC3 <- filter(filterednodes, leccluster==3)
filterednodes_LEC4 <- filter(filterednodes, leccluster==4)
filterednodes_LEC5 <- filter(filterednodes, leccluster==5)
filterednodes_EB1 <- filter(filterednodes, ebcluster==1)
filterednodes_EB2 <- filter(filterednodes, ebcluster==2)
filterednodes_EB3 <- filter(filterednodes, ebcluster==3)
filterednodes_EB4 <- filter(filterednodes, ebcluster==4)
filterednodes_EB5 <- filter(filterednodes, ebcluster==5)
filterednodes_EB6 <- filter(filterednodes, ebcluster==6)
filterednodes_FG1 <- filter(filterednodes, fgcluster==1)
filterednodes_FG2 <- filter(filterednodes, fgcluster==2)
filterednodes_FG3 <- filter(filterednodes, fgcluster==3)
filterednodes_FG4 <- filter(filterednodes, fgcluster==4)
filterednodes_WC1 <- filter(filterednodes, wccluster==1)
filterednodes_WC2 <- filter(filterednodes, wccluster==2)
filterednodes_WC3 <- filter(filterednodes, wccluster==3)
filterednodes_WC4 <- filter(filterednodes, wccluster==4)
filterednodes_LPC1 <- filter(filterednodes, lpccluster==1)
filterednodes_LPC2 <- filter(filterednodes, lpccluster==2)
filterednodes_LPC3 <- filter(filterednodes, lpccluster==3)
filterednodes_LPC4 <- filter(filterednodes, lpccluster==4)
filterednodes_CL1 <- filter(filterednodes, clcluster==1)
filterednodes_CL2 <- filter(filterednodes, clcluster==2)
filterednodes_CL3 <- filter(filterednodes, clcluster==3)
filterednodes_CL4 <- filter(filterednodes, clcluster==4)

# LEC filtered network plot #
ggplot() + 
  geom_segment(data=filterededges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y, size=weight), color="grey") +
  geom_point(data=filterednodes_LEC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=filterednodes_LEC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_LEC3, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=filterednodes_LEC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  geom_point(data=filterednodes_LEC5, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("PC_Dist2PQ2_19992021_Seen5_LEC_FilteredNetwork_Unweighted.jpg", dpi=600)

# LEC filtered network plot - weighted #
ggplot() +
  geom_segment(data=filterededgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=filterededgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=filterededgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=filterededgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=filterededgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=filterededgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=filterededgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=filterededgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=filterednodes_LEC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=filterednodes_LEC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_LEC3, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=filterednodes_LEC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  geom_point(data=filterednodes_LEC5, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("PC_Dist2PQ2_19992021_Seen5_LEC_FilteredNetwork_Weighted_v2.jpg", dpi=600)

# EB filtered network plot #
ggplot() + 
  geom_segment(data=filterededges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y),size=0.25, color="grey") +
  geom_point(data=filterednodes_EB1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=filterednodes_EB2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_EB3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_EB4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=filterednodes_EB5, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=filterednodes_EB6, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  theme_void() +
  theme(legend.position = "")
ggsave("PC_Dist2PQ2_19992021_Seen5_EB_FilteredNetwork_Unweighted.jpg", dpi=600)

# EB filtered network plot - weighted #
ggplot() +
  geom_segment(data=filterededgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=filterededgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=filterededgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=filterededgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=filterededgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=filterededgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=filterededgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=filterededgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=filterednodes_EB1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=filterednodes_EB2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_EB3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_EB4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=filterednodes_EB5, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=filterednodes_EB6, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  theme_void() +
  theme(legend.position = "")
ggsave("PC_Dist2PQ2_19992021_Seen5_EB_FilteredNetwork_Weighted_v2.jpg", dpi=600)

# FG filtered network plot #
ggplot() +
  geom_segment(data=filterededges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y),size=0.25, color="grey") +  
  geom_point(data=filterednodes_FG1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=filterednodes_FG2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_FG3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_FG4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_FG_FilteredNetwork_Unweighted.jpg", dpi=600)

# FG filtered network plot - weighted #
ggplot() +
  geom_segment(data=filterededgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=filterededgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=filterededgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=filterededgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=filterededgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=filterededgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=filterededgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=filterededgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=filterednodes_FG1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=filterednodes_FG2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_FG3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_FG4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_FG_FilteredNetwork_Weighted_v2.jpg", dpi=600)

# WC filtered network plot #
ggplot() +
  geom_segment(data=filterededges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +
  geom_point(data=filterednodes_WC1, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=filterednodes_WC2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_WC3, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_WC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_WC_FilteredNetwork_Unweighted.jpg", dpi=600)

# WC filtered network plot - weighted #
ggplot() +
  geom_segment(data=filterededgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=filterededgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=filterededgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=filterededgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=filterededgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=filterededgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=filterededgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=filterededgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=filterednodes_WC1, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=filterednodes_WC2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_WC3, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_WC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_WC_FilteredNetwork_Weighted.jpg", dpi=600)


# LPC filtered network plot #
ggplot() +
  geom_segment(data=filterededges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +  
  geom_point(data=filterednodes_LPC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=filterednodes_LPC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_LPC3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_LPC4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_LPC_FilteredNetwork_Unweighted.jpg", dpi=600)

# LPC filtered network plot - weighted #
ggplot() +
  geom_segment(data=filterededgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=filterededgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=filterededgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=filterededgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=filterededgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=filterededgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=filterededgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=filterededgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=filterednodes_LPC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=filterednodes_LPC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=filterednodes_LPC3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_LPC4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_LPC_FilteredNetwork_Weighted.jpg", dpi=600)

# CL filtered network plot #
ggplot() +
  geom_segment(data=filterededges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +  
  geom_point(data=filterednodes_CL1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=filterednodes_CL2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_CL3, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=filterednodes_CL4, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  theme_void() +
  theme(legend.position = "") 
ggsave("Pc_Dist2PQ2_19992021_Seen5_CL_FilteredNetwork_Unweighted.jpg", dpi=600)

# cL filtered network plot - weighted #
ggplot() +
  geom_segment(data=filterededgeslight_30, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.75, color="#808B96") +
  geom_segment(data=filterededgeslight_40, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1, color="#808B96") +
  geom_segment(data=filterededgeslight_50, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.25, color="#808B96") +
  geom_segment(data=filterededgeslight_60, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.5, color="#808B96") +
  geom_segment(data=filterededgeslight_70, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=1.75, color="#808B96") +
  geom_segment(data=filterededgeslight_80, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2, color="#808B96") +
  geom_segment(data=filterededgeslight_90, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.25, color="#808B96") +
  geom_segment(data=filterededgeslight_100, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=2.5, color="#808B96") +
  geom_point(data=filterednodes_CL1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=filterednodes_CL2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=filterednodes_CL3, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=filterednodes_CL4, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  theme_void() +
  theme(legend.position = "") +
  ggtitle("B") + theme(plot.title=element_text(hjust=0.01))
ggsave("Pc_Dist2PQ2_19992021_Seen5_CL_FilteredNetwork_Weighted_B.pdf", dpi=600)








## Probability graphs ##
# Each network is reproduced below, for ease of future reproduction #
# can copy color assignments from above graphs #
# C1= Yellow circles, shape 21 ,#fde725 
# C2=Green triangles, shape 24, #35b779 
# C3/C5=Blue squares, shape 22, #31688e
# C4=Purple diamonds, shape 23, #440154
# Additional cluster 1 = white circles, shape 21, #FFFFFF
# Additional cluster 2 = black circles, shape 21, #000000

# LEC #
# upload nodes and edges #
LECnodes <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_LEC_probabilitynodes.csv")
LECedges <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_LEC_probabilityedgelist.csv")
# Break down cluster members into individual items #
LECnodes_LEC1 <- filter(LECnodes, membership==1)
LECnodes_LEC2 <- filter(LECnodes, membership==2)
LECnodes_LEC3 <- filter(LECnodes, membership==3)
LECnodes_LEC4 <- filter(LECnodes, membership==4)
LECnodes_LEC5 <- filter(LECnodes, membership==5)
# probability graph #
ggplot() + 
  geom_segment(data=LECedges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y),size=0.25, color="grey") +
  geom_point(data=LECnodes_LEC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=LECnodes_LEC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=LECnodes_LEC3, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=LECnodes_LEC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  geom_point(data=LECnodes_LEC5, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_LEC_Probabilitygraph.jpg", dpi=600)

# EB #
# upload nodes and edges #
EBnodes <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_EB_probabilitynodes.csv")
EBedges <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_EB_probabilityedgelist.csv")
# Break down cluster members into individual items # 
EBnodes_EB1 <- filter(EBnodes, membership==1)
EBnodes_EB2 <- filter(EBnodes, membership==2)
EBnodes_EB3 <- filter(EBnodes, membership==3)
EBnodes_EB4 <- filter(EBnodes, membership==4)
EBnodes_EB5 <- filter(EBnodes, membership==5)
EBnodes_EB6 <- filter(EBnodes, membership==6)
# probability graph #
ggplot() +
  geom_segment(data=EBedges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y),size=0.25, color="grey") +
  geom_point(data=EBnodes_EB1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=EBnodes_EB2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=EBnodes_EB3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=EBnodes_EB4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=EBnodes_EB5, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="black") +
  geom_point(data=EBnodes_EB6, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="white") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_EB_Probabilitygraph.jpg", dpi=600)

# FG #
# upload nodes and edges #
FGnodes <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_FG_probabilitynodes.csv")
FGedges <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_FG_probabilityedgelist.csv")
# Break down cluster members into individual items # 
FGnodes_FG1 <- filter(FGnodes, membership==1)
FGnodes_FG2 <- filter(FGnodes, membership==2)
FGnodes_FG3 <- filter(FGnodes, membership==3)
FGnodes_FG4 <- filter(FGnodes, membership==4)
# probability graph #
ggplot() +
  geom_segment(data=FGedges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +  
  geom_point(data=FGnodes_FG1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=FGnodes_FG2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=FGnodes_FG3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=FGnodes_FG4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "") + coord_flip() + scale_y_reverse() +ggtitle("B") + theme(plot.title=element_text(hjust=0.01))
ggsave("Pc_Dist2PQ2_19992021_Seen5_FG_Probabilitygraph.pdf", dpi=600)

# WC #
# upload nodes and edges #
WCnodes <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_WC_probabilitynodes.csv")
WCedges <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_WC_probabilityedgelist.csv")
# adjust projection of graph #
WCnodes$X2 <- WCnodes$X2*-1
WCedges$from.y <- WCedges$from.y*-1
WCedges$to.y <- WCedges$to.y*-1
# Break down cluster members into individual items # 
WCnodes_WC1 <- filter(WCnodes, membership==1)
WCnodes_WC2 <- filter(WCnodes, membership==2)
WCnodes_WC3 <- filter(WCnodes, membership==3)
WCnodes_WC4 <- filter(WCnodes, membership==4)
# probability graph #
ggplot() +
  geom_segment(data=WCedges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +
  geom_point(data=WCnodes_WC1, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=WCnodes_WC2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=WCnodes_WC3, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=WCnodes_WC4, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  theme_void() +
  theme(legend.position = "") + ggtitle("C") + theme(plot.title=element_text(hjust=0.01))
ggsave("Pc_Dist2PQ2_19992021_Seen5_WC_Probabilitygraph_v2.pdf", dpi=600)

# LPC #
# upload nodes and edges #
LPCnodes <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_LPC_probabilitynodes.csv")
LPCedges <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_LPC_probabilityedgelist.csv")
# Break down cluster members into individual items #
LPCnodes_LPC1 <- filter(LPCnodes, membership==1)
LPCnodes_LPC2 <- filter(LPCnodes, membership==2)
LPCnodes_LPC3 <- filter(LPCnodes, membership==3)
LPCnodes_LPC4 <- filter(LPCnodes, membership==4)
# probability graph #
ggplot() +
  geom_segment(data=LPCedges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +
  geom_point(data=LPCnodes_LPC1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") + 
  geom_point(data=LPCnodes_LPC2, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  geom_point(data=LPCnodes_LPC3, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=LPCnodes_LPC4, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  theme_void() +
  theme(legend.position = "")
ggsave("Pc_Dist2PQ2_19992021_Seen5_LPC_Probabilitygraph.jpg", dpi=600)

# CL #
# upload nodes and edges #
CLnodes <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_CL_probabilitynodes.csv")
CLedges <- read.csv("Probability Graph CSVs/Pc_Dist2PQ2_19992021_Seen5_CL_probabilityedgelist.csv")
# adjust projection of graph #
CLnodes$X1 <- CLnodes$X1*-1
CLedges$from.x <- CLedges$from.x*-1
CLedges$to.x <- CLedges$to.x*-1
# break down cluster members into individual items #
CLnodes_CL1 <- filter(CLnodes, membership==1)
CLnodes_CL2 <- filter(CLnodes, membership==2)
CLnodes_CL3 <- filter(CLnodes, membership==3)
CLnodes_CL4 <- filter(CLnodes, membership==4)
# probability graph #
ggplot() +
  geom_segment(data=CLedges, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.25, color="grey") +  
  geom_point(data=CLnodes_CL1, aes(x=X1, y=X2), shape=21, size=4, color="black", fill="#fde725") +
  geom_point(data=CLnodes_CL2, aes(x=X1, y=X2), shape=24, size=4, color="black", fill="#35b779") +
  geom_point(data=CLnodes_CL3, aes(x=X1, y=X2), shape=22, size=4, color="black", fill="#31688e") +
  geom_point(data=CLnodes_CL4, aes(x=X1, y=X2), shape=23, size=4, color="black", fill="#440154") +
  theme_void() +
  theme(legend.position = "") + ggtitle("A") + theme(plot.title=element_text(hjust=0.01))
ggsave("Pc_Dist2PQ2_19992021_Seen5_CL_Probabilitygraph_v2.pdf", dpi=600)  
  
  


# Distribution of IDs by Island/Year/Cluster #
# load up data file #
sights <- read.csv("PseudorcaData_CorrectedNumSights_Raw.csv")
sights <- filter(sights, timesseen >4)
sights$Date <- as.POSIXct(sights$Date, format="%m/%d/%Y")
sights$Year <- year(sights$Date)
sights$Cluster <- nodes$clcluster[match(sights$ID, nodes$names)]
sightsunique <- sights %>%
  group_by(Year) %>%
  distinct(ID, .keep_all=TRUE) %>%
  ungroup() 
sightsunique <- sightsunique[order(sightsunique$Year, -sightsunique$Cluster),]
sightsunique$Cluster <- as.factor(sightsunique$Cluster) 
# make barplot #
ggplot() +
  geom_bar(data=sightsunique, aes(x=Year, fill=forcats::fct_rev(Cluster)), position="stack") + 
  scale_fill_manual(values=c("#440154", "#31688e", "#35b779", "#fde725")) + 
  theme(legend.position="top", legend.title = element_blank()) + theme_classic() +
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(breaks=c(1999:2021), expand=c(0,0)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  theme(axis.text=element_text(color="black", size=12), legend.position = "top", legend.title=element_blank()) +
  xlab("Year") + ylab("Number of distinctive individuals") + theme(text=element_text(color="black", size=12))
ggsave("Pc_Dist2PQ2_19992021_Seen5_HistogramDistIndsvsYear.jpg", dpi=600) 
# Histogram of number of identifications vs year, no sight restrictions (Figure 1) #
sights <- read.csv("PseudorcaData_CorrectedNumSights_Raw.csv")
sights$Date <- as.POSIXct(sights$Date, format="%m/%d/%Y")
sights$Year <- year(sights$Date)
ggplot() +
  geom_bar(data=sights, aes(x=Year), color="black", fill="grey") + 
  theme(legend.position="none") + theme_classic() +
  scale_y_continuous(expand=c(0,0), breaks=c(0, 50, 100, 150, 200, 250, 300)) + scale_x_continuous(breaks=c(1999:2021), expand=c(0,0)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  theme(axis.text=element_text(color="black", size=12)) +
  xlab("Year") + ylab("Number of identifications") + theme(text=element_text(color="black", size=12)) +
  ggtitle("B") + theme(plot.title=element_text(hjust=-.09))
ggsave("Pc_Dist2PQ2_19992021_nosightrestrictions_HistogramYearDistIdentifications_noclusters_B.pdf", dpi=600)




## make barplot of Island vs cluster #
# first we need to redo our data #
sights <- read.csv("PseudorcaData_CorrectedNumSights_Raw.csv")
sights <- filter(sights, timesseen >4)
sights$Cluster <- nodes$clcluster[match(sights$ID, nodes$names)]
sights$isl_area <- 
  dplyr::recode(sights$Island,
                Hawaii = "Hawai'i",
                Maui = "Maui Nui",
                Lanai = "Maui Nui",
                MauiKah = "Maui Nui",
                MauiNui = "Maui Nui",
                Kahoolawe = "Maui Nui",
                Kauai = "Kaua'i",
                Molokai = "Maui Nui",
                Oahu = "O'ahu")
sightsunique <- sights %>%
  group_by(isl_area) %>%
  distinct(ID, .keep_all=TRUE) %>%
  ungroup() 
sightsunique <- sightsunique[order(sightsunique$isl_area, -sightsunique$Cluster),]
sightsunique$Cluster <- as.factor(sightsunique$Cluster)
sightsunique$isl_area <- factor(sightsunique$isl_area, levels=c("Kaua'i", "O'ahu", "Maui Nui", "Hawai'i"))
ggplot() +
  geom_bar(data=sightsunique, aes(x=isl_area, fill=forcats::fct_rev(Cluster)), position="stack") + 
  scale_fill_manual(values=c("#440154", "#31688e", "#35b779", "#fde725")) + 
  theme(legend.position="top", legend.title = element_blank()) + theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text=element_text(color="black", size=12), legend.position = "top", legend.title=element_blank()) +
  xlab("Year") + ylab("Number of distinctive individuals (1999-2021)") + theme(text=element_text(color="black", size=12))
ggsave("Pc_Dist2PQ2_19992021_Seen5_HistogramDistIndsvsIsland.jpg", dpi=600)

# Histogram of identifications for island area vs cluster #
sights <- read.csv("PseudorcaData_CorrectedNumSights_Raw.csv")
sights <- filter(sights, timesseen >4)
sights$Cluster <- nodes$clcluster[match(sights$ID, nodes$names)]
sights$isl_area <- 
  dplyr::recode(sights$Island,
                Hawaii = "Hawai'i",
                Maui = "Maui Nui",
                Lanai = "Maui Nui",
                MauiKah = "Maui Nui",
                MauiNui = "Maui Nui",
                Kahoolawe = "Maui Nui",
                Kauai = "Kaua'i",
                Molokai = "Maui Nui",
                Oahu = "O'ahu")
sights <- sights[order(sights$isl_area, sights$Cluster),]
sights$Cluster <- as.factor(sights$Cluster)
sights$isl_area <- factor(sights$isl_area, levels=c("Kaua'i", "O'ahu", "Maui Nui", "Hawai'i"))
ggplot() +
  geom_bar_pattern(data=sights, aes(x=isl_area, fill=forcats::fct_rev(Cluster), pattern=Cluster),color="black", position="stack") + 
  scale_fill_manual(values=c("#440154", "#31688e", "#35b779", "#fde725")) + 
  theme(legend.position="top", legend.title = element_blank()) + theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text=element_text(color="black", size=12), legend.position = "top", legend.title=element_blank()) +
  guides(fill = guide_legend(reverse=TRUE, override.aes = list(pattern=c("none","stripe","none","circle")))) +
  guides(pattern = "none") +
  xlab("") + ylab("Number of identifications") + theme(text=element_text(color="black", size=12)) +
  scale_pattern_manual(values=c("none","stripe","none","circle")) + ggtitle("A") + theme(plot.title=element_text(hjust=-.08))
ggsave("Pc_Dist2PQ2_19992021_Seen5_HistogramDistIdentificationsvsIsland_v6.pdf", dpi=600)




# Histogram of identifications for island area vs cluster, for individuals Dist2+PQ2+, #
# WITHOUT restrictions on the number of times seen #
sights <- read.csv("PseudorcaData_CorrectedNumSights_Raw.csv")
fewerrestrictionnodes <- read.csv("Pc_Dist2PQ2_Nosightrestrictions_CLcommunityassignments.csv")
sights$Cluster <- fewerrestrictionnodes$Louvain[match(sights$ID, fewerrestrictionnodes$ID)]
sights$isl_area <- 
  dplyr::recode(sights$Island,
                Hawaii = "Hawai'i",
                Maui = "Maui Nui",
                Lanai = "Maui Nui",
                MauiKah = "Maui Nui",
                MauiNui = "Maui Nui",
                Kahoolawe = "Maui Nui",
                Kauai = "Kaua'i",
                Molokai = "Maui Nui",
                Oahu = "O'ahu")
sights <- sights[order(sights$isl_area, sights$Cluster),]
sights$Cluster <- as.factor(sights$Cluster)
sights$isl_area <- factor(sights$isl_area, levels=c("Kaua'i", "O'ahu", "Maui Nui", "Hawai'i"))
ggplot() +
  geom_bar_pattern(data=sights, aes(x=isl_area, fill=forcats::fct_rev(Cluster), pattern=Cluster), color="black", position="stack") + 
  scale_fill_manual(values=c("#440154", "#31688e", "#35b779", "#fde725")) + 
  theme(legend.position="top", legend.title = element_blank()) + theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text=element_text(color="black", size=12), legend.position = "none", legend.title=element_blank()) +
  guides(fill = guide_legend(reverse=TRUE, override.aes = list(pattern=c("none","stripe","none","circle")))) +
  guides(pattern = "none") +
  xlab("") + ylab("Number of identifications") + theme(text=element_text(color="black", size=12)) +
  scale_pattern_manual(values=c("none","stripe","none","circle")) + ggtitle("B") +
  theme(plot.title=element_text(hjust=-.08))
ggsave("Pc_Dist2PQ2_19992021_Nosightrestrictions_HistogramDistIdentificationsvsIsland_v3.pdf", dpi=600)
# Histogram of identifications for island area, without pattern/color for cluster - Figure 1 #
ggplot() +
  geom_bar(data=sights, aes(x=isl_area),color="black", fill="grey") + 
  theme(legend.position="none") + theme_classic() +
  ggtitle("A") +
  scale_y_continuous(expand=c(0,0), breaks=c(250, 500, 750, 1000, 1250)) +
  theme(axis.text=element_text(color="black", size=12)) +
  xlab("") + ylab("Number of identifications") + theme(text=element_text(color="black", size=12)) +
  theme(plot.title=element_text(hjust=-.1))
ggsave("Pc_Dist2PQ2_19992021_Nosightrestrictions_HistogramDistIdentifications_noclusters_A.pdf", dpi=600)
