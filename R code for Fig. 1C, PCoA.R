library(zCompositions)
library(ggplot2)
library(vegan)
library(ggrepel)

mydata <- read.csv('Table_S9.csv', skip=1, header=T)
otu <- mydata[,-(1:6)]
otu[1:10, 1:10]
dim(otu)
otu_t <- t(otu)
otu_t[1:5, 1:10]

otu_tr <- cmultRepl((otu_t), method="CZM", output="p-counts")
dim(otu_tr)
otu_tr[1:10, 1:8]

otu_hel <- decostand(otu_tr, "hellinger")
otu_bray <- vegdist(otu_hel, "bray")
pcoa <- cmdscale(otu_bray, k=3, eig=TRUE, add=TRUE)

points <- as.data.frame(pcoa$points)
eig = pcoa$eig

group_l <- substring(rownames(points), 1, 1)
points$group <- ifelse(group_l=="P", "JIA", "Control")

ggplot(points, aes(x=V1, y=V2, color=group)) + geom_point(alpha=.7, size=2) +
  scale_x_continuous(limits=c(-0.25, 0.23)) +
  scale_y_continuous(limits=c(-0.13, 0.15)) +
  guides(color=guide_legend("Groups", reverse=T)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=3), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=3), "%)", sep="")) +
  stat_ellipse(level=0.7) +
  theme(plot.margin=unit(c(5,5,5,5), "mm"), 
        panel.background=element_rect(fill="white"), 
        axis.line=element_line(), 
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15), 
        axis.text.x=element_text(color='black', size=12),
        axis.text.y=element_text(color='black', size=12),
        axis.ticks=element_line(color='black', size=1),
        legend.position=c(0.90, 0.93), 
        legend.title=element_text(color='black', size=15),
        legend.text=element_text(color='black', size=15))
