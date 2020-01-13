###Load the packages
library(ggplot2)
library(gridExtra)
library(ggsignif)  #This package is used for the comparisons between two groups.
library(vegan)
library(coin)  #This package is used for permutation tests.

###Load and pre-process the data
mydata <- read.csv('Table_S9.csv', skip=1, header=T)
mydata <- mydata[,-c(1:6)]  #Delete taxonomy
dim(mydata)
otu1 <- t(mydata)
otu1[1:10, 1:10]
"PBJ70" %in% rownames(otu1)
otu <- data.frame(trunc(otu1/rowSums(otu1)*1000000))
rowSums(otu)
otu[1:10, 1:10]

###Calculate diversity indices
S <- specnumber(otu)
estimateR(otu)
chao1 <- estimateR(otu)[2,]
shannon <- diversity(otu, "shannon")
simpson <- diversity(otu, "simpson")
pielou <- shannon/log(S)
alpha <- data.frame(chao1, shannon, simpson, pielou)

####Grouping factor
group_letter <- substring(rownames(otu), 1, 1)
group <- ifelse(group_letter=="P", "JIA", "Control")
alpha$group <- group
alpha$group <- as.factor(alpha$group)
str(alpha)
alpha$group <- factor(alpha$group, levels=c("JIA", "Control"))
#write.csv(alpha, file="Alpha_diversities.csv")

###Compare differences in α diversity indices between the two groups
#Chao1
wilcox.test(chao1~group, alpha)  #Wilcoxon rank sum test，P-value=0.002647
set.seed(100)
oneway_test(chao1~group, data=alpha, distribution=approximate(nresample=10000))  #Permutation test, P-value=0.0021
#shannon
wilcox.test(shannon~group, alpha)  #Wilcoxon rank sum test，P-value=0.03176
set.seed(100)
oneway_test(shannon~group, data=alpha, distribution=approximate(nresample=10000))  #Permutation test, P-value=0.0383
#simpson
wilcox.test(simpson~group, alpha)  #Wilcoxon rank sum test，P-value=0.2484
#pielou
wilcox.test(pielou~group, alpha)  #Wilcoxon rank sum test，P-value=0.1159

###Visualize results
#Chao1 Index
p1 <- ggplot(data=alpha, aes(x=group, y=chao1, fill=group)) + geom_boxplot(width=0.7) +
  scale_fill_brewer(palette='Set2', direction=-1) +
  scale_y_continuous(name='Chao1 Index', limits=c(100, 700), breaks=c(seq(100,700,200))) + 
  guides(fill=FALSE) +
  geom_signif(annotations="P=0.0026", y_position=600, xmin=1, xmax=2, tip_length=0.045, textsize=5) +
  geom_jitter(width=0.3, alpha=0.5, color="orange") +
  theme(panel.background=element_rect(fill="white"),
        axis.title.x=element_blank(), 
        axis.title.y=element_text(size=15), 
        axis.line=element_line(), 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(color='black', size=13),
        axis.ticks=element_line(color='black', size=1)) 
#Shannon Index
p2 <- ggplot(data=alpha, aes(x=group, y=shannon, fill=group)) + geom_boxplot(width=0.7) +
        scale_fill_brewer(palette='Set2', direction=-1) +
        scale_y_continuous(name='Shannon Index', limits=c(1.5, 4.8), breaks=c(seq(1.5,4.8,0.5))) + 
        guides(fill=FALSE) +
        geom_signif(annotations="P=0.0317", y_position=4.4, xmin=1, xmax=2, tip_length=0.045, textsize=5) +
        geom_jitter(width=0.3, alpha=0.5, color="orange") +
        theme(panel.background=element_rect(fill="white"), 
          axis.title.x=element_blank(), 
          axis.title.y=element_text(size=15),
          axis.line=element_line(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(color='black', size=13),
          axis.ticks=element_line(color='black', size=1))
#Simpson Index
p3 <- ggplot(data=alpha, aes(x=group, y=simpson, fill=group)) + geom_boxplot(width=0.7) +
        geom_jitter(width=0.3, alpha=0.5, color="orange") +
        scale_fill_brewer(palette='Set2', direction=-1) +
        guides(fill=guide_legend("Groups")) + 
        scale_y_continuous(name='Simpson Index', limits=c(0.8, 1.0), breaks=c(seq(0.8,1.0,0.05))) + 
        geom_signif(annotations="P=0.2484", y_position=0.975, xmin=1, xmax=2, tip_length=0.015, textsize=5) +
          theme(panel.background=element_rect(fill="white"),
          axis.title.x=element_blank(), 
          axis.title.y=element_text(size=15), 
          axis.line=element_line(), 
          axis.text.x=element_blank(), 
          axis.text.y=element_text(color='black', size=13), 
          axis.ticks=element_line(color='black', size=1),
          legend.title=element_text(color='black', size=15),
          legend.text=element_text(color='black', size=15))
#Pielou's Index
p4 <- ggplot(data=alpha, aes(x=group, y=pielou, fill=group)) + geom_boxplot(width=0.7) +
  geom_jitter(width=0.3, alpha=0.5, color="orange") +
  scale_fill_brewer(palette='Set2', direction=-1) +
  guides(fill=guide_legend("Groups")) + 
  scale_y_continuous(name="Pielou's Index", limits=c(0.4, 0.75), breaks=c(seq(0.4,0.75,0.1))) + 
  geom_signif(annotations="P=0.1159", y_position=0.72, xmin=1, xmax=2, tip_length=0.015, textsize=5) +
  theme(panel.background=element_rect(fill="white"),
        axis.title.x=element_blank(), 
        axis.title.y=element_text(size=15), 
        axis.line=element_line(), 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(color='black', size=13), 
        axis.ticks=element_line(color='black', size=1),
        legend.title=element_text(color='black', size=15),
        legend.text=element_text(color='black', size=15))


grid.arrange(p1, p2, p3, nrow=1, ncol=3, widths=c(1, 1, 1.6))



###Descriptive statistics and normal distribution assumptions of shannon and simpson
#Subsets
p_group <- subset(alpha_data, group=='JIA')
c_group <- subset(alpha_data, group=='Control')
#Descriptive statistics of shannon and simpson
summary(p_group$shannon)
summary(c_group$shannon)
summary(p_group$simpson)
summary(c_group$simpson)
#Histograms
par(mfrow=c(2,2))
hist(p_group$shannon, breaks=seq(1, 5, 0.5))
hist(c_group$shannon, breaks=seq(1, 5, 0.5))
hist(p_group$simpson, breaks=seq(0.4, 1, 0.1))
hist(c_group$simpson, breaks=seq(0.4, 1, 0.1))
dev.off()
#Assumptions of normal distribution
#shannon
shapiro.test(p_group$shannon)
shapiro.test(c_group$shannon)
par(mfrow=c(1,2))
qqnorm(p_group$shannon)
qqline(p_group$shannon)
qqnorm(c_group$shannon)
qqline(c_group$shannon)
dev.off()
#simpson
shapiro.test(p_group$simpson)
shapiro.test(c_group$simpson)
par(mfrow=c(1,2))
qqnorm(p_group$simpson)
qqline(p_group$simpson)
qqnorm(c_group$simpson)
qqline(c_group$simpson)
dev.off()
