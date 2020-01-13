library(randomForest)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggsignif)

###Load and pre-process data
mydata <- read.csv('Table_S9.csv', skip=1, header=T)
mydata[1:5, 1:10]
genera <- mydata[,-c(1,2,3,4,6)]  #Delete taxonomy except Genus
genera_sum <- aggregate(genera[ ,-1], by=list(genera[ ,1]), sum)
rownames(genera_sum) <- genera_sum[, 1]
otu_g <- genera_sum[,-1]
otu.t <- data.frame(t(otu_g))
df <- otu.t/rowSums(otu.t)
rowSums(df)
group_letter <- substring(rownames(df), 1, 1)
g <- ifelse(group_letter=="P", "JIA", "Control")
df$aa <- g
df <- df[,order(colnames(df), decreasing=F)]
str(df)
df$aa <- as.factor(df$aa)
str(df)
df[1:5, 1:10]

###Select an optimal mtry, and then perform random forest analysis using all microbiota.
#Select an optimal mtry
set.seed(100)
m.try <- tuneRF(df[-1], df$aa, ntreeTry=500, stepFactor=1.5, improve=0.01, trace=TRUE)
m.try[m.try[,2] == min(m.try[,2]), 1]  #The optimal mtry is 9 for genera; 5 for family; 4 for order; 2 for class; 2 for phylum

#Perform random forest analysis using all microbiota
set.seed(95)
rf <- randomForest(aa~., data=df, ntree=500, mtry=2, na.action=na.omit, importance=TRUE)
rf  #OOB: OTU(34.57%); Genera(32.1%); Family(39.51%); Order(40.74%); Class(38.27%); Phylum(53.09%)

###Pick out the first 12 genera, and then calculate the medians and IQRs of abundance in the two groups
im <- as.data.frame(rf$importance)
im <- im[order(im$MeanDecreaseGini, decreasing=T),]
im <- im[1:12,]
d <- df[,c("aa", rownames(im))]

sub_con <- subset(d, aa=="Control")
sub_jia <- subset(d, aa=="JIA")
con <- as.data.frame(apply(sub_con[,-1], 2, median)*100)
jia <- as.data.frame(apply(sub_jia[,-1], 2, median)*100)
per <- merge.data.frame(jia, con, by="row.names")
rownames(per) <- per[,1]
per <- per[,-1]
colnames(per) <- c("jia", "con")
subset(per, per$jia > per$con)
subset(per, per$jia < per$con)

d$aa <- factor(d$aa, levels=c("JIA", "Control"))
zw <- apply(d[,-1], 2, median)
sort(zw, decreasing=T)

md <- melt(d, id=c("aa"))  #reshape the data
fac <- names(sort(tapply(md$value, md$variable, IQR), decreasing=T))
md$variable <- factor(md$variable, levels=fac)

###visualization
ggplot(data=md, aes(x=variable, y=value*100, fill=aa)) + geom_boxplot(width=0.7, outlier.color="white") +
  scale_y_continuous(name='Relative abundance (%)', limits=c(0,8), breaks=c(seq(0,8,2))) +
  geom_signif(annotations = c("*","*","*","*"), y_position=c(6.1, 6.8, 4.7, 0.5),
              xmin=c(1.8, 2.8, 3.8, 11.8), xmax=c(2.2, 3.2, 4.2, 12.2), tip_length=0.005)+
  guides(fill=guide_legend('Groups')) +
  theme(panel.background=element_rect(fill="white"),
        plot.margin=unit(c(0.5,0.5,0.5,1.5), "cm"),
        axis.title.y=element_text(size=16, vjust=5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color='black', size=13, angle=1, hjust=-15),
        axis.text.x=element_text(color='black', size=13, angle=30, hjust=1, vjust=1),
        axis.line=element_line(),
        axis.ticks=element_line(color='black', size=1),
        legend.position=c(0.82, 0.82),
        legend.background=element_rect(color='black', linetype=2),
        legend.title=element_text(color='black', size=13),
        legend.text=element_text(color='black', size=13))


